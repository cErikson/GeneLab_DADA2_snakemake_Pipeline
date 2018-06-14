##	++======================================================++
##	|| DADA2 Analysis script for marker gene GLDS datasets  ||
##	++======================================================++
##
##	Christian Erikson: Bio.Erkson@gmail.com christian.b.erikson@nasa.gov
##

library(tidyverse)
library(dada2)
library(phyloseq)

# Log errors/ messages to the snakemake log
log_file=file(snakemake@log[[1]], open='wt')
sink(log_file)

# move to the workinging directory
setwd(snakemake@params[["wd"]])

# Forward and reverse fastq filenames have format: {ds}_{resource]_{sample}_{factor}_R1-trimmed.fastq
fnFs <- sort(snakemake@input[["fwd"]])
fnRs <- sort(snakemake@input[["rev"]])

print('File pairs to be processed in DADA2:')
print(data.frame(fwd_files=fnFs,rev_files=fnRs))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names = lapply(str_split(basename(fnFs), '_'), function(x){str_c(x[snakemake@params[['samp']]],collapse = "_")})
# Plot the quality of the reads
plotQualityProfile(fnFs)+
ggsave(snakemake@output[['fwd_err']])
plotQualityProfile(fnRs)+
ggsave(snakemake@output[['rev_err']])

print(paste('Quality plots saved to:', dirname(snakemake@output[['fwd_err']])))

# Place filtered files in filtered/ subdirectory
filtFs = file.path( 'data/filtered', paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs = file.path( 'data/filtered', paste0(sample.names, "_R2_filt.fastq.gz"))
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs,
					 maxN=0, maxEE=c(5,5), truncQ=5, rm.phix=TRUE,
					 compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

print('Reads that survived filtering:')
print(out)

# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)+
	ggsave(snakemake@output[['base_err']])

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# sample Infrence
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# merge pair reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# construct ASV table
seqtab <- makeSequenceTable(mergers)
print('sequence lengths:')
table(nchar(getSequences(seqtab)))

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
write.table(seqtab.nochim, file=snakemake@output[['asv_table']], sep="\t")

# Track Reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

print('Reads Surviving each step')
print(track)






