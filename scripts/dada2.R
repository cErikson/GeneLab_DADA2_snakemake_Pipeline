##	++======================================================++
##	|| DADA2 Analysis script for marker gene GLDS datasets  ||
##	++======================================================++
##
##	Christian Erikson: Bio.Erkson@gmail.com christian.b.erikson@nasa.gov
##

library(tidyverse)
library(dada2)


# Log errors/ messages to the snakemake log
log_file=file(snakemake@log[[1]], open='wt')
sink(log_file)

# move to the workinging directory
setwd(snakemake@params[["wd"]])

# Forward and reverse fastq filenames have format: {ds}_{resource]_{sample}_{factor}_R1-trimmed.fastq
fnFs = sort(snakemake@input[["fwd"]])
fnRs = sort(snakemake@input[["rev"]])
stopifnot(all(file.exists(fnFs)))
stopifnot(all(file.exists(fnRs)))

print('File pairs to be processed in DADA2:')
print(data.frame(fwd_files=fnFs, rev_files=fnRs))

# Extract sample names, assuming filenames have format: {ds}_{resource]_{sample}_{factor}_R1-trimmed.fastq
sample.names = lapply(str_split(basename(fnFs), '_'), function(x){str_c(x[snakemake@params[['samp']]],collapse = "_")})
study.name = str_split(basename(fnFs), '_')[[1]][1]

# Place filtered files in filtered/ subdirectory
print('Filtering reads')
filtFs = file.path( 'data/filtered', paste0(study.name,'_16s-amplicon-sequencing_',sample.names, "_R1-filtered.fastq.gz"))
filtRs = file.path( 'data/filtered', paste0(study.name,'_16s-amplicon-sequencing_',sample.names, "_R2-filtered.fastq.gz"))
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs,
					 maxN=0, maxEE=c(snakemake@params[['maxEE_fwd']],snakemake@params[['maxEE_rev']]), truncQ=snakemake@params[['truncQ']],
					 rm.phix=snakemake@params[['rmphix']], maxLen = snakemake@params[['maxLen']], minLen=snakemake@params[['minLen']],
					 compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]
sample.names.filt = lapply(str_split(basename(filtFs), '_'), function(x){str_c(x[snakemake@params[['samp']]],collapse = "_")})

if (any(exists==FALSE)){
	print('The fliter has removed all the reads from the following files, thus they will no longer be analyized')
	print(filtFs[exists==FALSE])
	print(filtRs[exists==FALSE])
}

# Learn error rates
print('Learning error rates')
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)+
	ggsave(snakemake@output[['base_err']])

# Dereplication
print('Dereplicating')
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names.filt
names(derepRs) <- sample.names.filt

# sample Infrence
print('DADA2 Core')
setDadaOpt(HOMOPOLYMER_GAP_PENALTY=snakemake@params[['HOMOPOLYMER_GAP_PENALTY']], BAND_SIZE=snakemake@params[['BAND_SIZE']])
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# merge pair reads
print('Merge pair end reads')
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# construct ASV table
print('Creating ASV table')
seqtab <- makeSequenceTable(mergers)
print('Sequence lengths of ASVs:')
table(nchar(getSequences(seqtab)))

# Remove Chimeras
print('Remove Chimeras')
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
write.table(seqtab.nochim, file=snakemake@output[['asv_table']], sep="\t")

# Track Reads
getN = function(x) sum(getUniques(x))
prefilt_reads = data.frame(samp=unlist(sample.names), reads_in = out[,1], reads_filt = out[,2]) 
filt_reads = data.frame(samp=unlist(sample.names.filt), denoisedF = sapply(dadaFs, getN), denoisedR = sapply(dadaRs, getN), merged = sapply(mergers, getN), nonchim = rowSums(seqtab.nochim))
track = full_join(prefilt_reads, filt_reads)
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
print('Reads Surviving each step')
print(track)

