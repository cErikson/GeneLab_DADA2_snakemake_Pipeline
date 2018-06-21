##	++==============================================++
##	|| DADA2 to BIOM for marker gene GLDS datasets  ||
##	++==============================================++
##
##	Christian Erikson: Bio.Erkson@gmail.com christian.b.erikson@nasa.gov
##

library(biom)
library(readr)
library(Risa)
library(stringr)


# Log errors/ messages to the snakemake log
log_file=file(snakemake@log[[1]], open='wt')
sink(log_file)

# move to the workinging directory
setwd(snakemake@params[["wd"]])

# read in all the files
isa=readISAtab(path = snakemake@params[['isa_dir']])
counts=read.table(snakemake@input[['counts']])
taxa=read_tsv(snakemake@input[['taxa']])

# HaCk away the duplicate columns, The isa standard needs to have unique columns. Debuging this error was cryptic and took forever
dedup_isa=isa['study.files'][[snakemake@params[['isa_samp_file']]]][, !duplicated(colnames(isa['study.files'][[snakemake@params[['isa_samp_file']]]]))]  ##!!!!! This is a hack, The isa file has duplicate col names for unit, ref protocol, and terms. The isa file standard needs to have unique column names 
dedup_isa=dedup_isa[order(dedup_isa[snakemake@params[['isa_samp_feild']]]),]

# If the orders do not match, stop. 
stopifnot(colnames(counts)==taxa$ASV) # Count sequences not in same order as Taxa seq 
stopifnot(str_split_fixed(rownames(counts), '_', 2)[,1]==dedup_isa[snakemake@params[['isa_samp_feild']]]) # count sample names not in same order as isa study sample names

# Create the biom object from data
the_biom=make_biom(t(counts),observation_metadata = taxa, sample_metadata = dedup_isa )

# And write it out
write_biom(the_biom, snakemake@output[['biom']])
