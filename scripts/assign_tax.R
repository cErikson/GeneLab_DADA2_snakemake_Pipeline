##	++=======================================================++
##	|| DADA2 taxnonomy script for marker gene GLDS datasets  ||
##	++=======================================================++
##
##	Christian Erikson: Bio.Erkson@gmail.com christian.b.erikson@nasa.gov
##

library(dada2)
library(dplyr)
library(readr)

# Log errors/ messages to the snakemake log
log_file=file(snakemake@log[[1]], open='wt')
sink(log_file)

# move to the workinging directory
setwd(snakemake@params[["wd"]])

# Read in seqtab file

seqtab=as.matrix(read.table(snakemake@input[['asv_table']]))
# assignTaxonomy
taxa = assignTaxonomy(seqtab, snakemake@params[['tax_train']], multithread=TRUE, verbose = TRUE) # main training set
taxa = addSpecies(taxa, snakemake@params[['tax_species']]) # Add species 
taxa = add_rownames(as.data.frame(taxa), var = 'ASV')
if (!is.null(snakemake@params[['add_taxa_ds']])){
	for (i in snakemake@params[['add_taxa_ds']]){
		add_taxa = assignTaxonomy(seqtab, i, multithread=TRUE, verbose = TRUE)
		add_taxa = add_rownames(as.data.frame(add_taxa), var = 'ASV')
		taxa = bind_rows(taxa, add_taxa) %>%
			distinct(., ASV)
	}
}

# Save
write_tsv(taxa, snakemake@output[['tax_table']])
