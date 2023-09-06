#!/usr/bin/env Rscript
################################################################################
# Title: batch_merge.R
# Discription: Merge DADA batches 
# Author: Christian Erikson
# GitHub: https://github.com/cErikson/DADA2_cli_and_galaxy_wrapper
# Date: sep 2023
################################################################################
packages = c("argparse", 'phyloseq','biomformat')
if (any(!(packages %in% installed.packages()[,"Package"]))) stop(sprintf('These packages are required: %s', packages[!(packages %in% installed.packages()[,"Package"])]))
library(argparse)
library(phyloseq)
library(biomformat)
parser <- ArgumentParser()
parser$add_argument('-i', "--input_phylo", type="character", help="Phyloseq object from taxonomy step")
parser$add_argument('-t', "--tree", type="character", help="Newick tree")
parser$add_argument('-p', "--output_phylo", default= F ,type="character",help="Name of (or path to) the output phylo.")
#parser$add_argument('-b', "--output_biom", default= F ,type="character",help="Name of (or path to) the output biom")
parser$add_argument("-s", "--meta",  help="Sample metadata, in which rows are samples", type="character", default=F)
parser$add_argument("-d", "--delim", help="The delimiter used in sample metadata.", type="character", default='\t')
#parser$add_argument("-c", "--samp_col", help="The column name, which has sample ids. Ids need to have a coraspnding entry in the ASV table, which may be a substring int the ASV.", type=str, default=False)
#parser$add_argument("-r", "--regex", help="Regex used to extract sample id names from the ASV ids", type=str, default='(.*)')
parser$add_argument("-n", "--study_id", help="Study ID to be used in the Biom ", type="character", default='DADA2 study')
args = parser$parse_args()
print(args)

data=readRDS(args$input_phylo)
tree=read_tree(args$tree)
if(args$meta != F){meta=read.delim(args$meta,sep = args$delim)}

data=merge_phyloseq(data, 
                    phy_tree(tree), 
                    if(args$meta != F){sample_data(meta)})


if(args$output_phylo != F){
  saveRDS(data, args$output_phylo)
}

# if(args$output_biom != F){
#   B=make_biom(
#     data = otu_table(data),
#     if(args$meta != F){sample_metadata = sample_data(data)},
#     observation_metadata = tax_table(data),
#     id = args$study_id
#   )
#   write_biom(B, args$output_biom)
# }