#!/usr/bin/env Rscript
################################################################################
# Title: batch_merge.R
# Discription: Merge DADA batches 
# Author: Christian Erikson
# GitHub: https://github.com/cErikson/DADA2_cli_and_galaxy_wrapper
# Date: sep 2023
################################################################################
packages = c("argparse", 'phyloseq')
if (any(!(packages %in% installed.packages()[,"Package"]))) stop(sprintf('These packages are required: %s', packages[!(packages %in% installed.packages()[,"Package"])]))
library(argparse)
library(phyloseq)

parser <- ArgumentParser()
parser$add_argument('-a', "--asv_tables", nargs='+', type="character", help="Name of (or path to) the ASV table output from DADA2. In R matrix format. aka first feild of col names is missing")
parser$add_argument('-o', "--output_table", default= F ,type="character",help="Name of (or path to) the output taxa table.")
args = parser$parse_args()

read_asv=function(x) phyloseq(otu_table(read.table(x, header=T),taxa_are_rows = T))
phylos=lapply(args.asv_tables, read_asv)
merged=do.call(merge_phyloseq, phylos)
write.table(merged, args.output_table)