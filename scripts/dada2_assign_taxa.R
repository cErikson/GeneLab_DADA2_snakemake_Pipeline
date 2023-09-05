#!/usr/bin/env Rscript
################################################################################
# Title: dada2_assign_taxa.R
# Discription: A wrapper for DADA2 taxonomic assignment
# Author: Christian Erikson
# GitHub: https://github.com/cErikson/DADA2_cli_and_galaxy_wrapper
# Date: 6/26/18
################################################################################
packages = c("dada2", "argparse", 'DECIPHER', 'phyloseq')
if (any(!(packages %in% installed.packages()[,"Package"]))) stop(sprintf('These packages are required: %s', packages[!(packages %in% installed.packages()[,"Package"])]))
library(argparse)

parser <- ArgumentParser()
  parser$add_argument('-p', "--phylo", default = F, type="character", help="Name of (or path to) the phyloseq output from DADA2")
  parser$add_argument('-a', "--asv_table", default = F, type="character", help="Name of (or path to) the ASV table output from DADA2. In R matrix format. aka first feild of col names is missing")
	parser$add_argument('-o', "--output_table", default= F ,type="character",help="Name of (or path to) the output taxa table.")
	parser$add_argument('-O', "--output_phylo", default= F ,type="character",help="Name of (or path to) the output phyloseq.")
	parser$add_argument('-m', "--method", default="idtaxa",type="character", help="There are two methods, idtaxa and bayes. The IDTAXA algo from DECIPHER, or naive bayes from dada2")
	parser$add_argument('--multithread', default= T, type="logical", help='(Optional). Default is TRUE. If TRUE, multithreading is enabled and the number of available threads is automatically determined. If an integer is provided, the number of threads to use is set by passing the argument on to setThreadOptions.')
	parser$add_argument('--verbose', default= T, type="logical", help='(Optional). Default FALSE. If TRUE, print status to standard output.')
	# Bayes
	parser$add_argument("-t", "--taxa_train", nargs='+', type="character", help="The path to one or more taxa trainning file(s)")
	parser$add_argument("-s","--species_train",default=FALSE,nargs='+', type="character", help="Bayes: The path to  species training file")
	parser$add_argument('--minBoot', default = 50, type="integer", help='(Optional). Default 50. The minimum bootstrap confidence for assigning a taxonomic level.')
	parser$add_argument('--tryRC', default= F, type="logical", help='(Optional). Default FALSE. If TRUE, the reverse-complement of each sequences will be used for classification if it is a better match to the reference sequences than the forward sequence.')
  parser$add_argument('--outputBootstraps', default= F, type="logical", help='(Optional). Default FALSE. If TRUE, bootstrap values will be retained in an integer matrix. A named list containing the assigned taxonomies (named "taxa") and the bootstrap values (named "boot") will be returned. Minimum bootstrap confidence filtering still takes place, to see full taxonomy set minBoot=0')
	parser$add_argument('--taxLevels', nargs='+', default=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") , help='(Optional). Default is <Kingdom Phylum Class Order Family Genus Species>. The taxonomic levels being assigned. Truncates if deeper levels not present in training fasta.')
	parser$add_argument('--allowMultiple', default= F, type="logical", help='(Optional). Default FALSE. Defines the behavior when multiple exact matches against different species are returned. By default only unambiguous identifications are return. If TRUE, a concatenated string of all exactly matched species is returned. If an integer is provided, multiple identifications up to that many are returned as a concatenated string.')
	#


args = parser$parse_args()
if (args$asv_table == F && args$phylo == F) stop("Need either a ASV_table or phyloseq object for input")
(args$output_table == F && args$output_phylo == F) stop("Need either a ASV_table or phyloseq object for output")
if (args$asv_table != F && args$phylo != F) stop("Need only one of either of the ASV_table or phyloseq object for input")
if (args$asv_table != F && any(!file.exists(args$asv_table))) stop(" The ASV_table does not exist")
if (args$phylo != F && any(!file.exists(args$phylo))) stop(" The phyloseq file does not exist")
if (!(args$method == 'idtaxa' || args$method == 'bayes')) stop("Invalid method")
if (args$method == 'bayes' && args$species_train != F && any(!file.exists(args$species_train))) stop("The species_trainning does not exist")
if (any(!file.exists(args$taxa_train))) stop("The taxa_trainning does not exist")


library(dada2)
print(args)
# Read in seqtab file
if (args$asv_table != F){
  phylo=phyloseq(otu_table(read.table(args$asv_table), taxa_are_rows = F))
  dna = Biostrings::DNAStringSet(taxa_names(phylo))
  names(dna) = taxa_names(phylo)
  phylo = merge_phyloseq(phylo, dna)
  taxa_names(phylo) = paste0("ASV", seq(ntaxa(phylo)))
}
if (args$phylo != F){
  phylo=readRDS(args$phylo)
}

if (args$method == 'bayes'){
  # assignTaxonomy
  taxaid = assignTaxonomy(refseq(phylo), args$taxa_train, minBoot = args$minBoot, tryRC = args$tryRC, outputBootstraps = args$outputBootstraps,
  					  taxLevels = args$taxLevels, multithread=args$multithread, verbose = args$verbose) # main training set
  if (args$species_train!=F) taxaid = addSpecies(taxaid, args$species_train, allowMultiple = args$allowMultiple, verbose = args$verbose) # Add species 
}

if (args$method == 'idtaxa'){
library("DECIPHER"); packageVersion("DECIPHER")

dna=refseq(phylo)
load(args$taxa_train) # CHANGE TO THE PATH OF YOUR TRAINING SET
ids = IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=T) # use all processors
ranks = args$taxLevels # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks
}


# Save
if (args.output_table != F){
  write.table(taxaid, args$output_table, sep="\t")
}
if (args.output_table != F){
  phylo=merge_phyloseq(phylo, tax_table(taxid))
  saveRDS(phylo, args.output_phylo)
}
