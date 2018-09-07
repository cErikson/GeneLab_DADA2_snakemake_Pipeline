```
                                                         _____
  ______________________________________________________/  O  \___/
 / DADA2 Snakemake pipeline for marker gene GLDS datasets ____/   \
<_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/__/ 
```

#Libraries and Programs 
```
Unix: fastqc
Python3: snakemake, cutadapt, multiqc, biom
R: tidyverse, dada2
Scripts in cli folder: https://github.com/cErikson/DADA2_cli_and_galaxy_wrapper
```

# Asumptions
0. Sequencing files are found under `data/sequencing`
1. Reads have been demuliplexed. Any barcode used for pcr/optical duplicates has been used.
2. Gzipped fastq file follow the `<Accession #>_<resource category>_<Sample Name>_(<Factor Level>_)+.<read>.fastq.gz` format. Where (...)+ indicates one or more factor levels, delimited by under_scores. Read feild is `R1` or `R2`

# Instructions 
0. Run `snakemake setup`
1. Add primers to /data/metadata/primers/[fwd_primer.fasta, rev_primers.fasta]
2. Unzip the metadata files into `data/metadata/`, there should not be any isa files in folders in this directory. 
3. Run `snakemake raw_multiqc -j <cores>` to get a sense of the read quality, look under `report/`
4. Adjust the parameters in this file to meet the needs of the experiment. Please read through each parameter discription.
Failure to do so may give you results that look good, but are bad. 
    Sections marked Optional probably don't need to be touched for a standard illumina run
 5. Run `snakemake all -j <num_cpu_cores>`, and go grab lunch, the pipeline should finish in two hours for a ~20gb dataset 
 on a 8 core machine using default parameters.

# Output
```
├── Config		# Configuration file for pipeline
├── data
│   ├── GLDS-146_metagenomics_dada2-asv-matrix.tsv		# A flat tsv containing counts of amplicon sequence variants(ASV) by sample
│   ├── GLDS-146_metagenomics_dada2-biom.biom			# A BIOM v1.0 (json) file containing the ASV counts, taxanomic assignment, asnd sample metadata
│   ├── GLDS-146_metagenomics_dada2-taxonomy-matrix.tsv	# A flat tsv whos order corasponds to the asv. Which contains taxanomic assignments to ASV.
│   ├── metadata
│   │   ├── isa_files
│   │   │   ├── a_{name}.txt		# Assay metadata
│   │   │   ├── i_{name}.txt		# Investigation metadata
│   │   │   └── s_{name}.txt		# Sample metadata
│   │   └── primers
│   │       ├── fwd_primers.fasta	# Forward primers that were used for trimming
│   │       └── rev_primers.fasta	# Reverse primers that were used for trimming
│   ├── sequencing
│   │   └── {study name}_16s-amplicon-sequencing_{sample id}_{factors}_{read #}-filtered.fastq.gz	# Raw sequencing reads 
│   ├── filtered
│   │   └── {study name}_16s-amplicon-sequencing_{sample id}_{factors}_{read #}-filtered.fastq.gz	# Reads that passed the DADA2 filter.
│   └── trimmed
│       └── {study name}_16s-amplicon-sequencing_{sample id}_{factors}_{read #}-filtered.fastq.gz	# Reads that passed cutadapt trimming.
├── logs
│   └── {pipeline step}.txt 	# Log files from stdout and stderr from the pipline.
├── notes
│   └── filename_changes.tsv 	# Document of file name changes
├── report
│   ├── cutadapt
│   │   ├── cutadapt_GLDS-146_16s-amplicon-sequencing_SAMN06277231_Mouse-1_day-10_gray-0.1
│   ├── fastqc
│   │   ├── filtered
│   │   │   └── {read id}_fastqc.html	# FastQC report
│   │   ├── raw
│   │   │   └── {read id}_fastqc.html	# FastQC report
│   │   └── trimmed
│   │   │   └── {read id}_fastqc.html	# FastQC report
│   ├── multiqc_data
│   │   └── *
│   ├── multiqc_report.html		# Aggregated report for FastQC and cutadapt.
│   ├── plots
│       └── base_error_rates.png 
├── scripts
│   ├── assign_tax.R	# script for taxanomic assignment
│   ├── dada2biom.py	# script for conversion to BIOM
│   └── dada2.R			# script for main DADA2 
└── Snakefile			# Snakemake pipeline file
```
