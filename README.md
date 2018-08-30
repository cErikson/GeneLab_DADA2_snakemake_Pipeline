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
