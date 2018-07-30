import os


## CONFIG ##
configfile: "Config"

## DEFINE WORKDIR ##
WORKDIR=workflow.basedir+'/'

##### DEFS #####

def yield_fasta(fasta, gz=False, rev_comp=False):
    '''Make a generator that yields (seq_header, seq) for each entry)'''
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X', 'M':'K', 'R':'Y', 'W':'W', 'S':'S', 'Y':'R', 'K':'M','V':'B', 'H':'D', 'D':'H', 'B':'V'}
    if gz==False:
        fhs=open(fasta, 'r')
    else:
        fhs=gzip.open(fasta, 'rt')
    l=fhs.readline().strip()
    while l != '':
        if l.startswith('>'):   #if header
            header=l  #save header
            if rev_comp is False:
                seq=fhs.readline().strip().upper()
            else:
                seq=''
                for base in fhs.readline().strip().upper()[::-1]:
                    seq+=complement[base]
            yield header, seq# yield data
        l=fhs.readline().strip()

##### SETUP #####
# Gather files to be processed.
file_ids = glob_wildcards(WORKDIR+"data/sequencing/{id}_{read}.fastq.gz")
# Grab the GLDS number.
DS_NUM=config['glds_num']

    
##### RULES #####

rule setup:
    shell:
        '''
        cd {WORKDIR}
        mkdir -p logs data/trim report/plots report/fastqc data/metadata/isa_files
        '''
        
rule help:
    run:
        print('''
        ++==============================================++
        || DADA2 pipeline for marker gene GLDS datasets ||
        ++==============================================++
        
        
        Please see the Config file for set up of the pipeline.
        ''')

rule fastqc_raw:
    input: 
        WORKDIR+'data/sequencing/{id}_{read}.fastq.gz'
    output:
        'report/fastqc/raw/{id}_{read}_fastqc.html'
    shell:
        '''
        fastqc {input} -o report/fastqc/raw/
        '''
        
rule fastqc_trimmed:
    input: 
        WORKDIR+'data/trimmed/{id}_{read}-trimmed.fastq'
    output:
        'report/fastqc/trimmed/{id}_{read}-trimmed_fastqc.html'
    shell:
        '''
        fastqc {input} -o report/fastqc/trimmed/
        '''
        
rule fastqc_filtered:
    input: 
        WORKDIR+'data/filtered/{id}_{read}-filtered.fastq.gz'
    output:
        'report/fastqc/filtered/{id}_{read}-filtered_fastqc.html'
    shell:
        '''
        fastqc {input} -o report/fastqc/filtered/
        '''
      
rule raw_multiqc:
    input:
        raw=expand('report/fastqc/raw/{IDS}_{READ}_fastqc.html', IDS=file_ids.id, READ=file_ids.read)
    shell:
        '''
        multiqc -d -f -ip -o report/ report/
        '''
        
rule full_multiqc:
    input:
        raw=expand('report/fastqc/raw/{IDS}_{READ}_fastqc.html', IDS=file_ids.id, READ=file_ids.read),
        trimmed=expand('report/fastqc/trimmed/{IDS}_{READ}-trimmed_fastqc.html', IDS=file_ids.id, READ=file_ids.read),
        filtered=expand('report/fastqc/filtered/{IDS}_{READ}-filtered_fastqc.html', IDS=file_ids.id, READ=file_ids.read)
    output:
        'report/multiqc_report.html'
    shell:
        '''
        multiqc -d -f -ip -o report/ report/
        '''

rule revcomp_primer:
    input:
        primers='data/metadata/primers/{dir}_primers.fasta'
    output:
        primers=temporary('data/metadata/primers/{dir}_revcomp_primers.fasta')
    run:
        with open(output[0], 'w') as fho:
                fho.writelines(['{}\n{}\n'.format(x[0],x[1]) for x in yield_fasta(input[0], rev_comp=True)])
        
rule cutadapt:
    input:
        fwd=WORKDIR+'data/sequencing/{ID}_R1.fastq.gz',
        rev=WORKDIR+'data/sequencing/{ID}_R2.fastq.gz',
        fwd_prime='data/metadata/primers/fwd_primers.fasta',
        rev_prime='data/metadata/primers/rev_primers.fasta',
        fwd_revcomp_prime='data/metadata/primers/fwd_revcomp_primers.fasta',
        rev_revcomp_prime='data/metadata/primers/rev_revcomp_primers.fasta'
    output:
        fwd_trim=WORKDIR+'data/trimmed/{ID}_R1-trimmed.fastq',
        rev_trim=WORKDIR+'data/trimmed/{ID}_R2-trimmed.fastq'
#    log:
#        lambda wc: WORKDIR+'report/cutadapt/cutadapt_{ID}'.format(ID=wc.ID)
    params:
        discard = '--discard-untrimmed ' if config['discard_untrimmed'] is True else '',
        indels = '--no-indels ' if config['no_indels'] is True else '',
        trim_n='--trim-n ' if config['trim_end_n'] is True else ''
    shell:
        '''
         cutadapt -g file:{input.fwd_prime} -G file:{input.rev_prime} -a file:{input.rev_revcomp_prime} -A file:{input.fwd_revcomp_prime} -o {output.fwd_trim} -p {output.rev_trim} --max-n={config[max_n]} -m {config[min_len_dis]} --error-rate {config[error_rate]} {params.discard}{params.indels}{params.trim_n}--pair-filter {config[pair_filter]} {input.fwd} {input.rev} > {log} 
        ''' 

rule dada2:
    input:
        fwd=expand(WORKDIR+'data/trimmed/{IDS}_R1-trimmed.fastq', IDS=set(file_ids.id)),
        rev=expand(WORKDIR+'data/trimmed/{IDS}_R2-trimmed.fastq', IDS=set(file_ids.id))
    output:
        base_err='report/plots/base_error_rates.png',
        asv_table='data/{GLDS}_metagenomics_dada2-asv-matrix.tsv'.format(GLDS=DS_NUM),
        fwd_filt=expand(WORKDIR+'data/filtered/{IDS}_R1-filtered.fastq.gz', IDS=set(file_ids.id)),
        rev_filt=expand(WORKDIR+'data/filtered/{IDS}_R2-filtered.fastq.gz', IDS=set(file_ids.id))
    params:
        wd=WORKDIR,
        samp=config['samp_fields'],
        maxEE_fwd=config['maxEE_fwd'],
        maxEE_rev=config['maxEE_rev'],
        maxLen=config['maxLen'],
        minLen=config['minLen'],
        rmphix=config['rmphix'],
        truncQ=config['truncQ'],
        BAND_SIZE=config['BAND_SIZE'],
        HOMOPOLYMER_GAP_PENALTY=config['HOMOPOLYMER_GAP_PENALTY']
    script:
        'scripts/dada2.R'
        
rule assign_tax:
    input:
        asv_table='data/{GLDS}_metagenomics_dada2-asv-matrix.tsv'
    output:
        tax_table='data/{GLDS}_metagenomics_dada2-taxonomy-matrix.tsv'
    params:
        wd=WORKDIR,
        tax_train=config['tax_train'],
        tax_species=config['tax_species'],
        add_taxa_ds=config['add_taxa_ds']
    script:
        'scripts/assign_tax.R'

rule DADA2BIOM:
    input:
        asv='data/{GLDS}_metagenomics_dada2-asv-matrix.tsv',
        taxa='data/{GLDS}_metagenomics_dada2-taxonomy-matrix.tsv'
    output:
        biom='data/{GLDS}_metagenomics_dada2-json.biom'
    params:
        isa_samp_file=config['isa_samp_file'],
        isa_samp_feild=config['isa_samp_feild']
    shell:
        '''
        scripts/dada2biom.py {input.asv} {output.biom} -t {input.taxa} -s data/metadata/isa_files/{config[isa_samp_file]} -c {config[isa_samp_feild]} -r {config[isa_regex]}
        '''

rule all:
    input: 
         biom='data/{GLDS}_metagenomics_dada2-json.biom'.format(GLDS=DS_NUM),
         report='report/multiqc_report.html'
