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
DS_NUM=file_ids[0][0].split('_')[0]

    
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
        ||             Paired End Edition               ||
        ++==============================================++
        
        
        Please see the Config file for set up of the pipeline.
        ''')

rule fastqc:
    shell:
        '''
        zcat data/sequencing/*R1.fastq.gz | fastqc stdin -o ./report/fastqc
        mv ./report/fastqc/stdin_fastqc.html ./report/fastqc/R1_fastqc.html
        mv ./report/fastqc/stdin_fastqc.zip ./report/fastqc/R1_fastqc.zip
        zcat data/sequencing/*R2.fastq.gz | fastqc stdin -o ./report/fastqc
        mv ./report/fastqc/stdin_fastqc.html ./report/fastqc/R2_fastqc.html
        mv ./report/fastqc/stdin_fastqc.zip ./report/fastqc/R2_fastqc.zip
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
        fwd_trim=WORKDIR+'data/trim/{ID}_R1-trimmed.fastq',
        rev_trim=WORKDIR+'data/trim/{ID}_R2-trimmed.fastq'
    log:
        lambda wc: WORKDIR+'logs/cutadapt_{ID}'.format(ID=wc.ID)
    params:
        discard = '--discard-untrimmed ' if config['discard_untrimmed'] is True else '',
        indels = '--no-indels ' if config['no_indels'] is True else ''
    shell:
        '''
         cutadapt -g file:{input.fwd_prime} -G file:{input.rev_revcomp_prime} -a file:{input.rev_prime} -A file:{input.fwd_revcomp_prime} -o {output.fwd_trim} -p {output.rev_trim} --max-n={config[max_n]} -m {config[min_len_dis]} --error-rate {config[error_rate]} {params.discard}{params.indels}--pair-filter {config[pair_filter]} {input.fwd} {input.rev} > {log} 
        ''' 

rule dada2:
    input:
        fwd=expand(WORKDIR+'data/trim/{IDS}_R1-trimmed.fastq', IDS=set(file_ids.id)),
        rev=expand(WORKDIR+'data/trim/{IDS}_R2-trimmed.fastq', IDS=set(file_ids.id))
    output:
        fwd_err='report/plots/fwd_quals.png',
        rev_err='report/plots/rev_quals.png',
        base_err='report/plots/base_error_rates.png',
        asv_table='data/{GLDS}_asvtable.tsv'
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
        
    
    log:
        'logs/dada2.txt'
    script:
        'scripts/dada2.R'
        
rule assign_tax:
    input:
        asv_table='data/{GLDS}_asvtable.tsv'
    output:
        tax_table='data/{GLDS}_taxtable.tsv'
    params:
        wd=WORKDIR,
        tax_train=config['tax_train'],
        tax_species=config['tax_species'],
        add_taxa_ds=config['add_taxa_ds']
    log:
        'logs/assign_tax.txt'
    script:
        'scripts/assign_tax.R'

rule DADA2BIOM:
    input:
        counts='data/{GLDS}_asvtable.tsv',
        taxa='data/{GLDS}_taxtable.tsv'
    output:
        biom='data/{GLDS}.biom'
    params:
        wd=WORKDIR,
        isa_dir=config['isa_dir'],
        isa_samp_file=config['isa_samp_file'],
        isa_samp_feild=config['isa_samp_feild']
    log:
        'logs/dada2biom.txt'
    script:
        'scripts/dada2biom.R'
    

rule all:
    input: 
         'data/{GLDS}.biom'.format(GLDS=DS_NUM)