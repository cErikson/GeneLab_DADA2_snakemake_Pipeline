import yaml
import os
import pdb
from snakemake.utils import R
from snakemake.utils import report

## DEFINE WORKDIR ##
WORKDIR=workflow.basedir+'/'

configfile: "Config"

# Gather files to be processed.
file_ids = glob_wildcards(WORKDIR+"data/sequencing/{id}_{read}.fastq.gz")
DS_NUM=file_ids[0][0].split('_')[0]

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

    
##### RULES #####

rule setup:
    shell:
        '''
        cd {WORKDIR}
        mkdir -p logs data/trim reports/plots
        
        '''
        
rule help:
    run:
        print('''
        ++==============================================++
        || DADA2 pipeline for marker gene GLDS datasets ||
        ++==============================================++
        ''')

rule revcomp_primer:
    input:
        primers='data/metadata/primers/{dir}_primers.txt'
    output:
        primers=temporary('data/metadata/primers/{dir}_revcomp_primers.txt')
    run:
        with open(output[0], 'w') as fho:
                fho.writelines(['{}\n{}\n'.format(x[0],x[1]) for x in yield_fasta(input[0], rev_comp=True)])
        
rule cutadapt:
    input:
        fwd=WORKDIR+'data/sequencing/{ID}_R1.fastq.gz',
        rev=WORKDIR+'data/sequencing/{ID}_R2.fastq.gz',
        fwd_prime='data/metadata/primers/fwd_primers.txt',
        rev_prime='data/metadata/primers/rev_primers.txt',
        fwd_revcomp_prime='data/metadata/primers/fwd_revcomp_primers.txt',
        rev_revcomp_prime='data/metadata/primers/rev_revcomp_primers.txt'
    output:
        fwd_trim=WORKDIR+'data/trim/{ID}_R1-trimmed.fastq',
        rev_trim=WORKDIR+'data/trim/{ID}_R2-trimmed.fastq'
    log:
        lambda wc: WORKDIR+'logs/cutadapt_{ID}'.format(ID=wc.ID)
        
    shell:
        '''
         cutadapt -g file:{input.fwd_prime} -G file:{input.rev_revcomp_prime} -a file:{input.rev_prime} -A file:{input.fwd_revcomp_prime} --max-n=0 -m {config[min_len_dis]} -o {output.fwd_trim} -p {output.rev_trim} {input.fwd} {input.rev} > {log} 
        '''

rule dada2:
    input:
        fwd=expand(WORKDIR+'data/trim/{IDS}_R1-trimmed.fastq', IDS=set(file_ids.id)),
        rev=expand(WORKDIR+'data/trim/{IDS}_R2-trimmed.fastq', IDS=set(file_ids.id))
    output:
        fwd_err='report/plots/fwd_quals.png',
        rev_err='report/plots/rev_quals.png',
        flt_fwd=expand("data/filtered/{IDS}_R1-trimmed.fastq.gz", IDS=set(file_ids.id)),
        flt_rev=expand("data/filtered/{IDS}_R2-trimmed.fastq.gz", IDS=set(file_ids.id)),
        base_err='report/plots/base_error_rates.png',
        asv_table='data/{GLDS}_asvtable.csv'
    params:
        wd=WORKDIR,
        samp=config['samp_fields']
    
    log:
        'logs/dada2.txt'
    script:
        'scripts/dada2.R'
        