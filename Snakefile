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

rule all:
    input: 
         biom='data/{GLDS}_metagenomics_dada2-json.biom'.format(GLDS=DS_NUM),
         report='report/multiqc_report.html'

rule setup:
    shell:
        '''
        mkdir -p logs data/trim report/plots report/fastqc data/metadata/
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
 
 
#### CUTADAPT #####
ruleorder: cutadapt_pe > cutadapt_se
        
rule cutadapt_pe:
    input:
        fwd='data/sequencing/{ID}_R1.fastq.gz',
        rev='data/sequencing/{ID}_R2.fastq.gz',
        fwd_prime='data/metadata/primers/fwd_primers.fasta',
        rev_prime='data/metadata/primers/rev_primers.fasta',
        fwd_revcomp_prime='data/metadata/primers/fwd_revcomp_primers.fasta',
        rev_revcomp_prime='data/metadata/primers/rev_revcomp_primers.fasta'
    output:
        fwd_trim='data/trimmed/{ID}_R1-trimmed.fastq',
        rev_trim='data/trimmed/{ID}_R2-trimmed.fastq'
    log:
        'report/cutadapt/cutadapt_{ID}'
    params:
        discard = '--discard-untrimmed ' if config['discard_untrimmed'] is True else '',
        indels = '--no-indels ' if config['no_indels'] is True else '',
        trim_n='--trim-n ' if config['trim_end_n'] is True else ''
    shell:
        '''
         cutadapt -g file:{input.fwd_prime} -G file:{input.rev_prime} -a file:{input.rev_revcomp_prime} -A file:{input.fwd_revcomp_prime} -o {output.fwd_trim} -p {output.rev_trim} --max-n={config[max_n]} -m {config[min_len_dis]} --error-rate {config[error_rate]} {params.discard}{params.indels}{params.trim_n}--pair-filter {config[pair_filter]} {input.fwd} {input.rev} > {log} 
        ''' 

rule cutadapt_se:
    input:
        fwd='data/sequencing/{ID}_R1.fastq.gz',
        fwd_prime='data/metadata/primers/fwd_primers.fasta',
        rev_prime='data/metadata/primers/rev_primers.fasta',
        fwd_revcomp_prime='data/metadata/primers/fwd_revcomp_primers.fasta',
        rev_revcomp_prime='data/metadata/primers/rev_revcomp_primers.fasta'
    output:
        fwd_trim='data/trimmed/{ID}_R1-trimmed.fastq'
    log:
        'report/cutadapt/cutadapt_{ID}'
    params:
        discard = '--discard-untrimmed ' if config['discard_untrimmed'] is True else '',
        indels = '--no-indels ' if config['no_indels'] is True else '',
        trim_n='--trim-n ' if config['trim_end_n'] is True else ''
    shell:
        '''
         cutadapt -g file:{input.fwd_prime} -G file:{input.rev_prime} -a file:{input.rev_revcomp_prime} -A file:{input.fwd_revcomp_prime} -o {output.fwd_trim} --max-n={config[max_n]} -m {config[min_len_dis]} --error-rate {config[error_rate]} {params.discard}{params.indels}{params.trim_n}--pair-filter {config[pair_filter]} {input.fwd} > {log} 
        ''' 

##### DADA FILTER #####
ruleorder: dada2_filter_pe > dada2_filter_se

rule dada2_filter_pe:
    input:
        fwd=expand('data/trimmed/{IDS}_R1-trimmed.fastq', IDS=set(file_ids.id)),
        rev=expand('data/trimmed/{IDS}_R2-trimmed.fastq', IDS=set(file_ids.id))
    output:
        fwd=expand('data/filtered/{IDS}_R1-filtered.fastq.gz', IDS=set(file_ids.id)),
        rev=expand('data/filtered/{IDS}_R2-filtered.fastq.gz', IDS=set(file_ids.id))
    threads: config['threads']
    shell:
        '''
        scripts/dada2_filter.R -f {input.fwd} -r {input.rev} -p data/filtered/ {config[samp_fields]} -q {config[truncq]} -l {config[trunclen]} -L {config[trimleft]} -x {config[maxlen]} -n {config[minlen]} -N {config[maxn]} -d {config[minq]} -E {config[maxee]} --rm_phix {config[rm_phix]} # --n_reads {config[n_reads]}
        '''
        
rule dada2_filter_se:
    input:
        fwd=expand('data/trimmed/{IDS}_R1-trimmed.fastq', IDS=set(file_ids.id))
    output:
        fwd=expand('data/filtered/{IDS}_R1-filtered.fastq.gz', IDS=set(file_ids.id))
    threads: config['threads']
    shell:
        '''
        scripts/dada2_filter.R -f {input.fwd} -p data/filtered/ {config[samp_fields]} -q {config[truncq]} -l {config[trunclen]} -L {config[trimleft]} -x {config[maxlen]} -n {config[minlen]} -N {config[maxn]} -d {config[minq]} -E {config[maxee]} --rm_phix {config[rm_phix]} # --n_reads {config[n_reads]}
        '''
 
##### DADA CORE #####
ruleorder: dada2_pe > dada2_se
rule dada2_pe:
    input:
        fwd=expand('data/filtered/{IDS}_R1-filtered.fastq.gz', IDS=set(file_ids.id)),
        rev=expand('data/filtered/{IDS}_R2-filtered.fastq.gz', IDS=set(file_ids.id))
    output:
        asv='data/{GLDS}_metagenomics_dada2-asv-matrix.tsv'.format(GLDS=DS_NUM),
    params:
        matrix = '' if config['score_matrix'] == 'nuc44' else '--SCORE_MATRIX '+config['score_matrix']
    threads: config['threads']
    shell:
        '''
        scripts/dada2.R -f {input.fwd} -r {input.rev} -o {output.asv} {config[samp_fields]} \
        -c {config[chimera_rm]} --LEARN_READS_N {config[learn_reads_n]} --OMEGA_A {config[omega_a]} --USE_QUALS {config[use_quals]}\
        --USE_KMERS {config[use_kmers]} --KDIST_CUTOFF {config[kdist_cutoff]} --BAND_SIZE {config[band_size]} {params.matrix}\
        --GAP_PENALTY {config[gap_penalty]} --HOMOPOLYMER_GAP_PENALTY {config[homopolymer_gap_penalty]} --MIN_FOLD {config[min_fold]}\
        --MIN_HAMMING {config[min_hamming]} --MAX_CLUST {config[max_clust]} --MAX_CONSIST {config[max_consist]}
        '''

rule dada2_se:
    input:
        fwd=expand('data/filtered/{IDS}_R1-filtered.fastq.gz', IDS=set(file_ids.id))
    output:
        asv='data/{GLDS}_metagenomics_dada2-asv-matrix.tsv'.format(GLDS=DS_NUM),
    params:
        matrix = '' if config['score_matrix'] == 'nuc44' else '--SCORE_MATRIX '+config['score_matrix']
    threads: config['threads']
    shell:
        '''
        scripts/dada2.R -f {input.fwd} -o {output.asv} {config[samp_fields]} \
        -c {config[chimera_rm]} --LEARN_READS_N {config[learn_reads_n]} --OMEGA_A {config[omega_a]} --USE_QUALS {config[use_quals]}\
        --USE_KMERS {config[use_kmers]} --KDIST_CUTOFF {config[kdist_cutoff]} --BAND_SIZE {config[band_size]} {params.matrix}\
        --GAP_PENALTY {config[gap_penalty]} --HOMOPOLYMER_GAP_PENALTY {config[homopolymer_gap_penalty]} --MIN_FOLD {config[min_fold]}\
        --MIN_HAMMING {config[min_hamming]} --MAX_CLUST {config[max_clust]} --MAX_CONSIST {config[max_consist]}
        '''

rule assign_tax:
    input:
        asv_table='data/{GLDS}_metagenomics_dada2-asv-matrix.tsv'
    output:
        tax_table='data/{GLDS}_metagenomics_dada2-taxonomy-matrix.tsv'
    threads: config['threads']
    shell:
        '''
        scripts/dada2_assign_taxa.R {input} {output}\
        -t {config[taxa_train]} -s {config[species_train]} --minBoot {config[minboot]} --tryRC {config[tryrc]}\
        --taxLevels {config[taxlevels]} --allowMultiple {config[allowmultiple]} --multithread {threads}
         '''



rule DADA2BIOM:
    input:
        asv='data/{GLDS}_metagenomics_dada2-asv-matrix.tsv',
        taxa='data/{GLDS}_metagenomics_dada2-taxonomy-matrix.tsv'
    output:
        biom='data/{GLDS}_metagenomics_dada2-json.biom'
    shell:
        '''
        scripts/dada2biom.py {input.asv} {output.biom} -t {input.taxa} -s data/metadata/{config[meta_samp_file]} -c {config[meta_samp_feild]} -r {config[meta_regex]} -n {DS_NUM}
        '''

