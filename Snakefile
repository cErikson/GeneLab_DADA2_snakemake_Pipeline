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
file_ids = glob_wildcards(WORKDIR+"data/seq/{id}_{read,R\d}{batch,.*}.fastq.gz")
# Grab the GLDS number.
DS_NUM=config['glds_num']

    
##### RULES #####

rule all:
    input: 
         biom=expand('data/{GLDS}{{batch}}_metagenomics_dada2-json.biom'.format(GLDS=DS_NUM), batch=file_ids.batch),
         report='report/multiqc_report.html'

rule setup:
    shell:
        '''
        mkdir -p logs data/trim report/plots report/fastqc data/metadata/ plots/
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
        'data/seq/{id}_{read}{batch}.fastq.gz'
    output:
        'report/fastqc/raw/{id}_{read,R.}{batch,_B.*}_fastqc.html'
    conda:
        'env/dada2.yaml'
    resources: cores = 1, runtime = config['fastqc_time'], part = config['fastqc_part']
    shell:
        '''
        fastqc {input} -o report/fastqc/raw/
        '''
        
rule fastqc_trimmed:
    input: 
        'data/trimmed/{id}_{read}{batch}-trimmed.fastq'
    output:
        'report/fastqc/trimmed/{id}_{read,R.}{batch,_B.*}-trimmed_fastqc.html'
    conda:
        "env/dada2.yaml"
    resources: cores = 1, runtime = config['fastqc_time'], part = config['fastqc_part'] 
    shell:
        '''
        fastqc {input} -o report/fastqc/trimmed/
        '''
        
rule fastqc_filtered:
    input: 
        'data/filtered/{id}_{read}{batch}-filtered.fastq.gz'
    output:
        'report/fastqc/filtered/{id}_{read,R.}{batch,_B.*}-filtered_fastqc.html'
    conda:
        "env/dada2.yaml"
    resources: cores = 1, runtime = config['fastqc_time'], part = config['fastqc_part'] 
    shell:
        '''
        fastqc {input} -o report/fastqc/filtered/
        '''
      
rule raw_multiqc:
    input:
        raw=['report/fastqc/raw/{IDS}_{READ}{BATCH}_fastqc.html'.format(IDS=I, READ=R, BATCH=B) for I,R,B in zip(file_ids.id, file_ids.read, file_ids.batch)]
    conda:
        "env/dada2.yaml"
    resources: cores = 1, runtime = config['fastqc_time'], part = config['fastqc_part'] 
    shell:
        '''
        multiqc -d -f -ip -o report/ report/
        '''
        
rule full_multiqc:
    input:
        raw=['report/fastqc/raw/{IDS}_{READ}{BATCH}_fastqc.html'.format(IDS=I, READ=R, BATCH=B) for I,R,B in zip(file_ids.id, file_ids.read, file_ids.batch)],
        trimmed=['report/fastqc/trimmed/{IDS}_{READ}{BATCH}-trimmed_fastqc.html'.format(IDS=I, READ=R, BATCH=B) for I,R,B in zip(file_ids.id, file_ids.read, file_ids.batch)],
        filtered=['report/fastqc/filtered/{IDS}_{READ}{BATCH}-filtered_fastqc.html'.format(IDS=I, READ=R, BATCH=B) for I,R,B in zip(file_ids.id, file_ids.read, file_ids.batch)]
    output:
        'report/multiqc_report.html'
    conda:
        "env/dada2.yaml"
    resources: cores = 1, runtime = config['fastqc_time'], part = config['fastqc_part'] 
    shell:
        '''
        multiqc -d -f -ip -o report/ report/
        '''

rule revcomp_primer:
    input:
        primers='data/metadata/primers/{dir}_primers.fasta'
    output:
        primers=temporary('data/metadata/primers/{dir}_revcomp_primers.fasta')
    resources: cores = 1, runtime = config['fastqc_time'], part = config['fastqc_part']
    run:
        with open(output[0], 'w') as fho:
                fho.writelines(['{}\n{}\n'.format(x[0],x[1]) for x in yield_fasta(input[0], rev_comp=True)])
 
 
#### CUTADAPT #####
ruleorder: cutadapt_pe > cutadapt_se
        
rule cutadapt_pe:
    input:
        fwd='data/seq/{ID}_R1{batch}.fastq.gz',
        rev='data/seq/{ID}_R2{batch}.fastq.gz',
        fwd_prime='data/metadata/primers/fwd_primers.fasta',
        rev_prime='data/metadata/primers/rev_primers.fasta',
        fwd_revcomp_prime='data/metadata/primers/fwd_revcomp_primers.fasta',
        rev_revcomp_prime='data/metadata/primers/rev_revcomp_primers.fasta'
    output:
        fwd_trim='data/trimmed/{ID}_R1{batch,_B.*}-trimmed.fastq',
        rev_trim='data/trimmed/{ID}_R2{batch,_B.*}-trimmed.fastq'
    log:
        'report/cutadapt/cutadapt_{ID}{batch}'
    params:
        discard = '--discard-untrimmed ' if config['discard_untrimmed'] is True else '',
        indels = '--no-indels ' if config['no_indels'] is True else '',
        trim_n='--trim-n ' if config['trim_end_n'] is True else ''
    conda:
        "env/dada2.yaml"
    resources: cores = 1, runtime = config['cutadapt_time'], part = config['cutadapt_part'] 
    shell:
        '''
         cutadapt -g file:{input.fwd_prime} -G file:{input.rev_prime} -a file:{input.rev_revcomp_prime} -A file:{input.fwd_revcomp_prime} -o {output.fwd_trim} -p {output.rev_trim} --max-n={config[max_n]} -m {config[min_len_dis]} --error-rate {config[error_rate]} {params.discard}{params.indels}{params.trim_n}--pair-filter {config[pair_filter]} {input.fwd} {input.rev} > {log} 
        ''' 

rule cutadapt_se:
    input:
        fwd='data/seq/{ID}_R1{batch}.fastq.gz',
        fwd_prime='data/metadata/primers/fwd_primers.fasta',
        rev_prime='data/metadata/primers/rev_primers.fasta',
        fwd_revcomp_prime='data/metadata/primers/fwd_revcomp_primers.fasta',
        rev_revcomp_prime='data/metadata/primers/rev_revcomp_primers.fasta'
    output:
        fwd_trim='data/trimmed/{ID}_R1{batch,_B.*}-trimmed.fastq'
    log:
        'report/cutadapt/cutadapt_{ID}_{batch}'
    params:
        discard = '--discard-untrimmed ' if config['discard_untrimmed'] is True else '',
        indels = '--no-indels ' if config['no_indels'] is True else '',
        trim_n='--trim-n ' if config['trim_end_n'] is True else ''
    conda:
        "env/dada2.yaml"
    resources: cores = 1, runtime = config['cutadapt_time'], part = config['cutadapt_part'] 
    shell:
        '''
         cutadapt -g file:{input.fwd_prime} -G file:{input.rev_prime} -a file:{input.rev_revcomp_prime} -A file:{input.fwd_revcomp_prime} -o {output.fwd_trim} --max-n={config[max_n]} -m {config[min_len_dis]} --error-rate {config[error_rate]} {params.discard}{params.indels}{params.trim_n}--pair-filter {config[pair_filter]} {input.fwd} > {log} 
        ''' 

##### DADA FILTER #####
ruleorder: dada2_filter_pe > dada2_filter_se

rule dada2_filter_pe:
    input:
        fwd='data/trimmed/{IDS}_R1{batch}-trimmed.fastq',
        rev='data/trimmed/{IDS}_R2{batch}-trimmed.fastq'
    output:
        fwd='data/filtered/{IDS}_R1{batch,_B.*}-filtered.fastq.gz',
        rev='data/filtered/{IDS}_R2{batch,_B.*}-filtered.fastq.gz'
    conda:
        "env/dada2.yaml"
    threads: config['filter_core']
    resources: runtime = config['filter_time'], part = config['filter_part'], mem_mb = config['filter_mem']
    shell:
        '''
        scripts/dada2_filter.R -f {input.fwd} -r {input.rev} --fwd_out {output.fwd} --rev_out {output.rev} -p data/filtered/ {config[samp_fields]} -q {config[truncq]} -l {config[trunclen]} -L {config[trimleft]} -x {config[maxlen]} -n {config[minlen]} -N {config[maxn]} -d {config[minq]} -E {config[maxee]} --rm_phix {config[rm_phix]} # --n_reads {config[n_reads]}
        '''
        
rule dada2_filter_se:
    input:
        fwd='data/trimmed/{IDS}_R1{batch}-trimmed.fastq'
    output:
        fwd='data/filtered/{IDS}_R1{batch,_B.*}-filtered.fastq.gz'
    conda:
        "env/dada2.yaml"
    threads: config['filter_core']
    resources: runtime = config['filter_time'], part = config['filter_part'], mem_mb = config['filter_mem']
    shell:
        '''
        scripts/dada2_filter.R -f {input.fwd} --fwd_out {output} -p data/filtered/ {config[samp_fields]} -q {config[truncq]} -l {config[trunclen]} -L {config[trimleft]} -x {config[maxlen]} -n {config[minlen]} -N {config[maxn]} -d {config[minq]} -E {config[maxee]} --rm_phix {config[rm_phix]} # --n_reads {config[n_reads]}
        '''
        
##### DADA CORE #####
ruleorder: dada2_pe > dada2_se
rule dada2_pe:
    input:
        #fwd=expand('data/filtered/{IDS}_R1{{batch}}-filtered.fastq.gz', IDS=sorted(set(file_ids.id))),
        #rev=expand('data/filtered/{IDS}_R2{{batch}}-filtered.fastq.gz', IDS=sorted(set(file_ids.id)))
        fwd=lambda wildcards: ['data/filtered/'+i+'_'+r+b+'-filtered.fastq.gz' for i,r,b in zip(file_ids.id, file_ids.read, file_ids.batch) if b == wildcards.batch and r == 'R1'],
        rev=lambda wildcards: ['data/filtered/'+i+'_'+r+b+'-filtered.fastq.gz' for i,r,b in zip(file_ids.id, file_ids.read, file_ids.batch) if b == wildcards.batch and r == 'R2']
    output:
        asv='data/core/{batch,_B.*}_metagenomics_dada2-asv-matrix.tsv'
    params:
        matrix = '' if config['score_matrix'] == 'nuc44' else '--SCORE_MATRIX '+config['score_matrix']
    threads: config['core_core']
    resources: runtime = config['core_time'], part = config['core_part'], mem_mb = config['core_mem']
    conda:
        "env/dada2.yaml"
    shell:
        '''
        scripts/dada2.R -f {input.fwd} -r {input.rev} -o {output.asv} -b {wildcards.batch} {config[samp_fields]} \
        -c {config[chimera_rm]} --LEARN_READS_N {config[learn_reads_n]} --OMEGA_A {config[omega_a]} --USE_QUALS {config[use_quals]}\
        --USE_KMERS {config[use_kmers]} --KDIST_CUTOFF {config[kdist_cutoff]} --BAND_SIZE {config[band_size]} {params.matrix}\
        --GAP_PENALTY {config[gap_penalty]} --HOMOPOLYMER_GAP_PENALTY {config[homopolymer_gap_penalty]} --MIN_FOLD {config[min_fold]}\
        --MIN_HAMMING {config[min_hamming]} --MAX_CLUST {config[max_clust]} --MAX_CONSIST {config[max_consist]}
        '''

rule dada2_se:
    input:
        fwd=lambda wildcards: ['data/filtered/'+i+'_'+r+b+'-filtered.fastq.gz' for i,r,b in zip(file_ids.id, file_ids.read, file_ids.batch) if b == wildcards.batch and r == 'R1']
    output:
        asv='data/core/{batch,_B.*}_metagenomics_dada2-asv-matrix.tsv'
    params:
        matrix = '' if config['score_matrix'] == 'nuc44' else '--SCORE_MATRIX '+config['score_matrix']
    threads: config['core_core']
    conda:
        "env/dada2.yaml"
    threads: config['core_core']
    resources: runtime = config['core_time'], part = config['core_part'], mem_mb = config['core_mem']
    shell:
        '''
        scripts/dada2.R -f {input.fwd} -o {output.asv} -b {wildcards.batch} {config[samp_fields]} \
        -c {config[chimera_rm]} --LEARN_READS_N {config[learn_reads_n]} --OMEGA_A {config[omega_a]} --USE_QUALS {config[use_quals]}\
        --USE_KMERS {config[use_kmers]} --KDIST_CUTOFF {config[kdist_cutoff]} --BAND_SIZE {config[band_size]} {params.matrix}\
        --GAP_PENALTY {config[gap_penalty]} --HOMOPOLYMER_GAP_PENALTY {config[homopolymer_gap_penalty]} --MIN_FOLD {config[min_fold]}\
        --MIN_HAMMING {config[min_hamming]} --MAX_CLUST {config[max_clust]} --MAX_CONSIST {config[max_consist]}
        '''

rule merge_batches:
    input:
        asv=expand('data/core/{BATCH}_metagenomics_dada2-asv-matrix.tsv', BATCH=sorted(set(file_ids.batch)))
    output:
        asv='data/{GLDS}_batch-merge_metagenomics_dada2-asv-matrix.tsv'

    resources: runtime = config['biom_time'], part = config['biom_part']
    threads: 1
    conda:
        "env/dada2.yaml"
    shell:
        '''
        scripts/batch_merge.R -a {input.asv} -o {output.asv}
        '''

rule assign_tax:
    input:
        asv='data/{GLDS}_batch-merge_metagenomics_dada2-asv-matrix.tsv'
    output:
        tax_table='data/{GLDS}_metagenomics_dada2-taxonomy-matrix.tsv',
        phylo='data/taxa/{GLDS}_metagenomics_dada2-taxonomy-phylo.Rds',
        fasta='data/{GLDS}_metagenomics_dada2-taxonomy-seqs.fasta'

    resources: runtime = config['core_time'], part = config['core_part'], mem_mb = config['core_mem']
    threads: config['core_core']
    conda:
        "env/dada2.yaml"
    shell:
        '''
        scripts/dada2_assign_taxa.R -a {input.asv} -o {output.tax_table} -O {output.phylo} -f {output.fasta}\
        -t {config[taxa_train]} -s {config[species_train]} -m {config[train_method]}  --minBoot {config[minboot]} --tryRC {config[tryrc]}\
        --taxLevels {config[taxlevels]} --allowMultiple {config[allowmultiple]} --multithread {threads}
         '''

rule align:
    input:
        fasta='data/{GLDS}_metagenomics_dada2-taxonomy-seqs.fasta'
    output:
        align='data/phy/{GLDS}_metagenomics_mafft-msa.aln'
    resources: runtime = config['aln_time'], part = config['aln_part'], mem_mb = config['aln_mem']
    threads: config['aln_core']
    conda:
        "env/phy.yaml"
    shell:
        '''
        mafft --thread {threads} --auto {input.fasta} > {output.align}
        '''

rule tree:
    input:
        align='data/phy/{GLDS}_metagenomics_mafft-msa.aln'
    output:
        tree='data/{GLDS}_metagenomics_fastree-tree-newick.tree'
    resources: runtime = config['tree_time'], part = config['tree_part'], mem_mb = config['tree_mem']
    threads: config['tree_core']
    log:
        'report/{GLDS}_fastree'
    conda:
        "env/phy.yaml"
    shell:
        '''
        export OMP_NUM_THREADS={threads}
        fasttree -nt -log {log} < {input.align} > {output.tree}
        '''

rule DADA2BIOM:
    input:
        phylo='data/taxa/{GLDS}_metagenomics_dada2-taxonomy-phylo.Rds',
        tree='data/{GLDS}_metagenomics_fastree-tree-newick.tree'
    output:
        #biom='data/{GLDS}_metagenomics_dada2-biom.json',
        phylo='data/{GLDS}_metagenomics_dada2-complete-phylo.Rds'
    conda:
        "env/dada2.yaml"
    resources: cores = 1, runtime = config['biom_time'], part = config['biom_part'],
    params:
        meta = '-s {config["meta_samp_file"]} ' if config['meta_samp_file'] is not False else ""
    shell:
        '''
        scripts/dada2biom_phylo.R -i {input.phylo} -t {input.tree} -p {output.phylo} -n {DS_NUM} {params.meta}
        '''

