import pandas as pd
import yaml
from pathlib import Path
import os
from scripts import getRunIds

BASE_DIR = Path(workflow.basedir)
configfile: str(BASE_DIR) + "/config/config.yaml"

INPUT = config['fastq_path']
OUTPUT = config['output_path']
REFERENCE = config['ref_path']
ANNOTATIONS = config['annotations_path']

runIds = getRunIds.getRids(INPUT)

rule all:
    input:
        OUTPUT + 'references/geneNames.csv',
        OUTPUT + 'references/reference.mmi',
        OUTPUT + 'references/reference.fa',
        OUTPUT + 'talon/database.db',
        OUTPUT + 'talon/config.csv',
        OUTPUT + 'talon/talon.done',
        OUTPUT + 'talon/talon_summarize.done',
        OUTPUT + 'talon/talon_abundance.done',
        expand(f"{OUTPUT}fastqc/{{runId}}.html", zip, runId=runIds),  
        expand(f"{OUTPUT}alignments/{{runId}}.sorted.bam.bai", zip, runId=runIds),   
        expand(f"{OUTPUT}bedfiles/{{runId}}.bed", zip, runId=runIds),   
        expand(f"{OUTPUT}stats/{{runId}}.bamstats", zip, runId=runIds),  
        expand(f"{OUTPUT}counts/{{runId}}.counts.txt", zip, runId=runIds),   
        expand(f"{OUTPUT}counts/{{runId}}.featureCounts", runId=runIds),
        expand(f"{OUTPUT}NanoCount/{{runId}}.tsv", runId=runIds),
        expand(f"{OUTPUT}talon/alignments/{{runId}}_labeled.sam", zip, runId=runIds),  


rule get_gene_names:
    input:
        annotations=ANNOTATIONS,
    output:
        OUTPUT + 'references/geneNames.csv',
    shell:
        "python scripts/getGeneNames.py {input.annotations} {output}"


rule fastqc:
    input:
        fastq=INPUT + '{runId}.fastq.gz',
    output:
        html=OUTPUT + "fastqc/{runId}.html",
        zip=OUTPUT + "fastqc/{runId}_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/{runId}.log"
    threads:
        8
    wrapper:
        "v1.29.0/bio/fastqc"

        
rule add_reference:
    input:
        refgenome=REFERENCE,
    output:
        OUTPUT + 'references/reference.fa.gz'
    shell:
        "cp {input} {output}"
        
        
rule unzip_reference:
    input:
        refgenome=OUTPUT + 'references/reference.fa.gz'
    output:
        OUTPUT + 'references/reference.fa'
    shell:
        "cat {input} | gzip -d > {output}"


rule minimap2_index:
    input:
        refgenome=OUTPUT + 'references/reference.fa.gz'
    output:
        OUTPUT + 'references/reference.mmi'
    shell:
        "minimap2 -d {output} {input.refgenome}"
        

rule minimap2_align:
   input:
       fastq=INPUT + '{runId}.fastq.gz',
       refgenome=OUTPUT + 'references/reference.fa.gz',
       refindex=OUTPUT + 'references/reference.mmi',
   output:        
       OUTPUT + 'alignments/{runId}.sam'
   params:
       args=config['minimap2_args'],
       threads=8
   wildcard_constraints:
        runId='|'.join([re.escape(x) for x in set(runIds)]),
   shell:
       "minimap2 {params.args} -t {params.threads} {input.refgenome} {input.fastq} > {output}"
 

rule make_bam:
    input:
        OUTPUT + 'alignments/{runId}.sam'
    output:       
        OUTPUT + 'alignments/{runId}.bam'
    wildcard_constraints:
        runId='|'.join([re.escape(x) for x in set(runIds)]),
    shell:
        "samtools view -Sb {input} > {output}"


rule bamtobed:
    input:
        OUTPUT + 'alignments/{runId}.bam'
    output:
        OUTPUT + 'bedfiles/{runId}.bed'
    log:
        "logs/bamtobed/{runId}.log",
    params:
        extra="-bedpe",  # optional parameters
    wrapper:
        "v2.2.1/bio/bedtools/bamtobed"


rule jellyfish_count:
    input:
        fastq=INPUT + '{runId}.fastq.gz',
    output:
        OUTPUT + 'jellyfish/{runId}.jf'
    log:
        OUTPUT + "logs/{runId}.jf.log",
    params:
        kmer_length=21,
        size="1G",
        extra="--canonical",
    threads:
        8
    wrapper:
        "v2.1.1/bio/jellyfish/count"


rule jellyfish_histo:
    input:
         OUTPUT + 'jellyfish/{runId}.jf'
    output:
         OUTPUT + 'jellyfish/{runId}.histo'
    log:
         OUTPUT + "logs/{runId}.histo.jf.log",
    threads:
        config['threads']
    wrapper:
        "v2.1.1/bio/jellyfish/histo"


rule bamtools_stats:
    input:
        OUTPUT + 'alignments/{runId}.bam'
    output:
        OUTPUT + 'stats/{runId}.bamstats'
    params:
        "-insert" # optional summarize insert size data
    log:
        "logs/bamtools/stats/{runId}.log"
    wrapper:
        "v2.1.1/bio/bamtools/stats"
        


rule feature_counts:
    input:
        samples=OUTPUT + "alignments/{runId}.bam",
        annotation=ANNOTATIONS,
    output:
        multiext(
            OUTPUT + "counts/{runId}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 
        config['threads']
    params:
        r_path="",  # implicitly sets the --Rpath flag
    log:
        "logs/{runId}.log",
    wrapper:
        "v1.29.0/bio/subread/featurecounts"


rule samtools_sort:
    input:
        OUTPUT + 'alignments/{runId}.bam'
    output:
        OUTPUT + 'alignments/{runId}.sorted.bam'
    wildcard_constraints:
        runId='|'.join([re.escape(x) for x in set(runIds)]),
    shell:
        "samtools sort -T {input} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        OUTPUT + 'alignments/{runId}.sorted.bam'
    output:
        OUTPUT + 'alignments/{runId}.sorted.bam.bai'
    wildcard_constraints:
        runId='|'.join([re.escape(x) for x in set(runIds)]),
    shell:
        "samtools index {input}"


rule nanocount:
    input:
        bam=OUTPUT + 'alignments/{runId}.sorted.bam',
    output:
        OUTPUT + 'NanoCount/{runId}.tsv',
    params:
    shell:
        "NanoCount -i {input.bam} -o {output}"
        
        
rule htseq_count:
    input:
        bam=OUTPUT + 'alignments/{runId}.sorted.bam',
        annotations=ANNOTATIONS,
    output:
        OUTPUT + "counts/{runId}.counts.txt"
    wildcard_constraints:
        runId='|'.join([re.escape(x) for x in set(runIds)]),
    params:
        minQual=int(config['minQual'])
    shell:
        "htseq-count -f bam -a {params.minQual} {input.bam} {input.annotations} > {output}"


rule talon_build:
    input:
        gtf=ANNOTATIONS,
    output:
        OUTPUT + 'talon/database.db'
    params:
        outpath_root=OUTPUT + 'talon/database',
        db_name=config['talon_db_name']
    shell:
        "talon_initialize_database --f {input.gtf} --g {params.db_name} --a {params.db_name} --o {params.outpath_root}"


rule talon_label_reads:
    input:
        sam=OUTPUT + 'alignments/{runId}.sam',
        ref=OUTPUT + 'references/reference.fa',
    output:
        OUTPUT + 'talon/alignments/{runId}_labeled.sam',
        OUTPUT + 'talon/alignments/{runId}_read_labels.tsv'
    threads:
        8
    params:
        outpath_root=OUTPUT + 'talon/alignments/{runId}',
        tmp=OUTPUT + 'talon/temp/{runId}'
    shell:
        "talon_label_reads --f {input.sam} --g {input.ref} --t {threads} --tmpDir {params.tmp} --deleteTmp --o {params.outpath_root}"

rule make_talon_config:
    input:
        expand(f"{OUTPUT}talon/alignments/{{runId}}_labeled.sam", zip, runId=runIds),  
    output:
        OUTPUT + 'talon/config.csv',
    params:
        platform=config['talon_platform'],
        sample=config['talon_sample_type'],
    shell:
        "python scripts/makeTalonConfig.py {params.platform} {params.sample} {output} {input}"

rule run_talon:
    input:
        db=OUTPUT + 'talon/database.db',
        cf=OUTPUT + 'talon/config.csv',
    output:
        touch(OUTPUT + 'talon/talon.done')
    params:
        db_name=config['talon_db_name'],
        prefix= OUTPUT + 'talon/talon_run',
        tmp=OUTPUT + 'talon/temp/'
    threads:
        16
    shell:
        "talon --f {input.cf} --db {input.db} --build {params.db_name} -t {threads} --tmpDir {params.tmp} --o {params.prefix}"

rule run_talon_summarize:
    input:
        db=OUTPUT + 'talon/database.db',
        flag=OUTPUT + 'talon/talon.done',
    output:
        touch(OUTPUT + 'talon/talon_summarize.done')
    params:
        prefix= OUTPUT + 'talon/talon_summarize',
    threads:
        16
    shell:
        "talon_summarize --db {input.db} --o {params.prefix}"


rule run_talon_abundance:
    input:
        db=OUTPUT + 'talon/database.db',
        flag=OUTPUT + 'talon/talon_summarize.done',
    output:
        touch(OUTPUT + 'talon/talon_abundance.done')
    params:
        db_name=config['talon_db_name'],
        prefix= OUTPUT + 'talon/talon_abundance',
    threads:
        16
    shell:
        "talon_abundance --db {input.db} --build {params.db_name} --a {params.db_name} --o {params.prefix}"
        