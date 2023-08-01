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
        expand(f"{OUTPUT}fastqc/{{runId}}.html", zip, runId=runIds),  
        expand(f"{OUTPUT}alignments/{{runId}}.sorted.bam.bai", zip, runId=runIds),   
        expand(f"{OUTPUT}stats/{{runId}}.bamstats", zip, runId=runIds),  
        expand(f"{OUTPUT}counts/{{runId}}.counts.txt", zip, runId=runIds),   
        expand(f"{OUTPUT}counts/{{runId}}.featureCounts", runId=runIds),
        expand(f"{OUTPUT}NanoCount/{{runId}}.tsv", runId=runIds),
        expand(f"{OUTPUT}jellyfish/{{runId}}.histo", zip, runId=runIds),  

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


rule jellyfish_count:
    input:
        fastq=INPUT + '{runId}.fastq.gz',
    output:
        OUTPUT + 'jellyfish/{runId}.jf'
    log:
        OUTPUT + "logs/{runId}.jf.log",
    params:
        kmer_length=21,
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