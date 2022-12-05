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

samples, runIds = getRunIds.getRids(INPUT)


rule all:
    input:
        OUTPUT + 'references/geneNames.csv',
        OUTPUT + 'references/minimap2_index.mmi',
        OUTPUT + 'references/reference.fa',
        expand(f"{OUTPUT}alignments/{{sample}}_{{runId}}.sorted.bam.bai", zip, sample=samples,runId=runIds),   
        expand(f"{OUTPUT}counts/{{sample}}_{{runId}}.counts.txt", zip, sample=samples,runId=runIds),   
        
        
        
rule get_gene_names:
    input:
        annotations=ANNOTATIONS,
    output:
        OUTPUT + 'references/geneNames.csv',
    shell:
        "python scripts/getGeneNames.py {input.annotations} {output}"
        
        
        
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
        OUTPUT + 'references/minimap2_index.mmi'
    shell:
        "minimap2 -d {output} {input.refgenome}"
        

rule minimap2_align:
   input:
       fastq=INPUT + '{sample}_{runId}.fastq.gz',
       refgenome=REFERENCE,
   output:        
       OUTPUT + 'alignments/{sample}_{runId}.sam'
   params:
       args=config['minimap2_args'],
       threads=config['threads']
   wildcard_constraints:
        sample='|'.join([re.escape(x) for x in set(samples)]),
        runId='|'.join([re.escape(x) for x in set(runIds)]),
   shell:
       "minimap2 {params.args} -t {params.threads} {input.refgenome} {input.fastq} > {output}"
 

rule make_bam:
    input:
        OUTPUT + 'alignments/{sample}_{runId}.sam'
    output:       
        OUTPUT + 'alignments/{sample}_{runId}.bam'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in set(samples)]),
        runId='|'.join([re.escape(x) for x in set(runIds)]),
    shell:
        "samtools view -Sb {input} > {output}"
        


rule samtools_sort:
    input:
        OUTPUT + 'alignments/{sample}_{runId}.bam'
    output:
        OUTPUT + 'alignments/{sample}_{runId}.sorted.bam'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in set(samples)]),
        runId='|'.join([re.escape(x) for x in set(runIds)]),
    shell:
        "samtools sort -T {input} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        OUTPUT + 'alignments/{sample}_{runId}.sorted.bam'
    output:
        OUTPUT + 'alignments/{sample}_{runId}.sorted.bam.bai'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in set(samples)]),
        runId='|'.join([re.escape(x) for x in set(runIds)]),
    shell:
        "samtools index {input}"
        
        
rule htseq_count:
    input:
        bam=OUTPUT + 'alignments/{sample}_{runId}.sorted.bam',
        annotations=ANNOTATIONS,
    output:
        OUTPUT + "counts/{sample}_{runId}.counts.txt"
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in set(samples)]),
        runId='|'.join([re.escape(x) for x in set(runIds)]),
    shell:
        "htseq-count -f bam {input.bam} {input.annotations} > {output}"