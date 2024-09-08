# DNA-Variant-NF-Focyte
A Nextflow pipeline for processing and calling of variants: FASTQ to VCF

Welcome to this **FOCYTE Pipeline** repository! This pipeline is designed to process FASTQ files, perform quality control, align reads to a reference genome, and call genetic variants. The workflow uses [Nextflow](https://www.nextflow.io/) to orchestrate the different steps of the analysis. Below you'll find details about the stages involved in this pipeline along with usage instructions.

---

## Table of Contents
1. [Overview](#overview)
2. [Pipeline Stages](#pipeline-stages)
    - [FASTQ Quality Control](#fastq-quality-control)
    - [Mapping Reads to Reference Genome](#mapping-reads-to-reference-genome)
    - [Variant Calling](#variant-calling)
3. [How to Run](#how-to-run)
4. [Requirements](#requirements)

---

## Overview

The pipeline processes paired-end FASTQ files and generates a VCF file containing the variants. It consists of two main phases:
1. **Quality Control & Trimming:** Raw reads are cleaned and trimmed using `fastp`.
2. **Variant Mapping & Calling:** Reads are mapped to a reference genome using `bwa` and variants are called using `bcftools`.

---

## Pipeline Stages

### 1. FASTQ Quality Control

This process uses the tool `fastp` to clean and trim raw FASTQ files, generating high-quality reads for downstream analysis. Each sample's paired-end reads are taken as input, and the output includes:
- Trimmed FASTQ files
- A quality control report (JSON and HTML formats)

#### Code Snippet:
```nextflow
process FASTP {
    tag "FASTP on $sample_id"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path('trim_*.fastq.gz'), emit: reads
    tuple val(sample_id), path("${sample_id}.json"), emit: json
    tuple val(sample_id), path("${sample_id}.html"), emit: html

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
          -o trim_${reads[0]} -O trim_${reads[1]} \
          --json ${sample_id}.json \
          --html ${sample_id}.html
    """
} 
```

### 2. Mapping Reads to Reference Genome

The cleaned FASTQ files are aligned to a reference genome using `bwa`, and the pipeline converts the SAM file into BAM format, sorts it, and indexes it for further analysis.

#### Code Snippet:
```nextflow
process INDEX {
    input:
    path transcriptome

    output:
    path "genome*"

    script:
    """
    bwa index $transcriptome
    """
}

process MAPPING {
    tag "BWA on $sample_id"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path index

    output:
    path "${sample_id}.sam"

    script:
    """
    bwa mem ${params.transcriptome_file} ${reads1} ${reads2} > ${sample_id}.sam
    """
}
```

The BAM file is generated, sorted, and indexed for variant calling:

#### Code Snippet:
```nextflow
process BAMCONVERT {
    input:
    path sam_file

    output:
    path "${sam_file.baseName}.bam"

    script:
    """
    samtools view -h -S -b -o ${sam_file.baseName}.bam ${sam_file}
    """
}

process BAMSORT {
    input:
    path bam_file

    output:
    path "${bam_file.baseName}_sorted.bam"

    script:
    """
    samtools sort ${bam_file} -o ${bam_file.baseName}_sorted.bam
    """
}

process BAMINDEX {
    input:
    path sorted_bam_file

    output:
    path "${sorted_bam_file.baseName}.bam.bai"

    script:
    """
    samtools index ${sorted_bam_file}
    """
}
```
