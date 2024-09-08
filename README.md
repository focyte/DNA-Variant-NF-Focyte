# DNA-Variant-NF-Focyte
A Nextflow pipeline for processing and calling of variants: FASTQ to VCF

Welcome to this **FOCYTE Pipeline** repository! This pipeline is designed to process FASTQ files, perform quality control, align reads to a reference genome, and call genetic variants. The workflow uses [Nextflow](https://www.nextflow.io/) to orchestrate the different steps of the analysis. Below you'll find details about the stages involved in this pipeline along with usage instructions.

---

## Table of Contents
1. [Overview](#overview)
2. [Pipeline Stages](#pipeline-stages)
    - [FASTQ Quality Control](#1-fastq-quality-control)
    - [Mapping Reads to Reference Genome](#2mapping-reads-to-reference-genome)
    - [Variant Calling](#3variant-calling)
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

