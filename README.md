# DNA Variant Analysis Pipeline
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
4. [Results](#results)

---

## Overview

The pipeline processes paired-end FASTQ files and generates a VCF file containing the variants. It consists of two main phases:
1. **Quality Control & Trimming:** Raw reads are cleaned and trimmed using `fastp`.
2. **Variant Mapping & Calling:** Reads are mapped to a reference genome using `bwa` and variants are called using `bcftools`.

---

## Pipeline Stages

### FASTQ Quality Control

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

### Mapping Reads to Reference Genome

The cleaned FASTQ files are aligned to a reference genome using `bwa`, and the pipeline converts the SAM file into BAM format, sorts it, and indexes it for further analysis.

The genome is indexed and this is used to align the trimmed reads from the Process pipeline to generate files in SAM format.

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

### Variant Calling

In this stage, variants are called using `bcftools`. The pipeline uses `mpileup` to create a BCF file, which is then processed to generate a VCF file with the called variants.

#### Code Snippet:
```nextflow
process BCFPILEUP {
    input:
    path flagged_bam_file
    path index

    output:
    path "${flagged_bam_file.baseName}.bcf"

    script:
    """
    bcftools mpileup -O b -o ${flagged_bam_file.baseName}.bcf -f ${params.transcriptome_file} ${flagged_bam_file}
    """
}

process BCFCALL {
    input:
    path flagged_bcf_file

    output:
    path "${flagged_bcf_file.baseName}.vcf"

    script:
    """
    bcftools call --ploidy 1 -m -v -o ${flagged_bcf_file.baseName}.vcf ${flagged_bcf_file}
    """
}
```

### How to Run

Clone this repository and navigate to the directory
Within your project directory create the following folders:

```console
mkdir ./data/ref-genome
```

In data, place your paired .fastq files for analysis
In ref-genome, place a FASTA file of your genome of interest

Create your conda environment containing the required tools and dependencies by loading the .yml file:

```console
conda env create -f variant.yml
```

Run the Procss.nf pipeline:

```console
nextflow Process.nf
```

Check the read quailty after trimming

Run the MapCall.nf pipeline:

```console
nextflow MapCall.nf
```

### Results

Clinical isolates from patients with SARS-CoV-2 infection collected by The COVID-19 Genomics UK (COG-UK) – Consortium were obtained from [The European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/PRJEB37886?show=analyses). Eight samples were randomly selected for variant calling compared to the [Wuhan-Hu-1 sequence GCA_009858895.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009858895.2/) using this pipeline.

VCF files were visualised against the Wuhan-Hu-1 genome assembly in IGV:

![IGV](https://github.com/focyte/nf-DNA-Variants/blob/main/SARS-CoV-2.png?raw=true)
