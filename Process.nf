#! /usr/bin/env nextflow

/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/*_{1,2}.fastq.gz"
params.index_files = "$projectDir/data/index"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"
log.info """\
    F O C Y T E   P I P E L I N E - V A R I A N T - P R O C E S S
    =============================================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

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

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    
    FASTP(read_pairs_ch)
}