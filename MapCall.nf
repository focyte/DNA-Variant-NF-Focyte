#! /usr/bin/env nextflow

/*
 * pipeline input parameters
 */
params.reads = "$projectDir/results/*_{1,2}.fastq.gz"
params.transcriptome_file = "$projectDir/data/ref_genome/genome.fasta"
params.outdir = "results"
log.info """\
    F O C Y T E   P I P E L I N E - V A R I A N T - M A P P I N G - C A L L
    =======================================================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * define the `INDEX` process
 */
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

/*
 * define the `MAPPING` process
 */
process MAPPING {
    tag "BWA on $sample_id"
    publishDir params.outdir, mode: 'copy'
    cpus 12

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

/*
 * BAMCONVERT Process
 */
process BAMCONVERT {
    tag "Sam to Bam on $sam_file"
    publishDir params.outdir, mode: 'copy'

    input:
    path sam_file

    output:
    path "${sam_file.baseName}.bam"

    script:
    """
    samtools view -h -S -b -o ${sam_file.baseName}.bam ${sam_file}
    """
}

/*
 * BAMSORT Process
 */
process BAMSORT {
    tag "Sorting $bam_file"
    publishDir params.outdir, mode: 'copy'

    input:
    path bam_file

    output:
    path "${bam_file.baseName}_sorted.bam"
    

    script:
    """
    samtools sort ${bam_file} -o ${bam_file.baseName}_sorted.bam 
    """
}

/*
 * BAMINDEX Process
 */
process BAMINDEX {
    tag "Indexing $sorted_bam_file"
    publishDir params.outdir, mode: 'copy'

    input:
    path sorted_bam_file

    output:
    path "${sorted_bam_file.baseName}.bam.bai"

    script:
    """
    samtools index ${sorted_bam_file}
    """
}

/*
 * BAMFLAG Process
 */
process BAMFLAG {
    tag "Flagging $sorted_bam_file"
    publishDir params.outdir, mode: 'copy'

    input:
    path sorted_bam_file

    output:
    path "${sorted_bam_file.baseName}.bam"

    script:
    """
    samtools flagstat ${sorted_bam_file}
    """
}

/*
 * BCFPILEUP Process
 */
process BCFPILEUP {
    tag "Piling $flagged_bam_file"
    publishDir params.outdir, mode: 'copy'

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

/*
 * BCFCALL Process
 */
process BCFCALL {
    tag "Calling $_file"
    publishDir params.outdir, mode: 'copy'

    input:
    path flagged_bcf_file

    output:
    path "${flagged_bcf_file.baseName}.vcf"

    script:
    """
    bcftools call --ploidy 1 -m -v -o ${flagged_bcf_file.baseName}.vcf ${flagged_bcf_file}
    """
}

/*
 * Workflow Definition
 */
workflow {
    read_files_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true, flat: true)
    
    index_ch = INDEX(params.transcriptome_file)
    
    mapping_ch = MAPPING(read_files_ch, index_ch)

    bam_ch = BAMCONVERT(mapping_ch)

    sorted_bam_ch = BAMSORT(bam_ch)

    BAMINDEX(sorted_bam_ch)

    flagged_bam_ch = BAMFLAG(sorted_bam_ch)

    piled_bam_ch = BCFPILEUP(flagged_bam_ch, index_ch)

    called_bam_ch = BCFCALL(piled_bam_ch)
}