#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process download_sra {
    tag { "sra:${sample_id}" }
    label 'sra_dwnld'
    publishDir "${params.fastq_dir}", mode: 'symlink', overwrite: false

    input:
      val sample_id

    output:
    	tuple val(sample_id),
          path("${sample_id}_*.fastq.gz"),
          emit: reads

    script:
    """
    # Download the .sra file into the current directory
    prefetch --output-directory . ${sample_id}

    # Find the first SRA file for this sample that prefetch downloaded
    sra_file=\$(find . -maxdepth 3 -type f -name "${sample_id}.sra" | head -n 1)

    # Convert to single-end FASTQ in the current directory
    fasterq-dump --split-files --threads ${task.cpus} --outdir . "\$sra_file"

    # Normalize single-end naming to *_1.fastq so downstream stays consistent
    if [[ -f "${sample_id}.fastq" && ! -f "${sample_id}_1.fastq" ]]; then
      mv "${sample_id}.fastq" "${sample_id}_1.fastq"
    fi

    # Compress the FASTQ file
    pigz -p ${task.cpus} -f ${sample_id}*.fastq
    """
}