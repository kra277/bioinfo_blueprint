#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process run_fastp {
  tag { "fastp:${sample_id}" }
  publishDir "${params.fastp_res}", mode: 'symlink', overwrite: true

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"),
                         path("${sample_id}_R2.trim.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}.fastp.json"),  emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"),  emit: html
    tuple val(sample_id), path("${sample_id}.fastp.log"),   emit: log

  script:
  """
  fastp \
    -i ${reads[0]} -I ${reads[1]} \
    -o ${sample_id}_R1.trim.fastq.gz \
    -O ${sample_id}_R2.trim.fastq.gz \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --length_required 30 \
    --correction \
    --qualified_quality_phred 20 \
    --cut_tail --cut_tail_mean_quality 20 \
    -w ${task.cpus} \
    -j ${sample_id}.fastp.json \
    -h ${sample_id}.fastp.html \
    2> ${sample_id}.fastp.log
  """
}