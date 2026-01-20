#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_fastqc {
  tag { "fastqc:${sample_id}" }
  publishDir "${params.fastqc_reports}", mode: 'symlink', overwrite: true

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("*_fastqc.zip"),  emit: zip
    tuple val(sample_id), path("*_fastqc.html"), emit: html

  script:
  """
  fastqc -t ${task.cpus} ${reads.join(' ')}
  """
}