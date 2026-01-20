#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_multiqc {
  tag { "Multiqc:QC:${params.gse}" }
  publishDir "${params.multiqc_res}", mode: 'copy', overwrite: true

  input:
    path "inputs/*"

  output:
    path "*.html"
    path "*data"

  script:
  """
  multiqc . --force --filename qc_report.html
  """
}