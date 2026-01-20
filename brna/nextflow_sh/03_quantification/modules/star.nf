#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_star {
  tag { "STAR:${sample_id}" }
  label 'big_mem'
  publishDir "${params.star_res}", mode: 'symlink', overwrite: true

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id), path("${sample_id}.bam"),       emit: bam
    tuple val(sample_id), path("${sample_id}.bam.bai"),   emit: bai
    tuple val(sample_id), path("*Log.final.out"),         emit: log

  script:
  """
  STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${params.index} \
    --readFilesIn ${r1} ${r2} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${sample_id}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes Standard \
    --twopassMode Basic \
    --outFilterMultimapNmax ${params.star_multimap} \
    --outFilterMismatchNmax ${params.star_mismatch}

  mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam
  samtools index -@ ${task.cpus} ${sample_id}.bam
  """
}