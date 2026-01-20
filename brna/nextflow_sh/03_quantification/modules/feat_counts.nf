process run_fc {
  tag { "FC:${sample_id}" }
  label 'counts'
  publishDir "${params.feat_counts}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(bam)

  output:
    tuple val(sample_id), path("${sample_id}_featureCounts.txt"),         emit: counts
    tuple val(sample_id), path("${sample_id}_featureCounts.txt.summary"), emit: log

  script:
  """
  # Create subset 
  samtools view -h ${bam} | head -n 200000 | \
  samtools view -bS - > ${sample_id}.subset.bam

  # Run Tests
  s1_rate=\$(featureCounts \
    -T ${task.cpus} \
    -a ${params.gtf_file} \
    -t exon -g gene_id \
    -p -B -C -s 1 ${sample_id}.subset.bam 2>&1 \
    | grep "Successfully assigned" | awk '{print \$6}' | sed 's/%//')

  s2_rate=\$(featureCounts \
    -T ${task.cpus} \
    -a ${params.gtf_file} \
    -t exon -g gene_id \
    -p -B -C -s 2 ${sample_id}.subset.bam 2>&1 \
    | grep "Successfully assigned" | awk '{print \$6}' | sed 's/%//')

  # Logic using AWK (Safer than 'bc' in containers)
  detected_strand=\$(awk -v s1="\$s1_rate" -v s2="\$s2_rate" 'BEGIN {
    if (s1 > s2 + 10) print "1";
    else if (s2 > s1 + 10) print "2";
    else print "0";
  }')

  # Final Run
  featureCounts \
    -T ${task.cpus} \
    -a ${params.gtf_file} \
    -o ${sample_id}_featureCounts.txt \
    -t exon \
    -g gene_id \
    --extraAttributes gene_name \
    -p -B -C \
    -s \$detected_strand \
    -Q 10 \
    ${bam}

  # Cleanup
  rm ${sample_id}.subset.bam
  """
}