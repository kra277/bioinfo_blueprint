#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {run_fastp}         from './modules/fastp.nf'
include {run_star}          from './modules/star.nf'
include {run_fc}            from './modules/feat_counts.nf'
include {run_multiqc}       from './modules/multiqc.nf'

// Define channels
// config file has reads = "${params.fastq_dir}/*_{1,2}.fastq.gz"
reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

// Set workflow
workflow {
  run_fastp(reads_ch)
  run_star(run_fastp.out.reads)
  run_fc(run_star.out.bam)

  Channel.empty()
        .mix(run_fastp.out.json.map {it[1]}) 
        .mix(run_star.out.log.map { it[1] })
        .mix(run_fc.out.log.map { it[1] })
        .collect()
        .set { log_files }
    run_multiqc(log_files)
}

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}
