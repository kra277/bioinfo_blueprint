#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include {get_sra_ids}       from './modules/get_sra_ids.nf'
include {download_sra}      from './modules/dwnld_sra.nf'
include {run_fastqc}        from './modules/fastqc.nf'
include {run_multiqc}       from './modules/multiqc.nf'

// Define channels
gse_ch = Channel.of(params.gse)

workflow {
    // Use pysradb to get SRR list from GSE
    get_sra_ids(gse_ch)

    // Prepare input for sratools
    sra_ids_ch = get_sra_ids.out.sra_ids
      .splitText { it.trim() } // Split AND trim each line
      .filter { it }           // Remove any empty lines

    // Download SRA and convert to FASTQ
    download_sra(sra_ids_ch)
    
    // Run FastQC
    run_fastqc(download_sra.out.reads)

    // Prepare input for MultiQC
    file_log = run_fastqc.out.zip
      .map { id, files -> files } 
      .flatten()                  // Ensures every file is an individual item
      .collect()                  // Puts them all into one big box

    // Run MultiQC
    run_multiqc(file_log)   
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