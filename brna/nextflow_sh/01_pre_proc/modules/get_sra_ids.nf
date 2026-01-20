#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process get_sra_ids {
    tag "pysradb:${gse_id}"
    label 'pysradb'
    publishDir "${params.data}", mode: 'symlink', overwrite: false

    input:
    val gse_id

    output:
    path("*.txt"), emit: sra_ids
    path "*.tsv", emit: metadata

    script:
    def metadata = "${gse_id}_metadata.tsv"
    def sra_list = "${gse_id}_sra_ids.txt"
    def gse_uc   = gse_id.toString().toUpperCase()
    """
    # store sra id in a variable, check if sra id is empty
    SRP=\$(pysradb gse-to-srp ${gse_uc} | awk 'NR==2 {print \$2}')
    [[ -z "\$SRP" ]] && exit 1

    # using psyradb get metadata
    pysradb metadata "\$SRP" --detailed > ${metadata}
    
    # get the first column of the metadata which is the run_accession id
    awk -F'\\t' 'NR>1 {print \$1}' ${metadata} > ${sra_list}
    """
}