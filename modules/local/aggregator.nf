process AGGREGATOR {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/aggregated", mode: 'copy'
    container 'docker.io/python:3.9'

    input:
    tuple val(meta), path(trim_log), path(fastqc_files), path(quast_dir), path(mlst_tsv), path(abricate_tabs), path(amrfinder_report), path(mash_sketch), path(spa_report), path(sccmec_report), path(agr_report), path(kma_res), path(metadata_json)

    output:
    tuple val(meta), path("${meta.id}_report.json"), emit: report
    tuple val(meta), path("${meta.id}_summary.csv"), emit: summary

    script:
    """
    python3 ${projectDir}/bin/staphit-aggregate --sample-id ${meta.id} --trim-log ${trim_log} --fastqc-dir . --quast-dir ${quast_dir} --mlst ${mlst_tsv} --abricate-dir . --amrfinder ${amrfinder_report} --mash ${mash_sketch} --spa ${spa_report} --sccmec ${sccmec_report} --agr ${agr_report} --kma ${kma_res} --metadata ${metadata_json} --outdir .
    """
}
