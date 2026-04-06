process AGGREGATOR {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/aggregated", mode: 'copy'
    stageInMode 'copy'
    container 'docker.io/python:3.9'

    input:
    tuple val(meta),
          path(trim_log),
          path(fastqc_files),
          path(quast_dir),
          path(mlst_tsv, stageAs: 'mlst_input.tsv'),
          path(abricate_tabs),
          path(amrfinder_report, stageAs: 'amrfinder_input.tsv'),
          path(mash_sketch),
          path(spa_report),
          path(sccmec_report),
          path(agr_report),
          path(kma_res, stageAs: 'kma_input.res'),
          path(metadata_json)

    output:
    tuple val(meta), path("${meta.id}_report.json"), emit: report
    tuple val(meta), path("${meta.id}_summary.csv"), emit: summary

    script:
    def trim_file = trim_log instanceof List ? trim_log[0] : trim_log
    """
    python3 ${projectDir}/bin/staphit-aggregate \
        --sample-id ${meta.id} \
        --trim-log ${trim_file} \
        --fastqc-dir . \
        --quast-dir ${quast_dir} \
        --mlst mlst_input.tsv \
        --abricate-dir . \
        --amrfinder amrfinder_input.tsv \
        --mash ${mash_sketch} \
        --spa ${spa_report} \
        --sccmec ${sccmec_report} \
        --agr ${agr_report} \
        --kma kma_input.res \
        --metadata ${metadata_json} \
        --outdir .
    """
}
