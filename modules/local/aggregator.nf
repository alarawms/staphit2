process AGGREGATOR {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/aggregated", mode: 'copy'
    container 'docker.io/python:3.9'

    input:
    tuple val(meta),
          path(trim_log),
          path(fastqc_files),
          path(quast_dir),
          path(mlst_tsv, stageAs: 'mlst_*.tsv'),
          path(abricate_tabs),
          path(amrfinder_report, stageAs: 'amrfinder_*.tsv'),
          path(mash_sketch),
          path(spa_report),
          path(sccmec_report),
          path(agr_report),
          path(kma_res, stageAs: 'kma_*.res'),
          path(metadata_json)

    output:
    tuple val(meta), path("${meta.id}_report.json"), emit: report
    tuple val(meta), path("${meta.id}_summary.csv"), emit: summary

    script:
    def trim_file = trim_log instanceof List ? trim_log[0] : trim_log
    def mlst_file = mlst_tsv instanceof List ? mlst_tsv[0] : mlst_tsv
    def amr_file = amrfinder_report instanceof List ? amrfinder_report[0] : amrfinder_report
    def kma_file = kma_res instanceof List ? kma_res[0] : kma_res
    """
    python3 ${projectDir}/bin/staphit-aggregate \
        --sample-id ${meta.id} \
        --trim-log ${trim_file} \
        --fastqc-dir . \
        --quast-dir ${quast_dir} \
        --mlst ${mlst_file} \
        --abricate-dir . \
        --amrfinder ${amr_file} \
        --mash ${mash_sketch} \
        --spa ${spa_report} \
        --sccmec ${sccmec_report} \
        --agr ${agr_report} \
        --kma ${kma_file} \
        --metadata ${metadata_json} \
        --outdir .
    """
}
