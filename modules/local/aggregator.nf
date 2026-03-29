process AGGREGATOR {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/aggregated", mode: 'copy'
    container 'docker.io/python:3.9'

    input:
    tuple val(meta),
          path(trim_log, stageAs: 'trim/*'),
          path(fastqc_files, stageAs: 'fastqc/*'),
          path(quast_dir, stageAs: 'quast/*'),
          path(mlst_tsv, stageAs: 'mlst/*'),
          path(abricate_tabs, stageAs: 'abricate/*'),
          path(amrfinder_report, stageAs: 'amrfinder/*'),
          path(mash_sketch, stageAs: 'mash/*'),
          path(spa_report, stageAs: 'spa/*'),
          path(sccmec_report, stageAs: 'sccmec/*'),
          path(agr_report, stageAs: 'agr/*'),
          path(kma_res, stageAs: 'kma/*'),
          path(metadata_json)

    output:
    tuple val(meta), path("${meta.id}_report.json"), emit: report
    tuple val(meta), path("${meta.id}_summary.csv"), emit: summary

    script:
    """
    python3 ${projectDir}/bin/staphit-aggregate \
        --sample-id ${meta.id} \
        --trim-log trim/${trim_log.name} \
        --fastqc-dir fastqc \
        --quast-dir quast \
        --mlst mlst/${mlst_tsv.name} \
        --abricate-dir abricate \
        --amrfinder amrfinder/${amrfinder_report.name} \
        --mash mash/${mash_sketch.name} \
        --spa spa/${spa_report.name} \
        --sccmec sccmec/${sccmec_report.name} \
        --agr agr/${agr_report.name} \
        --kma kma/${kma_res.name} \
        --metadata ${metadata_json} \
        --outdir .
    """
}
