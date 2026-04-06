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
          path(mlst_tsv),
          path(abricate_tabs),
          path(amrfinder_report),
          path(mash_sketch),
          path(spa_report),
          path(sccmec_report),
          path(agr_report),
          path(kma_res),
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
    # Stage files with unique names to avoid collisions
    mkdir -p _inputs
    cp -L ${trim_file} _inputs/trim.log 2>/dev/null || true
    cp -L ${mlst_file} _inputs/mlst.tsv 2>/dev/null || true
    cp -L ${amr_file} _inputs/amrfinder.tsv 2>/dev/null || true
    cp -L ${kma_file} _inputs/kma.res 2>/dev/null || true

    python3 ${projectDir}/bin/staphit-aggregate \
        --sample-id ${meta.id} \
        --trim-log _inputs/trim.log \
        --fastqc-dir . \
        --quast-dir ${quast_dir} \
        --mlst _inputs/mlst.tsv \
        --abricate-dir . \
        --amrfinder _inputs/amrfinder.tsv \
        --mash ${mash_sketch} \
        --spa ${spa_report} \
        --sccmec ${sccmec_report} \
        --agr ${agr_report} \
        --kma _inputs/kma.res \
        --metadata ${metadata_json} \
        --outdir .
    """
}
