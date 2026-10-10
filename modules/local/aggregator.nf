process AGGREGATOR {
    tag "$meta.id"
    label 'process_low'
    stageInMode 'copy'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

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
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def trim_file = trim_log instanceof List ? trim_log[0] : trim_log
    // Long-only samples have no fastp log; their KMA input was raw ONT reads
    def kma_ont   = trim_file.name == 'NO_TRIMLOG' ? '--kma-raw-ont' : ''
    """
    # v8: AMR consensus columns; v7: sccmec typing mode + cassette columns; v6: typer mec-locus proximity; v5: sccmec IWG columns; v4: harmonized year/host/source/sccmec_group columns; v3: fastp parser + metadata summary columns (bump to force re-aggregation after bin/staphit-aggregate changes)
    staphit-aggregate \
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
        --kma kma_input.res ${kma_ont} \
        --genome-size ${params.genome_size} --target-depth ${params.target_depth} \
        --metadata ${metadata_json} \
        --outdir .
    """

    stub:
    """
    echo '{}' > ${meta.id}_report.json
    printf "sample_id\n${meta.id}\n" > ${meta.id}_summary.csv
    """
}
