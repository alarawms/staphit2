process SCCMEC {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/alarawms/sccmec_typer:v2'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}/sccmec/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}_sccmec.tsv")              , emit: report
    tuple val(meta), path("${meta.id}_sccmec.json")             , emit: json
    tuple val(meta), path("${meta.id}_sccmec_elements.csv")     , emit: elements
    tuple val(meta), path("${meta.id}_sccmec_map.svg")          , optional: true, emit: svg
    tuple val(meta), path("${meta.id}_sccmec_report.html")      , optional: true, emit: html

    script:
    def viz_flag      = params.sccmec_viz      ? ''                                                   : '--no-viz'
    def bestfit_flag  = params.sccmec_best_fit ? "--best-fit --min-estimate-score ${params.sccmec_min_score}" : ''
    """
    python3 /app/bin/sccmec_typer.py \\
        --1 ${assembly} \\
        -d /app/db/sccmec_targets.fasta \\
        -o ${meta.id}_sccmec \\
        --threads ${task.cpus} \\
        ${viz_flag} ${bestfit_flag}
    """
}
