process SCCMEC {
    tag "$meta.id"
    label 'process_low'
    container "${params.sccmec_container}"
    containerOptions { workflow.containerEngine in ['docker', 'podman'] ? '--entrypoint ""' : '' }   // Docker-only flag; Apptainer rejects it

    input:
    tuple val(meta), path(assembly), path(reads, stageAs: 'reads/*')

    output:
    tuple val(meta), path("${meta.id}_sccmec.tsv")              , emit: report
    tuple val(meta), path("${meta.id}_sccmec.json")             , emit: json
    tuple val(meta), path("${meta.id}_sccmec_elements.csv")     , emit: elements
    tuple val(meta), path("${meta.id}_sccmec_map.svg")          , optional: true, emit: svg
    tuple val(meta), path("${meta.id}_sccmec_report.html")      , optional: true, emit: html
    tuple val("${task.process}"), val('sccmec_typer'), val(task.container), topic: versions, emit: versions_sccmec_typer

    script:
    def viz_flag      = params.sccmec_viz      ? ''                                                   : '--no-viz'
    def bestfit_flag  = params.sccmec_best_fit ? "--best-fit --min-estimate-score ${params.sccmec_min_score}" : ''
    // read fallback when the assembly splits the cassette (assembly_limited)
    def rl            = reads instanceof List ? reads : [ reads ]
    def fallback_flag = rl.size() >= 2 ? "--fallback-1 ${rl[0]} --fallback-2 ${rl[1]}"
                      : rl.size() == 1 ? "--fallback-1 ${rl[0]}" : ''
    """
    python3 /app/bin/sccmec_typer.py \\
        --1 ${assembly} \\
        -d /app/db/sccmec_targets.fasta \\
        -o ${meta.id}_sccmec \\
        --threads ${task.cpus} \\
        ${viz_flag} ${bestfit_flag} ${fallback_flag}
    """

    stub:
    """
    touch ${meta.id}_sccmec.tsv ${meta.id}_sccmec.json ${meta.id}_sccmec_elements.csv
    """
}
