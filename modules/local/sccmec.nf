process SCCMEC {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/alarawms/sccmec_typer:latest'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}/sccmec/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_sccmec.tsv"), emit: report
    tuple val(meta), path("*_sccmec.svg"), optional: true, emit: svg
    tuple val(meta), path("*_sccmec.html"), optional: true, emit: html

    script:
    def viz_flag = params.sccmec_viz ? '' : '--no-viz'
    """
    python3 /app/bin/sccmec_typer.py --1 ${assembly} -d /app/db -o ${meta.id}_sccmec ${viz_flag} 2>/dev/null || \
        echo -e "Sample\\tStatus\\tmecA_Present\\tSCCmec_Type\\n${meta.id}\\tND\\tND\\tND" > ${meta.id}_sccmec.tsv
    """
}
