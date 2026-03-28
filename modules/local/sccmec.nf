process SCCMEC {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/alarawms/sccmec_typer:latest'
    containerOptions '--entrypoint ""'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_sccmec.tsv"), emit: report

    script:
    """
    python3 /app/bin/sccmec_typer.py --1 ${assembly} -d /app/db -o ${meta.id}_sccmec --no-viz 2>/dev/null && \
        mv ${meta.id}_sccmec.tsv ${meta.id}_sccmec.tsv 2>/dev/null || \
        echo -e "Sample\\tStatus\\tmecA_Present\\tSCCmec_Type\\n${meta.id}\\tND\\tND\\tND" > ${meta.id}_sccmec.tsv
    """
}
