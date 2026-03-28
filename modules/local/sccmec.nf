process SCCMEC {
    tag "$meta.id"
    label 'process_low'
    errorStrategy 'ignore'
    container 'docker.io/alarawms/sccmec_typer:latest'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_sccmec.tsv"), emit: report

    script:
    """
    sccmec_typer.py --1 ${assembly} -d /opt/conda/share/sccmec_db -o ${meta.id}_sccmec.tsv 2>/dev/null || echo -e "sample_id\\tSCCmec_Type\\n${meta.id}\\tND" > ${meta.id}_sccmec.tsv
    """
}
