process SCCMEC {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/staphb/staphopia-sccmec:latest'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_sccmec.tsv"), emit: report

    script:
    """
    staphopia-sccmec --input ${assembly} --output ${meta.id}_sccmec.tsv 2>/dev/null || echo -e "sample_id\\tSCCmec_Type\\n${meta.id}\\tND" > ${meta.id}_sccmec.tsv
    """
}
