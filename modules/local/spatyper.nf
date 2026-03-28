process SPATYPER {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/biocontainers/spatyper:0.3.3--pyhdfd78af_3'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_spatyper.tsv"), emit: report

    script:
    """
    spaTyper -f ${assembly} --output ${meta.id}_spatyper.tsv
    """
}
