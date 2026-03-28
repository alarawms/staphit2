process SPATYPER {
    tag "$meta.id"
    label 'process_low'
    errorStrategy { task.exitStatus in [125,137] ? 'retry' : 'terminate' }
    maxRetries 2
    container 'docker.io/staphb/spatyper:latest'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_spatyper.tsv"), emit: report

    script:
    """
    spaTyper -f ${assembly} --output ${meta.id}_spatyper.tsv
    """
}
