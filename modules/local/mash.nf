process MASH {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/staphb/mash:latest'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.msh"), emit: sketch

    script:
    """
    mash sketch -o ${meta.id} ${assembly}
    """
}
