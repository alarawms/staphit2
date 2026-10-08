process MASH {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/staphb/mash:2.3'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.msh"), emit: sketch
    tuple val("${task.process}"), val('mash'), eval("mash --version"), topic: versions, emit: versions_mash

    script:
    """
    mash sketch -o ${meta.id} ${assembly}
    """

    stub:
    """
    touch ${meta.id}.msh
    """
}
