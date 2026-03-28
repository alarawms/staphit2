process SKESA {
    tag "$meta.id"
    label 'process_high'
    container 'docker.io/staphb/skesa:latest'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.scaffolds.fasta"), emit: scaffolds

    script:
    """
    skesa --reads ${reads[0]} ${reads[1]} --cores ${task.cpus} --memory ${task.memory.toGiga()} > ${meta.id}.scaffolds.fasta
    """
}
