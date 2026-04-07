process SNIPPY {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/snippy/${meta.id}", mode: 'link'
    errorStrategy 'ignore'
    container 'docker.io/staphb/snippy:4.6.0'

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("${meta.id}"), emit: results

    script:
    """
    snippy --cpus ${task.cpus} --ram ${task.memory.toGiga()} --outdir ${meta.id} --ref ${reference} --R1 ${reads[0]} --R2 ${reads[1]} --mincov ${params.snippy_mincov} --minqual ${params.snippy_minqual} --cleanup
    """
}

process SNIPPY_CORE {
    label 'process_medium'
    publishDir "${params.outdir}/snippy_core", mode: 'copy'
    container 'docker.io/staphb/snippy:4.6.0'

    input:
    path snippy_dirs
    path reference

    output:
    path "core.aln", emit: aln

    script:
    """
    snippy-core --ref ${reference} ${snippy_dirs}
    """
}
