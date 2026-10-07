process SKESA {
    tag "$meta.id"
    label 'process_high'
    container 'docker.io/staphb/skesa@sha256:5bffb46096deb74f4a7998109c07f2456c2160eb852e1238b71f1aad30e79843'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.scaffolds.fasta"), emit: scaffolds
    tuple val("${task.process}"), val('skesa'), eval("skesa --version 2>&1 | sed -n 's/^SKESA //p'"), topic: versions, emit: versions_skesa

    script:
    """
    skesa --reads ${reads[0]} ${reads[1]} --cores ${task.cpus} --memory ${task.memory.toGiga()} > ${meta.id}.scaffolds.fasta
    """

    stub:
    """
    { echo ">contig_1"; head -c 600000 /dev/zero | tr "\\0" A; echo; } > ${meta.id}.scaffolds.fasta   # passes the 500 kb assembly-size filter
    """
}
