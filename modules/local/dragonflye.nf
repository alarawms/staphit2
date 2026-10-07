process DRAGONFLYE {
    tag "$meta.id"
    label 'process_high'
    publishDir path: { "${params.outdir}/dragonflye/${meta.id}" }, mode: 'copy', pattern: '*.{log,gfa}'
    container 'docker.io/staphb/dragonflye:1.2.1'

    input:
    tuple val(meta), path(long_reads), path(short_reads)   // short_reads = [] for long-only

    output:
    tuple val(meta), path("${meta.id}.scaffolds.fasta"), emit: scaffolds
    tuple val(meta), path("${meta.id}.dragonflye.log"),  emit: log
    tuple val(meta), path("${meta.id}.gfa"),             emit: gfa

    script:
    // Long-only: Flye + Racon (+ Medaka if --medaka_model). Hybrid: additionally Polypolish with the short reads.
    def polish = meta.mode == 'hybrid' ? "--R1 ${short_reads[0]} --R2 ${short_reads[1]} --polypolish 1" : ''
    def medaka = params.medaka_model ? "--model ${params.medaka_model}" : ''
    def nanohq = params.flye_nanohq ? '--nanohq' : ''
    """
    dragonflye \\
        --reads ${long_reads} \\
        --gsize ${params.genome_size} \\
        --cpus ${task.cpus} \\
        --ram ${task.memory.toGiga()} \\
        ${nanohq} ${medaka} ${polish} \\
        --outdir out
    mv out/contigs.fa ${meta.id}.scaffolds.fasta
    mv out/dragonflye.log ${meta.id}.dragonflye.log
    mv out/flye*.gfa ${meta.id}.gfa   # flye-unpolished.gfa when polished, flye.gfa otherwise
    """
}
