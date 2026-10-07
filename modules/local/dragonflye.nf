process DRAGONFLYE {
    tag "$meta.id"
    label 'process_high'
    container 'docker.io/staphb/dragonflye:1.2.1'

    input:
    tuple val(meta), path(long_reads), path(short_reads)   // short_reads = [] for long-only

    output:
    tuple val(meta), path("${meta.id}.scaffolds.fasta"), emit: scaffolds
    tuple val(meta), path("${meta.id}.dragonflye.log"),  emit: log
    tuple val(meta), path("${meta.id}.gfa"),             emit: gfa
    tuple val("${task.process}"), val('dragonflye'), eval("dragonflye --version | sed 's/^dragonflye //'"), topic: versions, emit: versions_dragonflye
    tuple val("${task.process}"), val('flye'), eval("flye --version"), topic: versions, emit: versions_flye

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

    stub:
    """
    { echo ">contig_1"; head -c 600000 /dev/zero | tr "\\0" A; echo; } > ${meta.id}.scaffolds.fasta   # passes the 500 kb assembly-size filter
    touch ${meta.id}.dragonflye.log ${meta.id}.gfa
    """
}
