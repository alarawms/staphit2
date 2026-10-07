process NANOPLOT {
    tag "$meta.id"
    label 'process_low'
    publishDir path: { "${params.outdir}/nanoplot/${meta.id}" }, mode: 'copy'
    container 'quay.io/biocontainers/nanoplot:1.44.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(long_reads)

    output:
    tuple val(meta), path("${meta.id}_NanoStats.txt"), emit: stats
    tuple val(meta), path("*.html"),                   emit: html

    script:
    """
    NanoPlot --fastq ${long_reads} --threads ${task.cpus} --prefix ${meta.id}_ --outdir . --no_static
    """
}
