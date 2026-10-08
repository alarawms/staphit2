process NANOPLOT {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/biocontainers/nanoplot:1.44.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(long_reads)

    output:
    tuple val(meta), path("${meta.id}_NanoStats.txt"), emit: stats
    tuple val(meta), path("*.html"),                   emit: html
    tuple val("${task.process}"), val('nanoplot'), eval("NanoPlot --version | sed 's/^NanoPlot //'"), topic: versions, emit: versions_nanoplot

    script:
    """
    NanoPlot --fastq ${long_reads} --threads ${task.cpus} --prefix ${meta.id}_ --outdir . --no_static
    """

    stub:
    """
    touch ${meta.id}_NanoStats.txt ${meta.id}_NanoPlot-report.html
    """
}
