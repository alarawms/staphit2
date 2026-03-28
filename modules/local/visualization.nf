process VISUALIZATION {
    label 'process_low'
    publishDir "${params.outdir}/figures", mode: 'copy'
    container 'docker.io/python:3.9'

    input:
    path summary_tsv
    path clusters
    path trees

    output:
    path "*.png", optional: true, emit: plots
    path "*.svg", optional: true, emit: svg_plots

    script:
    def cluster_flag = clusters.name != 'NO_CLUSTERS' ? "--clusters ${clusters}" : ''
    def tree_flag = trees.name != 'NO_TREES' ? "--trees ${trees}" : ''
    """
    python3 ${projectDir}/bin/staphit-visualize --summary ${summary_tsv} ${cluster_flag} ${tree_flag} -o .
    """
}
