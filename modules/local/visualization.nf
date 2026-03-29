process VISUALIZATION {
    label 'process_low'
    publishDir "${params.outdir}/figures", mode: 'copy'
    errorStrategy 'ignore'
    container 'docker.io/python:3.9'

    input:
    path summary_tsv
    path clusters
    path tree_file

    output:
    path "plots/*.png", optional: true, emit: plots

    script:
    def tree_flag = tree_file.name != 'NO_TREE' ? "--tree ${tree_file}" : ''
    def cluster_flag = clusters.name != 'NO_CLUSTERS' ? "--clusters ${clusters}" : ''
    """
    pip install matplotlib seaborn pandas 2>&1 | tail -1
    mkdir -p plots
    python3 ${projectDir}/bin/staphit-visualize ${summary_tsv} plots ${tree_flag} ${cluster_flag} || echo "Visualization completed with warnings"
    """
}
