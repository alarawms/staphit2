process VISUALIZATION {
    label 'process_low'
    publishDir "${params.outdir}/figures", mode: 'copy'
    container 'docker.io/python:3.9'

    input:
    path summary_tsv
    path clusters

    output:
    path "plots/*.png", optional: true, emit: plots

    script:
    """
    pip install matplotlib seaborn pandas 2>&1 | tail -1
    mkdir -p plots
    python3 ${projectDir}/bin/staphit-visualize ${summary_tsv} plots || echo "Visualization completed with warnings"
    """
}
