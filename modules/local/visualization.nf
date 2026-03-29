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
    pip install matplotlib seaborn pandas > /dev/null 2>&1 || true
    mkdir -p plots
    python3 ${projectDir}/bin/staphit-visualize ${summary_tsv} plots 2>&1 || true
    ls -la plots/ || true
    """
}
