process VISUALIZATION {
    label 'process_low'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    input:
    path summary_tsv
    path clusters

    output:
    path "plots/*.png", optional: true, emit: plots
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    """
    mkdir -p plots
    staphit-visualize ${summary_tsv} plots
    """

    stub:
    """
    mkdir plots
    """
}
