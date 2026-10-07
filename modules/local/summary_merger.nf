process SUMMARY_MERGER {
    label 'process_low'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    input:
    path sample_summaries

    output:
    path "combined_summary.tsv", emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    """
    head -1 \$(ls *.csv | head -1) > combined_summary.tsv
    for f in *.csv; do
        tail -n +2 "\$f" >> combined_summary.tsv
    done
    """

    stub:
    """
    touch combined_summary.tsv
    """
}
