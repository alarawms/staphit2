process SUMMARY_MERGER {
    label 'process_low'
    publishDir "${params.outdir}/summary", mode: 'copy'
    container 'python:3.9-slim'

    input:
    path sample_summaries

    output:
    path "combined_summary.tsv", emit: summary

    script:
    """
    head -1 \$(ls *.csv | head -1) > combined_summary.tsv
    for f in *.csv; do
        tail -n +2 "\$f" >> combined_summary.tsv
    done
    """
}
