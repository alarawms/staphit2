process REPORT {
    label 'process_low'
    publishDir "${params.outdir}/report", mode: 'copy'
    container 'python:3.9-slim'

    input:
    path summary_tsv
    path clusters
    path trees

    output:
    path "staphit_report.html", emit: html
    path "staphit_report.json", emit: json

    script:
    def cluster_flag = clusters.name != 'NO_CLUSTERS' ? "--clusters ${clusters}" : ''
    def tree_flag = trees.name != 'NO_TREES' ? "--trees ${trees}" : ''
    """
    python3 ${projectDir}/bin/staphit-report --summary ${summary_tsv} ${cluster_flag} ${tree_flag} -o .
    """
}
