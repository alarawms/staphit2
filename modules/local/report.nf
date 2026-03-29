process REPORT {
    label 'process_low'
    publishDir "${params.outdir}/report", mode: 'copy'
    container 'docker.io/python:3.9'

    input:
    path summary_tsv
    path clusters
    path qc_report
    path plasmid_summary

    output:
    path "run_report.md", emit: report

    script:
    def cluster_flag = clusters.name != 'NO_CLUSTERS' ? "--clusters ${clusters}" : ''
    def qc_flag = qc_report.name != 'NO_QC' ? "--qc ${qc_report}" : ''
    def plasmid_flag = plasmid_summary.name != 'NO_PLASMIDS' ? "--plasmids ${plasmid_summary}" : ''
    """
    python3 ${projectDir}/bin/staphit-report --summary ${summary_tsv} ${cluster_flag} ${qc_flag} ${plasmid_flag} -o run_report.md
    """
}
