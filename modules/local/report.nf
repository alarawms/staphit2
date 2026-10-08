process REPORT {
    label 'process_low'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    input:
    path summary_tsv
    path clusters
    path qc_report
    path plasmid_summary
    path sample_status

    output:
    path "run_report.md", emit: report
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def cluster_flag = clusters.name != 'NO_CLUSTERS' ? "--clusters ${clusters}" : ''
    def qc_flag = qc_report.name != 'NO_QC' ? "--qc ${qc_report}" : ''
    def plasmid_flag = plasmid_summary.name != 'NO_PLASMIDS' ? "--plasmids ${plasmid_summary}" : ''
    """
    staphit-report --summary ${summary_tsv} ${cluster_flag} ${qc_flag} ${plasmid_flag} --status ${sample_status} -o run_report.md
    """

    stub:
    """
    touch run_report.md
    """
}
