process QC_GATE {
    label 'process_low'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    input:
    path "reports/*"

    output:
    path "qc_report.tsv", emit: report
    path "passed_samples.txt", emit: passed
    path "failed_samples.txt", emit: failed
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    """
    staphit-qc --checkm2-dir reports --min-completeness ${params.min_completeness} --max-contamination ${params.max_contamination} -o .
    """

    stub:
    """
    touch qc_report.tsv failed_samples.txt
    ls reports | sed 's/_quality_report.tsv//' > passed_samples.txt
    """
}
