process QC_GATE {
    label 'process_low'
    publishDir "${params.outdir}/qc_gate", mode: 'copy'
    container 'docker.io/python:3.9'

    input:
    path "reports/*"

    output:
    path "qc_report.tsv", emit: report
    path "passed_samples.txt", emit: passed
    path "failed_samples.txt", emit: failed

    script:
    """
    python3 ${projectDir}/bin/staphit-qc --checkm2-dir reports --min-completeness ${params.min_completeness} --max-contamination ${params.max_contamination} -o .
    """
}
