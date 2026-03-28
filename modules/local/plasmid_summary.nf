process PLASMID_SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}/plasmids", mode: 'copy'
    container 'python:3.9-slim'

    input:
    path mob_reports

    output:
    path "plasmid_summary.tsv", emit: summary
    path "plasmid_report.json", optional: true, emit: report

    script:
    """
    python3 ${projectDir}/bin/staphit-plasmids --mob-dir . -o .
    """
}
