process SAMPLE_STATUS {
    label 'process_single'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    input:
    path stages             // stage<TAB>sample_id lines collected by the workflow
    path qc_report          // qc_report.tsv or NO_QC
    path species_excluded   // species_excluded.tsv or NO_SPECIES_QC
    path summary            // combined_summary.tsv
    path tree               // tree file or NO_TREE

    output:
    path "sample_status.tsv", emit: tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    """
    staphit-sample-status \\
        --stages ${stages} \\
        --qc ${qc_report} \\
        --species-excluded ${species_excluded} \\
        --summary ${summary} \\
        --tree ${tree} \\
        -o sample_status.tsv
    """

    stub:
    """
    printf "sample_id\\tstatus\\tstage\\treason\\n" > sample_status.tsv
    """
}
