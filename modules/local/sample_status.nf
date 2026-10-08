process SAMPLE_STATUS {
    label 'process_single'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    input:
    path stages             // stage<TAB>sample_id lines collected by the workflow
    path qc_report          // qc_report.tsv or NO_QC
    path species_excluded   // species_excluded.tsv or NO_SPECIES_QC
    path summary, stageAs: 'summary_in.tsv'   // combined_summary.tsv from SUMMARY_MERGER
    path tree               // tree file or NO_TREE

    output:
    path "sample_status.tsv",    emit: tsv
    path "combined_summary.tsv", emit: summary   // the summary plus a qc_flags column
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    """
    staphit-sample-status \\
        --stages ${stages} \\
        --qc ${qc_report} \\
        --species-excluded ${species_excluded} \\
        --summary ${summary} \\
        --tree ${tree} \\
        --min-size ${params.min_assembly_size} --max-size ${params.max_assembly_size} \\
        --min-depth-short ${params.min_depth_short} --min-depth-long ${params.min_depth_long} \\
        --max-contigs ${params.max_contigs} --min-n50 ${params.min_n50} \\
        --summary-out combined_summary.tsv \\
        -o sample_status.tsv
    """

    stub:
    """
    printf "sample_id\\tstatus\\tstage\\treason\\tflags\\n" > sample_status.tsv
    cp ${summary} combined_summary.tsv
    """
}
