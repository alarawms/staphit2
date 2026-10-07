process VALIDATE_METADATA {
    label 'process_low'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    input:
    path metadata_csv
    path antibiogram_csv
    path samplesheet

    output:
    path "metadata.json", emit: json
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def abg_flag = antibiogram_csv.name != 'NO_ANTIBIOGRAM' ? "--antibiogram ${antibiogram_csv}" : ''
    """
    staphit-metadata validate \
        --metadata ${metadata_csv} \
        --samplesheet ${samplesheet} \
        ${abg_flag} || true
    staphit-metadata normalize \
        --metadata ${metadata_csv} \
        ${abg_flag} \
        -o metadata.json
    """

    stub:
    """
    echo '{}' > metadata.json
    """
}
