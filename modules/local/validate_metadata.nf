process VALIDATE_METADATA {
    label 'process_low'
    publishDir "${params.outdir}/metadata", mode: 'copy'
    container 'docker.io/python:3.9'

    input:
    path metadata_csv
    path antibiogram_csv
    path samplesheet

    output:
    path "metadata.json", emit: json

    script:
    def abg_flag = antibiogram_csv.name != 'NO_ANTIBIOGRAM' ? "--antibiogram ${antibiogram_csv}" : ''
    """
    python3 ${projectDir}/bin/staphit-metadata validate \
        --metadata ${metadata_csv} \
        --samplesheet ${samplesheet} \
        ${abg_flag} || true
    python3 ${projectDir}/bin/staphit-metadata normalize \
        --metadata ${metadata_csv} \
        ${abg_flag} \
        -o metadata.json
    """
}
