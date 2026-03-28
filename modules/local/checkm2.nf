process CHECKM2 {
    tag "$meta.id"
    label 'process_medium'
    container 'docker.io/staphb/checkm2:latest'

    input:
    tuple val(meta), path(assembly)
    path db_file

    output:
    tuple val(meta), path("${meta.id}_quality_report.tsv"), emit: report

    script:
    """
    mkdir -p input_dir
    # Handle both gzipped and plain FASTA
    if file ${assembly} | grep -q gzip; then
        zcat ${assembly} > input_dir/${meta.id}.fasta
    else
        cp ${assembly} input_dir/${meta.id}.fasta
    fi
    checkm2 predict --input input_dir --output-directory checkm2_out -x fasta --threads ${task.cpus} --force --remove_intermediates --database_path ${db_file}
    cp checkm2_out/quality_report.tsv ${meta.id}_quality_report.tsv
    """
}
