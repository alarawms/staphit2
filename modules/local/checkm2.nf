process CHECKM2 {
    tag "$meta.id"
    label 'process_medium'
    container 'docker.io/staphb/checkm2:1.1.0'

    input:
    tuple val(meta), path(assembly)
    path db_file

    output:
    tuple val(meta), path("${meta.id}_quality_report.tsv"), emit: report
    tuple val("${task.process}"), val('checkm2'), eval("checkm2 --version"), topic: versions, emit: versions_checkm2

    script:
    """
    mkdir -p input_dir
    if [[ "${assembly}" == *.gz ]]; then
        zcat ${assembly} > input_dir/${meta.id}.fasta
    else
        cp ${assembly} input_dir/${meta.id}.fasta
    fi
    checkm2 predict --input input_dir --output-directory checkm2_out -x fasta --threads ${task.cpus} --force --remove_intermediates --database_path ${db_file}
    cp checkm2_out/quality_report.tsv ${meta.id}_quality_report.tsv
    """

    stub:
    """
    printf "Name\tCompleteness\tContamination\n${meta.id}\t100\t0\n" > ${meta.id}_quality_report.tsv
    """
}
