process CHECKM2_DB {
    label 'process_low'
    container 'docker.io/staphb/checkm2:1.1.0'

    output:
    path "uniref100.KO.1.dmnd", emit: db
    tuple val("${task.process}"), val('checkm2'), eval("checkm2 --version"), topic: versions, emit: versions_checkm2
    tuple val("${task.process}"), val('checkm2_db_md5'), eval("md5sum uniref100.KO.1.dmnd | cut -d' ' -f1"), topic: versions, emit: versions_checkm2_db_md5

    script:
    """
    checkm2 database --download --path _tmp_db --no_write_json_db
    mv _tmp_db/CheckM2_database/uniref100.KO.1.dmnd .
    """

    stub:
    """
    touch uniref100.KO.1.dmnd
    """
}
