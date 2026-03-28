process CHECKM2_DB {
    label 'process_low'
    publishDir "${params.outdir}/databases", mode: 'copy'
    container 'staphb/checkm2:latest'

    output:
    path "uniref100.KO.1.dmnd", emit: db

    script:
    """
    checkm2 database --download --path _tmp_db
    mv _tmp_db/CheckM2_database/uniref100.KO.1.dmnd .
    """
}
