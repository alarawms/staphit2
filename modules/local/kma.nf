process FETCH_RESFINDER_DB {
    label 'process_low'
    container 'staphb/kma:1.4.14'
    storeDir "${params.outdir}/databases/resfinder"

    output:
    path "resfinder_db", emit: db

    script:
    """
    git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git
    """
}

process INDEX_DB {
    label 'process_low'
    container 'staphb/kma:1.4.14'

    input:
    path db

    output:
    path "indexed_db", emit: indexed_db

    script:
    """
    mkdir indexed_db
    for f in ${db}/*.fsa; do
        kma index -i "\$f" -o indexed_db/\$(basename "\$f" .fsa) 2>/dev/null || true
    done
    """
}

process KMA {
    tag "$meta.id"
    label 'process_low'
    container 'staphb/kma:1.4.14'

    input:
    tuple val(meta), path(reads)
    path indexed_db

    output:
    tuple val(meta), path("*.res"), emit: results

    script:
    """
    kma -i ${reads[0]} ${reads[1]} -o ${meta.id} -t_db indexed_db/\$(ls indexed_db/*.name | head -1 | sed 's/.name//') -1t1 || true
    touch ${meta.id}.res
    """
}
