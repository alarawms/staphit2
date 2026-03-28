process ABRICATE_MULTI {
    tag "$meta.id"
    label 'process_low'
    container 'staphb/abricate:latest'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tab"), emit: reports

    script:
    """
    abricate --db resfinder ${assembly} > ${meta.id}_resfinder.tab
    abricate --db vfdb ${assembly} > ${meta.id}_vfdb.tab
    abricate --db plasmidfinder ${assembly} > ${meta.id}_plasmidfinder.tab
    """
}
