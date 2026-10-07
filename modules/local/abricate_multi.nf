process ABRICATE_MULTI {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/staphb/abricate:1.2.0'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tab"), emit: reports
    tuple val("${task.process}"), val('abricate'), eval("abricate --version | sed 's/^abricate //'"), topic: versions, emit: versions_abricate

    script:
    """
    abricate --db resfinder ${assembly} > ${meta.id}_resfinder.tab
    abricate --db vfdb ${assembly} > ${meta.id}_vfdb.tab
    abricate --db plasmidfinder ${assembly} > ${meta.id}_plasmidfinder.tab
    """

    stub:
    """
    touch ${meta.id}_resfinder.tab ${meta.id}_vfdb.tab ${meta.id}_plasmidfinder.tab
    """
}
