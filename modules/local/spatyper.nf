process SPATYPER {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/biocontainers/spatyper:0.3.3--pyhdfd78af_3'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_spatyper.tsv"), emit: report
    tuple val("${task.process}"), val('spatyper'), eval("spaTyper --version 2>&1 | sed 's/^spaTyper //'"), topic: versions, emit: versions_spatyper

    script:
    """
    spaTyper -f ${assembly} --output ${meta.id}_spatyper.tsv 2>/dev/null || \
        echo -e "File\\tRepeats\\tType\\n${meta.id}\\tND\\tND" > ${meta.id}_spatyper.tsv
    """

    stub:
    """
    touch ${meta.id}_spatyper.tsv
    """
}
