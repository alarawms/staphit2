process ITOL_EXPORT {
    tag "itol"
    label 'process_low'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    input:
    path treefile
    path summary
    path clusters

    output:
    path "itol/*", emit: files
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def cl = clusters.name != 'NO_CLUSTERS' ? clusters : 'NONE'
    """
    # v2 — bump after editing bin/itol_annotations.py or bin/harmonize_metadata.py so -resume re-exports
    itol_annotations.py ${treefile} ${summary} ${cl} itol
    """

    stub:
    """
    mkdir itol
    touch itol/stub.txt
    """
}
