process ITOL_EXPORT {
    tag "itol"
    label 'process_low'
    container 'docker.io/python:3.9'
    publishDir "${params.outdir}", mode: 'copy'   // output path already starts with itol/

    input:
    path treefile
    path summary
    path clusters

    output:
    path "itol/*", emit: files

    script:
    def cl = clusters.name != 'NO_CLUSTERS' ? clusters : 'NONE'
    """
    # v2 — bump after editing bin/itol_annotations.py or bin/harmonize_metadata.py so -resume re-exports
    python3 ${projectDir}/bin/itol_annotations.py ${treefile} ${summary} ${cl} itol
    """
}
