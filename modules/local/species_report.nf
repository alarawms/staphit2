process SPECIES_REPORT {
    tag "species_qc"
    label 'process_low'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    input:
    path fastani_results
    path mash_results

    output:
    path "species_confirmed.tsv", emit: confirmed_tsv
    path "species_excluded.tsv" , emit: excluded_tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def mash_arg = mash_results.name != 'NO_MASH_RESULTS' ? "--mash-screen ${mash_results}" : ''
    """
    staphit-species-report \
        --fastani ${fastani_results} \
        ${mash_arg} \
        --threshold ${params.species_ani_threshold} \
        --outdir .
    """

    stub:
    """
    touch species_confirmed.tsv species_excluded.tsv
    """
}
