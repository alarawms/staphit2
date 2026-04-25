process SPECIES_REPORT {
    tag "species_qc"
    label 'process_low'
    publishDir "${params.outdir}/species_qc", mode: 'copy'
    container 'docker.io/python:3.9'

    input:
    path fastani_results
    path mash_results

    output:
    path "species_confirmed.tsv", emit: confirmed_tsv
    path "species_excluded.tsv" , emit: excluded_tsv

    script:
    def mash_arg = mash_results.name != 'NO_MASH_RESULTS' ? "--mash-screen ${mash_results}" : ''
    """
    python3 ${projectDir}/bin/staphit-species-report \
        --fastani ${fastani_results} \
        ${mash_arg} \
        --threshold ${params.species_ani_threshold} \
        --outdir .
    """
}
