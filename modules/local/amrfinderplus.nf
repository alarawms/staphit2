process AMRFINDERPLUS {
    tag "$meta.id"
    label 'process_medium'
    container 'docker.io/staphb/ncbi-amrfinderplus:latest'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tsv"), emit: report

    script:
    """
    # Handle gzipped input
    if file ${assembly} | grep -q gzip; then
        zcat ${assembly} > input.fasta
    else
        cp ${assembly} input.fasta
    fi

    amrfinder -n input.fasta --organism Staphylococcus_aureus --threads ${task.cpus} > ${meta.id}_amrfinder.tsv
    """
}
