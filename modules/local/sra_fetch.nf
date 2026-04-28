process SRA_FETCH {
    tag "$accession"
    label 'process_medium'
    container 'quay.io/biocontainers/sra-tools:3.1.0--h4304569_1'
    errorStrategy 'ignore'

    input:
    val accession

    output:
    tuple val(accession), path("${accession}_{1,2}.fastq.gz"), emit: reads

    script:
    """
    fasterq-dump ${accession} \\
        --split-3 \\
        --threads ${task.cpus} \\
        --outdir . \\
        --temp .

    # Require both mates — skip single-end or failed accessions
    if [ ! -f "${accession}_1.fastq" ] || [ ! -f "${accession}_2.fastq" ]; then
        echo "WARNING: ${accession} is not paired-end or download failed, skipping" >&2
        exit 1
    fi

    gzip ${accession}_1.fastq ${accession}_2.fastq
    """
}
