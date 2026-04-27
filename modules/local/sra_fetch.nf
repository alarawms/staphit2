process SRA_FETCH {
    tag "$accession"
    label 'process_medium'
    container 'quay.io/biocontainers/sra-tools:3.1.0--h4304569_1'

    input:
    val accession

    output:
    tuple val(accession), path("${accession}_{1,2}.fastq.gz"), emit: reads

    script:
    """
    fasterq-dump ${accession} \\
        --split-files \\
        --threads ${task.cpus} \\
        --outdir . \\
        --temp .
    gzip ${accession}_1.fastq
    gzip ${accession}_2.fastq
    """
}
