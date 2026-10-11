include { PUBMLST_FETCH } from '../../../modules/local/pubmlst_fetch'
include { SRA_FETCH     } from '../../../modules/local/sra_fetch'

workflow FETCH_PUBLIC {

    main:
    PUBMLST_FETCH(channel.value(file("${projectDir}/bin/pubmlst_fetch.py")))

    ch_accessions = PUBMLST_FETCH.out.accessions
        .splitText()
        .map { s -> s.trim() }
        .filter { s -> s }

    SRA_FETCH ( ch_accessions )

    ch_reads = SRA_FETCH.out.reads
        .map { acc, reads -> [ [id: acc, single_end: false, mode: 'short'], reads, [] ] }

    // Samplesheet of what was downloaded, pointing at the reads published to <outdir>/sra/:
    // rerun or extend the analysis with --input <outdir>/pubmlst/samplesheet.csv (no re-download)
    def sra_dir = file("${params.outdir}/sra").toUriString()
    SRA_FETCH.out.reads
        .map { acc, _reads -> "${acc},${sra_dir}/${acc}_1.fastq.gz,${sra_dir}/${acc}_2.fastq.gz," }
        .collectFile(
            name: 'samplesheet.csv',
            storeDir: "${params.outdir}/pubmlst",
            seed: 'sample,fastq_1,fastq_2,long_fastq',
            newLine: true,
            sort: true
        )

    emit:
    reads    = ch_reads
    metadata = PUBMLST_FETCH.out.metadata
    versions = channel.empty()
}
