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

    emit:
    reads    = ch_reads
    metadata = PUBMLST_FETCH.out.metadata
    versions = channel.empty()
}
