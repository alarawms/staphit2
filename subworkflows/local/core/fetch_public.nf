include { PUBMLST_FETCH } from '../../../modules/local/pubmlst_fetch'
include { SRA_FETCH     } from '../../../modules/local/sra_fetch'

workflow FETCH_PUBLIC {

    main:
    PUBMLST_FETCH()

    ch_accessions = PUBMLST_FETCH.out.accessions
        .splitText()
        .map { it.trim() }
        .filter { it }

    SRA_FETCH ( ch_accessions )

    ch_reads = SRA_FETCH.out.reads
        .map { acc, reads -> [ [id: acc, single_end: false], reads ] }

    emit:
    reads    = ch_reads
    metadata = PUBMLST_FETCH.out.metadata
    versions = Channel.empty()
}
