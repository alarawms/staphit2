/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AMR DETECTION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Assembly-based: AMRFinderPlus + ABRicate (resfinder, vfdb, plasmidfinder)
    Read-based: KMA against ResFinder DB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AMRFINDERPLUS_RUN           } from '../../../modules/nf-core/amrfinderplus/run/main'
include { ABRICATE_MULTI              } from '../../../modules/local/abricate_multi'
include { FETCH_RESFINDER_DB; INDEX_DB; KMA } from '../../../modules/local/kma'

workflow AMR_DETECTION {

    take:
    ch_assemblies     // channel: [ val(meta), path(assembly) ]
    ch_trimmed_reads  // channel: [ val(meta), path(reads) ]

    main:

    //
    // Assembly-based AMR detection with AMRFinderPlus
    //
    AMRFINDERPLUS_RUN ( ch_assemblies, [] )

    //
    // Assembly-based screening with ABRicate (resfinder + vfdb + plasmidfinder)
    //
    ABRICATE_MULTI ( ch_assemblies )

    //
    // Read-based AMR detection with KMA against ResFinder
    //
    FETCH_RESFINDER_DB ()
    INDEX_DB ( FETCH_RESFINDER_DB.out.db )
    KMA ( ch_trimmed_reads, INDEX_DB.out.indexed_db.collect() )

    emit:
    amrfinder = AMRFINDERPLUS_RUN.out.report  // [ val(meta), path(tsv) ]
    abricate  = ABRICATE_MULTI.out.reports    // [ val(meta), path(tabs) ]
    kma       = KMA.out.results               // [ val(meta), path(res) ]
}
