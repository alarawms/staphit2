/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    S. AUREUS TYPING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Species-specific typing: MLST, spa, SCCmec, agr, Mash
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MLST       } from '../../../../modules/nf-core/mlst/main'
include { SPATYPER   } from '../../../../modules/local/spatyper'
include { SCCMEC     } from '../../../../modules/local/sccmec'
include { AGR_TYPING } from '../../../../modules/local/agr_typing'
include { MASH       } from '../../../../modules/local/mash'

workflow SA_TYPING {

    take:
    ch_assemblies  // channel: [ val(meta), path(assembly) ]

    main:

    ch_versions = Channel.empty()

    MLST ( ch_assemblies )
    ch_versions = ch_versions.mix(MLST.out.versions.first())

    SPATYPER ( ch_assemblies )
    SCCMEC ( ch_assemblies )
    AGR_TYPING ( ch_assemblies )
    MASH ( ch_assemblies )

    emit:
    mlst     = MLST.out.tsv            // [ val(meta), path(tsv) ]
    spa      = SPATYPER.out.report     // [ val(meta), path(tsv) ]
    sccmec   = SCCMEC.out.report      // [ val(meta), path(tsv) ]
    agr      = AGR_TYPING.out.report   // [ val(meta), path(json) ]
    mash     = MASH.out.sketch         // [ val(meta), path(msh) ]
    versions = ch_versions             // channel: [ path(versions.yml) ]
}
