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
    ch_reads       // channel: [ val(meta), path(reads) ] — SCCmec read fallback (may be empty)

    main:

    ch_versions = channel.empty()

    MLST ( ch_assemblies )
    ch_versions = ch_versions.mix(MLST.out.versions.first())

    SPATYPER ( ch_assemblies )
    // SCCmec: assembly + the sample's reads (re-typed from reads when the assembly
    // splits the cassette); samples without reads get an empty list
    ch_sccmec_in = ch_assemblies.map { meta, fasta -> [ meta.id, meta, fasta ] }
        .join(ch_reads.map { meta, reads -> [ meta.id, reads ] }, remainder: true)
        // remainder gives [id, meta, fasta, null] (no reads) or [id, null, reads] (no assembly)
        .filter { row -> row[1] != null }
        .map { row -> [ row[1], row[2], (row.size() > 3 && row[3]) ? row[3] : [] ] }
    SCCMEC ( ch_sccmec_in )
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
