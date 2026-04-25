/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPECIES QC SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Confirms assemblies are S. aureus using fastANI (ANI >= 95% vs NCTC 8325).
    Non-aureus rejects are identified via Mash screen against RefSeq.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTANI         } from '../../../modules/nf-core/fastani/main'
include { MASH_SCREEN     } from '../../../modules/nf-core/mash/screen/main'
include { MASH_REFSEQ_DB  } from '../../../modules/local/mash_refseq_db'
include { SPECIES_REPORT  } from '../../../modules/local/species_report'

workflow SPECIES_QC {

    take:
    ch_assemblies   // channel: [ val(meta), path(fasta) ]

    main:

    ch_versions = Channel.empty()

    // Reference for fastANI (wrapped in tuple for nf-core module)
    ch_reference = Channel.of([ [id: 'NCTC8325'], file("${projectDir}/assets/references/NCTC8325.fasta") ])

    //
    // MODULE: Run fastANI — each assembly vs NCTC 8325
    //
    FASTANI (
        ch_assemblies,
        ch_reference.collect(),
        [],   // no query list
        []    // no reference list
    )

    //
    // Parse ANI results and split into confirmed / rejected
    //
    ch_ani_parsed = FASTANI.out.ani
        .map { meta, ani_file ->
            def lines = ani_file.text.trim().split('\n')
            def ani_val = 0.0f
            if (lines.size() > 0 && lines[0].split('\t').size() >= 3) {
                ani_val = lines[0].split('\t')[2] as Float
            }
            [ meta, ani_val ]
        }

    ch_ani_branched = ch_ani_parsed.branch {
        meta, ani ->
            pass: ani >= params.species_ani_threshold
            fail: true
    }

    // Get confirmed sample IDs as a set for filtering
    ch_confirmed_ids = ch_ani_branched.pass
        .map { meta, ani -> meta.id }
        .collect()
        .map { it.toSet() }

    ch_confirmed_assemblies = ch_assemblies
        .combine(ch_confirmed_ids)
        .filter { meta, fasta, confirmed_set -> meta.id in confirmed_set }
        .map { meta, fasta, confirmed_set -> [ meta, fasta ] }

    ch_rejected_assemblies = ch_assemblies
        .combine(ch_confirmed_ids)
        .filter { meta, fasta, confirmed_set -> !(meta.id in confirmed_set) }
        .map { meta, fasta, confirmed_set -> [ meta, fasta ] }

    //
    // MODULE: Download Mash RefSeq sketch (once, cached)
    //
    MASH_REFSEQ_DB ()

    //
    // MODULE: Mash screen on rejected samples to identify species
    //
    ch_mash_db = MASH_REFSEQ_DB.out.db.map { db -> [ [id: 'refseq'], db ] }

    MASH_SCREEN (
        ch_rejected_assemblies,
        ch_mash_db.collect()
    )
    ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())

    //
    // MODULE: Collate species report
    //
    ch_fastani_collected = FASTANI.out.ani
        .map { meta, f -> f }
        .collectFile(name: 'fastani_all.tsv')

    ch_mash_collected = MASH_SCREEN.out.screen
        .map { meta, screen -> screen }
        .collectFile(name: 'mash_screen_all.tsv')
        .ifEmpty(file('NO_MASH_RESULTS'))

    SPECIES_REPORT (
        ch_fastani_collected,
        ch_mash_collected
    )

    emit:
    confirmed_assemblies = ch_confirmed_assemblies  // [ val(meta), path(fasta) ]
    species_report       = SPECIES_REPORT.out.confirmed_tsv
    species_excluded     = SPECIES_REPORT.out.excluded_tsv
    versions             = ch_versions
}
