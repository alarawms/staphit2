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

    ch_versions = channel.empty()

    // Reference for fastANI (wrapped in tuple for nf-core module)
    ch_reference = channel.of([ [id: 'NCTC8325'], file("${projectDir}/assets/references/NCTC8325.fasta") ])

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
        _meta, ani ->
            pass: ani >= params.species_ani_threshold
            fail: true
    }

    // Get confirmed sample IDs as a set for filtering
    ch_confirmed_ids = ch_ani_branched.pass
        .map { meta, _ani -> meta.id }
        .collect()
        .map { ids -> ids.toSet() }

    ch_confirmed_assemblies = ch_assemblies
        .combine(ch_confirmed_ids)
        .filter { meta, _fasta, confirmed_set -> meta.id in confirmed_set }
        .map { meta, fasta, _confirmed_set -> [ meta, fasta ] }


    //
    // MODULE: Download Mash RefSeq sketch (once, cached)
    //
    // --mash_db reuses a local sketch (the download host is not always reachable from containers)
    if (params.mash_db) {
        ch_mash_raw = channel.value(file(params.mash_db, checkIfExists: true))
    } else {
        MASH_REFSEQ_DB ()
        ch_mash_raw = MASH_REFSEQ_DB.out.db
    }

    //
    // MODULE: Mash screen on every assembly: names the species of rejected samples and
    // reports a second bacterial species (possible contamination) for all samples
    //
    ch_mash_db = ch_mash_raw.map { db -> [ [id: 'refseq'], db ] }

    MASH_SCREEN (
        ch_assemblies,
        ch_mash_db.collect()
    )
    ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())

    //
    // MODULE: Collate species report
    //
    // fastANI writes nothing when a genome is too divergent (another species), so emit
    // an explicit ANI-0 line for those; otherwise they'd vanish from the report.
    ch_fastani_collected = FASTANI.out.ani
        .map { meta, f -> f.text.trim() ? f.text.trim() + '\n' : "${meta.id}.scaffolds.fasta\tnone\t0\t0\t0\n" }
        .collectFile(name: 'fastani_all.tsv')

    // Each sample's top 30 Mash hits, prefixed with its sample ID (phage hits are dropped
    // in the report, so keep enough lines for a second bacterial species to survive)
    ch_mash_collected = MASH_SCREEN.out.screen
        .map { meta, screen ->
            screen.text.readLines().findAll { l -> l.trim() }
                .sort { l -> -(l.split('\t')[0] as Double) }
                .take(30)
                .collect { l -> "${meta.id}\t${l}\n" }
                .join('')
        }
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
