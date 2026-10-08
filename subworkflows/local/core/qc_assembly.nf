/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC & ASSEMBLY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read mode is auto-detected per sample (meta.mode):
      short  : fastp -> Rasusa -> FastQC -> SKESA/SPAdes
      hybrid : short-read QC as above + NanoPlot -> Dragonflye (Flye + Polypolish)
      long   : NanoPlot -> Dragonflye (Flye + Racon[/Medaka])
    All assemblies -> QUAST -> CheckM2 -> QC Gate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                } from '../../../modules/nf-core/fastqc/main'
include { FASTP                 } from '../../../modules/local/fastp'
include { RASUSA                } from '../../../modules/nf-core/rasusa/main'
include { SPADES                } from '../../../modules/nf-core/spades/main'
include { QUAST                 } from '../../../modules/nf-core/quast/main'
include { SKESA                 } from '../../../modules/local/skesa'
include { CHECKM2_DB            } from '../../../modules/local/checkm2_db'
include { CHECKM2               } from '../../../modules/local/checkm2'
include { QC_GATE               } from '../../../modules/local/qc_gate'
include { NANOPLOT              } from '../../../modules/local/nanoplot'
include { DRAGONFLYE            } from '../../../modules/local/dragonflye'

workflow QC_ASSEMBLY {

    take:
    ch_reads  // channel: [ val(meta), [ short reads ], [ long reads ] ] — meta.mode: short | hybrid | long

    main:

    ch_versions = channel.empty()

    // ── Short reads (short + hybrid samples): trim → subsample → FastQC ─────
    ch_short = ch_reads.filter { _meta, sr, _lr -> sr }.map { meta, sr, _lr -> [ meta, sr ] }

    FASTP ( ch_short )
    ch_trim_log = FASTP.out.json

    ch_rasusa_input = FASTP.out.reads.map { meta, reads ->
        [ meta, reads, params.genome_size ?: 2800000 ]
    }
    RASUSA ( ch_rasusa_input, params.target_depth ?: 100 )
    ch_trimmed = RASUSA.out.reads

    FASTQC ( ch_trimmed )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // ── Long reads (long + hybrid samples): NanoPlot QC ─────────────────────
    ch_long = ch_reads.filter { _meta, _sr, lr -> lr }.map { meta, _sr, lr -> [ meta, lr ] }
    NANOPLOT ( ch_long )

    // ── Assembly, routed by mode ────────────────────────────────────────────
    ch_trimmed_short_only = ch_trimmed.filter { meta, _reads -> meta.mode == 'short' }

    if (params.assembler == 'spades') {
        SPADES (
            ch_trimmed_short_only.map { meta, reads -> [ meta, reads, [], [] ] },
            [],
            []
        )
        ch_short_scaffolds = SPADES.out.scaffolds
    } else {
        SKESA ( ch_trimmed_short_only )
        ch_short_scaffolds = SKESA.out.scaffolds
    }

    // Long-only: Flye; hybrid: Flye + Polypolish with the trimmed short reads
    ch_dragonflye_in = ch_long.filter { meta, _lr -> meta.mode == 'long' }
        .map { meta, lr -> [ meta, lr, [] ] }
        .mix(
            ch_long.filter { meta, _lr -> meta.mode == 'hybrid' }
                .map { meta, lr -> [ meta.id, meta, lr ] }
                .join(ch_trimmed.map { meta, sr -> [ meta.id, sr ] })
                .map { _id, meta, lr, sr -> [ meta, lr, sr ] }
        )
    DRAGONFLYE ( ch_dragonflye_in )

    ch_primary_scaffolds = ch_short_scaffolds.mix(DRAGONFLYE.out.scaffolds)

    ch_skesa_branched = ch_primary_scaffolds.branch {
        _meta, fasta ->
            pass: fasta.size() > 500000
            fail: true
    }
    ch_assemblies = ch_skesa_branched.pass

    // Log dropped samples (assembly too small)
    ch_assembly_dropped = ch_skesa_branched.fail
        .map { meta, fasta -> "${meta.id}\tassembly_too_small\t${fasta.size()}" }
        .collect()
        .map { lines -> lines.join('\n') }
        .ifEmpty('')

    //
    // MODULE: Assembly QC with QUAST
    //
    QUAST (
        ch_assemblies,
        [[:], []],
        [[:], []]
    )
    // ch_versions = ch_versions.mix(QUAST.out.versions.first()) // uses topic channels

    //
    // CheckM2 + QC gate. --skip_qc_gate skips both (and the ~3 GB CheckM2 database download).
    //
    if (params.skip_qc_gate) {
        ch_passed_assemblies = ch_assemblies
        ch_checkm2_reports   = channel.empty()
        ch_qc_report         = channel.value(file("${projectDir}/assets/NO_QC"))
        ch_qc_passed         = channel.empty()
        ch_qc_failed         = channel.empty()
    } else {
        // --checkm2_db reuses a local uniref100.KO.1.dmnd
        if (params.checkm2_db) {
            ch_checkm2_db = channel.value(file(params.checkm2_db, checkIfExists: true))
        } else {
            CHECKM2_DB ()
            ch_checkm2_db = CHECKM2_DB.out.db
        }

        CHECKM2 ( ch_assemblies, ch_checkm2_db )

        // Filter by completeness/contamination
        QC_GATE ( CHECKM2.out.report.map { _meta, report -> report }.collect() )

        ch_passed_ids = QC_GATE.out.passed
            .splitText()
            .map { s -> s.trim() }
            .filter { s -> s }
            .collect()
            .map { ids -> ids.toSet() }

        ch_passed_assemblies = ch_assemblies
            .combine(ch_passed_ids)
            .filter { meta, _fasta, passed_set -> meta.id in passed_set }
            .map { meta, fasta, _passed_set -> [ meta, fasta ] }

        ch_checkm2_reports = CHECKM2.out.report
        ch_qc_report       = QC_GATE.out.report
        ch_qc_passed       = QC_GATE.out.passed
        ch_qc_failed       = QC_GATE.out.failed
    }

    emit:
    trimmed_reads     = ch_trimmed              // channel: [ val(meta), [ path(reads) ] ]  — short + hybrid only
    long_reads        = ch_long                 // channel: [ val(meta), path(long fastq) ]  — long + hybrid only
    nanoplot_stats    = NANOPLOT.out.stats      // channel: [ val(meta), path(NanoStats.txt) ]
    trim_log          = ch_trim_log             // channel: [ val(meta), path(fastp.json) ]
    fastqc_zip        = FASTQC.out.zip          // channel: [ val(meta), path(zip) ]
    fastqc_html       = FASTQC.out.html         // channel: [ val(meta), path(html) ]
    assemblies        = ch_assemblies           // channel: [ val(meta), path(scaffolds) ]  — all assemblies
    passed_assemblies = ch_passed_assemblies    // channel: [ val(meta), path(scaffolds) ]  — QC-passed only
    primary_scaffolds = ch_primary_scaffolds    // channel: [ val(meta), path(scaffolds) ]
    quast_results     = QUAST.out.results       // channel: [ val(meta), path(results) ]
    checkm2_reports   = ch_checkm2_reports      // channel: [ val(meta), path(report) ]
    qc_report         = ch_qc_report            // channel: path(qc_report.tsv)
    qc_passed         = ch_qc_passed            // channel: path(passed_samples.txt)
    qc_failed         = ch_qc_failed            // channel: path(failed_samples.txt)
    assembly_dropped  = ch_assembly_dropped     // channel: val(string) — TSV of samples dropped by size filter
    versions          = ch_versions             // channel: [ path(versions.yml) ]
}
