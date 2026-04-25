/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC & ASSEMBLY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Reads -> Trim -> FastQC -> Assemble (SPAdes + SKESA) -> QUAST -> CheckM2 -> QC Gate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                } from '../../../modules/nf-core/fastqc/main'
include { TRIMGALORE            } from '../../../modules/nf-core/trimgalore/main'
include { RASUSA                } from '../../../modules/nf-core/rasusa/main'
include { SPADES                } from '../../../modules/nf-core/spades/main'
include { QUAST                 } from '../../../modules/nf-core/quast/main'
include { SKESA                 } from '../../../modules/local/skesa'
include { CHECKM2_DB            } from '../../../modules/local/checkm2_db'
include { CHECKM2               } from '../../../modules/local/checkm2'
include { QC_GATE               } from '../../../modules/local/qc_gate'

workflow QC_ASSEMBLY {

    take:
    ch_reads  // channel: [ val(meta), [ path(reads) ] ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Trim reads with Trim Galore
    //
    TRIMGALORE ( ch_reads )
    ch_trim_log = TRIMGALORE.out.log
    // ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first()) // uses topic channels

    //
    // MODULE: Subsample reads to target coverage with Rasusa
    // S. aureus genome ~2.8 Mb, target 100x
    //
    ch_rasusa_input = TRIMGALORE.out.reads.map { meta, reads ->
        [ meta, reads, params.genome_size ?: 2800000 ]
    }
    RASUSA ( ch_rasusa_input, params.target_depth ?: 100 )
    ch_trimmed = RASUSA.out.reads

    //
    // MODULE: FastQC on subsampled reads
    //
    FASTQC ( ch_trimmed )
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first()) // uses topic channels

    //
    // MODULE: Assemble with selected assembler (--assembler skesa|spades)
    //
    if (params.assembler == 'spades') {
        SPADES (
            ch_trimmed.map { meta, reads -> [ meta, reads, [], [] ] },
            [],
            []
        )
        ch_primary_scaffolds = SPADES.out.scaffolds
    } else {
        SKESA ( ch_trimmed )
        ch_primary_scaffolds = SKESA.out.scaffolds
    }

    //
    // Filter out junk assemblies (S. aureus ~2.8 Mb; anything under 500 KB is junk)
    //
    ch_skesa_branched = ch_primary_scaffolds.branch {
        meta, fasta ->
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
    // MODULE: Download CheckM2 database (runs once)
    //
    CHECKM2_DB ()

    //
    // MODULE: Assess assembly completeness with CheckM2
    //
    CHECKM2 (
        ch_assemblies,
        CHECKM2_DB.out.db
    )

    //
    // MODULE: QC Gate — filter samples by completeness/contamination
    //
    QC_GATE (
        CHECKM2.out.report.map { meta, report -> report }.collect()
    )

    //
    // Filter assemblies to only QC-passed samples (skip filter when --skip_qc_gate)
    //
    if (params.skip_qc_gate) {
        ch_passed_assemblies = ch_assemblies
    } else {
        ch_passed_ids = QC_GATE.out.passed
            .splitText()
            .map { it.trim() }
            .filter { it }
            .collect()
            .map { it.toSet() }

        ch_passed_assemblies = ch_assemblies
            .combine(ch_passed_ids)
            .filter { meta, fasta, passed_set -> meta.id in passed_set }
            .map { meta, fasta, passed_set -> [ meta, fasta ] }
    }

    emit:
    trimmed_reads     = ch_trimmed              // channel: [ val(meta), [ path(reads) ] ]
    trim_log          = ch_trim_log             // channel: [ val(meta), path(log) ]
    fastqc_zip        = FASTQC.out.zip          // channel: [ val(meta), path(zip) ]
    fastqc_html       = FASTQC.out.html         // channel: [ val(meta), path(html) ]
    assemblies        = ch_assemblies           // channel: [ val(meta), path(scaffolds) ]  — all assemblies
    passed_assemblies = ch_passed_assemblies    // channel: [ val(meta), path(scaffolds) ]  — QC-passed only
    primary_scaffolds = ch_primary_scaffolds    // channel: [ val(meta), path(scaffolds) ]
    quast_results     = QUAST.out.results       // channel: [ val(meta), path(results) ]
    checkm2_reports   = CHECKM2.out.report      // channel: [ val(meta), path(report) ]
    qc_report         = QC_GATE.out.report      // channel: path(qc_report.tsv)
    qc_passed         = QC_GATE.out.passed      // channel: path(passed_samples.txt)
    qc_failed         = QC_GATE.out.failed      // channel: path(failed_samples.txt)
    assembly_dropped  = ch_assembly_dropped     // channel: val(string) — TSV of samples dropped by size filter
    versions          = ch_versions             // channel: [ path(versions.yml) ]
}
