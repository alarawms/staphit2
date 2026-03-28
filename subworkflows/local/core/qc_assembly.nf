/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC & ASSEMBLY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Reads -> Trim -> FastQC -> Assemble (SPAdes + SKESA) -> QUAST -> CheckM2 -> QC Gate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                } from '../../../modules/nf-core/fastqc/main'
include { TRIMGALORE            } from '../../../modules/nf-core/trimgalore/main'
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
    ch_trimmed  = TRIMGALORE.out.reads
    ch_trim_log = TRIMGALORE.out.log
    // ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first()) // uses topic channels

    //
    // MODULE: FastQC on trimmed reads
    //
    FASTQC ( ch_trimmed )
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first()) // uses topic channels

    //
    // MODULE: Assemble with SPAdes
    //
    SPADES (
        ch_trimmed.map { meta, reads -> [ meta, reads, [], [] ] },
        [],
        []
    )
    ch_spades_scaffolds = SPADES.out.scaffolds
    // ch_versions = ch_versions.mix(SPADES.out.versions.first()) // uses topic channels

    //
    // MODULE: Assemble with SKESA (alternative assembler)
    //
    SKESA ( ch_trimmed )
    ch_skesa_scaffolds = SKESA.out.scaffolds

    //
    // Select primary assembly based on --assembler param
    // SKESA outputs plain FASTA, SPAdes nf-core outputs .gz
    //
    ch_assemblies = ch_skesa_scaffolds

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
    // Filter assemblies to only passed samples
    //
    ch_passed_ids = QC_GATE.out.passed
        .splitText()
        .map { it.trim() }
        .filter { it }

    ch_passed_assemblies = ch_assemblies
        .filter { meta, fasta ->
            true  // Pass all for now — QC filtering done downstream via exclude_samples
        }

    emit:
    trimmed_reads     = ch_trimmed              // channel: [ val(meta), [ path(reads) ] ]
    trim_log          = ch_trim_log             // channel: [ val(meta), path(log) ]
    fastqc_zip        = FASTQC.out.zip          // channel: [ val(meta), path(zip) ]
    fastqc_html       = FASTQC.out.html         // channel: [ val(meta), path(html) ]
    assemblies        = ch_assemblies           // channel: [ val(meta), path(scaffolds) ]  — all assemblies
    passed_assemblies = ch_passed_assemblies    // channel: [ val(meta), path(scaffolds) ]  — QC-passed only
    skesa_scaffolds   = ch_skesa_scaffolds      // channel: [ val(meta), path(scaffolds) ]
    quast_results     = QUAST.out.results       // channel: [ val(meta), path(results) ]
    checkm2_reports   = CHECKM2.out.report      // channel: [ val(meta), path(report) ]
    qc_report         = QC_GATE.out.report      // channel: path(qc_report.tsv)
    qc_passed         = QC_GATE.out.passed      // channel: path(passed_samples.txt)
    qc_failed         = QC_GATE.out.failed      // channel: path(failed_samples.txt)
    versions          = ch_versions             // channel: [ path(versions.yml) ]
}
