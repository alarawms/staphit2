/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STAPHIT2 — MRSA GENOMIC SURVEILLANCE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Phase 1: QC & Assembly
    Phase 2: S. aureus Typing (MLST, spa, SCCmec, agr, Mash)
    Phase 3: AMR Detection (AMRFinderPlus, ABRicate, KMA)
    Phase 4: Plasmid Analysis (MOB-suite)
    Phase 5: Phylogeny (Prokka → Panaroo/Snippy → FastTree/IQ-TREE → SNP dists)
    Phase 6: Per-sample Aggregation & Summary Merge
    Phase 7: Outbreak Clustering
    Phase 8: Reporting & Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_staphit2_pipeline'

// Subworkflows
include { QC_ASSEMBLY      } from '../subworkflows/local/core/qc_assembly'
include { SA_TYPING        } from '../subworkflows/local/species/s_aureus/typing'
include { AMR_DETECTION    } from '../subworkflows/local/core/amr'
include { PHYLOGENY        } from '../subworkflows/local/core/phylogeny'
include { PLASMID_ANALYSIS } from '../subworkflows/local/core/plasmids'
include { CLUSTERING       } from '../subworkflows/local/core/clustering'

// Reporting modules (wired directly for flexible channel joining)
include { AGGREGATOR       } from '../modules/local/aggregator'
include { SUMMARY_MERGER   } from '../modules/local/summary_merger'
include { REPORT           } from '../modules/local/report'
include { VISUALIZATION    } from '../modules/local/visualization'

// MultiQC
include { MULTIQC          } from '../modules/nf-core/multiqc/main'

workflow STAPHIT2 {

    take:
    ch_samplesheet  // channel: samplesheet read in from --input

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // ── Phase 1: QC & Assembly ──────────────────────────────────────────────
    QC_ASSEMBLY ( ch_samplesheet )
    ch_assemblies = QC_ASSEMBLY.out.passed_assemblies
    ch_trimmed    = QC_ASSEMBLY.out.trimmed_reads
    ch_versions   = ch_versions.mix(QC_ASSEMBLY.out.versions)

    // ── Phase 2: S. aureus Typing ───────────────────────────────────────────
    SA_TYPING ( ch_assemblies )
    ch_versions = ch_versions.mix(SA_TYPING.out.versions)

    // ── Phase 3: AMR Detection ──────────────────────────────────────────────
    AMR_DETECTION ( ch_assemblies, ch_trimmed )

    // ── Phase 4: Plasmid Analysis ───────────────────────────────────────────
    PLASMID_ANALYSIS ( ch_assemblies )

    // ── Phase 5: Phylogeny ──────────────────────────────────────────────────
    PHYLOGENY ( ch_assemblies, ch_trimmed )
    ch_versions = ch_versions.mix(PHYLOGENY.out.versions)

    // ── Phase 6: Per-sample Aggregation ─────────────────────────────────────
    // Normalize all channels to [sample_id, path] before joining
    // (nf-core modules may modify meta maps, breaking join on meta)
    def to_id = { meta, path -> [ meta.id, path ] }

    ch_agg_in = QC_ASSEMBLY.out.trim_log.map(to_id)
        .join(QC_ASSEMBLY.out.fastqc_zip.map(to_id))
        .join(QC_ASSEMBLY.out.quast_results.map(to_id))
        .join(SA_TYPING.out.mlst.map(to_id))
        .join(AMR_DETECTION.out.abricate.map(to_id))
        .join(AMR_DETECTION.out.amrfinder.map(to_id))
        .join(SA_TYPING.out.mash.map(to_id))
        .join(SA_TYPING.out.spa.map(to_id))
        .join(SA_TYPING.out.sccmec.map(to_id))
        .join(SA_TYPING.out.agr.map(to_id))
        .join(AMR_DETECTION.out.kma.map(to_id))
        .map { items ->
            // Reconstruct meta from sample_id
            def sid = items[0]
            def meta = [id: sid]
            [ meta ] + items[1..-1]
        }

    // Broadcast metadata JSON to all samples (or placeholder if not provided)
    ch_metadata = params.metadata
        ? Channel.fromPath(params.metadata, checkIfExists: true)
        : Channel.of(file('NO_METADATA'))

    ch_agg_final = ch_agg_in.combine(ch_metadata)

    AGGREGATOR ( ch_agg_final )

    // Merge per-sample summaries into a single run-level table
    SUMMARY_MERGER (
        AGGREGATOR.out.summary.map { meta, csv -> csv }.collect()
    )

    // ── Phase 7: Outbreak Clustering ────────────────────────────────────────
    ch_cgmlst_dists = Channel.of(file('NO_CGMLST_DISTS'))

    // Use Snippy whole-genome distances for clustering (more discriminatory)
    // Falls back to Panaroo core gene distances if Snippy not run
    CLUSTERING (
        PHYLOGENY.out.cluster_dists.map { meta, tsv -> tsv },
        ch_cgmlst_dists,
        SUMMARY_MERGER.out.summary
    )

    // ── Phase 8: Reporting & Visualization ──────────────────────────────────
    REPORT (
        SUMMARY_MERGER.out.summary,
        CLUSTERING.out.report,
        QC_ASSEMBLY.out.qc_report,
        PLASMID_ANALYSIS.out.plasmid_summary
    )

    VISUALIZATION (
        SUMMARY_MERGER.out.summary,
        CLUSTERING.out.clusters
    )

    // ── MultiQC ─────────────────────────────────────────────────────────────
    ch_multiqc_files = ch_multiqc_files.mix(
        QC_ASSEMBLY.out.fastqc_zip.collect { it[1] }
    )

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'staphit2_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files    = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList()  // path/to/multiqc_report.html
    versions       = ch_versions                   // [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
