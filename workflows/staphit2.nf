/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STAPHIT2 — MRSA GENOMIC SURVEILLANCE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Phase 1: QC & Assembly (auto short / hybrid / long-read per sample)
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
include { SPECIES_QC       } from '../subworkflows/local/core/species_qc'

// Reporting modules (wired directly for flexible channel joining)
include { AGGREGATOR       } from '../modules/local/aggregator'
include { VALIDATE_METADATA } from '../modules/local/validate_metadata'
include { SUMMARY_MERGER   } from '../modules/local/summary_merger'
include { REPORT           } from '../modules/local/report'
include { VISUALIZATION    } from '../modules/local/visualization'
include { PLOT_TREE        } from '../modules/local/plot_tree'
include { ITOL_EXPORT      } from '../modules/local/itol_export'

// MultiQC
include { MULTIQC          } from '../modules/nf-core/multiqc/main'
include { SAMPLE_STATUS    } from '../modules/local/sample_status'

workflow STAPHIT2 {

    take:
    ch_samplesheet  // channel: samplesheet read in from --input

    main:

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    // ── Phase 1: QC & Assembly ──────────────────────────────────────────────
    QC_ASSEMBLY ( ch_samplesheet )
    ch_trimmed    = QC_ASSEMBLY.out.trimmed_reads
    ch_versions   = ch_versions.mix(QC_ASSEMBLY.out.versions)

    // ── Phase 1b: Species Confirmation ──────────────────────────────────────
    if (!params.skip_species_qc) {
        SPECIES_QC ( QC_ASSEMBLY.out.passed_assemblies )
        ch_assemblies        = SPECIES_QC.out.confirmed_assemblies
        ch_species_excluded  = SPECIES_QC.out.species_excluded
        ch_versions          = ch_versions.mix(SPECIES_QC.out.versions)
    } else {
        ch_assemblies        = QC_ASSEMBLY.out.passed_assemblies
        ch_species_excluded  = channel.value(file("${projectDir}/assets/NO_SPECIES_QC"))
    }

    // ── Phase 2: S. aureus Typing ───────────────────────────────────────────
    // Reads per sample: trimmed short reads (short/hybrid) or ONT reads (long-only).
    // Used by KMA and as the SCCmec typer's read fallback for split cassettes.
    // Snippy cannot align raw ONT reads, so long-only samples use their assembly (--ctgs).
    def is_long = { meta, _x -> meta.mode == 'long' }
    ch_kma_reads    = ch_trimmed.mix(QC_ASSEMBLY.out.long_reads.filter(is_long))
    ch_snippy_reads = ch_trimmed.mix(ch_assemblies.filter(is_long))

    SA_TYPING ( ch_assemblies, ch_kma_reads )
    ch_versions = ch_versions.mix(SA_TYPING.out.versions)

    // ── Phase 3: AMR Detection ──────────────────────────────────────────────
    AMR_DETECTION ( ch_assemblies, ch_kma_reads )

    // ── Phase 4: Plasmid Analysis ───────────────────────────────────────────
    PLASMID_ANALYSIS ( ch_assemblies )

    // ── Phase 5: Phylogeny ──────────────────────────────────────────────────
    PHYLOGENY ( ch_assemblies, ch_snippy_reads )
    ch_versions = ch_versions.mix(PHYLOGENY.out.versions)

    // ── Phase 6: Per-sample Aggregation ─────────────────────────────────────
    // Normalize all channels to [sample_id, path] before joining
    // (nf-core modules may modify meta maps, breaking join on meta)
    def to_id = { meta, path -> [ meta.id, path ] }

    // Long-only samples have no fastp/FastQC output: give them empty placeholders
    ch_long_only_ids = QC_ASSEMBLY.out.long_reads.filter(is_long).map { meta, _lr -> meta.id }
    ch_trim_log  = QC_ASSEMBLY.out.trim_log.map(to_id)
        .mix(ch_long_only_ids.map { id -> [ id, file("${projectDir}/assets/NO_TRIMLOG") ] })
    ch_fastqc    = QC_ASSEMBLY.out.fastqc_zip.map(to_id)
        .mix(ch_long_only_ids.map { id -> [ id, file("${projectDir}/assets/NO_FASTQC") ] })

    ch_agg_in = ch_trim_log
        .join(ch_fastqc)
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

    // Broadcast metadata JSON to all samples (or placeholder if not provided).
    // A CSV (e.g. from bin/fetch_metadata.py) is validated and converted to JSON first.
    if (params.metadata && params.metadata.toString().endsWith('.csv')) {
        VALIDATE_METADATA (
            file(params.metadata, checkIfExists: true),
            params.antibiogram ? file(params.antibiogram, checkIfExists: true)
                               : file("${projectDir}/assets/NO_ANTIBIOGRAM"),
            file(params.input, checkIfExists: true)
        )
        ch_metadata = VALIDATE_METADATA.out.json
    } else {
        ch_metadata = params.metadata
            ? channel.fromPath(params.metadata, checkIfExists: true)
            : channel.of(file("${projectDir}/assets/NO_METADATA"))
    }

    ch_agg_final = ch_agg_in.combine(ch_metadata)

    AGGREGATOR ( ch_agg_final )

    // Merge per-sample summaries into a single run-level table
    SUMMARY_MERGER (
        AGGREGATOR.out.summary.map { _meta, csv -> csv }.collect()
    )

    // ── Phase 7: Outbreak Clustering ────────────────────────────────────────
    ch_cgmlst_dists = channel.of(file("${projectDir}/assets/NO_CGMLST_DISTS"))

    // Use Snippy whole-genome distances for clustering (more discriminatory)
    // Falls back to Panaroo core gene distances if Snippy not run
    CLUSTERING (
        PHYLOGENY.out.cluster_dists.map { _meta, tsv -> tsv },
        ch_cgmlst_dists,
        SUMMARY_MERGER.out.summary
    )

    // ── Phase 8: Reporting & Visualization ──────────────────────────────────
    // Which stages each sample reached; the first one missing is why it was dropped
    ch_stages = ch_samplesheet.map { meta, _sr, _lr -> "input\t${meta.id}" }
        .mix(
            QC_ASSEMBLY.out.primary_scaffolds.map { meta, _x -> "assembled\t${meta.id}" },
            QC_ASSEMBLY.out.primary_scaffolds.map { meta, fasta -> "size\t${meta.id}\t${fasta.size()}" },
            QC_ASSEMBLY.out.assemblies.map { meta, _x -> "size_pass\t${meta.id}" },
            QC_ASSEMBLY.out.passed_assemblies.map { meta, _x -> "qc_pass\t${meta.id}" },
            ch_assemblies.map { meta, _x -> "species_pass\t${meta.id}" }
        )
        .collectFile(name: 'stages.tsv', newLine: true, sort: true)

    SAMPLE_STATUS (
        ch_stages,
        QC_ASSEMBLY.out.qc_report,
        ch_species_excluded,
        SUMMARY_MERGER.out.summary,
        PHYLOGENY.out.tree.map { _meta, tree -> tree }.ifEmpty(file("${projectDir}/assets/NO_TREE")).first()
    )

    REPORT (
        SUMMARY_MERGER.out.summary,
        CLUSTERING.out.report,
        QC_ASSEMBLY.out.qc_report,
        PLASMID_ANALYSIS.out.plasmid_summary,
        SAMPLE_STATUS.out.tsv
    )

    VISUALIZATION (
        SUMMARY_MERGER.out.summary,
        CLUSTERING.out.clusters
    )

    if (params.plot_tree) {
        PLOT_TREE (
            PHYLOGENY.out.tree.map { _meta, tree -> tree }.first(),
            SUMMARY_MERGER.out.summary,
            CLUSTERING.out.clusters
        )
    }

    // iTOL annotation files (same harmonized labels/colours as the PDF tree)
    ITOL_EXPORT (
        PHYLOGENY.out.tree.map { _meta, tree -> tree }.first(),
        SUMMARY_MERGER.out.summary,
        CLUSTERING.out.clusters
    )

    // ── MultiQC ─────────────────────────────────────────────────────────────
    ch_multiqc_files = ch_multiqc_files.mix(
        QC_ASSEMBLY.out.fastqc_zip.collect { row -> row[1] },
        QC_ASSEMBLY.out.trim_log.collect { row -> row[1] },
        QC_ASSEMBLY.out.nanoplot_stats.collect { row -> row[1] }
    )

    // Collate and save software versions: versions.yml files plus [process, tool, version] topic tuples
    ch_topic_versions = channel.topic('versions')
        .unique()
        .map { process, tool, version -> [ process.tokenize(':')[-1], "    ${tool}: ${version.toString().trim()}" ] }
        .groupTuple()
        .map { process, tools -> "${process}:\n${tools.unique().sort().join('\n')}" }
    softwareVersionsToYAML(ch_versions)
        .mix(ch_topic_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'staphit2_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config
        ? channel.fromPath(params.multiqc_config, checkIfExists: true)
        : channel.empty()
    ch_multiqc_logo          = params.multiqc_logo
        ? channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files    = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
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
