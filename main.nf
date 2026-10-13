#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    alarawms/staphit2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/alarawms/staphit2
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAPHIT2  } from './workflows/staphit2'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_staphit2_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_staphit2_pipeline'
include { FETCH_PUBLIC            } from './subworkflows/local/core/fetch_public'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TYPED PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Command-line values arrive as strings unless the script declares a type, so
    `--use_beast` would fail schema validation as "true" (string). Only the
    boolean and numeric params are declared here; defaults stay in nextflow.config.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {
    // Public data fetch
    fetch_year_from: Integer?
    fetch_year_to: Integer?
    fetch_max: Integer?

    // References
    igenomes_ignore: Boolean

    // Phylogeny
    panaroo_threshold: Float
    panaroo_identity: Float
    panaroo_family_threshold: Float
    panaroo_entropy_filter: Float?
    snippy_mincov: Integer
    snippy_minqual: Integer
    plot_tree: Boolean
    skip_gubbins: Boolean
    use_beast: Boolean
    beast_chain_length: Integer
    beast_log_every: Integer

    // SCCmec typing
    sccmec_viz: Boolean
    sccmec_best_fit: Boolean
    sccmec_min_score: Float

    // Read subsampling
    genome_size: Integer
    target_depth: Integer

    // QC gate
    min_completeness: Integer
    max_contamination: Integer
    flye_nanohq: Boolean
    medaka_min_read_q: Float
    skip_qc_gate: Boolean
    min_assembly_size: Integer
    max_assembly_size: Integer
    min_depth_short: Float
    min_depth_long: Float
    max_contigs: Integer
    min_n50: Integer
    max_mixed_sites: Integer
    max_cgmlst_missing: Integer
    max_truncated_genes: Integer
    species_ani_threshold: Float
    skip_species_qc: Boolean

    // Outbreak clustering
    skip_cgmlst: Boolean
    cgmlst_pubmlst_scheme: Integer

    // Boilerplate (help stays undeclared: it takes a boolean or a parameter name)
    plaintext_email: Boolean
    monochrome_logs: Boolean
    help_full: Boolean
    show_hidden: Boolean
    version: Boolean
    validate_params: Boolean
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow ALARAWMS_STAPHIT2 {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    STAPHIT2 (
        samplesheet
    )
    emit:
    multiqc_report = STAPHIT2.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    def fetch_mode = params.cc || params.st
    def input_mode = params.input as Boolean

    if (!fetch_mode && !input_mode) {
        error "Provide --input samplesheet.csv, or use --cc / --st to fetch public data"
    }

    ch_reads = channel.empty()

    if (fetch_mode) {
        FETCH_PUBLIC ()
        ch_reads = ch_reads.mix(FETCH_PUBLIC.out.reads)
    }

    if (input_mode) {
        PIPELINE_INITIALISATION (
            params.version,
            params.validate_params,
            params.monochrome_logs,
            args,
            params.outdir,
            params.input,
            params.help,
            params.help_full,
            params.show_hidden
        )
        ch_reads = ch_reads.mix(PIPELINE_INITIALISATION.out.samplesheet)
    }

    ALARAWMS_STAPHIT2 ( ch_reads )

    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        ALARAWMS_STAPHIT2.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
