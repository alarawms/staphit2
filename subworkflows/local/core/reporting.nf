/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPORTING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Per-sample aggregation -> merged summary -> HTML/JSON report + visualizations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AGGREGATOR     } from '../../../modules/local/aggregator'
include { SUMMARY_MERGER } from '../../../modules/local/summary_merger'
include { REPORT         } from '../../../modules/local/report'
include { VISUALIZATION  } from '../../../modules/local/visualization'

workflow REPORTING {

    take:
    ch_aggregator_input   // channel: [ val(meta), path(trim_log), path(fastqc), path(quast), path(mlst), path(abricate), path(amrfinder), path(mash), path(spa), path(sccmec), path(agr), path(kma), path(metadata_json) ]
    ch_clusters           // path(clusters.tsv) or NO_CLUSTERS placeholder
    ch_qc_report          // path(qc_report.tsv) or NO_QC placeholder
    ch_plasmid_summary    // path(plasmid_summary.tsv) or NO_PLASMIDS placeholder
    ch_tree_file          // path(treefile) or NO_TREES placeholder

    main:

    //
    // Per-sample aggregation
    //
    AGGREGATOR ( ch_aggregator_input )

    //
    // Merge all per-sample summaries into one table
    //
    SUMMARY_MERGER (
        AGGREGATOR.out.summary.map { meta, csv -> csv }.collect()
    )

    //
    // Final report (HTML + JSON)
    //
    REPORT (
        SUMMARY_MERGER.out.summary,
        ch_clusters,
        ch_qc_report
    )

    //
    // Visualizations (plots)
    //
    VISUALIZATION (
        SUMMARY_MERGER.out.summary,
        ch_clusters,
        ch_tree_file
    )

    emit:
    per_sample_reports = AGGREGATOR.out.report       // [ val(meta), path(json) ]
    summary            = SUMMARY_MERGER.out.summary  // path(combined_summary.tsv)
    run_report         = REPORT.out.html             // path(staphit_report.html)
}
