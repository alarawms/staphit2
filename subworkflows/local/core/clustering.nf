/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLUSTERING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Outbreak detection from SNP/cgMLST distance matrices + summary data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { OUTBREAK_CLUSTER } from '../../../modules/local/outbreak_cluster'

workflow CLUSTERING {

    take:
    ch_snp_dists     // path(snp_dists.tsv)
    ch_cgmlst_dists  // path(cgmlst_dists.tsv) or NO_CGMLST_DISTS placeholder
    ch_summary       // path(combined_summary.tsv)

    main:

    OUTBREAK_CLUSTER ( ch_snp_dists, ch_cgmlst_dists, ch_summary )

    emit:
    clusters = OUTBREAK_CLUSTER.out.clusters  // path(clusters.tsv)
    report   = OUTBREAK_CLUSTER.out.report    // path(cluster_report.json)
    pairs    = OUTBREAK_CLUSTER.out.pairs     // path(transmission_pairs.tsv)
}
