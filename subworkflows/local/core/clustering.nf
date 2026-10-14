/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLUSTERING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cgMLST (PubMLST scheme 20, chewBBACA) allele distances -> tiered outbreak clusters.
    cgMLST distances do not depend on which other samples are in the run, unlike a
    core alignment across lineages, so they are the default (--cluster_by cgmlst).
    SNP distances are still used for transmission pairs and cluster diameters.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CGMLST_SCHEMA    } from '../../../modules/local/chewbbaca'
include { CHEWBBACA_PREP   } from '../../../modules/local/chewbbaca'
include { CHEWBBACA_ALLELE } from '../../../modules/local/chewbbaca'
include { CGMLST_DISTS     } from '../../../modules/local/chewbbaca'
include { OUTBREAK_CLUSTER } from '../../../modules/local/outbreak_cluster'
include { GRAPETREE        } from '../../../modules/local/grapetree'
include { MST_PLOT         } from '../../../modules/local/grapetree'

workflow CLUSTERING {

    take:
    ch_snp_dists     // path(snp_dists.tsv)
    ch_assemblies    // [ val(meta), path(assembly) ]
    ch_summary       // path(combined_summary.tsv)

    main:

    if (params.skip_cgmlst) {
        ch_cgmlst_dists = channel.value(file("${projectDir}/assets/NO_CGMLST_DISTS"))
        ch_cgmlst_stats = channel.empty()
    } else {
        // --cgmlst_schema reuses a prepared schema (<outdir>/databases/cgmlst_schema/prepared_schema)
        if (params.cgmlst_schema) {
            ch_schema = channel.value(file(params.cgmlst_schema, checkIfExists: true))
        } else {
            CGMLST_SCHEMA ()
            CHEWBBACA_PREP ( CGMLST_SCHEMA.out.fasta )
            ch_schema = CHEWBBACA_PREP.out.schema
        }
        ch_calls_in = ch_assemblies
            .map { meta, fasta -> [ meta.id, fasta ] }
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { pairs -> [ pairs.collect { p -> p[0] }, pairs.collect { p -> p[1] } ] }
        CHEWBBACA_ALLELE ( ch_calls_in, ch_schema )
        CGMLST_DISTS ( CHEWBBACA_ALLELE.out.profiles )
        ch_cgmlst_dists = CGMLST_DISTS.out.dists
        ch_cgmlst_stats = CGMLST_DISTS.out.stats
    }

    OUTBREAK_CLUSTER ( ch_snp_dists, ch_cgmlst_dists, ch_summary )

    // cgMLST minimum spanning tree, annotated with the summary metadata and clusters
    if (!params.skip_cgmlst) {
        GRAPETREE ( CGMLST_DISTS.out.grapetree )
        MST_PLOT ( GRAPETREE.out.tree, ch_summary, OUTBREAK_CLUSTER.out.clusters )
    }

    emit:
    clusters     = OUTBREAK_CLUSTER.out.clusters  // path(clusters.tsv)
    report       = OUTBREAK_CLUSTER.out.report    // path(cluster_report.json)
    pairs        = OUTBREAK_CLUSTER.out.pairs     // path(transmission_pairs.tsv)
    cgmlst_stats = ch_cgmlst_stats                // path(cgmlst_stats.tsv); empty with --skip_cgmlst
}
