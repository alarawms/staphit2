/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PHYLOGENY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Three modes:
      'panaroo' — core gene alignment for tree + SNP distances
      'snippy'  — reference-based SNPs for tree + distances
      'both'    — Panaroo tree + Snippy distances (recommended for surveillance)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROKKA             } from '../../../modules/nf-core/prokka/main'
include { PANAROO            } from '../../../modules/local/panaroo'
include { SNIPPY             } from '../../../modules/local/snippy'
include { SNIPPY_CORE        } from '../../../modules/local/snippy'
include { FASTTREE           } from '../../../modules/local/fasttree'
include { IQTREE             } from '../../../modules/nf-core/iqtree/main'
include { SNPDISTS as SNPDISTS_CORE   } from '../../../modules/nf-core/snpdists/main'
include { SNPDISTS as SNPDISTS_SNIPPY } from '../../../modules/nf-core/snpdists/main'
include { BEAST2                      } from '../../../modules/local/beast2'

workflow PHYLOGENY {

    take:
    ch_assemblies     // channel: [ val(meta), path(assembly) ]
    ch_trimmed_reads  // channel: [ val(meta), path(reads) ]

    main:

    ch_versions   = channel.empty()
    ch_snippy_dists = channel.empty()
    ch_beast2_tree  = channel.empty()

    //
    // MODULE: Annotate assemblies with Prokka (needed for Panaroo)
    //
    PROKKA ( ch_assemblies, [], [] )
    ch_versions = ch_versions.mix(PROKKA.out.versions.first())

    //
    // PANAROO: core gene alignment → tree (runs if 'panaroo' or 'both')
    //
    if (params.phylo_method == 'panaroo' || params.phylo_method == 'both') {
        PANAROO ( PROKKA.out.gff.map { _meta, gff -> gff }.collect() )
        ch_panaroo_aln = PANAROO.out.aln

        // Tree from Panaroo alignment
        if (params.tree_builder == 'iqtree') {
            IQTREE (
                ch_panaroo_aln.map { aln -> [ [id: 'core'], aln, [] ] },
                [], [], [], [], [], [], [], [], [], [], [], []
            )
            ch_tree = IQTREE.out.phylogeny
            ch_versions = ch_versions.mix(IQTREE.out.versions)
        } else {
            FASTTREE ( ch_panaroo_aln )
            ch_tree = FASTTREE.out.tree.map { tree -> [ [id: 'core'], tree ] }
        }

        // SNP distances from core genes (for tree visualization)
        SNPDISTS_CORE ( ch_panaroo_aln.map { aln -> [ [id: 'panaroo_core'], aln ] } )
        ch_versions = ch_versions.mix(SNPDISTS_CORE.out.versions)
    }

    //
    // SNIPPY: reference-based SNPs → distances for clustering (runs if 'snippy' or 'both')
    //
    if (params.phylo_method == 'snippy' || params.phylo_method == 'both') {
        ch_ref = channel.fromPath(params.reference, checkIfExists: true)
        SNIPPY ( ch_trimmed_reads, ch_ref.collect() )
        SNIPPY_CORE (
            SNIPPY.out.results.map { _meta, dir -> dir }.collect(),
            ch_ref.collect()
        )

        // SNP distances from whole-genome (for outbreak clustering)
        SNPDISTS_SNIPPY ( SNIPPY_CORE.out.aln.map { aln -> [ [id: 'snippy_wgs'], aln ] } )
        ch_snippy_dists = SNPDISTS_SNIPPY.out.tsv
        ch_versions = ch_versions.mix(SNPDISTS_SNIPPY.out.versions)

        // Optional Bayesian phylodynamics (--use_beast)
        if (params.use_beast) {
            ch_metadata = params.metadata
                ? channel.fromPath(params.metadata, checkIfExists: true)
                : channel.fromPath("${projectDir}/assets/NO_METADATA")
            BEAST2(
                SNIPPY_CORE.out.aln,
                ch_metadata,
                channel.value(file("${projectDir}/bin/beast2_prep.py"))
            )
            ch_beast2_tree = BEAST2.out.mcc_tree
        }

        // If snippy-only mode, also build tree from snippy alignment
        if (params.phylo_method == 'snippy') {
            if (params.tree_builder == 'iqtree') {
                IQTREE (
                    SNIPPY_CORE.out.aln.map { aln -> [ [id: 'core'], aln, [] ] },
                    [], [], [], [], [], [], [], [], [], [], [], []
                )
                ch_tree = IQTREE.out.phylogeny
                ch_versions = ch_versions.mix(IQTREE.out.versions)
            } else {
                FASTTREE ( SNIPPY_CORE.out.aln )
                ch_tree = FASTTREE.out.tree.map { tree -> [ [id: 'core'], tree ] }
            }
            SNPDISTS_CORE ( SNIPPY_CORE.out.aln.map { aln -> [ [id: 'snippy_core'], aln ] } )
            ch_versions = ch_versions.mix(SNPDISTS_CORE.out.versions)
        }
    }

    //
    // Select which SNP distances to use for clustering
    // 'both'/'panaroo' mode: Panaroo core gene distances (preserves cross-lineage
    //   resolution in diverse multi-ST collections; Snippy WGS core collapses when
    //   the reference-based core narrows across divergent lineages).
    // 'snippy' mode: Snippy WGS distances (only option when Panaroo not run).
    //
    ch_cluster_dists = (params.phylo_method == 'snippy')
        ? ch_snippy_dists
        : SNPDISTS_CORE.out.tsv

    emit:
    prokka_gff    = PROKKA.out.gff        // [ val(meta), path(gff) ]
    tree          = ch_tree               // [ val(meta), path(treefile) ]
    snp_dists     = SNPDISTS_CORE.out.tsv // [ val(meta), path(tsv) ] — for tree visualization
    cluster_dists = ch_cluster_dists      // [ val(meta), path(tsv) ] — for outbreak clustering
    beast2_tree   = ch_beast2_tree        // path(beast_mcc.tree) — MCC annotated tree (empty if --use_beast not set)
    versions      = ch_versions
}
