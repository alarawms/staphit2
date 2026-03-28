/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PHYLOGENY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Prokka annotation -> Panaroo (pan-genome) or Snippy (reference-based)
    -> core alignment -> tree (FastTree or IQ-TREE) + SNP distances
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROKKA       } from '../../../modules/nf-core/prokka/main'
include { PANAROO      } from '../../../modules/local/panaroo'
include { SNIPPY       } from '../../../modules/local/snippy'
include { SNIPPY_CORE  } from '../../../modules/local/snippy'
include { FASTTREE     } from '../../../modules/local/fasttree'
include { IQTREE       } from '../../../modules/nf-core/iqtree/main'
include { SNPDISTS     } from '../../../modules/nf-core/snpdists/main'

workflow PHYLOGENY {

    take:
    ch_assemblies     // channel: [ val(meta), path(assembly) ]
    ch_trimmed_reads  // channel: [ val(meta), path(reads) ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Annotate assemblies with Prokka (needed for Panaroo)
    //
    PROKKA ( ch_assemblies, [], [] )
    ch_versions = ch_versions.mix(PROKKA.out.versions.first())

    //
    // Core alignment: Snippy (reference-based) or Panaroo (pan-genome)
    //
    if (params.phylo_method == 'snippy') {
        ch_ref = Channel.fromPath(params.reference, checkIfExists: true)
        SNIPPY ( ch_trimmed_reads, ch_ref.collect() )
        SNIPPY_CORE (
            SNIPPY.out.results.map { meta, dir -> dir }.collect(),
            ch_ref.collect()
        )
        ch_core_aln = SNIPPY_CORE.out.aln
    } else {
        // Default: Panaroo pan-genome alignment
        PANAROO ( PROKKA.out.gff.map { meta, gff -> gff }.collect() )
        ch_core_aln = PANAROO.out.aln
    }

    //
    // Tree building: FastTree (default) or IQ-TREE
    //
    if (params.tree_builder == 'iqtree') {
        // IQTREE nf-core module: tuple(meta, alignment, tree) + 11 optional path args
        IQTREE (
            ch_core_aln.map { aln -> [ [id: 'core'], aln, [] ] },
            [], [], [], [], [], [], [], [], [], []
        )
        ch_tree = IQTREE.out.phylogeny
        ch_versions = ch_versions.mix(IQTREE.out.versions.first())
    } else {
        // Default: FastTree
        FASTTREE ( ch_core_aln )
        ch_tree = FASTTREE.out
    }

    //
    // SNP distance matrix
    //
    SNPDISTS ( ch_core_aln.map { aln -> [ [id: 'core'], aln ] } )
    ch_versions = ch_versions.mix(SNPDISTS.out.versions.first())

    emit:
    prokka_gff = PROKKA.out.gff       // [ val(meta), path(gff) ]
    core_aln   = ch_core_aln          // path(alignment)
    tree       = ch_tree               // [ val(meta), path(treefile) ] or path(treefile)
    snp_dists  = SNPDISTS.out.tsv      // [ val(meta), path(tsv) ]
    versions   = ch_versions           // channel: [ path(versions.yml) ]
}
