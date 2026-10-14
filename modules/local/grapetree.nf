process GRAPETREE {
    tag "cgmlst"
    label 'process_low'
    container 'quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0'

    input:
    path profile

    output:
    path "cgmlst_mst.nwk", emit: tree
    tuple val("${task.process}"), val('grapetree'), eval("python -c 'import importlib.metadata as m; print(m.version(\"grapetree\"))'"), topic: versions, emit: versions_grapetree

    script:
    """
    # MSTreeV2: minimum spanning tree that tolerates missing loci (Zhou et al. 2018)
    grapetree -p ${profile} -m MSTreeV2 -n ${task.cpus} > cgmlst_mst.nwk
    """

    stub:
    """
    echo '(A:0,B:1);' > cgmlst_mst.nwk
    """
}

process MST_PLOT {
    tag "cgmlst"
    label 'process_low'
    container 'ghcr.io/alarawms/staphit2-r:1.0.0'

    input:
    path tree
    path summary
    path clusters

    output:
    path "cgmlst_mst.pdf",         emit: pdf, optional: true
    path "grapetree_metadata.tsv", emit: metadata
    tuple val("${task.process}"), val('ggtree'), eval("Rscript -e 'cat(as.character(packageVersion(\"ggtree\")))'"), topic: versions, emit: versions_ggtree

    script:
    """
    plot_mst.R ${tree} ${summary} ${clusters} cgmlst_mst
    """

    stub:
    """
    touch cgmlst_mst.pdf grapetree_metadata.tsv
    """
}
