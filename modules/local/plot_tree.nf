process PLOT_TREE {
    tag "phylogeny"
    label 'process_medium'
    container 'ghcr.io/alarawms/staphit2-r:1.0.0'   // all R packages pre-installed (docker/staphit2-r)

    input:
    path treefile
    path summary
    path clusters

    output:
    path "*.pdf", emit: pdf, optional: true
    path "*.svg", emit: svg, optional: true
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(R.version[[\"major\"]], R.version[[\"minor\"]], sep=\".\")'"), topic: versions, emit: versions_r_base
    tuple val("${task.process}"), val('ggtree'), eval("Rscript -e 'cat(as.character(packageVersion(\"ggtree\")))'"), topic: versions, emit: versions_ggtree

    script:
    def run_name = params.outdir.split('/').last()
    """
    # v3: sccmec_group from summary; v2: harmonized year/host/source rings, explicit ring width, no svglite/readr
    #     (bump this line after editing bin/plot_tree.R so -resume re-plots)

    # Recreate the directory layout plot_tree.R expects
    mkdir -p iqtree summary clusters
    ln -sf \$(realpath ${treefile}) iqtree/core.treefile
    ln -sf \$(realpath ${summary})  summary/combined_summary.tsv
    ln -sf \$(realpath ${clusters}) clusters/clusters.tsv

    plot_tree.R . ${run_name}
    """

    stub:
    """
    touch stub_tree.pdf
    """
}
