process PLOT_TREE {
    tag "phylogeny"
    label 'process_medium'
    container 'quay.io/biocontainers/bioconductor-ggtreeextra:1.20.1--r45hdfd78af_0'
    errorStrategy 'ignore'

    publishDir "${params.outdir}/phylogeny", mode: 'copy', pattern: '*.{pdf,svg}'

    input:
    path treefile
    path summary
    path clusters

    output:
    path "*.pdf", emit: pdf, optional: true
    path "*.svg", emit: svg, optional: true

    script:
    def run_name = params.outdir.split('/').last()
    """
    # v3: sccmec_group from summary; v2: harmonized year/host/source rings, explicit ring width, no svglite/readr
    #     (bump this line after editing bin/plot_tree.R so -resume re-plots)
    # Install missing CRAN packages into local lib (cached by work dir on -resume)
    mkdir -p rlibs
    Rscript -e "
        .libPaths(c('rlibs', .libPaths()))
        pkgs <- c('ggnewscale', 'RColorBrewer', 'dplyr', 'stringr')
        new  <- pkgs[!pkgs %in% installed.packages()[,'Package']]
        if (length(new)) install.packages(new, repos='https://cloud.r-project.org', lib='rlibs', quiet=TRUE)
    "

    # Recreate the directory layout plot_tree.R expects
    mkdir -p iqtree summary clusters
    ln -sf \$(realpath ${treefile}) iqtree/core.treefile
    ln -sf \$(realpath ${summary})  summary/combined_summary.tsv
    ln -sf \$(realpath ${clusters}) clusters/clusters.tsv

    R_LIBS=\$PWD/rlibs Rscript ${projectDir}/bin/plot_tree.R . ${run_name}
    """
}
