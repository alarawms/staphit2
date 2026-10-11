process GUBBINS {
    label 'process_high'
    container 'quay.io/biocontainers/gubbins:3.4.3--py310hfc0ef84_1'
    // Gubbins' numba workers need a large /dev/shm (SIGBUS, exit 135, with Docker's 64 MB default)
    containerOptions { workflow.containerEngine in ['docker', 'podman'] ? '--shm-size=16g' : '' }

    input:
    path full_aln           // snippy-core core.full.aln (Reference + samples)

    output:
    path "snps.masked.fasta",                       emit: snps       // variable sites after masking recombination
    path "constant_sites.txt",                      emit: constant   // invariant A C G T counts (BEAST2 ascertainment)
    path "gubbins.recombination_predictions.gff",   emit: gff
    path "gubbins.final_tree.tre",                  emit: tree
    path "gubbins.per_branch_statistics.csv",       emit: stats
    tuple val("${task.process}"), val('gubbins'), eval('run_gubbins.py --version'), topic: versions, emit: versions_gubbins

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache
    run_gubbins.py --prefix gubbins --threads ${task.cpus} ${task.ext.args ?: ''} ${full_aln}
    mask_gubbins_aln.py --aln ${full_aln} --gff gubbins.recombination_predictions.gff --out masked.full.aln
    beast_alignment.py masked.full.aln snps.masked.fasta constant_sites.txt --drop Reference
    """

    stub:
    """
    printf ">s1\\nACGA\\n>s2\\nACTT\\n" > snps.masked.fasta
    echo "700000 700000 700000 700000" > constant_sites.txt
    touch gubbins.recombination_predictions.gff gubbins.final_tree.tre gubbins.per_branch_statistics.csv
    """
}
