process BEAST2 {
    tag "beast2"
    label 'process_high'
    container 'quay.io/biocontainers/beast2:2.7.7--he65b2d3_0'

    input:
    path alignment       // core.aln from snippy-core
    path metadata        // metadata TSV (may be NO_METADATA sentinel)
    path prep_script     // bin/beast2_prep.py staged in

    output:
    path "beast_mcc.tree",   emit: mcc_tree,  optional: true
    path "beast_*.log",      emit: log,        optional: true
    path "beast_*.trees",    emit: trees,      optional: true
    tuple val("${task.process}"), val('beast2'), eval("beast -version 2>&1 | grep -o 'v[0-9.]*' | head -1 | tr -d v"), topic: versions, emit: versions_beast2

    script:
    def meta_arg = metadata.name != 'NO_METADATA' ? "--metadata ${metadata}" : ""
    def chain    = params.beast_chain_length ?: 10000000
    def log_ev   = params.beast_log_every    ?: 1000
    """
    python ${prep_script} \\
        --alignment ${alignment} \\
        ${meta_arg} \\
        --output beast_run.xml \\
        --chain-length ${chain} \\
        --log-every ${log_ev}

    beast -overwrite beast_run.xml

    treeannotator -burnin 10 -heights mean beast_run.trees beast_mcc.tree
    """

    stub:
    """
    touch beast_mcc.tree
    """
}
