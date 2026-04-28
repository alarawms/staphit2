process BEAST2 {
    tag "beast2"
    label 'process_high'
    publishDir "${params.outdir}/beast2", mode: 'copy'
    container 'quay.io/biocontainers/beast2:2.7.7--hdfd78af_1'

    input:
    path alignment       // core.aln from snippy-core
    path metadata        // metadata TSV (may be NO_METADATA sentinel)
    path prep_script     // bin/beast2_prep.py staged in

    output:
    path "beast_mcc.tree",   emit: mcc_tree,  optional: true
    path "beast_*.log",      emit: log,        optional: true
    path "beast_*.trees",    emit: trees,      optional: true

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
}
