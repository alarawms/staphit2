process BEAST2_PREP {
    label 'process_single'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    input:
    path snps               // recombination-masked SNP alignment (GUBBINS)
    path constant_sites     // invariant A C G T counts (GUBBINS)
    path metadata           // metadata TSV (may be NO_METADATA sentinel)

    output:
    path "beast_run.xml", emit: xml
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def meta_arg = metadata.name != 'NO_METADATA' ? "--metadata ${metadata}" : ""
    def chain    = params.beast_chain_length ?: 10000000
    def log_ev   = params.beast_log_every    ?: 1000
    """
    beast2_prep.py \\
        --alignment ${snps} \\
        --constant-sites ${constant_sites} \\
        ${meta_arg} \\
        --output beast_run.xml \\
        --chain-length ${chain} \\
        --log-every ${log_ev}
    """

    stub:
    """
    touch beast_run.xml
    """
}

process BEAST2 {
    tag "beast2"
    label 'process_high'
    container 'quay.io/biocontainers/beast2:2.7.7--he65b2d3_0'

    input:
    path xml                // from BEAST2_PREP

    output:
    // not optional: BEAST exits 0 even when the XML fails to parse, so a missing file must fail the task
    path "beast_mcc.tree",   emit: mcc_tree
    path "beast_run.log",    emit: log
    path "beast_run.trees",  emit: trees
    tuple val("${task.process}"), val('beast2'), eval("beast -version 2>&1 | grep -o 'v[0-9.]*' | head -1 | tr -d v"), topic: versions, emit: versions_beast2

    script:
    """
    # Java takes user.home from /etc/passwd, which has no entry for the container's uid under
    # Docker -u; BEAST then cannot set up its package directory and fails to load any class
    export JAVA_TOOL_OPTIONS="-Duser.home=\$PWD/.beast_home"
    beast -overwrite -threads ${task.cpus} ${xml}
    treeannotator -burnin 10 -height mean beast_run.trees beast_mcc.tree
    """

    stub:
    """
    touch beast_mcc.tree beast_run.log beast_run.trees
    """
}
