process FASTTREE {
    label 'process_medium'
    container 'quay.io/biocontainers/fasttree:2.2.0--h7b50bb2_1'

    input:
    path alignment

    output:
    path "*.treefile", optional: true, emit: tree
    tuple val("${task.process}"), val('fasttree'), eval("FastTree 2>&1 | head -1 | cut -d' ' -f3"), topic: versions, emit: versions_fasttree

    script:
    """
    seq_count=\$(grep -c "^>" $alignment)
    if [ "\$seq_count" -lt 3 ]; then
        echo "Too few sequences. Skipping."
        exit 0
    fi
    FastTreeMP -gtr -nt $alignment > core.treefile
    """

    stub:
    """
    echo "(a,b,c);" > core.treefile
    """
}
