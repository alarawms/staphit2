process FASTTREE {
    label 'process_medium'
    publishDir "${params.outdir}/fasttree", mode: 'copy'
    container 'staphb/lyveset:2.0.1'

    input:
    path alignment

    output:
    path "*.treefile", optional: true

    script:
    """
    seq_count=\$(grep -c "^>" $alignment)
    if [ "\$seq_count" -lt 3 ]; then
        echo "Too few sequences. Skipping."
        exit 0
    fi
    FastTreeMP -gtr -nt $alignment > core.treefile
    """
}
