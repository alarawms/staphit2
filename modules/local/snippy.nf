process SNIPPY {
    tag "$meta.id"
    label 'process_medium'
    container 'docker.io/staphb/snippy:4.6.0'

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    tuple val("${task.process}"), val('snippy'), eval("snippy --version 2>&1 | sed 's/^snippy //'"), topic: versions, emit: versions_snippy

    script:
    // Long-only samples arrive as their assembly (Snippy cannot align raw ONT reads)
    def input = meta.mode == 'long' ? "--ctgs ${reads}" : "--R1 ${reads[0]} --R2 ${reads[1]}"
    """
    snippy --cpus ${task.cpus} --ram ${task.memory.toGiga()} --outdir ${meta.id} --ref ${reference} ${input} --mincov ${params.snippy_mincov} --minqual ${params.snippy_minqual} --cleanup
    # Remove broken symlinks left by --cleanup so publishDir doesn't fail on them
    find ${meta.id} -xtype l -delete
    """

    stub:
    """
    mkdir ${meta.id}
    """
}

process SNIPPY_CORE {
    label 'process_high'
    container 'docker.io/staphb/snippy:4.6.0'

    input:
    path snippy_dirs
    path reference

    output:
    path "core.aln", emit: aln
    tuple val("${task.process}"), val('snippy'), eval("snippy --version 2>&1 | sed 's/^snippy //'"), topic: versions, emit: versions_snippy

    script:
    """
    # Filter out samples with >80% unaligned genome (bad data or non-S.aureus)
    GOOD_DIRS=""
    REMOVED=0
    for d in */; do
        d=\${d%/}
        if [ -f "\$d/snps.txt" ]; then
            UNALIGNED=\$(grep "^Unaligned" "\$d/snps.txt" 2>/dev/null | awk '{print \$2}' || echo "0")
            TOTAL=\$(grep "^Length" "\$d/snps.txt" 2>/dev/null | awk '{print \$2}' || echo "2800000")
            if [ "\$TOTAL" -gt 0 ] 2>/dev/null; then
                PCT=\$((UNALIGNED * 100 / TOTAL))
                if [ "\$PCT" -lt 80 ]; then
                    GOOD_DIRS="\$GOOD_DIRS \$d"
                else
                    echo "FILTERED: \$d (\$PCT% unaligned)" >&2
                    REMOVED=\$((REMOVED + 1))
                fi
            else
                GOOD_DIRS="\$GOOD_DIRS \$d"
            fi
        else
            GOOD_DIRS="\$GOOD_DIRS \$d"
        fi
    done

    echo "Filtered \$REMOVED samples with >80% unaligned genome" >&2
    echo "Running snippy-core on \$(echo \$GOOD_DIRS | wc -w) samples" >&2

    if [ \$(echo \$GOOD_DIRS | wc -w) -lt 3 ]; then
        echo "ERROR: Too few samples after filtering" >&2
        exit 1
    fi

    snippy-core --ref ${reference} \$GOOD_DIRS
    """

    stub:
    """
    printf ">Reference\\nACGTACGT\\n>s1\\nACGTACGA\\n>s2\\nACTTACGT\\n" > core.aln
    """
}
