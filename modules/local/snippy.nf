process SNIPPY {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/snippy/${meta.id}", mode: 'link'
    errorStrategy 'ignore'
    container 'docker.io/staphb/snippy:4.6.0'

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("${meta.id}"), emit: results

    script:
    """
    snippy --cpus ${task.cpus} --ram ${task.memory.toGiga()} --outdir ${meta.id} --ref ${reference} --R1 ${reads[0]} --R2 ${reads[1]} --mincov ${params.snippy_mincov} --minqual ${params.snippy_minqual} --cleanup
    """
}

process SNIPPY_CORE {
    label 'process_high'
    publishDir "${params.outdir}/snippy_core", mode: 'copy'
    container 'docker.io/staphb/snippy:4.6.0'

    input:
    path snippy_dirs
    path reference

    output:
    path "core.aln", emit: aln

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
}
