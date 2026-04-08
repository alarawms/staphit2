process PANAROO {
    label 'process_high'
    publishDir "${params.outdir}/panaroo", mode: 'copy'
    errorStrategy 'ignore'
    container 'docker.io/staphb/panaroo:1.3.4'

    input:
    path gffs

    output:
    path "core_gene_alignment.aln", optional: true, emit: aln
    path "gene_presence_absence.csv", optional: true
    path "pan_genome_reference.fa", optional: true

    script:
    """
    for f in *.gff; do
        cds_count=\$(grep -c "CDS" "\$f" 2>/dev/null || echo 0)
        if [ ! -s "\$f" ] || [ "\$cds_count" -lt 500 ]; then
            rm -f "\$f"
        fi
    done
    count=\$(ls *.gff 2>/dev/null | wc -l)
    if [ "\$count" -lt 2 ]; then
        echo "Not enough valid samples for Panaroo. Skipping."
        exit 0
    fi
    panaroo -i *.gff -o . \
        --clean-mode ${params.panaroo_clean} \
        --remove-invalid-genes \
        -a core \
        --core_threshold ${params.panaroo_threshold} \
        --aligner ${params.panaroo_aligner} \
        -c ${params.panaroo_identity} \
        -f ${params.panaroo_family_threshold} \
        ${params.panaroo_entropy_filter ? "--core_entropy_filter ${params.panaroo_entropy_filter}" : ''} \
        -t ${task.cpus}
    """
}
