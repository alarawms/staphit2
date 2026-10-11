process MIXED_STRAINS {
    tag "$meta.id"
    label 'process_medium'
    container 'docker.io/staphb/snippy:4.6.0'   // minimap2, samtools, bcftools, python3

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.mixed.tsv"), emit: tsv
    tuple val("${task.process}"), val('minimap2'), eval("minimap2 --version"), topic: versions, emit: versions_minimap2
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -1 | sed 's/bcftools //'"), topic: versions, emit: versions_bcftools

    script:
    // Reads mapped back to the sample's own assembly: a single strain agrees with itself,
    // a mix of two S. aureus strains shows a second base wherever they differ (B2)
    def unpack = assembly.name.endsWith('.gz') ? "gzip -cd ${assembly} > asm.fa" : "ln -s ${assembly} asm.fa"
    """
    ${unpack}
    samtools faidx asm.fa
    minimap2 -t ${task.cpus} -ax sr asm.fa ${reads[0]} ${reads[1]} | samtools sort -@ ${task.cpus} -o aln.bam -
    bcftools mpileup -f asm.fa -a AD -q 30 -Q 20 -d 1000 -Ou aln.bam \\
        | bcftools query -f '[%AD]\\n' \\
        | staphit-mixed ${meta.id} > ${meta.id}.mixed.tsv
    rm aln.bam
    """

    stub:
    """
    printf "sample_id\\tsites_assessed\\tmixed_sites\\tminor_fraction\\n${meta.id}\\t0\\t0\\t0.000\\n" > ${meta.id}.mixed.tsv
    """
}
