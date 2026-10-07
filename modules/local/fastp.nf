process FASTP {
    tag "$meta.id"
    label 'process_medium'
    publishDir path: { "${params.outdir}/fastp" }, mode: 'copy', pattern: '*.{json,html}'
    container 'quay.io/biocontainers/fastp:0.24.0--heae3180_1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed_R{1,2}.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"),             emit: json
    tuple val(meta), path("${meta.id}.fastp.html"),             emit: html

    script:
    // Single process (no internal pipes, unlike Trim Galore + cutadapt).
    // Mirrors the old Trim Galore settings: adapter auto-detect, 3' Q20 trimming, min length 20.
    """
    fastp \\
        --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --detect_adapter_for_pe \\
        --cut_tail --cut_tail_mean_quality 20 \\
        --length_required 20 \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
}
