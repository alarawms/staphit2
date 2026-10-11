process TRUNCATED_GENES {
    tag "$meta.id"
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl5321hdfd78af_5' :
        'quay.io/biocontainers/prokka:1.14.6--pl5321hdfd78af_5' }"   // blastp + python3

    input:
    tuple val(meta), path(faa)
    path reference_gbk

    output:
    tuple val(meta), path("${meta.id}.truncated.tsv"), emit: tsv
    tuple val("${task.process}"), val('blastp'), eval("blastp -version | head -1 | sed 's/blastp: //'"), topic: versions, emit: versions_blastp

    script:
    // Proteins shorter than 90% of the same NCTC 8325 protein: frameshifts from ONT indel errors (B6)
    """
    staphit-truncated ref ${reference_gbk} > ref.faa
    makeblastdb -in ref.faa -dbtype prot -out ref > /dev/null
    blastp -query ${faa} -db ref -outfmt '6 qseqid sseqid pident qlen slen' -max_target_seqs 1 -evalue 1e-10 -num_threads ${task.cpus} > hits.tsv
    staphit-truncated count ${meta.id} hits.tsv > ${meta.id}.truncated.tsv
    """

    stub:
    """
    printf "sample_id\\tproteins_with_hit\\ttruncated\\ttruncated_pct\\n${meta.id}\\t0\\t0\\t0.0\\n" > ${meta.id}.truncated.tsv
    """
}
