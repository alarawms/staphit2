process AGR_TYPING {
    tag "$meta.id"
    label 'process_low'
    errorStrategy { task.exitStatus in [125,137] ? 'retry' : 'terminate' }
    maxRetries 2
    container 'docker.io/alarawms/staph_agr_typer:latest'
    containerOptions '--entrypoint ""'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_agr.json"), emit: report

    script:
    """
    staph_agr_typer run --fasta ${assembly} -o agr_out 2>/dev/null && \
        cp agr_out/*.json ${meta.id}_agr.json 2>/dev/null || \
        echo '{"agr_group": "ND", "confidence": 0.0}' > ${meta.id}_agr.json
    """
}
