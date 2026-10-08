process AGR_TYPING {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/alarawms/staph_agr_typer@sha256:49f5a4177b8f91656a4a900d8d12d43148da40874dd8b604aaaeb9db7f533e6d'
    containerOptions { workflow.containerEngine in ['docker', 'podman'] ? '--entrypoint ""' : '' }   // Docker-only flag; Apptainer rejects it

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_agr.json"), emit: report
    tuple val("${task.process}"), val('staph_agr_typer'), val(task.container), topic: versions, emit: versions_staph_agr_typer

    script:
    """
    staph_agr_typer run --fasta ${assembly} -o agr_out 2>/dev/null && \
        cp agr_out/*.json ${meta.id}_agr.json 2>/dev/null || \
        echo '{"agr_group": "ND", "confidence": 0.0}' > ${meta.id}_agr.json
    """

    stub:
    """
    echo '{"agr_group": "ND", "confidence": 0.0}' > ${meta.id}_agr.json
    """
}
