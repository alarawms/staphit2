process AGR_TYPING {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/alarawms/staph_agr_typer:latest'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_agr.json"), emit: report

    script:
    """
    agr_typer ${assembly} > ${meta.id}_agr.json 2>/dev/null || echo '{"agr_group": "ND", "confidence": 0.0}' > ${meta.id}_agr.json
    """
}
