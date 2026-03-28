process MOB_RECON {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'
    container 'docker.io/staphb/mob-suite:latest'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}_mob/contig_report.txt"), emit: contig_report
    tuple val(meta), path("${meta.id}_mob/mobtyper_results.txt"), emit: mobtyper
    tuple val(meta), path("${meta.id}_mob"), emit: full_output

    script:
    """
    # Initialize MOB-suite databases if needed
    mob_init 2>/dev/null || true
    mob_recon -i ${assembly} -o ${meta.id}_mob -n ${task.cpus} -s ${meta.id} --force
    if [ ! -f ${meta.id}_mob/mobtyper_results.txt ]; then
        echo -e "sample_id\\tnum_contigs\\tsize\\tgc\\tmd5\\trep_type(s)\\trep_type_accession(s)\\trelaxase_type(s)\\trelaxase_type_accession(s)\\tmpf_type\\tmpf_type_accession(s)\\torit_type(s)\\torit_accession(s)\\tpredicted_mobility\\tmash_nearest_neighbor\\tmash_neighbor_distance\\tmash_neighbor_identification\\tprimary_cluster_id\\tsecondary_cluster_id\\tpredicted_host_range_overall_rank\\tpredicted_host_range_overall_name\\tobserved_host_range_ncbi_rank\\tobserved_host_range_ncbi_name\\treported_host_range_lit_rank\\treported_host_range_lit_name\\tassociated_pmid(s)" > ${meta.id}_mob/mobtyper_results.txt
    fi
    """
}
