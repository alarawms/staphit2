process OUTBREAK_CLUSTER {
    label 'process_low'
    publishDir "${params.outdir}/clusters", mode: 'copy'
    container 'docker.io/python:3.9'

    input:
    path snp_dists
    path cgmlst_dists
    path summary_tsv

    output:
    path "clusters.tsv", emit: clusters
    path "cluster_report.json", emit: report
    path "transmission_pairs.tsv", emit: pairs

    script:
    def snp_flag = snp_dists.name != 'NO_SNP_DISTS' ? "--snp-dists ${snp_dists}" : ''
    def cgmlst_flag = cgmlst_dists.name != 'NO_CGMLST_DISTS' ? "--cgmlst-dists ${cgmlst_dists}" : ''
    def summary_flag = summary_tsv.name != 'NO_SUMMARY' ? "--summary ${summary_tsv}" : ''
    """
    python3 ${projectDir}/bin/staphit-cluster ${snp_flag} ${cgmlst_flag} ${summary_flag} --snp-tiers "${params.cluster_snp_tiers}" --cgmlst-tiers "${params.cluster_cgmlst_tiers}" -o .
    """
}
