process PUBMLST_FETCH {
    tag "${params.cc ? "CC${params.cc}" : "ST${params.st}"}"
    label 'process_single'
    container 'docker.io/python:3.11'

    input:
    path fetch_script

    output:
    path "accessions.txt", emit: accessions
    path "metadata.tsv",   emit: metadata

    script:
    def cc_arg  = params.cc               ? "--cc ${params.cc}"                     : ""
    def st_arg  = params.st               ? "--st '${params.st}'"                   : ""
    def src_arg = params.fetch_source     ? "--source ${params.fetch_source}"       : "--source all"
    def ctr_arg = params.fetch_country    ? "--country '${params.fetch_country}'"   : ""
    def con_arg = params.fetch_continent  ? "--continent ${params.fetch_continent}" : ""
    def hst_arg = params.fetch_host       ? "--host '${params.fetch_host}'"         : ""
    def yf_arg  = params.fetch_year_from  ? "--year-from ${params.fetch_year_from}" : ""
    def yt_arg  = params.fetch_year_to    ? "--year-to ${params.fetch_year_to}"     : ""
    def max_arg = params.fetch_max        ? "--max-downloads ${params.fetch_max}"   : ""
    """
    pip install -q --user requests oauthlib requests-oauthlib
    python ${fetch_script} \\
        ${cc_arg} ${st_arg} ${src_arg} \\
        ${ctr_arg} ${con_arg} ${hst_arg} \\
        ${yf_arg} ${yt_arg} ${max_arg} \\
        --sra-only \\
        --out metadata.tsv \\
        --accessions accessions.txt
    """
}
