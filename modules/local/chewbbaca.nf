// cgMLST: PubMLST S. aureus scheme 20 (1,716 loci) -> chewBBACA allele calls -> allele distances.
// Reuse the prepared schema with --cgmlst_schema <outdir>/databases/cgmlst_schema/prepared_schema.

process CGMLST_SCHEMA {
    label 'process_low'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    output:
    path "schema_fasta", emit: fasta

    script:
    """
    staphit-cgmlst-schema -o schema_fasta --db ${params.cgmlst_pubmlst_db} --scheme ${params.cgmlst_pubmlst_scheme}
    """

    stub:
    """
    mkdir schema_fasta
    """
}

process CHEWBBACA_PREP {
    label 'process_medium'
    container 'docker.io/ummidock/chewbbaca:v3.3.10'

    input:
    path schema_fasta_dir

    output:
    path "prepared_schema", emit: schema
    tuple val("${task.process}"), val('chewbbaca'), eval("chewBBACA.py --version | sed 's/.*: //'"), topic: versions, emit: versions_chewbbaca

    script:
    """
    chewBBACA.py PrepExternalSchema -g ${schema_fasta_dir} -o prepared_schema --cpu ${task.cpus}
    """

    stub:
    """
    mkdir prepared_schema
    """
}

// All samples in one AlleleCall run: with --no-inferred, new (INF-) allele numbers are only
// consistent within a run, so per-sample runs would give unrelated new alleles the same id.
process CHEWBBACA_ALLELE {
    label 'process_high'
    container 'docker.io/ummidock/chewbbaca:v3.3.10'

    input:
    tuple val(ids), path(assemblies, stageAs: 'asm/*')
    path schema

    output:
    path "results_alleles.tsv",    emit: profiles
    path "results_statistics.tsv", emit: stats
    tuple val("${task.process}"), val('chewbbaca'), eval("chewBBACA.py --version | sed 's/.*: //'"), topic: versions, emit: versions_chewbbaca

    script:
    // chewBBACA names each profile after its file: give every assembly <sample_id>.fasta
    def files = assemblies instanceof List ? assemblies : [assemblies]   // one sample: a single path
    def links = [ids, files].transpose().collect { id, f ->
        f.name.endsWith('.gz') ? "gzip -cd ${f} > input_dir/${id}.fasta" : "ln -s ../${f} input_dir/${id}.fasta"
    }.join('\n    ')
    """
    mkdir input_dir
    ${links}
    chewBBACA.py AlleleCall -i input_dir -g ${schema} -o results --no-inferred --cpu ${task.cpus} --mode 4
    mv results/results_alleles.tsv results/results_statistics.tsv .
    """

    stub:
    """
    touch results_alleles.tsv results_statistics.tsv
    """
}

process CGMLST_DISTS {
    label 'process_low'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    input:
    path profiles

    output:
    path "cgmlst_distances.tsv", emit: dists
    path "cgmlst_stats.tsv",     emit: stats
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    """
    staphit-cgmlst-dists ${profiles} cgmlst_distances.tsv cgmlst_stats.tsv
    """

    stub:
    """
    touch cgmlst_distances.tsv
    printf "sample_id\\tcalled_loci\\tmissing_loci\\ttotal_loci\\n" > cgmlst_stats.tsv
    """
}
