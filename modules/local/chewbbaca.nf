process CHEWBBACA_PREP {
    label 'process_medium'
    publishDir "${params.outdir}/cgmlst/schema", mode: 'copy'
    container 'ummidock/chewbbaca:latest'

    input:
    path schema_fasta_dir

    output:
    path "prepared_schema", emit: schema

    script:
    """
    chewBBACA.py PrepExternalSchema -g ${schema_fasta_dir} -o prepared_schema --cpu ${task.cpus}
    """
}

process CHEWBBACA_ALLELE {
    tag "$meta.id"
    label 'process_medium'
    container 'ummidock/chewbbaca:latest'

    input:
    tuple val(meta), path(assembly)
    path schema

    output:
    tuple val(meta), path("results/results_alleles.tsv"), emit: profile
    tuple val(meta), path("results/results_statistics.tsv"), emit: stats

    script:
    """
    mkdir -p input_dir
    cp ${assembly} input_dir/${meta.id}.fasta
    chewBBACA.py AlleleCall -i input_dir -g ${schema} -o results --no-inferred --cpu ${task.cpus} --mode 4
    if ls results/results_*/results_alleles.tsv 1>/dev/null 2>&1; then
        mv results/results_*/* results/
    fi
    """
}

process CHEWBBACA_JOIN {
    label 'process_low'
    publishDir "${params.outdir}/cgmlst", mode: 'copy'
    container 'ummidock/chewbbaca:latest'

    input:
    path profiles

    output:
    path "cgmlst_profiles.tsv", emit: profiles
    path "cgmlst_stats.tsv", emit: stats

    script:
    """
    mkdir -p profile_files
    i=0
    for f in ${profiles}; do
        cp "\$f" "profile_files/profile_\${i}.tsv"
        i=\$((i+1))
    done
    profile_list=\$(ls profile_files/*.tsv | tr '\\n' ' ')
    chewBBACA.py JoinProfiles -p \$profile_list -o cgmlst_profiles.tsv
    echo -e "sample_id\\tcalled_loci\\tmissing_loci\\ttotal_loci" > cgmlst_stats.tsv
    python3 -c "
import csv
with open('cgmlst_profiles.tsv') as f:
    reader = csv.DictReader(f, delimiter='\\t')
    loci = [c for c in reader.fieldnames if c != 'FILE']
    total = len(loci)
    for row in reader:
        sid = row['FILE'].replace('.fasta', '')
        called = sum(1 for l in loci if row[l] not in ('0', 'LNF', 'PLOT3', 'PLOT5', 'NIPH', 'NIPHEM', 'ASM', 'ALM', 'LOTSC'))
        print(f'{sid}\\t{called}\\t{total - called}\\t{total}')
" >> cgmlst_stats.tsv
    """
}

process CGMLST_DISTS {
    label 'process_low'
    publishDir "${params.outdir}/cgmlst", mode: 'copy'
    container 'python:3.9-slim'

    input:
    path profiles

    output:
    path "cgmlst_distances.tsv"

    script:
    """
    python3 ${projectDir}/bin/staphit-cgmlst-stats ${profiles} cgmlst_distances.tsv 2>/dev/null || python3 -c "
import csv
with open('${profiles}') as f:
    reader = csv.DictReader(f, delimiter='\\t')
    loci = [c for c in reader.fieldnames if c != 'FILE']
    samples = {}
    for row in reader:
        sid = row['FILE'].replace('.fasta', '')
        samples[sid] = [row[l] for l in loci]
sample_ids = sorted(samples.keys())
missing = {'0','LNF','PLOT3','PLOT5','NIPH','NIPHEM','ASM','ALM','LOTSC'}
with open('cgmlst_distances.tsv', 'w') as out:
    out.write('\\t' + '\\t'.join(sample_ids) + '\\n')
    for s1 in sample_ids:
        dists = []
        for s2 in sample_ids:
            if s1 == s2: dists.append('0')
            else:
                diff = sum(1 for a, b in zip(samples[s1], samples[s2]) if a not in missing and b not in missing and a != b)
                dists.append(str(diff))
        out.write(s1 + '\\t' + '\\t'.join(dists) + '\\n')
"
    """
}
