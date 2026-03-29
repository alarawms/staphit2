process PLASMID_SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}/plasmids", mode: 'copy'
    container 'docker.io/python:3.9'

    input:
    path mob_dirs

    output:
    path "plasmid_summary.tsv", emit: summary

    script:
    '''
    python3 << 'PYEOF'
import csv, os

out_rows = []
for d in os.listdir('.'):
    cr_path = os.path.join(d, 'contig_report.txt')
    mt_path = os.path.join(d, 'mobtyper_results.txt')
    if not os.path.isdir(d) or not os.path.exists(cr_path):
        continue
    plasmid_clusters = set()
    sample_id = ''
    with open(cr_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if not sample_id:
                sample_id = row.get('sample_id', '')
            if row.get('molecule_type') == 'plasmid':
                cid = row.get('primary_cluster_id', '')
                if cid and cid != '-':
                    plasmid_clusters.add(cid)
    mobility_parts = []
    if os.path.exists(mt_path):
        with open(mt_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                cid = row.get('primary_cluster_id', '-')
                mob = row.get('predicted_mobility', '-')
                mobility_parts.append(f"{cid}:{mob}")
    out_rows.append({
        'sample_id': sample_id,
        'plasmid_count': str(len(plasmid_clusters)),
        'plasmid_mobility': ';'.join(mobility_parts) if mobility_parts else '',
    })

with open('plasmid_summary.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['sample_id', 'plasmid_count', 'plasmid_mobility'], delimiter='\t')
    w.writeheader()
    w.writerows(sorted(out_rows, key=lambda x: x['sample_id']))
print(f"Summarized {len(out_rows)} samples")
PYEOF
    '''
}
