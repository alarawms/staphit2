#!/usr/bin/env python3
"""Write iTOL annotation files for a staphit2 tree.

Inputs : core.treefile, combined_summary.tsv, [clusters.tsv]
Output : <outdir>/  numbered iTOL dataset files + the tree itself. Upload the tree to
         https://itol.embl.de and drag the .txt files onto it.

Labels and colours for year / host / isolation source come from the summary's
harmonized columns and bin/harmonize_metadata.py — the same as the PDF tree plot.

Usage: itol_annotations.py <treefile> <combined_summary.tsv> <clusters.tsv|NONE> <outdir>
"""
import csv
import os
import re
import shutil
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harmonize_metadata as hm  # noqa: E402

QUAL = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d',
        '#1f78b4', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6', '#6a3d9a', '#b15928', '#33a02c']
AGR_COLORS = {'gp1': '#1b9e77', 'gp2': '#d95f02', 'gp3': '#7570b3', 'gp4': '#e7298a'}


def tree_tips(newick):
    """Leaf names from Newick text (internal node labels such as bootstraps are skipped)."""
    return re.findall(r'[(,]\s*([^(),:;\s]+)\s*(?=[:,);])', newick)


def top_n(values, n, other='Other'):
    keep = {v for v, _ in Counter(v for v in values if v not in ('', '-', hm.UNKNOWN)).most_common(n)}
    return [v if v in keep else (hm.UNKNOWN if v in ('', '-', hm.UNKNOWN) else other) for v in values]


def colorstrip(label, ids, cats, colors, order=None):
    present = [c for c in (order or sorted(set(cats), key=lambda c: (c in ('Other', hm.UNKNOWN), c)))
               if c in set(cats)]
    lines = ['DATASET_COLORSTRIP', 'SEPARATOR TAB', f'DATASET_LABEL\t{label}', 'COLOR\t#333333',
             f'LEGEND_TITLE\t{label}', 'LEGEND_SHAPES\t' + '\t'.join('1' for _ in present),
             'LEGEND_COLORS\t' + '\t'.join(colors[c] for c in present),
             'LEGEND_LABELS\t' + '\t'.join(present),
             'STRIP_WIDTH\t28', 'MARGIN\t4', 'BORDER_WIDTH\t0', 'SHOW_INTERNAL\t0', 'DATA']
    lines += [f'{i}\t{colors[c]}\t{c}' for i, c in zip(ids, cats)]
    return lines


def palette_for(cats, fixed=None):
    fixed = dict(fixed or {})
    free = [c for c in sorted(set(cats)) if c not in fixed and c not in ('Other', hm.UNKNOWN)]
    fixed.update({c: QUAL[i % len(QUAL)] for i, c in enumerate(free)})
    fixed.setdefault('Other', '#bdbdbd')
    fixed.setdefault(hm.UNKNOWN, '#f0f0f0')
    for c in set(cats) - set(fixed):
        fixed[c] = hm.OTHER_COLOR
    return fixed


def build(tips, rows, clusters):
    by_id = {r['sample_id']: r for r in rows}
    ids = [t for t in tips if t in by_id]
    R = [by_id[i] for i in ids]
    files = {}

    # Tip labels: ID | ST | year — readable without opening any dataset
    files['01_tip_labels.txt'] = ['LABELS', 'SEPARATOR TAB', 'DATA'] + [
        f"{i}\t{i} | ST{r['mlst_st'] if r['mlst_st'] not in ('', '-') else '?'} | {r.get('year') or hm.UNKNOWN}"
        for i, r in zip(ids, R)]

    st = top_n([r['mlst_st'] for r in R], 15)
    st = [f'ST{s}' if s not in ('Other', hm.UNKNOWN) else s for s in st]
    files['02_ST_strip.txt'] = colorstrip('ST (top 15)', ids, st, palette_for(st))

    scc = [r.get('sccmec_group') or hm.sccmec_group(r.get('sccmec_type')) for r in R]
    files['03_SCCmec_strip.txt'] = colorstrip('SCCmec', ids, scc, palette_for(scc, hm.SCCMEC_COLORS),
                                              order=list(hm.SCCMEC_COLORS))

    agr = [r['agr_group'] if r['agr_group'] not in ('', '-') else hm.UNKNOWN for r in R]
    files['04_agr_strip.txt'] = colorstrip('agr group', ids, agr, palette_for(agr, AGR_COLORS))

    yrs = [r.get('year') or hm.UNKNOWN for r in R]
    yc = hm.year_colors(yrs)
    files['05_year_strip.txt'] = colorstrip('Collection year', ids, yrs, yc,
                                            order=sorted(c for c in yc if c != hm.UNKNOWN) + [hm.UNKNOWN])

    host = [r.get('host_group') or hm.UNKNOWN for r in R]
    files['06_host_strip.txt'] = colorstrip('Host', ids, host, palette_for(host, hm.HOST_COLORS),
                                            order=list(hm.HOST_COLORS))

    src = [r.get('source_group') or hm.UNKNOWN for r in R]
    files['07_source_strip.txt'] = colorstrip('Isolation source', ids, src,
                                              palette_for(src, hm.SOURCE_COLORS), order=list(hm.SOURCE_COLORS))

    cl = [clusters.get(i, '') for i in ids]
    multi = {c for c, n in Counter(c for c in cl if c).items() if n > 1}
    top = {c for c, _ in Counter(c for c in cl if c in multi).most_common(15)}
    cl = [c if c in top else ('Other cluster' if c in multi else 'Not clustered') for c in cl]
    ccol = palette_for([c for c in cl if c in top])
    ccol.update({'Other cluster': '#bdbdbd', 'Not clustered': '#ffffff'})
    files['08_outbreak_cluster_strip.txt'] = colorstrip('Outbreak cluster (top 15)', ids, cl, ccol,
                                                        order=sorted(top) + ['Other cluster', 'Not clustered'])

    def flag(r, marker):
        vs = r.get('virulence_summary') or ''
        return '1' if f'{marker}+' in vs else ('0' if f'{marker}-' in vs else '-1')
    files['09_PVL_TSST_mecA_binary.txt'] = [
        'DATASET_BINARY', 'SEPARATOR TAB', 'DATASET_LABEL\tPVL / TSST / mecA', 'COLOR\t#b2182b',
        'FIELD_SHAPES\t2\t2\t2', 'FIELD_LABELS\tPVL\tTSST\tmecA', 'FIELD_COLORS\t#d73027\t#4575b4\t#b2182b',
        'LEGEND_TITLE\tToxins / mecA', 'LEGEND_SHAPES\t2\t2\t2',
        'LEGEND_COLORS\t#d73027\t#4575b4\t#b2182b', 'LEGEND_LABELS\tPVL+\tTSST+\tmecA+',
        'HEIGHT_FACTOR\t0.9', 'SHOW_INTERNAL\t0', 'DATA'] + [
        f"{i}\t{flag(r, 'PVL')}\t{flag(r, 'TSST')}\t{'1' if 'mecA' in (r.get('amrfinder_genes') or '') else '0'}"
        for i, r in zip(ids, R)]

    files['10_AMR_gene_count_bar.txt'] = [
        'DATASET_SIMPLEBAR', 'SEPARATOR TAB', 'DATASET_LABEL\tAMR genes (AMRFinderPlus)',
        'COLOR\t#d45500', 'WIDTH\t120', 'SHOW_INTERNAL\t0', 'DATA'] + [
        f"{i}\t{len([g for g in (r.get('amrfinder_genes') or '').split(';') if g])}" for i, r in zip(ids, R)]

    # Hover/search metadata (iTOL shows it in the leaf popup and allows searching)
    fields = ['mlst_st', 'spa_type', 'sccmec_type', 'agr_group', 'year', 'host_group',
              'source_group', 'geo_loc_country', 'collection_date', 'isolation_source', 'biosample_accession']
    files['11_popup_metadata.txt'] = ['METADATA', 'SEPARATOR TAB',
                                      'FIELD_LABELS\t' + '\t'.join(fields + ['outbreak_cluster']), 'DATA'] + [
        i + '\t' + '\t'.join((r.get(f) or '') for f in fields) + '\t' + clusters.get(i, '')
        for i, r in zip(ids, R)]
    return files, len(ids)


def main(treefile, summary, cluster_file, outdir):
    rows = list(csv.DictReader(open(summary), delimiter='\t'))
    clusters = {}
    if cluster_file != 'NONE' and os.path.exists(cluster_file):
        for r in csv.DictReader(open(cluster_file), delimiter='\t'):
            clusters[r['sample_id']] = r.get('outbreak_cluster') or ''
    tips = tree_tips(open(treefile).read())
    files, n = build(tips, rows, clusters)
    os.makedirs(outdir, exist_ok=True)
    for name, lines in files.items():
        with open(os.path.join(outdir, name), 'w') as f:
            f.write('\n'.join(lines) + '\n')
    shutil.copy(treefile, os.path.join(outdir, '00_tree.newick'))
    missing = len(tips) - n
    print(f'{len(files)} iTOL files for {n} of {len(tips)} tree tips -> {outdir}'
          + (f' ({missing} tips without summary rows)' if missing else ''))


def _selftest():
    tips = tree_tips('((A:0.1,B:0.2)90:0.3,C:0.4);')
    assert tips == ['A', 'B', 'C'], tips
    rows = [dict(sample_id='A', mlst_st='5', sccmec_type='Type IV', agr_group='gp2', year='2020',
                 host_group='Human', source_group='Blood', virulence_summary='PVL+;TSST-',
                 amrfinder_genes='mecA;blaZ'),
            dict(sample_id='B', mlst_st='-', sccmec_type='Negative', agr_group='', year='Unknown',
                 host_group='Camel', source_group='Food: meat', virulence_summary='', amrfinder_genes='')]
    files, n = build(tips, rows, {'A': 'OB-1', 'B': 'OB-1'})
    assert n == 2
    assert files['01_tip_labels.txt'][-2:] == ['A\tA | ST5 | 2020', 'B\tB | ST? | Unknown']
    assert files['03_SCCmec_strip.txt'][-1] == 'B\t#d9d9d9\tMSSA', files['03_SCCmec_strip.txt'][-1]
    assert files['06_host_strip.txt'][-1] == 'B\t#a6761d\tCamel'
    assert files['09_PVL_TSST_mecA_binary.txt'][-2:] == ['A\t1\t0\t1', 'B\t-1\t-1\t0']
    assert files['10_AMR_gene_count_bar.txt'][-2:] == ['A\t2', 'B\t0']
    legend = [l for l in files['07_source_strip.txt'] if l.startswith('LEGEND_LABELS')][0]
    assert legend == 'LEGEND_LABELS\tBlood\tFood: meat', legend
    print('selftest ok')


if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        _selftest()
    elif len(sys.argv) == 5:
        main(*sys.argv[1:])
    else:
        sys.exit(__doc__)
