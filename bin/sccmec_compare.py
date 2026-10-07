#!/usr/bin/env python3
"""Extract SCCmec regions from staphit2 results and draw a stacked, annotated comparison.

  extract  <results_dir> <samples.tsv> <outdir>
      samples.tsv: sample_id<TAB>label  (one per line; drawn top to bottom in this order)
      Writes <outdir>/regions/<sample>.gbk — the SCCmec region of each assembly, cut from the
      Prokka GenBank: the span of SCCmec typer hits (mec, ccr, IS) on the mecA contig, extended
      to orfX (rlmH, the integration site) when it is on the same contig, padded, and oriented
      so orfX is on the left. Gene categories are stored in a /sccmec_class qualifier.

  scaffold <results_dir> <sample> <ref_region.gbk> <contigs_vs_ref.paf> <out.gbk> <label>
      For fragmented (short-read) assemblies: place the sample's contigs onto a complete
      long-read SCCmec region of the same lineage (minimap2 PAF of the sample's Prokka .fna
      vs the reference region). Aligned pieces keep their Prokka genes and typer hits;
      uncovered reference stretches become N gaps (/sccmec_class=gap); contig ends are
      marked (misc_feature /note=contig_break).

  coverage <outdir> <groups.tsv>
      Read-level check: groups.tsv rows = sample<TAB>reference_sample<TAB>label, with
      <outdir>/coverage/<sample>.depth.tsv (samtools depth -a over the reference's element)
      and <sample>.meandepth (chromosome mean). Draws <outdir>/sccmec_read_coverage.{pdf,png}
      (reference genes on top, depth/chromosome-mean per sample below) and writes
      <outdir>/sccmec_gene_coverage.tsv (per labelled gene: normalized depth, % bases covered).

  plot     <outdir> <clinker_links.csv> [title]
      Reads the regions + clinker's gene-to-gene identities (clinker -o, -dl ',') and draws
      <outdir>/sccmec_comparison.{pdf,png,svg}: one track per sample, genes coloured by class,
      grey links between homologous genes of adjacent tracks (darker = higher identity).
"""
import csv
import glob
import os
import re
import sys

from Bio import SeqIO

PAD = 3000           # bp added around the element
MAX_SPAN = 90000     # cap on region length

CLASS_COLORS = {
    'orfX (rlmH)':          '#000000',
    'mecA':                 '#d7191c',
    'mec regulators':       '#fdae61',
    'ccr recombinase':      '#2c7bb6',
    'IS / transposase':     '#7b3294',
    'Metal resistance':     '#1a9641',
    'Other AMR':            '#e66101',
    'Other':                '#d9d9d9',
}


# ── Classification ────────────────────────────────────────────────────────────
def classify(feat, typer_hits):
    """(class, short label or '') for one CDS."""
    q = feat.qualifiers
    gene = (q.get('gene', [''])[0]).split('_')[0]
    prod = q.get('product', [''])[0].lower()
    s, e = int(feat.location.start), int(feat.location.end)
    ov = [h for h in typer_hits if min(e, h['end']) - max(s, h['start']) > 0.5 * (e - s)]
    hit = max(ov, key=lambda h: float(h['id_pct'])) if ov else None
    hname = hit['gene'] if hit else ''
    if gene == 'rlmH' or 'methyltransferase h' in prod:
        return 'orfX (rlmH)', 'orfX'
    if 'pbp2a' in prod or 'penicillin-binding protein 2' in prod or hname in ('mecA', 'mecC'):
        return 'mecA', hname or 'mecA'
    if hname.startswith(('mecR', 'mecI')) or gene in ('mecR1', 'mecI', 'blaR1', 'blaI') \
            or 'methicillin resistance' in prod:
        return 'mec regulators', hname or gene
    if hname.startswith('ccr'):
        return 'ccr recombinase', hname
    if hname.startswith('IS') or 'transposase' in prod or 'insertion element' in prod:
        return 'IS / transposase', hname.split('_')[0] if hname.startswith('IS') else ''
    if re.search(r'mercur|merA|merB|cadmium|cad[ACD]|arsen|copper|zinc', prod + ' ' + gene, re.I):
        return 'Metal resistance', gene if gene and not gene.startswith('hypothetical') else ''
    if re.search(r'aminoglycoside|tetracycline|erythromycin|streptomycin|kanamycin|bleomycin|'
                 r'fosfomycin|lincosamide|chloramphenicol|trimethoprim|beta-lactamase', prod):
        return 'Other AMR', gene
    return 'Other', ''


# ── Extraction ────────────────────────────────────────────────────────────────
def load_hits(path):
    if not os.path.exists(path):
        return []
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        r['start'], r['end'] = int(r['start']), int(r['end'])
    return rows


def extract_region(gbk, hits):
    """(SeqRecord of the oriented region, note) or (None, reason)."""
    records = {r.id: r for r in SeqIO.parse(gbk, 'genbank')}
    pbp = [(cid, f) for cid, r in records.items() for f in r.features if f.type == 'CDS'
           and ('pbp2a' in f.qualifiers.get('product', [''])[0].lower()
                or 'penicillin-binding protein 2' in f.qualifiers.get('product', [''])[0].lower())]
    if not pbp:
        return None, 'no mecA/PBP2a CDS'
    cid, mec = pbp[0]
    rec = records[cid]
    on = [h for h in hits if h['contig'] == cid] or [{'start': int(mec.location.start), 'end': int(mec.location.end)}]
    lo, hi = min(h['start'] for h in on), max(h['end'] for h in on)
    orfx = [f for f in rec.features if f.type == 'CDS' and f.qualifiers.get('gene', [''])[0] == 'rlmH']
    note = 'orfX on contig'
    if orfx and abs(int(orfx[0].location.start) - lo) < MAX_SPAN:
        o = orfx[0]
        lo, hi = min(lo, int(o.location.start)), max(hi, int(o.location.end))
        flip = o.location.strand == -1
    else:
        note = 'orfX not on mecA contig (fragmented)' if not orfx else 'orfX far from mec'
        flip = mec.location.strand == -1
    lo, hi = max(0, lo - PAD), min(len(rec), hi + PAD)
    if hi - lo > MAX_SPAN:
        hi = lo + MAX_SPAN
    sub = rec[lo:hi]
    if flip:
        sub = sub.reverse_complement(id=True, name=True, description=True, annotations=True)
    sub.annotations['molecule_type'] = 'DNA'
    sub.annotations['contig'] = cid
    return sub, f'{cid}:{lo}-{hi} ({hi - lo} bp, {"rev" if flip else "fwd"}; {note})'


def cmd_extract(results, samples_tsv, outdir):
    os.makedirs(f'{outdir}/regions', exist_ok=True)
    rows = [l.rstrip('\n').split('\t') for l in open(samples_tsv) if l.strip() and not l.startswith('#')]
    with open(f'{outdir}/regions.tsv', 'w') as log:
        log.write('order\tsample\tlabel\tregion\n')
        for i, (sid, label) in enumerate(rows, 1):
            hits = load_hits(f'{results}/sccmec/{sid}/{sid}_sccmec_elements.csv')
            sub, note = extract_region(f'{results}/prokka/{sid}/{sid}.gbk', hits)
            print(f'{sid}: {note}')
            if sub is None:
                continue
            # coordinates of typer hits inside the region (for classification)
            off = int(re.search(r':(\d+)-', note).group(1))
            end = int(re.search(r'-(\d+) ', note).group(1))
            flip = 'rev;' in note
            local = []
            for h in hits:
                if h['contig'] != sub.annotations['contig']:
                    continue
                s, e = h['start'] - off, h['end'] - off
                if flip:
                    s, e = (end - off) - e, (end - off) - s
                local.append(dict(h, start=s, end=e))
            for f in sub.features:
                if f.type == 'CDS':
                    cls, lab = classify(f, local)
                    f.qualifiers['sccmec_class'] = [cls]
                    if lab:
                        f.qualifiers['sccmec_label'] = [lab]
                    f.qualifiers.setdefault('locus_tag', [f'{sid}_{int(f.location.start)}'])
            sub.id = sub.name = sid[-10:]
            sub.description = label
            SeqIO.write(sub, f'{outdir}/regions/{i:02d}_{sid}.gbk', 'genbank')
            log.write(f'{i}\t{sid}\t{label}\t{note}\n')


# ── Reference-guided scaffold for fragmented assemblies ───────────────────────
def cmd_scaffold(results, sid, ref_gbk, paf, out_gbk, label, min_len=500):
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.SeqRecord import SeqRecord
    ref = SeqIO.read(ref_gbk, 'genbank')
    tlen = len(ref)
    contigs = {r.id: r for r in SeqIO.parse(f'{results}/prokka/{sid}/{sid}.gbk', 'genbank')}
    hits = load_hits(f'{results}/sccmec/{sid}/{sid}_sccmec_elements.csv')
    blocks = []
    for line in open(paf):
        p = line.split('\t')
        qn, qs, qe, strand, ts, te, alen = p[0], int(p[2]), int(p[3]), p[4], int(p[7]), int(p[8]), int(p[10])
        if alen >= min_len and qn in contigs:
            blocks.append((alen, qn, qs, qe, strand, ts, te))
    placed, used = [], []                        # greedy: longest alignments first, no overlap
    for b in sorted(blocks, reverse=True):
        ts, te = b[5], b[6]
        if all(te <= u0 or ts >= u1 for u0, u1 in used):
            placed.append(b); used.append((ts, te))
    seq = ['N'] * tlen
    feats = []
    for alen, qn, qs, qe, strand, ts, te in sorted(placed, key=lambda b: b[5]):
        c = contigs[qn]
        piece = c.seq[qs:qe] if strand == '+' else c.seq[qs:qe].reverse_complement()
        n = min(len(piece), tlen - ts)
        seq[ts:ts + n] = str(piece[:n])
        def to_ref(x):
            return ts + (x - qs) if strand == '+' else ts + (qe - x)
        local_hits = [dict(h, start=min(to_ref(h['start']), to_ref(h['end'])),
                           end=max(to_ref(h['start']), to_ref(h['end'])))
                      for h in hits if h['contig'] == qn and h['start'] >= qs and h['end'] <= qe]
        for f in c.features:
            if f.type != 'CDS' or int(f.location.start) < qs or int(f.location.end) > qe:
                continue
            a, b = sorted((to_ref(int(f.location.start)), to_ref(int(f.location.end))))
            st = f.location.strand if strand == '+' else -f.location.strand
            nf = SeqFeature(FeatureLocation(max(0, a), min(tlen, b), strand=st), type='CDS',
                            qualifiers=dict(f.qualifiers))
            nf.qualifiers.pop('translation', None)
            cls, lab = classify(nf, local_hits)
            nf.qualifiers['sccmec_class'] = [cls]
            if lab:
                nf.qualifiers['sccmec_label'] = [lab]
            nf.qualifiers['translation'] = f.qualifiers.get('translation', [''])
            feats.append(nf)
        for x in (ts, min(tlen, ts + n)):        # contig ends inside the region
            feats.append(SeqFeature(FeatureLocation(x, x), type='misc_feature',
                                    qualifiers={'note': ['contig_break'], 'contig': [qn]}))
    # explicit gaps
    i = 0
    while i < tlen:
        if seq[i] == 'N':
            j = i
            while j < tlen and seq[j] == 'N':
                j += 1
            feats.append(SeqFeature(FeatureLocation(i, j), type='misc_feature',
                                    qualifiers={'sccmec_class': ['gap'], 'note': ['not in assembly / not placed']}))
            i = j
        else:
            i += 1
    rec = SeqRecord(Seq(''.join(seq)), id=sid[-10:], name=sid[-10:], description=label,
                    features=feats, annotations={'molecule_type': 'DNA'})
    SeqIO.write(rec, out_gbk, 'genbank')
    cov = sum(1 for x in seq if x != 'N') / tlen
    print(f'{sid}: {len(placed)} contig pieces placed on {ref.id} ({cov:.0%} of {tlen} bp covered)')


# ── Read coverage over a reference element ────────────────────────────────────
def _region_info(outdir):
    """{sample: (gbk_path, contig, lo, hi, flipped)} from regions.tsv written by extract."""
    info = {}
    for r in csv.DictReader(open(f'{outdir}/regions.tsv'), delimiter='\t'):
        m = re.match(r'(\S+):(\d+)-(\d+) \(\d+ bp, (rev|fwd)', r['region'])
        info[r['sample']] = (glob.glob(f"{outdir}/regions/*_{r['sample']}.gbk")[0],
                             m.group(1), int(m.group(2)), int(m.group(3)), m.group(4) == 'rev')
    return info


def cmd_coverage(outdir, groups_tsv):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrow
    info = _region_info(outdir)
    rows = [l.rstrip('\n').split('\t') for l in open(groups_tsv) if l.strip() and not l.startswith('#')]
    groups = {}
    for sid, ref, label in rows:
        groups.setdefault(ref, []).append((sid, label))
    table = ['sample\treference\tgene\tref_start\tref_end\tnorm_depth\tpct_covered\tcall']
    fig, axes = plt.subplots(len(groups) * 2, 1, figsize=(16, 4.2 * len(groups)),
                             gridspec_kw={'height_ratios': [1, 2.4] * len(groups)})
    for gi, (ref, members) in enumerate(groups.items()):
        gbk, contig, lo, hi, flip = info[ref]
        rec = SeqIO.read(gbk, 'genbank')
        L = hi - lo
        axg, axd = axes[2 * gi], axes[2 * gi + 1]
        labelled = []
        for f in rec.features:
            if f.type != 'CDS':
                continue
            s0, e0 = int(f.location.start), int(f.location.end)
            col = CLASS_COLORS[f.qualifiers.get('sccmec_class', ['Other'])[0]]
            hl = min(0.45, (e0 - s0) / 1000 * 0.35)
            x, dx = (s0 / 1000, (e0 - s0) / 1000) if f.location.strand != -1 else (e0 / 1000, -(e0 - s0) / 1000)
            axg.add_patch(FancyArrow(x, 0, dx, 0, width=0.42, head_width=0.62, head_length=hl,
                                     length_includes_head=True, color=col, ec='#333333', lw=0.3))
            lab = f.qualifiers.get('sccmec_label', [''])[0]
            if lab:
                axg.text((s0 + e0) / 2000, 0.45, lab, ha='center', va='bottom', fontsize=7, rotation=35,
                         fontstyle='italic' if lab[0].islower() else 'normal')
                labelled.append((lab, s0, e0, f.qualifiers['sccmec_class'][0]))
        axg.plot([0, L / 1000], [0, 0], color='#555555', lw=0.8, zorder=0)
        axg.set_xlim(0, L / 1000); axg.set_ylim(-0.6, 1.5); axg.axis('off')
        axg.set_title(f'Reference: {rec.description}  ({contig}:{lo}-{hi})', fontsize=10, loc='left')
        for lab, s0, e0, cls in labelled:          # shade ccr / mecA positions in the depth panel
            if cls in ('ccr recombinase', 'mecA'):
                axd.axvspan(s0 / 1000, e0 / 1000, color=CLASS_COLORS[cls], alpha=0.12, lw=0)
        for k, (sid, label) in enumerate(members):
            mean = float(open(f'{outdir}/coverage/{sid}.meandepth').read().strip() or 0) or 1.0
            depth = [0.0] * L
            for line in open(f'{outdir}/coverage/{sid}.depth.tsv'):
                _, pos, d = line.split('\t')
                x = int(pos) - 1 - lo
                if flip:
                    x = L - 1 - x
                if 0 <= x < L:
                    depth[x] = int(d) / mean
            w = 200                                 # 200-bp running mean for the plot
            xs = list(range(0, L, w))
            ys = [sum(depth[i:i + w]) / len(depth[i:i + w]) for i in xs]
            axd.plot([x / 1000 for x in xs], ys, lw=1.1, label=label, color=QUAL_LINES[k % len(QUAL_LINES)])
            for lab, s0, e0, cls in labelled:
                seg = depth[s0:e0]
                nd = sum(seg) / len(seg)
                cov = sum(1 for v in seg if v > 0.05) / len(seg)
                call = 'present' if nd >= 0.4 and cov >= 0.9 else ('partial' if cov >= 0.3 else 'absent')
                table.append(f'{sid}\t{ref}\t{lab}\t{s0}\t{e0}\t{nd:.2f}\t{100 * cov:.0f}\t{call}')
        axd.axhline(1, color='#999999', lw=0.6, ls='--')
        axd.set_xlim(0, L / 1000); axd.set_ylim(0, 3.2)
        axd.set_ylabel('depth / chromosome mean'); axd.set_xlabel('kb along reference element (orfX left)')
        axd.legend(fontsize=8, loc='upper right', frameon=False)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(f'{outdir}/sccmec_read_coverage.{ext}', dpi=180 if ext == 'png' else None)
    open(f'{outdir}/sccmec_gene_coverage.tsv', 'w').write('\n'.join(table) + '\n')
    print(f'wrote {outdir}/sccmec_read_coverage.pdf/.png and sccmec_gene_coverage.tsv')


QUAL_LINES = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']


# ── Plot ──────────────────────────────────────────────────────────────────────
def load_links(path):
    """{(locus_a, locus_b): identity} from clinker -o output with ',' delimiter."""
    links = {}
    for line in open(path):
        p = [x.strip() for x in line.split(',')]
        if len(p) >= 3:
            try:
                ident = float(p[2])
            except ValueError:
                continue
            links[(p[0], p[1])] = links[(p[1], p[0])] = ident
    return links


def cmd_plot(outdir, links_csv, title='SCCmec region comparison'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrow, Polygon, Patch

    files = sorted(glob.glob(f'{outdir}/regions/*.gbk'))
    recs = [SeqIO.read(f, 'genbank') for f in files]
    links = load_links(links_csv) if links_csv and os.path.exists(links_csv) else {}
    maxlen = max(len(r) for r in recs)
    gap = 1.6
    fig, ax = plt.subplots(figsize=(18, 1.05 * len(recs) + 1.8))
    genes = []  # per track: [(locus, x0, x1)]
    for i, r in enumerate(recs):
        y = -i * gap
        ax.plot([0, len(r) / 1000], [y, y], color='#555555', lw=0.8, zorder=1)
        for f in r.features:
            if f.type != 'misc_feature':
                continue
            if f.qualifiers.get('sccmec_class') == ['gap']:
                g0, g1 = int(f.location.start) / 1000, int(f.location.end) / 1000
                ax.plot([g0, g1], [y, y], color='white', lw=2.2, zorder=2)
                ax.plot([g0, g1], [y, y], color='#bbbbbb', lw=0.8, ls=(0, (2, 2)), zorder=2)
            elif f.qualifiers.get('note') == ['contig_break']:
                x = int(f.location.start) / 1000
                ax.plot([x, x], [y - 0.38, y + 0.38], color='#e31a1c', lw=1.1, zorder=4)
        ax.text(-0.8, y, r.description, ha='right', va='center', fontsize=8.5)
        track = []
        for f in r.features:
            if f.type != 'CDS':
                continue
            s, e = int(f.location.start) / 1000, int(f.location.end) / 1000
            cls = f.qualifiers.get('sccmec_class', ['Other'])[0]
            col = CLASS_COLORS[cls]
            L = e - s
            hl = min(0.45, L * 0.35)
            if f.location.strand == -1:
                ax.add_patch(FancyArrow(e, y, -L, 0, width=0.42, head_width=0.62, head_length=hl,
                                        length_includes_head=True, color=col, ec='#333333', lw=0.3, zorder=3))
            else:
                ax.add_patch(FancyArrow(s, y, L, 0, width=0.42, head_width=0.62, head_length=hl,
                                        length_includes_head=True, color=col, ec='#333333', lw=0.3, zorder=3))
            lab = f.qualifiers.get('sccmec_label', [''])[0]
            if lab:
                ax.text((s + e) / 2, y + 0.42, lab, ha='center', va='bottom', fontsize=6.5,
                        rotation=35, fontstyle='italic' if lab[0].islower() else 'normal', zorder=5)
            track.append((f.qualifiers['locus_tag'][0], s, e))
        genes.append(track)
    # homology links between adjacent tracks
    for i in range(len(recs) - 1):
        ya, yb = -i * gap - 0.32, -(i + 1) * gap + 0.32
        for la, a0, a1 in genes[i]:
            for lb, b0, b1 in genes[i + 1]:
                ident = links.get((la, lb))
                if ident is None or ident < 0.3:
                    continue
                shade = str(round(0.92 - 0.6 * (ident - 0.3) / 0.7, 2))
                ax.add_patch(Polygon([(a0, ya), (a1, ya), (b1, yb), (b0, yb)], closed=True,
                                     color=shade, alpha=0.75, lw=0, zorder=0))
    ax.set_xlim(-0.5, maxlen / 1000 + 1)
    ax.set_ylim(-(len(recs) - 1) * gap - 1.1, 1.6)
    ax.set_yticks([])
    for sp in ('left', 'right', 'top'):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel('kb (orfX on the left; region oriented 5\'→3\')')
    ax.set_title(title, fontsize=12, fontweight='bold', loc='left')
    handles = [Patch(color=c, label=k) for k, c in CLASS_COLORS.items()]
    from matplotlib.lines import Line2D
    handles += [Line2D([0], [0], color='#e31a1c', lw=1.5, label='contig break (short reads)'),
                Line2D([0], [0], color='#bbbbbb', lw=1, ls='--', label='not in assembly'),
                Patch(color='0.3', alpha=0.75, label='link ≥ 99% aa identity'),
                Patch(color='0.85', alpha=0.75, label='link ~40% identity')]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.06 if len(recs) > 6 else -0.12),
              ncol=5, fontsize=8, frameon=False)
    fig.tight_layout()
    for ext in ('pdf', 'png', 'svg'):
        fig.savefig(f'{outdir}/sccmec_comparison.{ext}', dpi=200 if ext == 'png' else None,
                    bbox_inches='tight')
    print(f'wrote {outdir}/sccmec_comparison.pdf/.png/.svg ({len(recs)} tracks, {len(links) // 2} links)')


if __name__ == '__main__':
    if len(sys.argv) >= 4 and sys.argv[1] == 'extract':
        cmd_extract(*sys.argv[2:5])
    elif len(sys.argv) == 4 and sys.argv[1] == 'coverage':
        cmd_coverage(*sys.argv[2:4])
    elif len(sys.argv) == 8 and sys.argv[1] == 'scaffold':
        cmd_scaffold(*sys.argv[2:8])
    elif len(sys.argv) >= 3 and sys.argv[1] == 'plot':
        cmd_plot(*sys.argv[2:5])
    else:
        sys.exit(__doc__)
