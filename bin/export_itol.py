#!/usr/bin/env python3
"""
Generate iTOL-compatible annotation files from staphit2 pipeline outputs.

Usage:
    python bin/export_itol.py <results_subdir> [run_name] [pubmlst_metadata.tsv]

Output:
    itol/<run_name>/  — numbered .txt files, drag all into the iTOL annotation panel
"""
import sys, csv, re
from pathlib import Path
from collections import Counter

if len(sys.argv) < 2:
    sys.exit("Usage: export_itol.py <results_subdir> [run_name] [metadata.tsv]")

outdir    = Path(sys.argv[1])
run_name  = sys.argv[2] if len(sys.argv) > 2 else outdir.name
meta_path = Path(sys.argv[3]) if len(sys.argv) > 3 else None

summary_file = outdir / "summary" / "combined_summary.tsv"
cluster_file = outdir / "clusters" / "clusters.tsv"

if not summary_file.exists():
    sys.exit(f"Not found: {summary_file}")

# ── Helpers ────────────────────────────────────────────────────────────────────
def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

PALETTE_20 = [
    "#4daf4a","#377eb8","#ff7f00","#984ea3","#e41a1c",
    "#a65628","#f781bf","#999999","#1b9e77","#d95f02",
    "#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d",
    "#8dd3c7","#fb8072","#80b1d3","#fdb462","#b3de69",
]

def palette(n, colors=None):
    base = colors or PALETTE_20
    return [base[i % len(base)] for i in range(n)]

def header(dtype, label, color="#333333", extra_lines=""):
    lines = [dtype, "SEPARATOR\tTAB", f"DATASET_LABEL\t{label}", f"COLOR\t{color}"]
    if extra_lines:
        lines.append(extra_lines.rstrip())
    return "\n".join(lines)

def has_gene(gene_str, gene):
    if not gene_str:
        return 0
    return 1 if re.search(r'\b' + re.escape(gene) + r'\b', gene_str) else 0

def count_genes(gene_str):
    if not gene_str or gene_str.strip() in ("-", ""):
        return 0
    return len([g for g in re.split(r'[;,\s]+', gene_str.strip()) if g])

def parse_binary_virulence(vs, marker):
    """Return 1/0/-1 (present/absent/unknown). iTOL shows -1 as empty."""
    if not vs:
        return -1
    return 1 if f"{marker}+" in vs else 0

# ── Load data ──────────────────────────────────────────────────────────────────
summary = read_tsv(summary_file)
by_id   = {r["sample_id"]: r for r in summary}
ids     = [r["sample_id"] for r in summary]

clusters = {}
if cluster_file.exists():
    for r in read_tsv(cluster_file):
        clusters[r["sample_id"]] = r

pub_meta = {}
if meta_path and meta_path.exists():
    for r in read_tsv(meta_path):
        acc = r.get("run_accession", "").strip()
        if acc:
            pub_meta[acc] = r

if not pub_meta:
    for candidate in [
        outdir.parent / f"{run_name}_metadata.tsv",
        outdir.parent.parent / f"{run_name}_metadata.tsv",
        Path(f"{run_name}_metadata.tsv"),
        Path("cc97_metadata.tsv"),
    ]:
        if candidate.exists():
            for r in read_tsv(candidate):
                acc = r.get("run_accession", "").strip()
                if acc:
                    pub_meta[acc] = r
            if pub_meta:
                break

# ── Output directory ───────────────────────────────────────────────────────────
itol_dir = Path("itol") / run_name
itol_dir.mkdir(parents=True, exist_ok=True)
written = []

def write_file(name, lines):
    p = itol_dir / name
    p.write_text("\n".join(lines) + "\n")
    written.append(name)

# ── 1. Country — DATASET_SYMBOL ────────────────────────────────────────────────
country_map = {}
for sid in ids:
    pm = pub_meta.get(sid, {})
    c  = pm.get("country", "").split(":")[0].strip()
    if c:
        country_map[sid] = c

unique_countries = sorted(set(country_map.values()))
country_colors   = dict(zip(unique_countries, palette(len(unique_countries))))

write_file("01_country_symbols.txt", [
    header("DATASET_SYMBOL", "Country", "#333333",
        "LEGEND_TITLE\tCountry\n"
        "LEGEND_SHAPES\t" + "\t".join(["1"] * len(unique_countries)) + "\n"
        "LEGEND_COLORS\t" + "\t".join(country_colors[c] for c in unique_countries) + "\n"
        "LEGEND_LABELS\t" + "\t".join(unique_countries)),
    "DATA",
    "#node_id\tsymbol\tsize\tcolor\tfill\tposition",
    *[f"{sid}\t1\t10\t{country_colors[country_map[sid]]}\t1\t1"
      for sid in ids if sid in country_map],
])

# ── 2. Year — DATASET_TEXT ─────────────────────────────────────────────────────
year_rows = []
for sid in ids:
    year = pub_meta.get(sid, {}).get("year", "").strip()
    if year:
        year_rows.append(f"{sid}\t{year}\t1\t#555555\t0.8\tnormal\t0")

write_file("02_year_labels.txt", [
    header("DATASET_TEXT", "Year", "#555555", "SHOW_INTERNAL\t0"),
    "DATA",
    "#node_id\tlabel\tposition\tcolor\tsize_factor\tstyle\trotation",
    *year_rows,
])

# ── 3. ST — DATASET_TEXT ───────────────────────────────────────────────────────
st_rows = []
for sid in ids:
    st = by_id[sid].get("mlst_st", "").strip()
    if st and st not in ("-", ""):
        st_rows.append(f"{sid}\tST{st}\t1\t#222222\t0.8\tnormal\t0")

write_file("03_st_labels.txt", [
    header("DATASET_TEXT", "ST", "#222222", "SHOW_INTERNAL\t0"),
    "DATA",
    "#node_id\tlabel\tposition\tcolor\tsize_factor\tstyle\trotation",
    *st_rows,
])

# ── 4. spa type — DATASET_TEXT ─────────────────────────────────────────────────
spa_skip = {"", "-", "not_determined", "undetermined", "not determined"}
spa_rows = []
for sid in ids:
    spa = by_id[sid].get("spa_type", "").strip().lower()
    raw = by_id[sid].get("spa_type", "").strip()
    if raw and spa not in spa_skip:
        spa_rows.append(f"{sid}\t{raw}\t1\t#444444\t0.8\tnormal\t0")

write_file("04_spa_labels.txt", [
    header("DATASET_TEXT", "spa type", "#444444", "SHOW_INTERNAL\t0"),
    "DATA",
    "#node_id\tlabel\tposition\tcolor\tsize_factor\tstyle\trotation",
    *spa_rows,
])

# ── 5. SCCmec — DATASET_COLORSTRIP ────────────────────────────────────────────
SCCMEC_COLORS = {
    "MSSA":                "#d9d9d9",
    "Type IV":             "#4daf4a",
    "Type IVa":            "#4daf4a",
    "Type IVb":            "#6dbf6d",
    "Type IVc":            "#8fd08f",
    "Type IVd":            "#b2e0b2",
    "Type V":              "#377eb8",
    "Type VII":            "#ff7f00",
    "Type XI":             "#984ea3",
    "Composite (Type IV)": "#a8d978",
}
SCCMEC_DEFAULT = "#eeeeee"
SCCMEC_MSSA    = {"", "-", "ND", "not detected", "Negative", "MSSA"}

def norm_sccmec(val):
    v = (val or "").strip()
    return "MSSA" if v in SCCMEC_MSSA else v

sccmec_vals   = [norm_sccmec(by_id[sid].get("sccmec_type", "")) for sid in ids]
unique_sccmec = sorted(set(sccmec_vals))
sccmec_leg_c  = "\t".join(SCCMEC_COLORS.get(v, SCCMEC_DEFAULT) for v in unique_sccmec)

write_file("05_sccmec_strip.txt", [
    header("DATASET_COLORSTRIP", "SCCmec", "#333333",
        "COLOR_BRANCHES\t0\nSHOW_LABELS\t1\nLABEL_SIZE\t0.8\n"
        "LEGEND_TITLE\tSCCmec\n"
        "LEGEND_SHAPES\t" + "\t".join(["1"] * len(unique_sccmec)) + "\n"
        "LEGEND_COLORS\t" + sccmec_leg_c + "\n"
        "LEGEND_LABELS\t" + "\t".join(unique_sccmec)),
    "DATA",
    "#node_id\tcolor\tlabel",
    *[f"{sid}\t{SCCMEC_COLORS.get(v, SCCMEC_DEFAULT)}\t{v}"
      for sid, v in zip(ids, sccmec_vals)],
])

# ── 6. agr group — DATASET_COLORSTRIP ─────────────────────────────────────────
AGR_COLORS  = {"gp1": "#1b9e77", "gp2": "#d95f02", "gp3": "#7570b3", "gp4": "#e7298a"}
AGR_DEFAULT = "#cccccc"

agr_vals   = [by_id[sid].get("agr_group", "").strip() for sid in ids]
unique_agr = sorted(set(v for v in agr_vals if v))

write_file("06_agr_strip.txt", [
    header("DATASET_COLORSTRIP", "agr group", "#333333",
        "COLOR_BRANCHES\t0\nSHOW_LABELS\t1\nLABEL_SIZE\t0.8\n"
        "LEGEND_TITLE\tagr group\n"
        "LEGEND_SHAPES\t" + "\t".join(["1"] * len(unique_agr)) + "\n"
        "LEGEND_COLORS\t" + "\t".join(AGR_COLORS.get(v, AGR_DEFAULT) for v in unique_agr) + "\n"
        "LEGEND_LABELS\t" + "\t".join(unique_agr)),
    "DATA",
    "#node_id\tcolor\tlabel",
    *[f"{sid}\t{AGR_COLORS.get(v, AGR_DEFAULT)}\t{v if v else 'ND'}"
      for sid, v in zip(ids, agr_vals)],
])

# ── 7. PVL / TSST / mecA — DATASET_BINARY ─────────────────────────────────────
write_file("07_virulence_binary.txt", [
    header("DATASET_BINARY", "PVL / TSST / mecA", "#d73027",
        "FIELD_SHAPES\t1\t1\t1\n"
        "FIELD_LABELS\tPVL\tTSST\tmecA\n"
        "FIELD_COLORS\t#d73027\t#4575b4\t#b2182b"),
    "DATA",
    "#node_id\tPVL\tTSST\tmecA",
    *["{}\t{}\t{}\t{}".format(
        sid,
        parse_binary_virulence(by_id[sid].get("virulence_summary", ""), "PVL"),
        parse_binary_virulence(by_id[sid].get("virulence_summary", ""), "TSST"),
        has_gene(by_id[sid].get("amrfinder_genes", ""), "mecA"),
    ) for sid in ids],
])

# ── 8. AMR gene count — DATASET_SIMPLEBAR ─────────────────────────────────────
counts    = {sid: count_genes(by_id[sid].get("amrfinder_genes", "")) for sid in ids}
max_count = max(counts.values(), default=1)

write_file("08_amr_count_bar.txt", [
    header("DATASET_SIMPLEBAR", "AMR gene count", "#d45500",
        f"WIDTH\t300\nHEIGHT_FACTOR\t0.5\nMAX_VALUE\t{max_count}"),
    "DATA",
    "#node_id\tvalue",
    *[f"{sid}\t{counts[sid]}" for sid in ids],
])

# ── 9. AMR top-5 genes — DATASET_BINARY ───────────────────────────────────────
all_genes = []
for r in summary:
    gs = r.get("amrfinder_genes", "")
    if gs and gs.strip() not in ("-", ""):
        all_genes.extend(g.strip() for g in re.split(r'[;,\s]+', gs) if g.strip())

top5 = [g for g, _ in Counter(all_genes).most_common(5)]
top5_colors = ["#e41a1c", "#ff7f00", "#4daf4a", "#377eb8", "#984ea3"]

write_file("09_amr_top5_binary.txt", [
    header("DATASET_BINARY", "Top 5 AMR genes", "#e41a1c",
        "FIELD_SHAPES\t" + "\t".join(["1"] * len(top5)) + "\n"
        "FIELD_LABELS\t" + "\t".join(top5) + "\n"
        "FIELD_COLORS\t" + "\t".join(top5_colors[:len(top5)])),
    "DATA",
    "#node_id\t" + "\t".join(top5),
    *["{}\t{}".format(sid, "\t".join(str(has_gene(by_id[sid].get("amrfinder_genes", ""), g)) for g in top5))
      for sid in ids],
])

# ── 10. Outbreak cluster — DATASET_COLORSTRIP ─────────────────────────────────
ob_counts = Counter(
    clusters[sid].get("outbreak_cluster", "")
    for sid in ids
    if sid in clusters and clusters[sid].get("outbreak_cluster", "")
)
ob_top = [c for c, _ in ob_counts.most_common(12) if c]
ob_colors = dict(zip(ob_top, palette(len(ob_top), [
    "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c",
    "#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928",
])))

ob_data = [
    f"{sid}\t{ob_colors[ob]}\t{ob}"
    for sid in ids
    for ob in [clusters.get(sid, {}).get("outbreak_cluster", "")]
    if ob in ob_colors
]

write_file("10_outbreak_strip.txt", [
    header("DATASET_COLORSTRIP", "Outbreak cluster", "#333333",
        "COLOR_BRANCHES\t0\nSHOW_LABELS\t1\nLABEL_SIZE\t0.8\n"
        "LEGEND_TITLE\tOutbreak cluster\n"
        "LEGEND_SHAPES\t" + "\t".join(["1"] * len(ob_top)) + "\n"
        "LEGEND_COLORS\t" + "\t".join(ob_colors[c] for c in ob_top) + "\n"
        "LEGEND_LABELS\t" + "\t".join(ob_top)),
    "DATA",
    "#node_id\tcolor\tlabel",
    *ob_data,
])

# ── Summary ────────────────────────────────────────────────────────────────────
print(f"\niTOL files written to: {itol_dir}/")
for f in written:
    print(f"  {f}")
if top5:
    print(f"\nTop 5 AMR genes detected: {', '.join(top5)}")
print(f"\nDrag all {len(written)} files into the iTOL tree annotation panel.")
