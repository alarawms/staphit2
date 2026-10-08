#!/usr/bin/env python3
"""
Export a self-contained .microreact project file from staphit2 outputs.
The file can be dragged into microreact.org or loaded into a self-hosted instance.

Usage:
    python bin/export_microreact.py <results_subdir> [run_name] [pubmlst_metadata.tsv]

Output:
    <run_name>.microreact  — drag into microreact.org to view
"""
import sys, json, csv, re
from pathlib import Path

# ── Arguments ──────────────────────────────────────────────────────────────────
if len(sys.argv) < 2:
    sys.exit("Usage: export_microreact.py <results_subdir> [run_name] [metadata.tsv]")

outdir    = Path(sys.argv[1])
run_name  = sys.argv[2] if len(sys.argv) > 2 else outdir.name
meta_path = Path(sys.argv[3]) if len(sys.argv) > 3 else None

tree_file    = outdir / "iqtree"   / "core.treefile"
summary_file = outdir / "summary"  / "combined_summary.tsv"
cluster_file = outdir / "clusters" / "clusters.tsv"

for f in [tree_file, summary_file]:
    if not f.exists():
        sys.exit(f"Not found: {f}")

# ── Load data ──────────────────────────────────────────────────────────────────
tree_str = tree_file.read_text().strip()

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

summary  = {r["sample_id"]: r for r in read_tsv(summary_file)}
clusters = {}
if cluster_file.exists():
    for r in read_tsv(cluster_file):
        clusters[r["sample_id"]] = r

# PubMLST/fetch metadata (country, year, source)
pub_meta = {}
if meta_path and meta_path.exists():
    for r in read_tsv(meta_path):
        acc = r.get("run_accession", "").strip()
        if acc:
            pub_meta[acc] = r

# Infer metadata path from project dir if not supplied
if not pub_meta:
    for candidate in [
        outdir.parent / "cc97_metadata.tsv",
        outdir.parent.parent / "cc97_metadata.tsv",
        Path("cc97_metadata.tsv"),
    ]:
        if candidate.exists():
            for r in read_tsv(candidate):
                acc = r.get("run_accession", "").strip()
                if acc:
                    pub_meta[acc] = r
            if pub_meta:
                break

# ── Build metadata CSV ────────────────────────────────────────────────────────
def parse_pvl(vs):
    if not vs:
        return ""
    return "PVL+" if "PVL+" in vs else ("PVL-" if "PVL-" in vs else "")

def parse_tsst(vs):
    if not vs:
        return ""
    return "TSST+" if "TSST+" in vs else ("TSST-" if "TSST-" in vs else "")

def has_gene(genes, gene):
    if not genes:
        return ""
    return "+" if re.search(r'\b' + gene + r'\b', genes) else "-"

rows = []
for sid, s in summary.items():
    pm  = pub_meta.get(sid, {})
    cl  = clusters.get(sid, {})

    country   = pm.get("country", "")
    # Strip sub-national qualifier: "Australia: Queensland" → "Australia"
    country   = country.split(":")[0].strip() if ":" in country else country

    rows.append({
        "__id":             sid,
        "Country":          country,
        "Continent":        pm.get("continent", ""),
        "Year":             pm.get("year", ""),
        "Source":           pm.get("source_category", pm.get("disease", "")),
        "Host":             pm.get("host", ""),
        "ST":               s.get("mlst_st", ""),
        "spa_type":         s.get("spa_type", ""),
        "SCCmec":           s.get("sccmec_type", ""),
        "agr_group":        s.get("agr_group", ""),
        "PVL":              parse_pvl(s.get("virulence_summary", "")),
        "TSST":             parse_tsst(s.get("virulence_summary", "")),
        "mecA":             has_gene(s.get("amrfinder_genes", ""), "mecA"),
        "Outbreak_cluster": cl.get("outbreak_cluster", ""),
        "Related_cluster":  cl.get("related_cluster", ""),
        "AMR_genes":        s.get("amrfinder_genes", ""),
    })

# Serialise CSV string
import io
buf = io.StringIO()
if rows:
    writer = csv.DictWriter(buf, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
csv_str = buf.getvalue()

# ── Build .microreact project JSON ────────────────────────────────────────────
project = {
    "version": 2,
    "info": {
        "name":        f"{run_name} — CC97 MRSA Core SNP Phylogeny",
        "description": f"staphit2 pipeline output · {len(rows)} isolates · IQ-TREE GTR+F+I+R3",
    },
    "files": {
        "tree-1": {
            "id":      "tree-1",
            "name":    "core.treefile",
            "format":  "newick",
            "content": tree_str,
        },
        "data-1": {
            "id":      "data-1",
            "name":    f"{run_name}_metadata.csv",
            "format":  "text/csv",
            "content": csv_str,
        },
    },
    "trees": {
        "tree-1": {
            "id":       "tree-1",
            "title":    "Core SNP tree",
            "file":     "tree-1",
            "labelField": "__id",
        }
    },
    "datasets": {
        "dataset-1": {
            "id":       "dataset-1",
            "file":     "data-1",
            "idField":  "__id",
        }
    },
    "views": {
        "default": {
            "tree":    "tree-1",
            "dataset": "dataset-1",
        }
    },
}

out_file = Path(f"{run_name}.microreact")
out_file.write_text(json.dumps(project, indent=2))

n_with_country = sum(1 for r in rows if r["Country"])
n_with_year    = sum(1 for r in rows if r["Year"])
print(f"Saved: {out_file}  ({len(rows)} isolates, {n_with_country} with country, {n_with_year} with year)")
print(f"Upload to:  https://microreact.org/upload")
print(f"Or drag the file directly into https://microreact.org")
