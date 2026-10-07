#!/usr/bin/env python3
"""Build a staphit2 samplesheet from any file containing ENA/SRA accessions.

Every line of the input (CSV, TSV, manifest, plain list — any column) is scanned for
accessions. Per line only the most specific level found is used, so a manifest row
holding project + sample + run resolves to just that run:
    run        SRR/ERR/DRR
    experiment SRX/ERX/DRX
    sample     SRS/ERS/DRS, SAMN/SAMEA/SAMD
    study      SRP/ERP/DRP, PRJNA/PRJEB/PRJDB

Runs are merged into one sample when they share a BioSample OR appear on the same line
(so "MRSA_94  ERR_illumina  ERR_nanopore" pairs two BioSamples into one hybrid sample).
Each sample's read mode is detected from instrument_platform:
    short  — paired-end short reads only          (fastq_1, fastq_2)
    hybrid — paired-end short + ONT long reads    (fastq_1, fastq_2, long_fastq)
    long   — ONT long reads only                  (long_fastq)

FASTQ columns are ENA https URLs (Nextflow downloads them at run time), or, with
--download DIR, local paths: files are fetched once into DIR (parallel, retried,
size + MD5 checked; files already present with the right size are skipped).
Prefer --download for large sets: EBI refuses bursts of connections and Nextflow
treats a refused remote file as fatal for the whole run.
If a sample has several runs of one kind, the run with the most bases is kept.
PacBio and single-end short-read runs are skipped.

Usage: accessions_to_samplesheet.py <any file with accessions> <samplesheet.csv> [--download DIR]
"""
import argparse
import hashlib
import os
import time
import csv
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

# ── 1. Find accessions ───────────────────────────────────────────────────────
LEVELS = {  # most specific first
    "run":        r"[EDS]RR\d{6,}",
    "experiment": r"[EDS]RX\d{6,}",
    "sample":     r"[EDS]RS\d{6,}|SAM(?:N|EA|D)\d+",
    "study":      r"[EDS]RP\d{6,}|PRJ(?:NA|EB|DB)\d+",
}
PATTERNS = {lvl: re.compile(rf"(?<![A-Za-z0-9])(?:{p})(?![0-9])") for lvl, p in LEVELS.items()}


def find_accessions(line):
    """Accessions of the most specific level present in a line (empty list if none)."""
    for lvl, pat in PATTERNS.items():
        hits = pat.findall(line)
        if hits:
            return sorted(set(hits))
    return []


# ── 2. Resolve accessions to runs (ENA) ──────────────────────────────────────
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = "run_accession,sample_accession,instrument_platform,base_count,fastq_ftp,fastq_bytes,fastq_md5"


def ena_runs(accession, tries=6):
    """ENA run rows for an accession. Retries on network errors and on ENA's
    'ERROR occurred' body, which it returns with HTTP 200 when overloaded."""
    q = urllib.parse.urlencode({"accession": accession, "result": "read_run", "fields": FIELDS, "format": "tsv"})
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(f"{ENA}?{q}", timeout=120) as r:
                text = r.read().decode()
            if "ERROR occurred" in text:
                raise IOError("ENA returned a partial/error response")
            return list(csv.DictReader(text.splitlines(), delimiter="\t"))
        except urllib.error.HTTPError as e:
            if e.code in (204, 404):  # no runs for this accession
                return []
            err = e
        except Exception as e:
            err = e
        if attempt < tries:
            time.sleep(min(120, 5 * 2 ** attempt))
    raise RuntimeError(f"ENA lookup failed for {accession}: {err}")


def resolve(accessions, fetch=ena_runs, threads=4):
    """{accession: [run rows]}"""
    with ThreadPoolExecutor(threads) as pool:
        return dict(zip(accessions, pool.map(fetch, accessions)))


# ── 3. Group runs into samples ───────────────────────────────────────────────
def group(line_accessions, runs_by_acc):
    """Union runs that share a BioSample or a line. Returns list of run-row lists."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        parent[find(a)] = find(b)

    runs = {}
    for accs in line_accessions:
        line_runs = [r for a in accs for r in runs_by_acc.get(a, [])]
        for r in line_runs:
            runs[r["run_accession"]] = r
            union("run:" + r["run_accession"], "bs:" + r["sample_accession"])
            union("run:" + r["run_accession"], "run:" + line_runs[0]["run_accession"])
    groups = defaultdict(list)
    for acc, r in runs.items():
        groups[find("run:" + acc)].append(r)
    return list(groups.values())


# ── 4. Classify runs and pick files ──────────────────────────────────────────
LONG_PLATFORMS = {"OXFORD_NANOPORE"}
SKIP_PLATFORMS = {"PACBIO_SMRT"}  # ponytail: Dragonflye is ONT-only; add a Flye --pacbio-hifi path if PacBio is needed


def classify(run):
    """(kind, [urls]) for one ENA run, or (None, reason)."""
    platform = run["instrument_platform"]
    urls = ["https://" + u for u in run["fastq_ftp"].split(";") if u]
    sizes = [int(b or 0) for b in run["fastq_bytes"].split(";")] if run["fastq_bytes"] else [0] * len(urls)
    if platform in SKIP_PLATFORMS:
        return None, f"{platform} not supported"
    if not urls:
        return None, "no FASTQ in ENA"
    if platform in LONG_PLATFORMS:
        return "long", [max(zip(sizes, urls))[1]]  # largest file (skips ENA stub files)
    r1 = [u for u in urls if u.endswith("_1.fastq.gz")]
    r2 = [u for u in urls if u.endswith("_2.fastq.gz")]
    if r1 and r2:
        return "short", [r1[0], r2[0]]
    return None, f"{platform} single-end short reads not supported"


def to_row(runs, skipped):
    """One samplesheet row from a group of runs (None if nothing usable)."""
    best = {}
    for r in runs:
        kind, files = classify(r)
        if kind is None:
            skipped.append(f"{r['run_accession']}\t{r['sample_accession']}\t{files}")
        elif int(r["base_count"] or 0) >= best.get(kind, (-1,))[0]:
            best[kind] = (int(r["base_count"] or 0), files)
    if not best:
        return None
    s, l = best.get("short"), best.get("long")
    return {
        "sample": "_".join(sorted({r["sample_accession"] for r in runs})),
        "fastq_1": s[1][0] if s else "",
        "fastq_2": s[1][1] if s else "",
        "long_fastq": l[1][0] if l else "",
        "mode": "hybrid" if s and l else ("long" if l else "short"),
    }


# ── 5. Optional local download ───────────────────────────────────────────────
def file_info(runs_by_acc):
    """{url: (bytes, md5)} for every FASTQ of the resolved runs."""
    info = {}
    for runs in runs_by_acc.values():
        for r in runs:
            urls = [u for u in r["fastq_ftp"].split(";") if u]
            sizes = r["fastq_bytes"].split(";") if r["fastq_bytes"] else [""] * len(urls)
            md5s = r.get("fastq_md5", "").split(";") if r.get("fastq_md5") else [""] * len(urls)
            for u, b, m in zip(urls, sizes, md5s):
                info["https://" + u] = (int(b or 0), m)
    return info


def fetch_file(url, dest, size, md5, tries=6):
    """Download url -> dest unless dest already has the expected size. Returns dest."""
    if os.path.exists(dest) and (not size or os.path.getsize(dest) == size):
        return dest
    part = dest + ".part"
    for attempt in range(1, tries + 1):
        try:
            h = hashlib.md5()
            with urllib.request.urlopen(url, timeout=300) as r, open(part, "wb") as out:
                for chunk in iter(lambda: r.read(1 << 20), b""):
                    out.write(chunk)
                    h.update(chunk)
            if size and os.path.getsize(part) != size:
                raise IOError(f"size {os.path.getsize(part)} != {size}")
            if md5 and h.hexdigest() != md5:
                raise IOError(f"md5 {h.hexdigest()} != {md5}")
            os.replace(part, dest)
            return dest
        except Exception as e:
            print(f"retry {attempt}/{tries} {os.path.basename(dest)}: {e}", file=sys.stderr)
            time.sleep(min(300, 10 * 2 ** attempt))
    raise RuntimeError(f"failed to download {url}")


def download(rows, info, outdir, threads=4):
    """Fetch every URL in rows into outdir and rewrite rows to local absolute paths."""
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    cols = ("fastq_1", "fastq_2", "long_fastq")
    urls = sorted({r[c] for r in rows for c in cols if r[c].startswith("https://")})
    local = {u: os.path.join(outdir, os.path.basename(u)) for u in urls}
    done = 0
    with ThreadPoolExecutor(threads) as pool:
        jobs = [pool.submit(fetch_file, u, local[u], *info.get(u, (0, ""))) for u in urls]
        for j in jobs:
            j.result()
            done += 1
            print(f"\rdownloaded {done}/{len(urls)}", end="", file=sys.stderr, flush=True)
    print(file=sys.stderr)
    for r in rows:
        for c in cols:
            r[c] = local.get(r[c], r[c])


# ── 6. Pipeline ──────────────────────────────────────────────────────────────
TABLE_COLS = {"run_accession", "sample_accession", "instrument_platform", "fastq_ftp"}


def runs_from_table(lines):
    """If lines are an ENA run table (TSV/CSV with TABLE_COLS), return its rows; else None.
    Lets ENA metadata exports be used directly, without the (sometimes down) ENA API."""
    if not lines:
        return None
    delim = "\t" if "\t" in lines[0] else ","
    header = set(lines[0].rstrip("\r\n").split(delim))
    if not TABLE_COLS <= header:
        return None
    rows = list(csv.DictReader(lines, delimiter=delim))
    for r in rows:
        r.setdefault("base_count", "")
        r.setdefault("fastq_bytes", "")
        r.setdefault("fastq_md5", "")
    return rows


def build(lines, fetch=ena_runs):
    lines = list(lines)
    table = runs_from_table(lines)
    if table is not None:  # ENA run table: one line per run, no API lookups
        line_accs = [[r["run_accession"]] for r in table]
        runs_by_acc = {r["run_accession"]: [r] for r in table}
    else:
        line_accs = [a for a in (find_accessions(l) for l in lines) if a]
        runs_by_acc = resolve(sorted({a for accs in line_accs for a in accs}), fetch)
    skipped = [f"{a}\t-\tno runs in ENA" for a, runs in runs_by_acc.items() if not runs]
    rows = [r for r in (to_row(g, skipped) for g in group(line_accs, runs_by_acc)) if r]
    return sorted(rows, key=lambda r: r["sample"]), skipped, file_info(runs_by_acc)


def main(src, dst, download_dir=None):
    with open(src) as f:
        rows, skipped, info = build(f)
    if download_dir:
        download(rows, info, download_dir)
    with open(dst, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "fastq_1", "fastq_2", "long_fastq"],
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for line in skipped:
        print("SKIPPED\t" + line, file=sys.stderr)
    n = Counter(r["mode"] for r in rows)
    print(f"{len(rows)} samples -> {dst}: short={n['short']} hybrid={n['hybrid']} long={n['long']}")


def _selftest():
    # find_accessions: most specific level wins; URLs and "_1" suffixes handled
    assert find_accessions("PRJEB1\tSAMEA9\tERR1043239\tftp://x/ERR1043239/ERR1043239_1.fastq.gz") == ["ERR1043239"]
    assert find_accessions("PRJEB83124,SAMEA111505969,SHA01M2") == ["SAMEA111505969"]
    assert find_accessions("strain s32 GZ261") == []

    def run(acc, bs, platform, ftp, sizes="", bases=0):
        return dict(run_accession=acc, sample_accession=bs, instrument_platform=platform,
                    fastq_ftp=ftp, fastq_bytes=sizes, base_count=str(bases))
    ena = {
        "ERR1000001": [run("ERR1000001", "SAMEA1", "ILLUMINA", "h/a_1.fastq.gz;h/a_2.fastq.gz")],
        "SAMEA2": [run("ERR1000002", "SAMEA2", "DNBSEQ", "h/b.fastq.gz;h/b_1.fastq.gz;h/b_2.fastq.gz", "169;9;9"),
                   run("ERR1000003", "SAMEA2", "OXFORD_NANOPORE", "h/small.fastq.gz", bases=10),
                   run("ERR1000004", "SAMEA2", "OXFORD_NANOPORE", "h/big.fastq.gz", bases=99)],
        "ERR1000005": [run("ERR1000005", "SAMEA5", "DNBSEQ", "h/c_1.fastq.gz;h/c_2.fastq.gz")],
        "ERR1000006": [run("ERR1000006", "SAMEA6", "OXFORD_NANOPORE", "h/d.fastq.gz")],
        "ERR1000007": [run("ERR1000007", "SAMEA7", "OXFORD_NANOPORE", "h/e_1.fastq.gz")],
        "ERR1000008": [run("ERR1000008", "SAMEA8", "PACBIO_SMRT", "h/f.fastq.gz")],
    }
    rows, skipped, _ = build([
        "sample,run\n",                              # header: no accessions
        "x,ERR1000001\n",                            # short
        "PRJEB1;SAMEA2\n",                           # hybrid via BioSample (sample beats project)
        "MRSA_94\tERR1000005\tERR1000006\n",         # hybrid via same line, two BioSamples
        "ERR1000007\n", "ERR1000007 again\n",        # long, duplicate line
        "ERR1000008\n",                              # PacBio -> skipped
    ], fetch=lambda a: ena.get(a, []))
    by = {r["sample"]: r for r in rows}
    assert set(by) == {"SAMEA1", "SAMEA2", "SAMEA5_SAMEA6", "SAMEA7"}, by
    assert by["SAMEA1"]["mode"] == "short"
    assert by["SAMEA2"]["mode"] == "hybrid" and by["SAMEA2"]["long_fastq"] == "https://h/big.fastq.gz"
    assert by["SAMEA2"]["fastq_1"] == "https://h/b_1.fastq.gz"
    assert by["SAMEA5_SAMEA6"]["mode"] == "hybrid"
    assert by["SAMEA7"]["mode"] == "long"
    assert len(skipped) == 1

    # ENA run table input: used directly, the fetch function must not be called
    def no_fetch(a):
        raise AssertionError("API called for table input")
    table = ["run_accession\tsample_accession\tinstrument_platform\tfastq_ftp\tfastq_bytes\n",
             "ERR2000001\tSAMEA20\tOXFORD_NANOPORE\th/x.fastq.gz\t5\n",
             "ERR2000002\tSAMEA20\tILLUMINA\th/y_1.fastq.gz;h/y_2.fastq.gz\t5;5\n"]
    rows, skipped, _ = build(table, fetch=no_fetch)
    assert [r["mode"] for r in rows] == ["hybrid"], rows
    print("selftest ok")


if __name__ == "__main__":
    if sys.argv[1:] == ["--selftest"]:
        _selftest()
    else:
        ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        ap.add_argument("input", help="any file containing accessions")
        ap.add_argument("samplesheet", help="output staphit2 samplesheet.csv")
        ap.add_argument("--download", metavar="DIR", help="download FASTQs into DIR and write local paths")
        a = ap.parse_args()
        main(a.input, a.samplesheet, a.download)
