#!/usr/bin/env python3
"""
Merge an existing samplesheet with newly downloaded ENA/SRA FASTQs from
nf-core/fetchngs output, producing an expanded samplesheet for the next
pipeline run.

fetchngs names files as {EXPERIMENT}_{RUN}_1.fastq.gz. This script
extracts the run accession from the filename and uses it as the sample ID,
skipping any accessions already present in the existing samplesheet.

Usage:
    python3 bin/build_cc97_samplesheet.py \
        --existing  cc97/samplesheet_cc97_correct.csv \
        --fastq_dir /path/to/cc97_fetchngs/fastq \
        --output    cc97/samplesheet_cc97_v2.csv

Options:
    --existing   Existing samplesheet CSV to extend (sample,fastq_1,fastq_2).
                 If omitted, only the newly downloaded samples are written.
    --fastq_dir  Directory containing fetchngs FASTQ output files.
    --output     Path for the new combined samplesheet CSV.
    --dry_run    Print the result without writing the output file.
"""

import sys
import csv
import re
import argparse
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Merge existing samplesheet with new fetchngs FASTQs."
    )
    p.add_argument("--existing",  type=Path, default=None,
                   help="Existing samplesheet CSV to extend.")
    p.add_argument("--fastq_dir", type=Path, required=True,
                   help="Directory containing fetchngs FASTQ files.")
    p.add_argument("--output",    type=Path, required=True,
                   help="Output samplesheet CSV path.")
    p.add_argument("--dry_run",   action="store_true",
                   help="Print result without writing.")
    return p.parse_args()


def load_existing(path):
    if path is None or not path.exists():
        return [], set()
    with open(path) as f:
        rows = list(csv.DictReader(f))
    ids = {r["sample"] for r in rows}
    return rows, ids


def find_new_fastqs(fastq_dir, existing_ids):
    new_rows = []
    missing_r2 = []

    for fq1 in sorted(fastq_dir.glob("*_1.fastq.gz")):
        m = re.match(r".+_([A-Z]{3}\d+)_1\.fastq\.gz$", fq1.name)
        if not m:
            continue
        run_acc = m.group(1)
        if run_acc in existing_ids:
            continue
        fq2 = fq1.parent / fq1.name.replace("_1.fastq.gz", "_2.fastq.gz")
        if not fq2.exists():
            missing_r2.append(run_acc)
            continue
        new_rows.append({
            "sample":  run_acc,
            "fastq_1": str(fq1.resolve()),
            "fastq_2": str(fq2.resolve()),
        })

    return new_rows, missing_r2


def main():
    args = parse_args()

    if not args.fastq_dir.is_dir():
        sys.exit(f"ERROR: --fastq_dir {args.fastq_dir} does not exist")

    existing_rows, existing_ids = load_existing(args.existing)
    new_rows, missing_r2 = find_new_fastqs(args.fastq_dir, existing_ids)

    all_rows = existing_rows + new_rows

    if not all_rows:
        sys.exit("ERROR: No samples found — check --existing and --fastq_dir paths")

    print(f"\n── BUILD CC97 SAMPLESHEET ────────────────────────────────")
    print(f"  Existing samples : {len(existing_rows)}")
    print(f"  New ENA samples  : {len(new_rows)}")
    if missing_r2:
        print(f"  Skipped (no R2)  : {len(missing_r2)}")
        for acc in missing_r2:
            print(f"    {acc}")
    print(f"  Total            : {len(all_rows)}")

    if args.dry_run:
        print("\n  [dry run — not written]")
        for r in all_rows:
            print(f"  {r['sample']}")
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "fastq_1", "fastq_2"])
        w.writeheader()
        w.writerows(all_rows)

    print(f"\n  Written: {args.output}")


if __name__ == "__main__":
    main()
