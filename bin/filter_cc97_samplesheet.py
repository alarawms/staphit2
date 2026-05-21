#!/usr/bin/env python3
"""
After a pipeline run, read combined_summary.tsv and filter the samplesheet
to only samples confirmed as CC97 (ST97/ST1153 and known SLVs) by MLST.
Writes a new samplesheet and reports what was kept/dropped.

Usage:
    python3 bin/filter_cc97_samplesheet.py \
        <results_dir> \
        <input_samplesheet.csv> \
        <output_samplesheet.csv>
"""

import sys, csv
from pathlib import Path
from collections import Counter

CC97_STS = {'97', '1153', '1465', '3187', '8064', '2996', '5435',
            '3528', '2458', '3009', '7570', '3636', '5264'}

def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: filter_cc97_samplesheet.py <results_dir> <input.csv> <output.csv>")

    results_dir  = Path(sys.argv[1])
    input_sheet  = Path(sys.argv[2])
    output_sheet = Path(sys.argv[3])

    summary_file = results_dir / "summary" / "combined_summary.tsv"
    if not summary_file.exists():
        sys.exit(f"ERROR: {summary_file} not found — run pipeline first")

    # Load MLST results
    mlst = {}
    with open(summary_file) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            st = r.get('mlst_st', '').strip().lstrip('ST')
            mlst[r['sample_id']] = st

    # Load samplesheet
    with open(input_sheet) as f:
        rows = list(csv.DictReader(f))

    kept, dropped_non_cc97, dropped_no_mlst = [], [], []
    for r in rows:
        sid = r['sample']
        st  = mlst.get(sid, None)
        if st is None:
            # Not in summary — failed QC before MLST; keep local IDs (already confirmed)
            if sid.startswith('ID'):
                kept.append(r)
            else:
                dropped_no_mlst.append((sid, 'failed QC / no MLST'))
        elif st in CC97_STS:
            kept.append(r)
        elif st in ('-', ''):
            # MLST ran but couldn't assign — keep local, drop public
            if sid.startswith('ID'):
                kept.append(r)
            else:
                dropped_no_mlst.append((sid, f'MLST unresolved (ST-)'))
        else:
            dropped_non_cc97.append((sid, f'ST{st}'))

    # Write filtered samplesheet
    with open(output_sheet, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(kept)

    print(f"\n── CC97 FILTER RESULT ────────────────────────────────────")
    print(f"  Input samples    : {len(rows)}")
    print(f"  Kept (CC97)      : {len(kept)}")
    print(f"  Dropped non-CC97 : {len(dropped_non_cc97)}")
    print(f"  Dropped no MLST  : {len(dropped_no_mlst)}")

    if dropped_non_cc97:
        print(f"\n  Non-CC97 excluded:")
        for sid, st in dropped_non_cc97:
            print(f"    {sid}  {st}")

    if dropped_no_mlst:
        print(f"\n  No MLST / QC-failed excluded:")
        for sid, reason in dropped_no_mlst:
            print(f"    {sid}  ({reason})")

    st_dist = Counter(mlst.get(r['sample'], '?').lstrip('ST') for r in kept)
    print(f"\n  ST distribution in kept samples:")
    for st, n in sorted(st_dist.items(), key=lambda x: -x[1]):
        print(f"    ST{st}: {n}")

    print(f"\n  Written: {output_sheet}")

if __name__ == '__main__':
    main()
