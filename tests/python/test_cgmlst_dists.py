"""Tests for bin/staphit-cgmlst-dists."""
import importlib.machinery
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / 'bin' / 'staphit-cgmlst-dists'
_loader = importlib.machinery.SourceFileLoader('cgmlst_dists', str(BIN))
_spec = importlib.util.spec_from_loader('cgmlst_dists', _loader)
cd = importlib.util.module_from_spec(_spec)
_loader.exec_module(cd)


def test_distances_ignore_missing_and_read_inferred(tmp_path):
    profiles = tmp_path / 'results_alleles.tsv'
    profiles.write_text('FILE\tL1\tL2\tL3\tL4\n'
                        'A.fasta\t1\t2\t3\tINF-9\n'
                        'B.fasta\t1\t5\tLNF\tINF-9\n'      # L2 differs, L3 missing, same new allele at L4
                        'C\t2\tPLOT3\tASM\t7\n')            # L1 and L4 differ
    dists, stats = tmp_path / 'd.tsv', tmp_path / 's.tsv'
    cd.main([str(profiles), str(dists), str(stats)])
    rows = [line.split('\t') for line in dists.read_text().splitlines()]
    assert rows[0][1:] == ['A', 'B', 'C']
    m = {r[0]: [int(x) for x in r[1:]] for r in rows[1:]}
    assert m['A'] == [0, 1, 2] and m['B'] == [1, 0, 2] and m['C'] == [2, 2, 0]
    assert stats.read_text().splitlines()[1:] == ['A\t4\t0\t4', 'B\t3\t1\t4', 'C\t2\t2\t4']
