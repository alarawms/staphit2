"""Tests for bin/staphit-truncated (B6)."""
import importlib.machinery
import importlib.util
import io
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_loader = importlib.machinery.SourceFileLoader('truncated', str(ROOT / 'bin' / 'staphit-truncated'))
_spec = importlib.util.spec_from_loader('truncated', _loader)
tr = importlib.util.module_from_spec(_spec)
_loader.exec_module(tr)


def test_count_truncated_best_hit_only():
    hits = [
        'p1\tSAOUHSC_1\t99.5\t300\t310',   # full length
        'p2\tSAOUHSC_2\t99.0\t120\t400',   # truncated
        'p2\tSAOUHSC_9\t40.0\t120\t130',   # not the best hit: ignored
        'p3\tSAOUHSC_3\t60.0\t50\t400',    # other gene family: not counted
    ]
    assert tr.count_truncated(hits) == (2, 1)


def test_reference_proteins_from_shipped_gbk():
    out = io.StringIO()
    n = tr.reference_proteins(str(ROOT / 'assets' / 'references' / 'NCTC8325.gbk'), out)
    assert n > 2500 and out.getvalue().startswith('>SAOUHSC_')
