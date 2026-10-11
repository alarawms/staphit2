"""Tests for bin/staphit-mixed (B2)."""
import importlib.machinery
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / 'bin' / 'staphit-mixed'
_loader = importlib.machinery.SourceFileLoader('mixed', str(BIN))
_spec = importlib.util.spec_from_loader('mixed', _loader)
mx = importlib.util.module_from_spec(_spec)
_loader.exec_module(mx)


def test_count_mixed():
    lines = [
        '50,0',      # pure
        '40,10',     # 20% minor: mixed
        '95,5',      # 5%: below --min-maf
        '10,5',      # depth 15: not assessed
        '30,0,3',    # 3 reads: below --min-reads
        '.',         # no coverage
        '20,0,10',   # minor allele is the 2nd ALT: mixed (33%)
    ]
    assessed, mixed, frac = mx.count_mixed(lines)
    assert (assessed, mixed) == (5, 2)
    assert abs(frac - (0.2 + 1 / 3) / 2) < 1e-9
