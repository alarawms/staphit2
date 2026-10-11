"""Tests for bin/staphit-calibrate-thresholds (B11)."""
import importlib.machinery
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / 'bin' / 'staphit-calibrate-thresholds'
_loader = importlib.machinery.SourceFileLoader('calibrate', str(BIN))
_spec = importlib.util.spec_from_loader('calibrate', _loader)
cal = importlib.util.module_from_spec(_spec)
_loader.exec_module(cal)


def test_best_threshold_separates_outbreaks(tmp_path, capsys):
    ids = ['A1', 'A2', 'B1', 'B2', 'S']
    d = {('A1', 'A2'): 3, ('B1', 'B2'): 8}   # within outbreaks; everything else 40
    m = tmp_path / 'd.tsv'
    m.write_text('x\t' + '\t'.join(ids) + '\n' + ''.join(
        a + '\t' + '\t'.join('0' if a == b else str(d.get((a, b), d.get((b, a), 40))) for b in ids) + '\n' for a in ids))
    lab = tmp_path / 'o.tsv'
    lab.write_text('sample_id\toutbreak\nA1\tob1\nA2\tob1\nB1\tob2\nB2\tob2\nS\tsporadic\nNOT_IN_MATRIX\tob1\n')
    linked, unlinked = cal.labelled_pairs(cal.read_matrix(str(m)), {'A1': 'o1', 'A2': 'o1', 'S': ''})
    assert linked == [3] and unlinked == [40, 40]
    cal.main([str(m), str(lab), '--max', '50'])
    assert 'best threshold 8 (sensitivity 1.00, specificity 1.00)' in capsys.readouterr().err
