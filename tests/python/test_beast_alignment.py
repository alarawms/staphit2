"""Tests for bin/beast_alignment.py."""
import importlib.machinery
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / 'bin' / 'beast_alignment.py'
_loader = importlib.machinery.SourceFileLoader('beast_alignment', str(BIN))
_spec = importlib.util.spec_from_loader('beast_alignment', _loader)
ba = importlib.util.module_from_spec(_spec)
_loader.exec_module(ba)


def test_split_masked_alignment(tmp_path):
    aln = tmp_path / 'masked.full.aln'
    # cols: A-const, C-const, SNP(G/T), N in one sample, recombination-masked (-), G-const, T-const
    aln.write_text('>Reference\nACGNAGT\n>s1\nACGAAGT\n>s2\nACTA-GT\n>s3\nACGNAGT\n')
    snps, const = tmp_path / 'snps.fasta', tmp_path / 'constant.txt'
    ba.main([str(aln), str(snps), str(const), '--drop', 'Reference'])
    assert snps.read_text() == '>s1\nG\n>s2\nT\n>s3\nG\n'
    assert const.read_text() == '1 1 1 1\n'          # A C G T
