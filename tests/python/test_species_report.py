"""Tests for staphit-species-report."""
import csv
import os
import subprocess
import sys

import pytest

TOOL_PATH = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'staphit-species-report')


def _write_fastani(tmp_path, rows):
    path = tmp_path / 'fastani.tsv'
    with open(path, 'w') as fh:
        for row in rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    return str(path)


def _write_mash_screen(tmp_path, rows):
    path = tmp_path / 'mash_screen.tsv'
    with open(path, 'w') as fh:
        for row in rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    return str(path)


class TestSpeciesReport:

    def test_all_pass(self, tmp_path):
        fastani = _write_fastani(tmp_path, [
            ['S1.fasta', 'ref.fasta', 98.5, 800, 900],
            ['S2.fasta', 'ref.fasta', 96.2, 750, 900],
        ])
        outdir = str(tmp_path / 'out')
        result = subprocess.run(
            [sys.executable, TOOL_PATH, '--fastani', fastani,
             '--threshold', '95.0', '--outdir', outdir],
            capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        confirmed = os.path.join(outdir, 'species_confirmed.tsv')
        assert os.path.exists(confirmed)
        with open(confirmed) as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        assert len(rows) == 2
        assert all(r['status'] == 'PASS' for r in rows)

    def test_one_fail(self, tmp_path):
        fastani = _write_fastani(tmp_path, [
            ['S1.fasta', 'ref.fasta', 98.5, 800, 900],
            ['S2.fasta', 'ref.fasta', 80.1, 400, 900],
        ])
        mash = _write_mash_screen(tmp_path, [
            [0.95, '800/1000', 50, 'GCF_000007645.1', 'S. epidermidis strain ATCC 12228'],
        ])
        outdir = str(tmp_path / 'out')
        result = subprocess.run(
            [sys.executable, TOOL_PATH, '--fastani', fastani,
             '--mash-screen', mash, '--threshold', '95.0', '--outdir', outdir],
            capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        excluded = os.path.join(outdir, 'species_excluded.tsv')
        assert os.path.exists(excluded)
        with open(excluded) as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        assert len(rows) == 1
        assert 'epidermidis' in rows[0]['identified_species']

    def test_sample_id_extraction(self, tmp_path):
        fastani = _write_fastani(tmp_path, [
            ['ID00001.scaffolds.fasta', 'ref.fasta', 99.1, 850, 900],
        ])
        outdir = str(tmp_path / 'out')
        result = subprocess.run(
            [sys.executable, TOOL_PATH, '--fastani', fastani,
             '--threshold', '95.0', '--outdir', outdir],
            capture_output=True, text=True)
        assert result.returncode == 0
        with open(os.path.join(outdir, 'species_confirmed.tsv')) as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        assert rows[0]['sample_id'] == 'ID00001'
