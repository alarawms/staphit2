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
            # sample_id, identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment
            ['S2', 0.95, '800/1000', 50, 0, 'GCF_000007645.1_ASM764v1_genomic.fna.gz',
             '[1 seqs] NC_004461.1 Staphylococcus epidermidis ATCC 12228'],
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

    def test_second_species_and_contamination_flag(self, tmp_path):
        fastani = _write_fastani(tmp_path, [
            ['CLEAN.scaffolds.fasta', 'ref.fasta', 99.1, 900, 950],
            ['MIXED.scaffolds.fasta', 'ref.fasta', 99.0, 880, 950],
        ])
        aureus = '[2 seqs] NZ_CP000253.1 Staphylococcus aureus subsp. aureus NCTC 8325'
        mash = _write_mash_screen(tmp_path, [
            # winner-take-all can hand the top hit to an unnamed genome; it must not hide the species
            ['CLEAN', 0.9996, '991/1000', 60, 0, 'GCF_sp_genomic.fna.gz', '[1 seqs] NZ_JA1.1 Staphylococcus sp. isolate 7'],
            ['CLEAN', 0.9995, '990/1000', 60, 0, 'GCF_000013425.1_genomic.fna.gz', aureus],
            # phages match every S. aureus genome and must never count as a second species
            ['CLEAN', 0.9900, '813/1000', 5, 0, 'GCF_x_genomic.fna.gz', '[1 seqs] NC_1.1 Staphylococcus phage 80alpha'],
            ['CLEAN', 0.9800, '700/1000', 3, 0, 'GCF_y_genomic.fna.gz', '[1 seqs] NC_2.1 Staphylococcus prophage phiPV83'],
            ['MIXED', 0.9994, '987/1000', 60, 0, 'GCF_000013425.1_genomic.fna.gz', aureus],
            ['MIXED', 0.9273, '205/1000', 2, 0, 'GCF_z_genomic.fna.gz', '[1 seqs] NZ_CP3.1 Lysinibacillus fusiformis strain X'],
        ])
        outdir = str(tmp_path / 'out')
        result = subprocess.run(
            [sys.executable, TOOL_PATH, '--fastani', fastani, '--mash-screen', mash,
             '--threshold', '95.0', '--outdir', outdir],
            capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        with open(os.path.join(outdir, 'species_confirmed.tsv')) as fh:
            rows = {r['sample_id']: r for r in csv.DictReader(fh, delimiter='\t')}
        assert rows['CLEAN']['mash_species'] == 'Staphylococcus aureus'
        assert rows['CLEAN']['second_species'] == ''
        assert rows['CLEAN']['flag'] == ''
        assert rows['MIXED']['second_species'] == 'Lysinibacillus fusiformis'
        assert rows['MIXED']['second_shared_hashes'] == '205/1000'
        assert rows['MIXED']['flag'].startswith('possible_contamination: Lysinibacillus fusiformis')
