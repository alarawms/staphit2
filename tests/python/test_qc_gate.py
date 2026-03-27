"""Tests for bin/staphit-qc."""
import csv
import importlib.machinery
import importlib.util
import os
import tempfile

import pytest

# ---------------------------------------------------------------------------
# Import the script as a module (no .py extension)
# ---------------------------------------------------------------------------

TOOL_PATH = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'staphit-qc')

loader = importlib.machinery.SourceFileLoader('staphit_qc', TOOL_PATH)
spec = importlib.util.spec_from_loader('staphit_qc', loader)
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)


# ===========================================================================
# Helpers
# ===========================================================================

CHECKM2_HEADER = (
    "Name\tCompleteness\tContamination\tGenome_Size\tGC_Content\tN50\n"
)


def _write_report(path, rows):
    """Write a minimal quality_report.tsv."""
    with open(path, 'w') as fh:
        fh.write(CHECKM2_HEADER)
        for row in rows:
            fh.write('\t'.join(str(v) for v in row) + '\n')


# ===========================================================================
# TestParseCheckm2
# ===========================================================================

class TestParseCheckm2:

    def test_parses_report(self, tmp_path):
        rpt = tmp_path / 'quality_report.tsv'
        _write_report(str(rpt), [
            ['SAMPLE1', '99.5', '0.3', '2800000', '32.8', '150000'],
        ])
        results = qc.parse_checkm2_report(str(rpt))
        assert len(results) == 1
        assert results[0]['sample_id'] == 'SAMPLE1'
        assert results[0]['completeness'] == 99.5
        assert results[0]['contamination'] == 0.3
        assert results[0]['genome_size'] == 2800000
        assert results[0]['gc'] == 32.8
        assert results[0]['n50'] == 150000

    def test_strips_scaffolds_suffix(self, tmp_path):
        rpt = tmp_path / 'quality_report.tsv'
        _write_report(str(rpt), [
            ['SAMPLE1.scaffolds', '98.0', '1.0', '2800000', '32.8', '100000'],
        ])
        results = qc.parse_checkm2_report(str(rpt))
        assert results[0]['sample_id'] == 'SAMPLE1'

    def test_strips_fasta_suffix(self, tmp_path):
        rpt = tmp_path / 'quality_report.tsv'
        _write_report(str(rpt), [
            ['SAMPLE2.fasta', '97.0', '2.0', '2700000', '33.0', '90000'],
        ])
        results = qc.parse_checkm2_report(str(rpt))
        assert results[0]['sample_id'] == 'SAMPLE2'

    def test_multiple_samples(self, tmp_path):
        rpt = tmp_path / 'quality_report.tsv'
        _write_report(str(rpt), [
            ['S1', '99.0', '0.1', '2800000', '32.8', '150000'],
            ['S2.scaffolds', '85.0', '10.0', '2600000', '33.0', '50000'],
        ])
        results = qc.parse_checkm2_report(str(rpt))
        assert len(results) == 2
        assert results[0]['sample_id'] == 'S1'
        assert results[1]['sample_id'] == 'S2'


# ===========================================================================
# TestAssessQuality
# ===========================================================================

class TestAssessQuality:

    def _sample(self, completeness=99.0, contamination=0.5):
        return {
            'sample_id': 'TEST',
            'completeness': completeness,
            'contamination': contamination,
            'genome_size': 2800000,
            'gc': 32.8,
            'n50': 150000,
        }

    def test_pass(self):
        result = qc.assess_quality(self._sample())
        assert result['pass'] is True
        assert result['reason'] == 'OK'

    def test_fail_low_completeness(self):
        result = qc.assess_quality(self._sample(completeness=80.0))
        assert result['pass'] is False
        assert 'completeness' in result['reason']

    def test_fail_high_contamination(self):
        result = qc.assess_quality(self._sample(contamination=12.0))
        assert result['pass'] is False
        assert 'contamination' in result['reason']

    def test_fail_both(self):
        result = qc.assess_quality(
            self._sample(completeness=50.0, contamination=20.0)
        )
        assert result['pass'] is False
        assert 'completeness' in result['reason']
        assert 'contamination' in result['reason']

    def test_custom_thresholds(self):
        # Should pass with relaxed thresholds
        result = qc.assess_quality(
            self._sample(completeness=80.0, contamination=8.0),
            min_completeness=70,
            max_contamination=10,
        )
        assert result['pass'] is True


# ===========================================================================
# TestCLI
# ===========================================================================

class TestCLI:

    def test_produces_all_output_files(self, tmp_path):
        # Create a checkm2 report structure
        rpt_dir = tmp_path / 'checkm2' / 'SAMPLE1' / 'checkm2_out'
        rpt_dir.mkdir(parents=True)
        _write_report(str(rpt_dir / 'quality_report.tsv'), [
            ['SAMPLE1.scaffolds', '99.0', '0.2', '2800000', '32.8', '150000'],
            ['SAMPLE2.fasta', '85.0', '1.0', '2700000', '33.0', '90000'],
        ])

        outdir = tmp_path / 'out'
        outdir.mkdir()

        qc.main([
            '--checkm2-dir', str(tmp_path / 'checkm2'),
            '--min-completeness', '90',
            '--max-contamination', '5',
            '-o', str(outdir),
        ])

        assert (outdir / 'qc_report.tsv').exists()
        assert (outdir / 'passed_samples.txt').exists()
        assert (outdir / 'failed_samples.txt').exists()

    def test_correctly_separates_passed_failed(self, tmp_path):
        rpt_dir = tmp_path / 'checkm2_out'
        rpt_dir.mkdir(parents=True)
        _write_report(str(rpt_dir / 'quality_report.tsv'), [
            ['GOOD', '99.0', '0.2', '2800000', '32.8', '150000'],
            ['BAD', '70.0', '15.0', '2600000', '33.0', '40000'],
        ])

        outdir = tmp_path / 'out'
        outdir.mkdir()

        qc.main([
            '--checkm2-dir', str(tmp_path),
            '--min-completeness', '90',
            '--max-contamination', '5',
            '-o', str(outdir),
        ])

        passed = (outdir / 'passed_samples.txt').read_text().strip().splitlines()
        failed = (outdir / 'failed_samples.txt').read_text().strip().splitlines()

        assert passed == ['GOOD']
        assert failed == ['BAD']

        # Verify qc_report.tsv content
        with open(outdir / 'qc_report.tsv') as fh:
            reader = list(csv.DictReader(fh, delimiter='\t'))
        good_row = [r for r in reader if r['sample_id'] == 'GOOD'][0]
        bad_row = [r for r in reader if r['sample_id'] == 'BAD'][0]
        assert good_row['pass'] == 'PASS'
        assert bad_row['pass'] == 'FAIL'
