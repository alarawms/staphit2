"""Tests for bin/staphit-watch."""
import csv
import importlib.machinery
import importlib.util
import os
import tempfile

import pytest

# ---------------------------------------------------------------------------
# Import the script as a module (no .py extension)
# ---------------------------------------------------------------------------

TOOL_PATH = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'staphit-watch')

loader = importlib.machinery.SourceFileLoader('staphit_watch', TOOL_PATH)
spec = importlib.util.spec_from_loader('staphit_watch', loader)
watch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(watch)


# ===========================================================================
# Helpers
# ===========================================================================

def _touch(path):
    """Create an empty file, making parent dirs as needed."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        pass


def _write_samplesheet(path, rows):
    """Write a samplesheet CSV with columns sample,fastq_1,fastq_2."""
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['sample', 'fastq_1', 'fastq_2'])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_transmission_pairs(path, rows):
    """Write a transmission_pairs TSV with columns sample_1,sample_2,snp_distance."""
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['sample_1', 'sample_2', 'snp_distance'],
                                delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ===========================================================================
# TestScanFastqs
# ===========================================================================

class TestScanFastqs:
    """Tests for scan_for_new_fastqs()."""

    def test_finds_paired_fastqs(self, tmp_path):
        """Two complete pairs and one orphan R1 → returns only the 2 pairs."""
        _touch(str(tmp_path / 'SAMPLE1_R1.fastq.gz'))
        _touch(str(tmp_path / 'SAMPLE1_R2.fastq.gz'))
        _touch(str(tmp_path / 'SAMPLE2_R1.fastq.gz'))
        _touch(str(tmp_path / 'SAMPLE2_R2.fastq.gz'))
        _touch(str(tmp_path / 'ORPHAN_R1.fastq.gz'))  # no R2

        result = watch.scan_for_new_fastqs(str(tmp_path))

        assert len(result) == 2
        sample_ids = [sid for sid, _, _ in result]
        assert 'SAMPLE1' in sample_ids
        assert 'SAMPLE2' in sample_ids
        assert 'ORPHAN' not in sample_ids

    def test_empty_dir(self, tmp_path):
        """Empty directory returns empty list."""
        result = watch.scan_for_new_fastqs(str(tmp_path))
        assert result == []

    def test_nested_dirs(self, tmp_path):
        """FASTQs in subdirectories are still found."""
        subdir = tmp_path / 'run1' / 'fastqs'
        _touch(str(subdir / 'DEEP_R1.fastq.gz'))
        _touch(str(subdir / 'DEEP_R2.fastq.gz'))

        result = watch.scan_for_new_fastqs(str(tmp_path))

        assert len(result) == 1
        assert result[0][0] == 'DEEP'


# ===========================================================================
# TestExtractSampleId
# ===========================================================================

class TestExtractSampleId:
    """Tests for extract_sample_id()."""

    def test_standard_format(self):
        assert watch.extract_sample_id('ID00001_R1.fastq.gz') == 'ID00001'

    def test_standard_r2(self):
        assert watch.extract_sample_id('ID00001_R2.fastq.gz') == 'ID00001'

    def test_underscore_1_format(self):
        assert watch.extract_sample_id('SAMPLE1_1.fastq.gz') == 'SAMPLE1'

    def test_complex_illumina_name(self):
        fname = 'M-21-3999_1_UDI001-UDI001_L001_R1_001.fastq.gz'
        sid = watch.extract_sample_id(fname)
        assert sid == 'M-21-3999_1_UDI001-UDI001_L001'

    def test_custom_pattern(self):
        fname = 'M-21-3999_1_UDI001-UDI001_L001_R1_001.fastq.gz'
        sid = watch.extract_sample_id(fname, pattern=r'^(M-\d+-\d+)')
        assert sid == 'M-21-3999'


# ===========================================================================
# TestUpdateSamplesheet
# ===========================================================================

class TestUpdateSamplesheet:
    """Tests for update_samplesheet()."""

    def test_appends_new_samples(self, tmp_path):
        """One existing + one new → samplesheet has 2 rows."""
        sheet = str(tmp_path / 'samplesheet.csv')
        _write_samplesheet(sheet, [
            {'sample': 'EXISTING', 'fastq_1': '/r1', 'fastq_2': '/r2'},
        ])

        new_pairs = [('NEW_SAMPLE', '/new_r1', '/new_r2')]
        added = watch.update_samplesheet(sheet, new_pairs)

        assert added == 1
        with open(sheet) as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)
        assert len(rows) == 2
        assert rows[0]['sample'] == 'EXISTING'
        assert rows[1]['sample'] == 'NEW_SAMPLE'

    def test_skips_existing(self, tmp_path):
        """Same sample_id already present → no change, returns 0."""
        sheet = str(tmp_path / 'samplesheet.csv')
        _write_samplesheet(sheet, [
            {'sample': 'EXISTING', 'fastq_1': '/r1', 'fastq_2': '/r2'},
        ])

        new_pairs = [('EXISTING', '/r1', '/r2')]
        added = watch.update_samplesheet(sheet, new_pairs)

        assert added == 0

    def test_creates_new_samplesheet(self, tmp_path):
        """If samplesheet does not exist, create it."""
        sheet = str(tmp_path / 'samplesheet.csv')
        new_pairs = [('FIRST', '/r1', '/r2')]
        added = watch.update_samplesheet(sheet, new_pairs)

        assert added == 1
        assert os.path.exists(sheet)
        with open(sheet) as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)
        assert len(rows) == 1
        assert rows[0]['sample'] == 'FIRST'


# ===========================================================================
# TestDetectNewClusters
# ===========================================================================

class TestDetectNewClusters:
    """Tests for detect_new_clusters()."""

    def test_finds_new_pairs(self, tmp_path):
        """Previous empty, current has 1 pair → returns 1."""
        prev = str(tmp_path / 'prev.tsv')
        curr = str(tmp_path / 'curr.tsv')

        _write_transmission_pairs(prev, [])
        _write_transmission_pairs(curr, [
            {'sample_1': 'A', 'sample_2': 'B', 'snp_distance': '5'},
        ])

        result = watch.detect_new_clusters(prev, curr)
        assert len(result) == 1
        assert result[0]['sample_1'] == 'A'
        assert result[0]['sample_2'] == 'B'

    def test_no_new(self, tmp_path):
        """Both have the same pair → returns empty list."""
        prev = str(tmp_path / 'prev.tsv')
        curr = str(tmp_path / 'curr.tsv')

        pair = {'sample_1': 'A', 'sample_2': 'B', 'snp_distance': '5'}
        _write_transmission_pairs(prev, [pair])
        _write_transmission_pairs(curr, [pair])

        result = watch.detect_new_clusters(prev, curr)
        assert result == []

    def test_prev_does_not_exist(self, tmp_path):
        """If prev file doesn't exist, all curr pairs are new."""
        curr = str(tmp_path / 'curr.tsv')
        _write_transmission_pairs(curr, [
            {'sample_1': 'X', 'sample_2': 'Y', 'snp_distance': '3'},
        ])

        result = watch.detect_new_clusters(None, curr)
        assert len(result) == 1
