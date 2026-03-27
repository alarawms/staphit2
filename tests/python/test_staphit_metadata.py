"""Tests for staphit-metadata CLI tool."""
import subprocess
import csv
import json
import os
import tempfile
import pytest

TOOL = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'staphit-metadata')

@pytest.fixture
def tmpdir():
    with tempfile.TemporaryDirectory() as d:
        yield d

@pytest.fixture
def samplesheet(tmpdir):
    path = os.path.join(tmpdir, 'samplesheet.csv')
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['sample', 'fastq_1', 'fastq_2'])
        writer.writerow(['ID00001', '/data/ID00001_R1.fastq.gz', '/data/ID00001_R2.fastq.gz'])
        writer.writerow(['ID00002', '/data/ID00002_R1.fastq.gz', '/data/ID00002_R2.fastq.gz'])
        writer.writerow(['ID00003', '/data/ID00003_R1.fastq.gz', '/data/ID00003_R2.fastq.gz'])
    return path


class TestTemplate:
    def test_blank_template(self, tmpdir):
        """Generate blank metadata template with header only."""
        out = os.path.join(tmpdir, 'sample_metadata.csv')
        result = subprocess.run(
            ['python', TOOL, 'template', '--output-dir', tmpdir],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        assert os.path.exists(out)
        with open(out) as f:
            reader = csv.reader(f)
            header = next(reader)
            assert header[0] == 'sample_id'
            assert 'organism' in header
            assert 'collection_date' in header
            rows = list(reader)
            assert len(rows) == 0  # blank, no data rows

    def test_prefilled_from_samplesheet(self, tmpdir, samplesheet):
        """Generate template with sample_id column pre-filled from samplesheet."""
        out = os.path.join(tmpdir, 'sample_metadata.csv')
        result = subprocess.run(
            ['python', TOOL, 'template', '--samplesheet', samplesheet, '--output-dir', tmpdir],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        with open(out) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) == 3
            assert rows[0]['sample_id'] == 'ID00001'
            assert rows[1]['sample_id'] == 'ID00002'
            assert rows[0]['organism'] == 'Staphylococcus aureus'

    def test_antibiogram_template(self, tmpdir, samplesheet):
        """Generate antibiogram template alongside metadata template."""
        result = subprocess.run(
            ['python', TOOL, 'template', '--samplesheet', samplesheet,
             '--antibiogram', '--output-dir', tmpdir],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        abg_path = os.path.join(tmpdir, 'antibiogram.csv')
        assert os.path.exists(abg_path)
        with open(abg_path) as f:
            reader = csv.reader(f)
            header = next(reader)
            assert header[0] == 'sample_id'
            assert 'antibiotic' in header
            assert 'resistance_phenotype' in header


class TestValidate:
    def _write_metadata(self, tmpdir, rows):
        path = os.path.join(tmpdir, 'metadata.csv')
        with open(path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        return path

    def _write_antibiogram(self, tmpdir, rows):
        path = os.path.join(tmpdir, 'antibiogram.csv')
        with open(path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_valid_metadata(self, tmpdir):
        meta = self._write_metadata(tmpdir, [{'sample_id': 'ID00001', 'organism': 'Staphylococcus aureus', 'collection_date': '2024-03-15', 'geo_loc_country': 'Saudi Arabia', 'geo_loc_region': 'Riyadh', 'host': 'Homo sapiens', 'isolation_source': 'blood'}])
        result = subprocess.run(['python', TOOL, 'validate', '--metadata', meta], capture_output=True, text=True)
        assert result.returncode == 0

    def test_missing_required_warns(self, tmpdir):
        meta = self._write_metadata(tmpdir, [{'sample_id': 'ID00001', 'organism': 'Staphylococcus aureus', 'collection_date': '', 'geo_loc_country': 'Saudi Arabia', 'geo_loc_region': 'Riyadh', 'host': 'Homo sapiens', 'isolation_source': 'blood'}])
        result = subprocess.run(['python', TOOL, 'validate', '--metadata', meta], capture_output=True, text=True)
        assert result.returncode == 0
        assert 'WARNING' in result.stderr
        assert 'collection_date' in result.stderr

    def test_invalid_date_warns(self, tmpdir):
        meta = self._write_metadata(tmpdir, [{'sample_id': 'ID00001', 'organism': 'Staphylococcus aureus', 'collection_date': '15/03/2024', 'geo_loc_country': 'Saudi Arabia', 'geo_loc_region': 'Riyadh', 'host': 'Homo sapiens', 'isolation_source': 'blood'}])
        result = subprocess.run(['python', TOOL, 'validate', '--metadata', meta], capture_output=True, text=True)
        assert result.returncode == 0
        assert 'WARNING' in result.stderr
        assert 'date' in result.stderr.lower()

    def test_invalid_sex_warns(self, tmpdir):
        meta = self._write_metadata(tmpdir, [{'sample_id': 'ID00001', 'organism': 'Staphylococcus aureus', 'collection_date': '2024-03-15', 'geo_loc_country': 'Saudi Arabia', 'geo_loc_region': 'Riyadh', 'host': 'Homo sapiens', 'isolation_source': 'blood', 'host_sex': 'M'}])
        result = subprocess.run(['python', TOOL, 'validate', '--metadata', meta], capture_output=True, text=True)
        assert result.returncode == 0
        assert 'WARNING' in result.stderr
        assert 'host_sex' in result.stderr

    def test_samplesheet_crossref(self, tmpdir, samplesheet):
        meta = self._write_metadata(tmpdir, [{'sample_id': 'ID99999', 'organism': 'Staphylococcus aureus', 'collection_date': '2024-03-15', 'geo_loc_country': 'Saudi Arabia', 'geo_loc_region': 'Riyadh', 'host': 'Homo sapiens', 'isolation_source': 'blood'}])
        result = subprocess.run(['python', TOOL, 'validate', '--metadata', meta, '--samplesheet', samplesheet], capture_output=True, text=True)
        assert result.returncode == 0
        assert 'WARNING' in result.stderr
        assert 'ID99999' in result.stderr

    def test_antibiogram_validation(self, tmpdir):
        meta = self._write_metadata(tmpdir, [{'sample_id': 'ID00001', 'organism': 'Staphylococcus aureus', 'collection_date': '2024-03-15', 'geo_loc_country': 'Saudi Arabia', 'geo_loc_region': 'Riyadh', 'host': 'Homo sapiens', 'isolation_source': 'blood'}])
        abg = self._write_antibiogram(tmpdir, [{'sample_id': 'ID00001', 'antibiotic': 'oxacillin', 'resistance_phenotype': 'R', 'measurement': '4', 'measurement_sign': '=', 'measurement_units': 'mg/L', 'laboratory_typing_method': 'MIC', 'testing_standard': 'CLSI'}])
        result = subprocess.run(['python', TOOL, 'validate', '--metadata', meta, '--antibiogram', abg], capture_output=True, text=True)
        assert result.returncode == 0

    def test_antibiogram_orphan_warns(self, tmpdir):
        meta = self._write_metadata(tmpdir, [{'sample_id': 'ID00001', 'organism': 'Staphylococcus aureus', 'collection_date': '2024-03-15', 'geo_loc_country': 'Saudi Arabia', 'geo_loc_region': 'Riyadh', 'host': 'Homo sapiens', 'isolation_source': 'blood'}])
        abg = self._write_antibiogram(tmpdir, [{'sample_id': 'ID00099', 'antibiotic': 'oxacillin', 'resistance_phenotype': 'R', 'measurement': '4', 'measurement_sign': '=', 'measurement_units': 'mg/L', 'laboratory_typing_method': 'MIC', 'testing_standard': 'CLSI'}])
        result = subprocess.run(['python', TOOL, 'validate', '--metadata', meta, '--antibiogram', abg], capture_output=True, text=True)
        assert result.returncode == 0
        assert 'WARNING' in result.stderr
        assert 'ID00099' in result.stderr

    def test_invalid_sir_warns(self, tmpdir):
        meta = self._write_metadata(tmpdir, [{'sample_id': 'ID00001', 'organism': 'Staphylococcus aureus', 'collection_date': '2024-03-15', 'geo_loc_country': 'Saudi Arabia', 'geo_loc_region': 'Riyadh', 'host': 'Homo sapiens', 'isolation_source': 'blood'}])
        abg = self._write_antibiogram(tmpdir, [{'sample_id': 'ID00001', 'antibiotic': 'oxacillin', 'resistance_phenotype': 'X', 'measurement': '4', 'measurement_sign': '=', 'measurement_units': 'mg/L', 'laboratory_typing_method': 'MIC', 'testing_standard': 'CLSI'}])
        result = subprocess.run(['python', TOOL, 'validate', '--metadata', meta, '--antibiogram', abg], capture_output=True, text=True)
        assert result.returncode == 0
        assert 'WARNING' in result.stderr
        assert 'resistance_phenotype' in result.stderr

    def test_malformed_csv_fails(self, tmpdir):
        bad_path = os.path.join(tmpdir, 'bad.csv')
        with open(bad_path, 'w') as f:
            f.write('not,a,valid\x00csv\nwith\x00nulls')
        result = subprocess.run(['python', TOOL, 'validate', '--metadata', bad_path], capture_output=True, text=True)
        assert result.returncode == 1


class TestNormalize:
    def _write_metadata(self, tmpdir, rows):
        path = os.path.join(tmpdir, 'metadata.csv')
        with open(path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        return path

    def _write_antibiogram(self, tmpdir, rows):
        path = os.path.join(tmpdir, 'antibiogram.csv')
        with open(path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_normalize_metadata_only(self, tmpdir):
        meta = self._write_metadata(tmpdir, [{'sample_id': 'ID00001', 'organism': 'Staphylococcus aureus', 'collection_date': '2024-03-15', 'geo_loc_country': 'Saudi Arabia', 'geo_loc_region': 'Riyadh', 'host': 'Homo sapiens', 'isolation_source': 'blood', 'infection_origin': 'hospital'}])
        out = os.path.join(tmpdir, 'metadata.json')
        result = subprocess.run(['python', TOOL, 'normalize', '--metadata', meta, '-o', out], capture_output=True, text=True)
        assert result.returncode == 0
        with open(out) as f:
            data = json.load(f)
        assert isinstance(data, list)
        assert len(data) == 1
        assert data[0]['sample_id'] == 'ID00001'
        assert data[0]['collection_date'] == '2024-03-15'
        assert data[0]['geo_loc_country'] == 'Saudi Arabia'
        assert 'antibiogram' not in data[0] or data[0]['antibiogram'] == []

    def test_normalize_with_antibiogram(self, tmpdir):
        meta = self._write_metadata(tmpdir, [{'sample_id': 'ID00001', 'organism': 'Staphylococcus aureus', 'collection_date': '2024-03-15', 'geo_loc_country': 'Saudi Arabia', 'geo_loc_region': 'Riyadh', 'host': 'Homo sapiens', 'isolation_source': 'blood'}])
        abg = self._write_antibiogram(tmpdir, [
            {'sample_id': 'ID00001', 'antibiotic': 'oxacillin', 'resistance_phenotype': 'R', 'measurement': '4', 'measurement_sign': '=', 'measurement_units': 'mg/L', 'laboratory_typing_method': 'MIC', 'testing_standard': 'CLSI', 'testing_standard_version': '', 'platform': ''},
            {'sample_id': 'ID00001', 'antibiotic': 'vancomycin', 'resistance_phenotype': 'S', 'measurement': '1', 'measurement_sign': '=', 'measurement_units': 'mg/L', 'laboratory_typing_method': 'MIC', 'testing_standard': 'CLSI', 'testing_standard_version': '', 'platform': ''},
        ])
        out = os.path.join(tmpdir, 'metadata.json')
        result = subprocess.run(['python', TOOL, 'normalize', '--metadata', meta, '--antibiogram', abg, '-o', out], capture_output=True, text=True)
        assert result.returncode == 0
        with open(out) as f:
            data = json.load(f)
        assert len(data[0]['antibiogram']) == 2
        assert data[0]['antibiogram'][0]['antibiotic'] == 'oxacillin'
        assert data[0]['antibiogram'][0]['sir'] == 'R'
        assert data[0]['antibiogram'][1]['antibiotic'] == 'vancomycin'

    def test_normalize_uses_run_id_key(self, tmpdir):
        meta = self._write_metadata(tmpdir, [{'sample_id': 'ID00001', 'organism': 'Staphylococcus aureus', 'collection_date': '2024-03-15', 'geo_loc_country': 'Saudi Arabia', 'geo_loc_region': 'Riyadh', 'host': 'Homo sapiens', 'isolation_source': 'blood'}])
        out = os.path.join(tmpdir, 'metadata.json')
        subprocess.run(['python', TOOL, 'normalize', '--metadata', meta, '-o', out], capture_output=True, text=True)
        with open(out) as f:
            data = json.load(f)
        assert data[0]['run_id'] == 'ID00001'
