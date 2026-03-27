"""Tests for staphit-metadata convert subcommand (Vitek 2 CSV conversion)."""
import csv
import os
import subprocess
import tempfile

import pytest

TOOL = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'staphit-metadata')

# Import the helper directly for unit tests (file has no .py extension)
import importlib.util
import importlib.machinery
_tool_path = os.path.abspath(TOOL)
_loader = importlib.machinery.SourceFileLoader('staphit_metadata', _tool_path)
spec = importlib.util.spec_from_loader('staphit_metadata', _loader)
mod = importlib.util.module_from_spec(spec)
# Prevent main() from running during import
import unittest.mock as _mock
with _mock.patch('sys.argv', ['staphit-metadata', '--help']):
    try:
        spec.loader.exec_module(mod)
    except SystemExit:
        pass
_parse_vitek_cell = mod._parse_vitek_cell

CONVERT_OUTPUT_COLUMNS = [
    'sample_id', 'antibiotic', 'resistance_phenotype', 'measurement',
    'measurement_sign', 'measurement_units', 'laboratory_typing_method',
    'testing_standard', 'testing_standard_version', 'platform',
    'aes_modified', 'deduced',
]


@pytest.fixture
def tmpdir():
    with tempfile.TemporaryDirectory() as d:
        yield d


def _write_vitek_tsv(tmpdir, rows, filename='vitek.tsv'):
    """Write a wide-format Vitek TSV. rows is a list of dicts keyed by column."""
    path = os.path.join(tmpdir, filename)
    fieldnames = list(rows[0].keys())
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    return path


# --- Unit tests for _parse_vitek_cell ---

class TestParseVitekCell:
    def test_resistant_with_gte_sign(self):
        result = _parse_vitek_cell('R:>= 4')
        assert result['sir'] == 'R'
        assert result['measurement'] == '4'
        assert result['sign'] == '>='
        assert result['aes_modified'] is False
        assert result['deduced'] is False

    def test_susceptible_with_lte_sign(self):
        result = _parse_vitek_cell('S:<= 0.25')
        assert result['sir'] == 'S'
        assert result['measurement'] == '0.25'
        assert result['sign'] == '<='

    def test_resistant_pos_screen(self):
        result = _parse_vitek_cell('R:POS')
        assert result['sir'] == 'R'
        assert result['measurement'] == 'POS'
        assert result['sign'] == '='
        assert result['aes_modified'] is False

    def test_susceptible_neg_screen(self):
        result = _parse_vitek_cell('S:NEG')
        assert result['sir'] == 'S'
        assert result['measurement'] == 'NEG'
        assert result['sign'] == '='

    def test_intermediate_with_aes_asterisk(self):
        result = _parse_vitek_cell('I:64*')
        assert result['sir'] == 'I'
        assert result['measurement'] == '64'
        assert result['sign'] == '='
        assert result['aes_modified'] is True

    def test_resistant_with_aes_asterisk(self):
        result = _parse_vitek_cell('R:0.5*')
        assert result['sir'] == 'R'
        assert result['measurement'] == '0.5'
        assert result['sign'] == '='
        assert result['aes_modified'] is True

    def test_plain_mic_no_sign(self):
        result = _parse_vitek_cell('S:1')
        assert result['sir'] == 'S'
        assert result['measurement'] == '1'
        assert result['sign'] == '='

    def test_resistant_with_sign_and_aes(self):
        """R:>= 4* — sign + AES flag (theoretical edge case)."""
        result = _parse_vitek_cell('R:>= 4*')
        assert result['sir'] == 'R'
        assert result['measurement'] == '4'
        assert result['sign'] == '>='
        assert result['aes_modified'] is True

    def test_empty_cell_returns_none(self):
        assert _parse_vitek_cell('') is None
        assert _parse_vitek_cell('  ') is None

    def test_intermediate_plain(self):
        result = _parse_vitek_cell('I:1')
        assert result['sir'] == 'I'
        assert result['measurement'] == '1'
        assert result['sign'] == '='
        assert result['aes_modified'] is False


# --- Integration tests for convert subcommand ---

class TestConvertCSV:
    def test_basic_conversion(self, tmpdir):
        """Convert a minimal wide-format TSV to long-format CSV."""
        tsv = _write_vitek_tsv(tmpdir, [{
            'UID': 'ID00001',
            'Oxacillin': 'R:>= 4',
            'Vancomycin': 'S:<= 0.5',
        }])
        out = os.path.join(tmpdir, 'out.csv')
        result = subprocess.run(
            ['python', TOOL, 'convert', '--from-vitek-csv', tsv, '-o', out],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        with open(out) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 2
        oxa = [r for r in rows if r['antibiotic'] == 'Oxacillin'][0]
        assert oxa['sample_id'] == 'ID00001'
        assert oxa['resistance_phenotype'] == 'R'
        assert oxa['measurement'] == '4'
        assert oxa['measurement_sign'] == '>='
        assert oxa['measurement_units'] == 'mg/L'
        assert oxa['laboratory_typing_method'] == 'MIC'
        assert oxa['testing_standard'] == 'CLSI'
        assert oxa['platform'] == 'Vitek 2'

    def test_cefoxitin_screen(self, tmpdir):
        """Cefoxitin POS/NEG uses Screen method and screen units."""
        tsv = _write_vitek_tsv(tmpdir, [{
            'UID': 'ID00001',
            'Cefoxitin': 'R:POS',
        }])
        out = os.path.join(tmpdir, 'out.csv')
        subprocess.run(
            ['python', TOOL, 'convert', '--from-vitek-csv', tsv, '-o', out],
            capture_output=True, text=True
        )
        with open(out) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 1
        assert rows[0]['measurement'] == 'POS'
        assert rows[0]['measurement_units'] == 'screen'
        assert rows[0]['laboratory_typing_method'] == 'Screen'

    def test_aes_modified_flag(self, tmpdir):
        """AES-modified values (asterisk) set aes_modified=true."""
        tsv = _write_vitek_tsv(tmpdir, [{
            'UID': 'ID00001',
            'Nitrofurantoin': 'I:64*',
            'Oxacillin': 'R:1*',
        }])
        out = os.path.join(tmpdir, 'out.csv')
        subprocess.run(
            ['python', TOOL, 'convert', '--from-vitek-csv', tsv, '-o', out],
            capture_output=True, text=True
        )
        with open(out) as f:
            rows = list(csv.DictReader(f))
        nitro = [r for r in rows if r['antibiotic'] == 'Nitrofurantoin'][0]
        assert nitro['resistance_phenotype'] == 'I'
        assert nitro['measurement'] == '64'
        assert nitro['aes_modified'] == 'true'
        oxa = [r for r in rows if r['antibiotic'] == 'Oxacillin'][0]
        assert oxa['aes_modified'] == 'true'
        assert oxa['measurement'] == '1'

    def test_empty_cells_skipped(self, tmpdir):
        """Empty antibiotic columns are not emitted as rows."""
        tsv = _write_vitek_tsv(tmpdir, [{
            'UID': 'ID00001',
            'Oxacillin': 'R:>= 4',
            'Fosfomycin': '',
            'Mupirocin': '',
        }])
        out = os.path.join(tmpdir, 'out.csv')
        subprocess.run(
            ['python', TOOL, 'convert', '--from-vitek-csv', tsv, '-o', out],
            capture_output=True, text=True
        )
        with open(out) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 1
        assert rows[0]['antibiotic'] == 'Oxacillin'

    def test_multiple_samples(self, tmpdir):
        """Multiple samples produce correct row count."""
        tsv = _write_vitek_tsv(tmpdir, [
            {'UID': 'ID00001', 'Oxacillin': 'R:>= 4', 'Vancomycin': 'S:<= 0.5'},
            {'UID': 'ID00002', 'Oxacillin': 'S:0.5', 'Vancomycin': 'S:1'},
        ])
        out = os.path.join(tmpdir, 'out.csv')
        subprocess.run(
            ['python', TOOL, 'convert', '--from-vitek-csv', tsv, '-o', out],
            capture_output=True, text=True
        )
        with open(out) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 4  # 2 samples × 2 antibiotics

    def test_output_column_structure(self, tmpdir):
        """Output CSV has all expected columns."""
        tsv = _write_vitek_tsv(tmpdir, [{
            'UID': 'ID00001',
            'Oxacillin': 'R:>= 4',
        }])
        out = os.path.join(tmpdir, 'out.csv')
        subprocess.run(
            ['python', TOOL, 'convert', '--from-vitek-csv', tsv, '-o', out],
            capture_output=True, text=True
        )
        with open(out) as f:
            reader = csv.DictReader(f)
            _ = list(reader)
            header = reader.fieldnames
        for col in CONVERT_OUTPUT_COLUMNS:
            assert col in header, f"Missing column: {col}"

    def test_stderr_summary(self, tmpdir):
        """Convert prints a summary to stderr."""
        tsv = _write_vitek_tsv(tmpdir, [
            {'UID': 'ID00001', 'Oxacillin': 'R:>= 4'},
            {'UID': 'ID00002', 'Oxacillin': 'S:0.5'},
        ])
        out = os.path.join(tmpdir, 'out.csv')
        result = subprocess.run(
            ['python', TOOL, 'convert', '--from-vitek-csv', tsv, '-o', out],
            capture_output=True, text=True
        )
        assert '2 samples' in result.stderr
        assert '2 rows' in result.stderr
