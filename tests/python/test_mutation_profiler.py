"""Tests for staphit-mutations CLI tool."""
import importlib.machinery
import importlib.util
import json
import os
import subprocess
import sys
import tempfile

import pytest

# Load the bin/staphit-mutations script (no .py extension) as a module.
# Patch sys.argv so main() doesn't run during import if guarded by __name__.
TOOL_PATH = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'staphit-mutations')

loader = importlib.machinery.SourceFileLoader('staphit_mutations', TOOL_PATH)
spec = importlib.util.spec_from_loader('staphit_mutations', loader)
mut = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mut)

# ---------------------------------------------------------------------------
# TSV helpers
# ---------------------------------------------------------------------------

HEADER = (
    "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\t"
    "Element name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\t"
    "Target length\tReference sequence length\t% Coverage of reference\t"
    "% Identity to reference\tAlignment length\tClosest reference accession\t"
    "Closest reference name\tHMM accession\tHMM description"
)


def _make_row(symbol, subtype, drug_class, **overrides):
    """Build a single AMRFinderPlus TSV row with sensible defaults."""
    fields = {
        'Protein id': 'prot_001',
        'Contig id': 'contig_1',
        'Start': '100',
        'Stop': '200',
        'Strand': '+',
        'Element symbol': symbol,
        'Element name': 'some element',
        'Scope': 'core',
        'Type': 'AMR',
        'Subtype': subtype,
        'Class': drug_class,
        'Subclass': drug_class,
        'Method': 'POINTX',
        'Target length': '300',
        'Reference sequence length': '300',
        '% Coverage of reference': '100.00',
        '% Identity to reference': '99.50',
        'Alignment length': '300',
        'Closest reference accession': 'WP_000',
        'Closest reference name': 'ref',
        'HMM accession': '',
        'HMM description': '',
    }
    fields.update(overrides)
    return '\t'.join(fields[col] for col in HEADER.split('\t'))


def _write_tsv(rows_text, tmpdir):
    """Write header + rows to a temp TSV file and return path."""
    path = os.path.join(tmpdir, 'amrfinder.tsv')
    with open(path, 'w') as fh:
        fh.write(HEADER + '\n')
        for row in rows_text:
            fh.write(row + '\n')
    return path


# ---------------------------------------------------------------------------
# TestParseAmrfinderMutations
# ---------------------------------------------------------------------------

class TestParseAmrfinderMutations:

    def test_extracts_point_rows_only(self, tmp_path):
        """Only rows with Subtype==POINT should be extracted."""
        rows = [
            _make_row('gyrA_S84L', 'POINT', 'QUINOLONE'),
            _make_row('mecA', 'AMR', 'BETA-LACTAM'),
            _make_row('rpoB_H481N', 'POINT', 'RIFAMYCIN'),
        ]
        tsv = _write_tsv(rows, str(tmp_path))
        result = mut.parse_amrfinder_mutations(tsv)
        assert len(result) == 2
        genes = {r['gene'] for r in result}
        assert genes == {'gyrA', 'rpoB'}

    def test_parses_symbol_format(self, tmp_path):
        """Element symbol 'gene_mutation' is split correctly."""
        rows = [_make_row('parC_S80F', 'POINT', 'QUINOLONE')]
        tsv = _write_tsv(rows, str(tmp_path))
        result = mut.parse_amrfinder_mutations(tsv)
        assert len(result) == 1
        assert result[0]['gene'] == 'parC'
        assert result[0]['mutation'] == 'S80F'

    def test_empty_input(self, tmp_path):
        """Empty TSV (header only) returns empty list."""
        tsv = _write_tsv([], str(tmp_path))
        result = mut.parse_amrfinder_mutations(tsv)
        assert result == []


# ---------------------------------------------------------------------------
# TestClassifyMutations
# ---------------------------------------------------------------------------

class TestClassifyMutations:

    def test_quinolone_mutations_grouped(self):
        mutations = [
            {'gene': 'gyrA', 'mutation': 'S84L', 'drug_class': 'Fluoroquinolones'},
            {'gene': 'parC', 'mutation': 'S80F', 'drug_class': 'Fluoroquinolones'},
        ]
        classified = mut.classify_mutations(mutations)
        assert 'Fluoroquinolones' in classified
        assert len(classified['Fluoroquinolones']) == 2

    def test_rifampicin(self):
        mutations = [
            {'gene': 'rpoB', 'mutation': 'H481N', 'drug_class': 'Rifampicin'},
        ]
        classified = mut.classify_mutations(mutations)
        assert 'Rifampicin' in classified
        assert len(classified['Rifampicin']) == 1

    def test_empty(self):
        assert mut.classify_mutations([]) == {}


# ---------------------------------------------------------------------------
# TestPredictPhenotype
# ---------------------------------------------------------------------------

class TestPredictPhenotype:

    def test_predicts_r_for_fluoroquinolones(self):
        mutations = [
            {'gene': 'gyrA', 'mutation': 'S84L', 'drug_class': 'Fluoroquinolones'},
        ]
        predictions = mut.predict_phenotype(mutations)
        assert predictions == {'Fluoroquinolones': 'R'}

    def test_empty_returns_empty(self):
        assert mut.predict_phenotype([]) == {}


# ---------------------------------------------------------------------------
# TestGenoPheno
# ---------------------------------------------------------------------------

class TestGenoPheno:

    def test_concordant(self):
        """Predicted R matches observed R -> concordant."""
        predictions = {'Fluoroquinolones': 'R'}
        antibiogram = [{'antibiotic': 'Levofloxacin', 'sir': 'R'}]
        result = mut.compare_genotype_phenotype(predictions, antibiogram)
        assert result['Fluoroquinolones']['concordant'] is True

    def test_discordant(self):
        """Predicted R but observed S -> discordant."""
        predictions = {'Fluoroquinolones': 'R'}
        antibiogram = [{'antibiotic': 'Levofloxacin', 'sir': 'S'}]
        result = mut.compare_genotype_phenotype(predictions, antibiogram)
        assert result['Fluoroquinolones']['concordant'] is False

    def test_no_antibiogram_returns_empty(self):
        predictions = {'Fluoroquinolones': 'R'}
        result = mut.compare_genotype_phenotype(predictions, {})
        assert result == {}

    def test_no_predictions_returns_empty(self):
        antibiogram = [{'antibiotic': 'Levofloxacin', 'sir': 'R'}]
        result = mut.compare_genotype_phenotype({}, antibiogram)
        assert result == {}


# ---------------------------------------------------------------------------
# TestMutationSummary
# ---------------------------------------------------------------------------

class TestMutationSummary:

    def test_format(self):
        mutations = [
            {'gene': 'gyrA', 'mutation': 'S84L', 'drug_class': 'Fluoroquinolones'},
            {'gene': 'parC', 'mutation': 'S80F', 'drug_class': 'Fluoroquinolones'},
            {'gene': 'rpoB', 'mutation': 'H481N', 'drug_class': 'Rifampicin'},
        ]
        summary = mut.format_mutation_summary(mutations)
        assert summary == 'gyrA:S84L;parC:S80F;rpoB:H481N'

    def test_empty_returns_empty_string(self):
        assert mut.format_mutation_summary([]) == ''


# ---------------------------------------------------------------------------
# TestCLI
# ---------------------------------------------------------------------------

class TestCLI:

    def test_basic_invocation(self, tmp_path):
        """Run CLI on a mock TSV and check JSON output."""
        rows = [
            _make_row('gyrA_S84L', 'POINT', 'QUINOLONE'),
            _make_row('rpoB_H481N', 'POINT', 'RIFAMYCIN'),
        ]
        tsv = _write_tsv(rows, str(tmp_path))
        out_path = os.path.join(str(tmp_path), 'result.json')

        result = subprocess.run(
            [sys.executable, TOOL_PATH, tsv, '-o', out_path],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr

        with open(out_path) as fh:
            data = json.load(fh)

        assert 'mutations' in data
        assert 'classified' in data
        assert 'predictions' in data
        assert 'summary' in data

    def test_with_antibiogram(self, tmp_path):
        """Pass antibiogram JSON and check concordance section."""
        rows = [_make_row('gyrA_S84L', 'POINT', 'QUINOLONE')]
        tsv = _write_tsv(rows, str(tmp_path))

        abg = [{'antibiotic': 'Levofloxacin', 'sir': 'R'}]
        abg_path = os.path.join(str(tmp_path), 'antibiogram.json')
        with open(abg_path, 'w') as fh:
            json.dump(abg, fh)

        out_path = os.path.join(str(tmp_path), 'result.json')
        result = subprocess.run(
            [sys.executable, TOOL_PATH, tsv, '-o', out_path, '--antibiogram', abg_path],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr

        with open(out_path) as fh:
            data = json.load(fh)

        assert 'concordance' in data
        assert data['concordance']['Fluoroquinolones']['concordant'] is True


# ---------------------------------------------------------------------------
# TestRealData - validate against real AMRFinderPlus pipeline output
# ---------------------------------------------------------------------------

RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results', 'amrfinderplus')


@pytest.mark.skipif(
    not os.path.exists(os.path.join(os.path.dirname(__file__), '..', 'results', 'amrfinderplus', 'ID00080_amrfinder.tsv')),
    reason="Pipeline results not available",
)
class TestRealData:

    def test_id00080_has_quinolone_mutations(self):
        result = subprocess.run(
            [sys.executable, TOOL_PATH, os.path.join(RESULTS_DIR, 'ID00080_amrfinder.tsv')],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        data = json.loads(result.stdout)
        genes = [m['gene'] for m in data['mutations']]
        assert 'gyrA' in genes
        assert 'parC' in genes
        assert 'Fluoroquinolones' in data['predictions']

    def test_id00080_has_trimethoprim_mutation(self):
        result = subprocess.run(
            [sys.executable, TOOL_PATH, os.path.join(RESULTS_DIR, 'ID00080_amrfinder.tsv')],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        data = json.loads(result.stdout)
        genes = [m['gene'] for m in data['mutations']]
        assert 'dfrB' in genes
        assert 'Trimethoprim' in data['predictions']

    def test_id00001_has_quinolone_and_fosfomycin_mutations(self):
        """ID00001 has gyrA_S84L, parC_S80F (quinolone) and glpT_A100V (fosfomycin)."""
        tsv = os.path.join(RESULTS_DIR, 'ID00001_amrfinder.tsv')
        if not os.path.exists(tsv):
            pytest.skip("ID00001 not available")
        result = subprocess.run(
            [sys.executable, TOOL_PATH, tsv], capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        data = json.loads(result.stdout)
        genes = [m['gene'] for m in data['mutations']]
        # Should have exactly 3 POINT mutations, no AMR-type rows
        assert len(data['mutations']) == 3
        assert 'gyrA' in genes
        assert 'parC' in genes
        assert 'glpT' in genes
        assert 'Fluoroquinolones' in data['predictions']
        assert 'Fosfomycin' in data['predictions']
