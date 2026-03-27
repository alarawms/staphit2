"""Tests for staphit-plasmids CLI tool."""
import importlib.machinery
import importlib.util
import os
import json
import subprocess
import tempfile
import csv

import pytest

# Load the bin/staphit-plasmids script (no .py extension) as a module.
TOOL_PATH = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'staphit-plasmids')

loader = importlib.machinery.SourceFileLoader('staphit_plasmids', TOOL_PATH)
spec = importlib.util.spec_from_loader('staphit_plasmids', loader)
plasmid = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plasmid)


# ---------------------------------------------------------------------------
# Helpers to create mock MOB-suite files
# ---------------------------------------------------------------------------

CONTIG_REPORT_HEADER = [
    'sample_id', 'molecule_type', 'primary_cluster_id', 'contig_id',
    'size', 'rep_type(s)', 'predicted_mobility',
]

MOBTYPER_HEADER = [
    'sample_id', 'size', 'predicted_mobility', 'rep_type(s)',
    'host_range', 'num_contigs',
]


def _write_contig_report(tmpdir, rows, filename='contig_report.txt'):
    """Write a mock contig_report.txt and return its path.

    rows: list of dicts with keys matching CONTIG_REPORT_HEADER.
    """
    path = os.path.join(tmpdir, filename)
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CONTIG_REPORT_HEADER, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def _write_mobtyper(tmpdir, rows, filename='mobtyper_results.txt'):
    """Write a mock mobtyper_results.txt and return its path.

    rows: list of dicts with keys matching MOBTYPER_HEADER.
    """
    path = os.path.join(tmpdir, filename)
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=MOBTYPER_HEADER, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestParseContigReport:
    """Test parse_contig_report separates chromosome and plasmid contigs."""

    def test_separates_chromosome_and_plasmid(self, tmp_path):
        rows = [
            {
                'sample_id': 'SAMPLE1',
                'molecule_type': 'chromosome',
                'primary_cluster_id': '-',
                'contig_id': 'contig_1',
                'size': '2800000',
                'rep_type(s)': '-',
                'predicted_mobility': '-',
            },
            {
                'sample_id': 'SAMPLE1',
                'molecule_type': 'plasmid',
                'primary_cluster_id': 'AA411',
                'contig_id': 'contig_2',
                'size': '24000',
                'rep_type(s)': 'rep16',
                'predicted_mobility': 'mobilizable',
            },
        ]
        path = _write_contig_report(str(tmp_path), rows)
        chrom, plas = plasmid.parse_contig_report(path)

        assert len(chrom) == 1
        assert len(plas) == 1
        assert chrom[0]['contig_id'] == 'contig_1'
        assert chrom[0]['molecule_type'] == 'chromosome'
        assert plas[0]['contig_id'] == 'contig_2'
        assert plas[0]['primary_cluster_id'] == 'AA411'


class TestParseMobtyper:
    """Test parse_mobtyper reads plasmid typing data."""

    def test_parses_plasmid_typing(self, tmp_path):
        rows = [
            {
                'sample_id': 'SAMPLE1:AA411',
                'size': '24000',
                'predicted_mobility': 'mobilizable',
                'rep_type(s)': 'rep16',
                'host_range': 'Staphylococcus',
                'num_contigs': '5',
            },
        ]
        path = _write_mobtyper(str(tmp_path), rows)
        result = plasmid.parse_mobtyper(path)

        assert 'AA411' in result
        assert result['AA411']['size'] == 24000
        assert result['AA411']['mobility'] == 'mobilizable'
        assert result['AA411']['rep_types'] == ['rep16']
        assert result['AA411']['host_range'] == 'Staphylococcus'
        assert result['AA411']['num_contigs'] == 5


class TestBuildProfile:
    """Test build_plasmid_profile combines both inputs."""

    def test_full_profile(self, tmp_path):
        contig_rows = [
            {
                'sample_id': 'SAMPLE1',
                'molecule_type': 'chromosome',
                'primary_cluster_id': '-',
                'contig_id': 'contig_1',
                'size': '2800000',
                'rep_type(s)': '-',
                'predicted_mobility': '-',
            },
            {
                'sample_id': 'SAMPLE1',
                'molecule_type': 'plasmid',
                'primary_cluster_id': 'AA411',
                'contig_id': 'contig_2',
                'size': '12000',
                'rep_type(s)': 'rep16',
                'predicted_mobility': 'mobilizable',
            },
            {
                'sample_id': 'SAMPLE1',
                'molecule_type': 'plasmid',
                'primary_cluster_id': 'AA411',
                'contig_id': 'contig_3',
                'size': '12000',
                'rep_type(s)': 'rep16',
                'predicted_mobility': 'mobilizable',
            },
            {
                'sample_id': 'SAMPLE1',
                'molecule_type': 'plasmid',
                'primary_cluster_id': 'AE018',
                'contig_id': 'contig_4',
                'size': '4500',
                'rep_type(s)': '-',
                'predicted_mobility': 'non-mobilizable',
            },
        ]
        mobtyper_rows = [
            {
                'sample_id': 'SAMPLE1:AA411',
                'size': '24000',
                'predicted_mobility': 'mobilizable',
                'rep_type(s)': 'rep16',
                'host_range': 'Staphylococcus',
                'num_contigs': '2',
            },
            {
                'sample_id': 'SAMPLE1:AE018',
                'size': '4500',
                'predicted_mobility': 'non-mobilizable',
                'rep_type(s)': '-',
                'host_range': 'Enterococcus',
                'num_contigs': '1',
            },
        ]
        contig_path = _write_contig_report(str(tmp_path), contig_rows)
        mobtyper_path = _write_mobtyper(str(tmp_path), mobtyper_rows)

        profile = plasmid.build_plasmid_profile(contig_path, mobtyper_path)

        assert profile['plasmid_count'] == 2

        # Verify plasmid details
        plasmid_ids = [p['cluster_id'] for p in profile['plasmids']]
        assert 'AA411' in plasmid_ids
        assert 'AE018' in plasmid_ids

        aa411 = next(p for p in profile['plasmids'] if p['cluster_id'] == 'AA411')
        assert aa411['size'] == 24000
        assert aa411['mobility'] == 'mobilizable'
        assert aa411['rep_types'] == ['rep16']

        # Verify contig_map
        assert profile['contig_map']['contig_1'] == 'chromosome'
        assert profile['contig_map']['contig_2'] == 'AA411'
        assert profile['contig_map']['contig_3'] == 'AA411'
        assert profile['contig_map']['contig_4'] == 'AE018'

        # Verify _summary
        assert '_summary' in profile


class TestAssignGenes:
    """Test assign_genes_to_locations maps genes to chromosome or plasmid."""

    def test_assigns_genes_to_plasmid_or_chromosome(self):
        contig_map = {
            'contig_1': 'chromosome',
            'contig_2': 'AA411',
            'contig_3': 'AE018',
        }
        gene_list = [
            {'gene': 'mecA', 'contig': 'contig_1'},
            {'gene': 'blaZ', 'contig': 'contig_2'},
        ]
        result = plasmid.assign_genes_to_locations(gene_list, contig_map)

        assert result['mecA'] == 'chromosome'
        assert result['blaZ'] == 'AA411'


class TestSummary:
    """Test format_plasmid_summary output format."""

    def test_format(self):
        profile = {
            'plasmids': [
                {
                    'cluster_id': 'AA411',
                    'mobility': 'mobilizable',
                    'rep_types': ['rep16'],
                    'size': 24000,
                    'host_range': 'Staphylococcus',
                    'num_contigs': 5,
                },
                {
                    'cluster_id': 'AE018',
                    'mobility': 'non-mobilizable',
                    'rep_types': [],
                    'size': 4500,
                    'host_range': 'Enterococcus',
                    'num_contigs': 1,
                },
            ],
        }
        summary = plasmid.format_plasmid_summary(profile)
        assert summary == 'AA411:mobilizable(rep16);AE018:non-mobilizable'


class TestCLI:
    """Test command-line interface end-to-end."""

    def test_basic_invocation(self, tmp_path):
        contig_rows = [
            {
                'sample_id': 'SAMPLE1',
                'molecule_type': 'chromosome',
                'primary_cluster_id': '-',
                'contig_id': 'contig_1',
                'size': '2800000',
                'rep_type(s)': '-',
                'predicted_mobility': '-',
            },
            {
                'sample_id': 'SAMPLE1',
                'molecule_type': 'plasmid',
                'primary_cluster_id': 'AA411',
                'contig_id': 'contig_2',
                'size': '24000',
                'rep_type(s)': 'rep16',
                'predicted_mobility': 'mobilizable',
            },
        ]
        mobtyper_rows = [
            {
                'sample_id': 'SAMPLE1:AA411',
                'size': '24000',
                'predicted_mobility': 'mobilizable',
                'rep_type(s)': 'rep16',
                'host_range': 'Staphylococcus',
                'num_contigs': '1',
            },
        ]
        contig_path = _write_contig_report(str(tmp_path), contig_rows)
        mobtyper_path = _write_mobtyper(str(tmp_path), mobtyper_rows)
        out_path = os.path.join(str(tmp_path), 'output.json')

        result = subprocess.run(
            ['python', TOOL_PATH,
             '--contig-report', contig_path,
             '--mobtyper', mobtyper_path,
             '-o', out_path],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"

        with open(out_path) as fh:
            data = json.load(fh)

        assert data['plasmid_count'] == 1
        assert data['plasmids'][0]['cluster_id'] == 'AA411'
        assert data['contig_map']['contig_1'] == 'chromosome'
        assert data['contig_map']['contig_2'] == 'AA411'
        assert 'AA411:mobilizable(rep16)' in data['_summary']
