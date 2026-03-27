"""Tests for staphit-cluster outbreak clustering tool."""
import csv
import importlib.machinery
import importlib.util
import json
import os
import subprocess
import sys

import pytest

# Load the bin/staphit-cluster script (no .py extension) as a module.
TOOL_PATH = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'staphit-cluster')

loader = importlib.machinery.SourceFileLoader('staphit_cluster', TOOL_PATH)
spec = importlib.util.spec_from_loader('staphit_cluster', loader)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_dist_matrix(tmpdir, samples, distances_dict, filename='snp_dists.tsv'):
    """Write a square distance matrix TSV.

    Parameters
    ----------
    tmpdir : str
        Directory to write the file in.
    samples : list of str
        Sample names.
    distances_dict : dict
        Mapping of (s1, s2) -> int distance. Missing pairs default to 0 for
        self-comparisons. Symmetric entries are filled automatically.
    filename : str
        Filename for the matrix.

    Returns
    -------
    str
        Path to the written file.
    """
    # Build full symmetric dict
    full = {}
    for (s1, s2), d in distances_dict.items():
        full[(s1, s2)] = d
        full[(s2, s1)] = d

    path = os.path.join(str(tmpdir), filename)
    with open(path, 'w', newline='') as fh:
        # Header: empty first cell + sample names
        fh.write('\t'.join([''] + samples) + '\n')
        for s1 in samples:
            row = [s1]
            for s2 in samples:
                if s1 == s2:
                    row.append('0')
                else:
                    row.append(str(full.get((s1, s2), 999)))
            fh.write('\t'.join(row) + '\n')
    return path


def _write_summary(tmpdir, rows_list_of_dicts, filename='final_summary.tsv'):
    """Write a mock final_summary.tsv.

    Parameters
    ----------
    tmpdir : str
        Directory to write the file in.
    rows_list_of_dicts : list of dict
        Each dict is a row. All unique keys become fieldnames.
    filename : str
        Filename for the summary.

    Returns
    -------
    str
        Path to the written file.
    """
    path = os.path.join(str(tmpdir), filename)
    if not rows_list_of_dicts:
        with open(path, 'w', newline='') as fh:
            fh.write('sample_id\n')
        return path

    # Collect all fieldnames preserving order
    fieldnames = []
    seen = set()
    for row in rows_list_of_dicts:
        for k in row:
            if k not in seen:
                fieldnames.append(k)
                seen.add(k)

    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames,
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows_list_of_dicts)
    return path


# ---------------------------------------------------------------------------
# TestParseDistanceMatrix
# ---------------------------------------------------------------------------


class TestParseDistanceMatrix:

    def test_parse_square_matrix(self, tmp_path):
        """3 samples, verify symmetric access."""
        samples = ['S1', 'S2', 'S3']
        dists = {('S1', 'S2'): 5, ('S1', 'S3'): 10, ('S2', 'S3'): 15}
        path = _write_dist_matrix(tmp_path, samples, dists)

        matrix, parsed_samples = mod.parse_distance_matrix(path)

        assert set(parsed_samples) == {'S1', 'S2', 'S3'}
        assert matrix[('S1', 'S2')] == 5
        assert matrix[('S2', 'S1')] == 5
        assert matrix[('S1', 'S3')] == 10
        assert matrix[('S3', 'S1')] == 10
        assert matrix[('S2', 'S3')] == 15
        assert matrix[('S3', 'S2')] == 15

    def test_snp_dists_format(self, tmp_path):
        """Handle 'snp-dists 0.8.2' in top-left cell."""
        samples = ['A', 'B']
        path = os.path.join(str(tmp_path), 'snpdists.tsv')
        with open(path, 'w') as fh:
            fh.write('snp-dists 0.8.2\tA\tB\n')
            fh.write('A\t0\t7\n')
            fh.write('B\t7\t0\n')

        matrix, parsed_samples = mod.parse_distance_matrix(path)

        assert set(parsed_samples) == {'A', 'B'}
        assert matrix[('A', 'B')] == 7
        assert matrix[('B', 'A')] == 7


# ---------------------------------------------------------------------------
# TestClusterAtThreshold
# ---------------------------------------------------------------------------


class TestClusterAtThreshold:

    def test_two_clusters(self):
        """A-B=3, C-D=4, all others=50, threshold=5 -> {A,B} and {C,D}."""
        samples = ['A', 'B', 'C', 'D']
        matrix = {
            ('A', 'B'): 3, ('B', 'A'): 3,
            ('C', 'D'): 4, ('D', 'C'): 4,
            ('A', 'C'): 50, ('C', 'A'): 50,
            ('A', 'D'): 50, ('D', 'A'): 50,
            ('B', 'C'): 50, ('C', 'B'): 50,
            ('B', 'D'): 50, ('D', 'B'): 50,
        }
        result = mod.cluster_at_threshold(matrix, samples, 5)

        # A and B same cluster
        assert result['A'] == result['B']
        # C and D same cluster
        assert result['C'] == result['D']
        # Different clusters
        assert result['A'] != result['C']

    def test_single_linkage_chain(self):
        """A-B=3, B-C=4, A-C=8, threshold=5 -> all same cluster via chaining."""
        samples = ['A', 'B', 'C']
        matrix = {
            ('A', 'B'): 3, ('B', 'A'): 3,
            ('B', 'C'): 4, ('C', 'B'): 4,
            ('A', 'C'): 8, ('C', 'A'): 8,
        }
        result = mod.cluster_at_threshold(matrix, samples, 5)

        # All in same cluster via single-linkage
        assert result['A'] == result['B'] == result['C']

    def test_no_clusters_at_low_threshold(self):
        """A-B=50, threshold=5 -> different clusters."""
        samples = ['A', 'B']
        matrix = {('A', 'B'): 50, ('B', 'A'): 50}
        result = mod.cluster_at_threshold(matrix, samples, 5)

        assert result['A'] != result['B']

    def test_singletons(self):
        """Single sample gets a cluster ID."""
        samples = ['A']
        matrix = {}
        result = mod.cluster_at_threshold(matrix, samples, 5)

        assert 'A' in result
        assert result['A'].startswith('C')


# ---------------------------------------------------------------------------
# TestTieredClusters
# ---------------------------------------------------------------------------


class TestTieredClusters:

    def test_three_tiers(self):
        """A-B=3, C-D=12, all cross=30.

        At t=5: {A,B} separate from {C,D}.
        At t=15: both pairs clustered.
        At t=40: all together.
        """
        samples = ['A', 'B', 'C', 'D']
        matrix = {
            ('A', 'B'): 3, ('B', 'A'): 3,
            ('C', 'D'): 12, ('D', 'C'): 12,
            ('A', 'C'): 30, ('C', 'A'): 30,
            ('A', 'D'): 30, ('D', 'A'): 30,
            ('B', 'C'): 30, ('C', 'B'): 30,
            ('B', 'D'): 30, ('D', 'B'): 30,
        }
        thresholds = [5, 15, 40]
        tiers = mod.assign_tiered_clusters(matrix, samples, thresholds)

        assert len(tiers) == 3

        # Tier 0 (t=5): A,B together; C,D separate from each other
        assert tiers[0]['A'] == tiers[0]['B']
        assert tiers[0]['C'] != tiers[0]['A']
        assert tiers[0]['C'] != tiers[0]['D']  # C-D=12 > 5

        # Tier 1 (t=15): A,B together; C,D together; but the two groups separate
        assert tiers[1]['A'] == tiers[1]['B']
        assert tiers[1]['C'] == tiers[1]['D']
        assert tiers[1]['A'] != tiers[1]['C']

        # Tier 2 (t=40): all together
        assert tiers[2]['A'] == tiers[2]['B'] == tiers[2]['C'] == tiers[2]['D']


# ---------------------------------------------------------------------------
# TestCLI
# ---------------------------------------------------------------------------


class TestCLI:

    def test_basic_invocation(self, tmp_path):
        """3 samples, check output files exist and have correct content."""
        samples = ['S1', 'S2', 'S3']
        dists = {('S1', 'S2'): 3, ('S1', 'S3'): 50, ('S2', 'S3'): 50}
        snp_path = _write_dist_matrix(tmp_path, samples, dists, 'snp.tsv')

        summary_rows = [
            {'sample_id': 'S1', 'st': 'ST8', 'spa_type': 't008',
             'location': 'ICU', 'collection_date': '2025-01-01',
             'amrfinder_genes': 'mecA;blaZ'},
            {'sample_id': 'S2', 'st': 'ST8', 'spa_type': 't008',
             'location': 'ICU', 'collection_date': '2025-01-05',
             'amrfinder_genes': 'mecA;blaZ;ermC'},
            {'sample_id': 'S3', 'st': 'ST5', 'spa_type': 't002',
             'location': 'ER', 'collection_date': '2025-02-01',
             'amrfinder_genes': 'mecA'},
        ]
        summary_path = _write_summary(tmp_path, summary_rows, 'summary.tsv')

        outdir = os.path.join(str(tmp_path), 'output')
        result = subprocess.run(
            [sys.executable, TOOL_PATH,
             '--snp-dists', snp_path,
             '--summary', summary_path,
             '-o', outdir],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr

        # Check files exist
        clusters_path = os.path.join(outdir, 'clusters.tsv')
        report_path = os.path.join(outdir, 'cluster_report.json')
        pairs_path = os.path.join(outdir, 'transmission_pairs.tsv')

        assert os.path.exists(clusters_path)
        assert os.path.exists(report_path)
        assert os.path.exists(pairs_path)

        # Check clusters.tsv content
        with open(clusters_path, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            rows = list(reader)
        assert len(rows) == 3
        assert 'sample_id' in reader.fieldnames
        assert 'direct_transmission_cluster' in reader.fieldnames
        assert 'outbreak_cluster' in reader.fieldnames
        assert 'related_cluster' in reader.fieldnames

        # S1 and S2 should be in same direct_transmission cluster (dist=3 <= 5)
        s1_dt = [r for r in rows if r['sample_id'] == 'S1'][0]['direct_transmission_cluster']
        s2_dt = [r for r in rows if r['sample_id'] == 'S2'][0]['direct_transmission_cluster']
        s3_dt = [r for r in rows if r['sample_id'] == 'S3'][0]['direct_transmission_cluster']
        assert s1_dt == s2_dt
        assert s1_dt != s3_dt

        # Check cluster_report.json
        with open(report_path) as fh:
            report = json.load(fh)
        assert isinstance(report, list)
        # Should have at least one multi-sample cluster
        dt_clusters = [c for c in report if c['tier'] == 'direct_transmission']
        assert len(dt_clusters) >= 1
        dt = dt_clusters[0]
        assert set(dt['samples']) == {'S1', 'S2'}
        assert dt['shared_st'] == 'ST8'
        assert dt['shared_spa'] == 't008'
        assert 'mecA' in dt['shared_resistance']
        assert 'blaZ' in dt['shared_resistance']
        assert dt['max_snp_distance'] == 3

        # Check transmission_pairs.tsv
        with open(pairs_path, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            pairs = list(reader)
        assert len(pairs) == 1
        assert pairs[0]['sample_1'] == 'S1'
        assert pairs[0]['sample_2'] == 'S2'
        assert pairs[0]['snp_distance'] == '3'
        assert pairs[0]['same_st'] == 'True'

    def test_without_summary(self, tmp_path):
        """Works without --summary flag."""
        samples = ['A', 'B', 'C']
        dists = {('A', 'B'): 3, ('A', 'C'): 50, ('B', 'C'): 50}
        snp_path = _write_dist_matrix(tmp_path, samples, dists, 'snp.tsv')

        outdir = os.path.join(str(tmp_path), 'output')
        result = subprocess.run(
            [sys.executable, TOOL_PATH,
             '--snp-dists', snp_path,
             '-o', outdir],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr

        # All three output files should exist
        assert os.path.exists(os.path.join(outdir, 'clusters.tsv'))
        assert os.path.exists(os.path.join(outdir, 'cluster_report.json'))
        assert os.path.exists(os.path.join(outdir, 'transmission_pairs.tsv'))

        # Clusters file should have all samples
        with open(os.path.join(outdir, 'clusters.tsv'), newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            rows = list(reader)
        sample_ids = {r['sample_id'] for r in rows}
        assert sample_ids == {'A', 'B', 'C'}
