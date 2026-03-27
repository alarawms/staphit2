"""Tests for bin/staphit-report."""
import csv
import importlib.machinery
import importlib.util
import json
import os

import pytest

# ---------------------------------------------------------------------------
# Import the script as a module (no .py extension)
# ---------------------------------------------------------------------------

TOOL_PATH = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'staphit-report')

loader = importlib.machinery.SourceFileLoader('staphit_report', TOOL_PATH)
spec = importlib.util.spec_from_loader('staphit_report', loader)
report_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report_mod)


# ===========================================================================
# Helpers — mock data writers
# ===========================================================================

SUMMARY_HEADERS = [
    'sample_id', 'assembly_length', 'n50', 'contigs',
    'mlst_st', 'spa_type', 'sccmec_type',
    'amrfinder_genes', 'point_mutations', 'geno_pheno_discordance',
    'virulence_summary',
]


def _write_summary(path, n=10):
    """Write a mock final_summary.tsv with n samples."""
    sts = ['ST8', 'ST8', 'ST5', 'ST5', 'ST5', 'ST22', 'ST22', 'ST30', 'ST239', 'ST398']
    spa_types = ['t008', 't008', 't002', 't002', 't002', 't032', 't032', 't019', 't037', 't011']
    sccmec = ['IV', 'IV', 'II', 'II', 'II', 'IV', 'IV', 'IV', 'III', 'V']
    amr = [
        'mecA;blaZ', 'mecA;blaZ;ermC', 'mecA;blaZ', 'mecA;aac(6\')-Ie',
        'mecA;blaZ;tetM', 'mecA;blaZ', 'mecA;ermA', 'mecA;blaZ;ermC',
        'mecA;blaZ;aac(6\')-Ie;ermA', 'mecA;tetM',
    ]
    mutations = [
        'gyrA:S84L', 'gyrA:S84L;grlA:S80F', 'gyrA:S84L', '',
        'rpoB:H481N', '', 'gyrA:S84L', '', 'gyrA:S84L;rpoB:H481N', '',
    ]
    discord = ['0', '1', '0', '', '0', '', '0', '', '1', '']
    virulence = [
        'PVL+;TSST-;IEC-B;ica:complete;hlg:complete',
        'PVL-;TSST-;IEC-A;ica:complete;hlg:complete',
        'PVL+;TSST+;IEC-B;ica:partial;hlg:complete',
        'PVL-;TSST-;IEC-negative;ica:absent;hlg:absent',
        'PVL-;TSST-;IEC-B;ica:complete;hlg:complete',
        'PVL+;TSST-;IEC-E;ica:complete;hlg:complete',
        'PVL-;TSST-;IEC-B;ica:complete;hlg:complete',
        'PVL-;TSST+;IEC-A;ica:partial;hlg:complete',
        'PVL+;TSST-;IEC-B;ica:complete;hlg:complete',
        'PVL-;TSST-;IEC-negative;ica:absent;hlg:absent',
    ]

    rows_data = []
    for i in range(n):
        rows_data.append({
            'sample_id': f'SAMPLE{i+1:02d}',
            'assembly_length': str(2700000 + i * 10000),
            'n50': str(100000 + i * 5000),
            'contigs': str(30 + i),
            'mlst_st': sts[i % len(sts)],
            'spa_type': spa_types[i % len(spa_types)],
            'sccmec_type': sccmec[i % len(sccmec)],
            'amrfinder_genes': amr[i % len(amr)],
            'point_mutations': mutations[i % len(mutations)],
            'geno_pheno_discordance': discord[i % len(discord)],
            'virulence_summary': virulence[i % len(virulence)],
        })

    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_HEADERS, delimiter='\t')
        writer.writeheader()
        for row in rows_data:
            writer.writerow(row)


def _write_clusters(path):
    """Write a mock cluster_report.json with one outbreak cluster."""
    clusters = [
        {
            'cluster_id': 'C001',
            'tier': 'high',
            'size': 4,
            'shared_st': 'ST8',
            'locations': ['ICU', 'Ward-A'],
            'date_range': '2026-01-15 to 2026-02-10',
        }
    ]
    with open(path, 'w') as fh:
        json.dump(clusters, fh, indent=2)


def _write_qc(path):
    """Write a mock qc_report.tsv with one pass and one fail."""
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=['sample_id', 'pass', 'reason', 'completeness', 'contamination'],
            delimiter='\t',
        )
        writer.writeheader()
        writer.writerow({
            'sample_id': 'SAMPLE01',
            'pass': 'PASS',
            'reason': 'OK',
            'completeness': '99.5',
            'contamination': '0.3',
        })
        writer.writerow({
            'sample_id': 'SAMPLE02',
            'pass': 'FAIL',
            'reason': 'completeness 82.0% < 90%',
            'completeness': '82.0',
            'contamination': '1.2',
        })


# ===========================================================================
# TestReport
# ===========================================================================

class TestReport:

    def test_generates_markdown(self, tmp_path):
        """Full report with all optional inputs produces valid Markdown."""
        summary = tmp_path / 'final_summary.tsv'
        clusters = tmp_path / 'cluster_report.json'
        qc = tmp_path / 'qc_report.tsv'
        out = tmp_path / 'run_report.md'

        _write_summary(str(summary))
        _write_clusters(str(clusters))
        _write_qc(str(qc))

        report_mod.main([
            '--summary', str(summary),
            '--clusters', str(clusters),
            '--qc', str(qc),
            '-o', str(out),
        ])

        text = out.read_text()
        assert '# Staphit Run Report' in text
        assert '**Samples:** 10' in text

        # Assembly QC section
        assert '## Assembly QC' in text
        assert 'Assembly length' in text
        assert 'N50' in text
        assert 'Contigs' in text

        # Typing overview
        assert '## Typing Overview' in text
        assert 'ST8' in text
        assert 't008' in text

        # Resistance
        assert '## Resistance' in text
        assert 'mecA' in text
        assert 'gyrA:S84L' in text

        # Virulence
        assert '## Virulence' in text
        assert 'PVL+ rate' in text
        assert 'TSST+ rate' in text

        # Outbreak clusters
        assert '## Outbreak Clusters' in text
        assert 'C001' in text
        assert 'ICU' in text

        # QC flags
        assert '## QC Flags' in text
        assert 'Passed' in text
        assert 'SAMPLE02' in text

    def test_with_qc(self, tmp_path):
        """Report with --qc includes QC flags section."""
        summary = tmp_path / 'final_summary.tsv'
        qc = tmp_path / 'qc_report.tsv'
        out = tmp_path / 'run_report.md'

        _write_summary(str(summary))
        _write_qc(str(qc))

        report_mod.main([
            '--summary', str(summary),
            '--qc', str(qc),
            '-o', str(out),
        ])

        text = out.read_text()
        assert '## QC Flags' in text
        assert '**Passed:** 1' in text
        assert '**Failed:** 1' in text
        assert 'SAMPLE02' in text
        assert 'completeness 82.0% < 90%' in text

        # No outbreak or plasmid sections
        assert '## Outbreak Clusters' not in text
        assert '## Plasmids' not in text

    def test_summary_only(self, tmp_path):
        """Report with only --summary produces core sections, no optional ones."""
        summary = tmp_path / 'final_summary.tsv'
        out = tmp_path / 'run_report.md'

        _write_summary(str(summary))

        report_mod.main([
            '--summary', str(summary),
            '-o', str(out),
        ])

        text = out.read_text()
        assert '# Staphit Run Report' in text
        assert '**Samples:** 10' in text
        assert '## Assembly QC' in text
        assert '## Typing Overview' in text
        assert '## Resistance' in text
        assert '## Virulence' in text

        # Optional sections are absent
        assert '## Outbreak Clusters' not in text
        assert '## QC Flags' not in text
        assert '## Plasmids' not in text

    def test_empty_summary(self, tmp_path):
        """Report with an empty summary TSV produces header only."""
        summary = tmp_path / 'final_summary.tsv'
        out = tmp_path / 'run_report.md'

        # Write headers only, no data rows
        with open(str(summary), 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=SUMMARY_HEADERS, delimiter='\t')
            writer.writeheader()

        report_mod.main([
            '--summary', str(summary),
            '-o', str(out),
        ])

        text = out.read_text()
        assert '# Staphit Run Report' in text
        assert '**Samples:** 0' in text

        # No data sections
        assert '## Assembly QC' not in text
        assert '## Typing Overview' not in text
        assert '## Resistance' not in text
        assert '## Virulence' not in text


# ===========================================================================
# TestCountFrequencies
# ===========================================================================

class TestCountFrequencies:

    def test_basic(self):
        series = ['mecA;blaZ', 'mecA;ermC', 'mecA']
        result = report_mod._count_frequencies(series)
        assert result['mecA'] == 3
        assert result['blaZ'] == 1
        assert result['ermC'] == 1

    def test_skips_empty_and_dash(self):
        series = ['', '-', 'mecA', '']
        result = report_mod._count_frequencies(series)
        assert result['mecA'] == 1
        assert len(result) == 1

    def test_empty_series(self):
        result = report_mod._count_frequencies([])
        assert len(result) == 0
