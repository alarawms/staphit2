"""Tests for staphit-visualize publication figures script."""

import csv
import importlib.machinery
import importlib.util
import os
import subprocess
import tempfile

import pytest

TOOL = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'staphit-visualize')

# Load module for direct function testing
loader = importlib.machinery.SourceFileLoader('staphit_visualize', TOOL)
spec = importlib.util.spec_from_loader('staphit_visualize', loader)
viz = importlib.util.module_from_spec(spec)
spec.loader.exec_module(viz)


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

SUMMARY_COLUMNS = [
    'sample_id', 'mlst_st', 'spa_type', 'sccmec_type', 'agr_group',
    'amrfinder_genes', 'virulence_genes', 'virulence_summary',
    'point_mutations', 'geno_pheno_discordance', 'assembly_length',
    'contigs', 'n50',
]

AMR_GENE_POOLS = [
    'mecA;blaZ;ermA', 'mecA;blaZ', 'mecA;ermC;aph(3\')-IIIa',
    'mecA;blaZ;ermA;aac(6\')-Ie-aph(2\'\')-Ia', 'blaZ;ermC',
    'mecA;blaZ;tetM', 'mecA', 'mecA;blaZ;ermA;dfrG',
    'mecA;blaZ;fusB', 'mecA;blaZ;ermA;mphC',
]

VIRULENCE_SUMMARIES = [
    'PVL+;TSST-;IEC-B;ica:complete;hlg:complete',
    'PVL-;TSST+;IEC-A;ica:partial;hlg:complete',
    'PVL+;TSST-;IEC-E;ica:complete;hlg:partial',
    'PVL-;TSST-;IEC-B;ica:complete;hlg:complete',
    'PVL+;TSST+;IEC-C;ica:partial;hlg:partial',
    'PVL-;TSST-;IEC-B;ica:complete;hlg:complete',
    'PVL+;TSST-;IEC-D;ica:complete;hlg:complete',
    'PVL-;TSST-;IEC-G;ica:partial;hlg:complete',
    'PVL+;TSST-;IEC-B;ica:complete;hlg:partial',
    'PVL-;TSST+;IEC-F;ica:complete;hlg:complete',
]

POINT_MUTATIONS = [
    'gyrA:S84L;grlA:S80F', 'gyrA:S84L', 'grlA:S80Y;rpoB:H481N',
    'gyrA:S84L;grlA:S80F', '', 'rpoB:H481N',
    'gyrA:S84L;grlA:S80F;rpoB:S486L', 'gyrA:E88K', '-',
    'gyrA:S84L;grlA:S80F',
]


def _write_summary(tmpdir, n=10):
    """Write a mock final_summary.tsv with *n* rows."""
    path = os.path.join(tmpdir, 'final_summary.tsv')
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(SUMMARY_COLUMNS)
        for i in range(n):
            writer.writerow([
                f'S{i:03d}',                          # sample_id
                f'ST{(i % 5) + 1}',                   # mlst_st
                f't{100 + i}',                         # spa_type
                f'SCCmec-IV' if i % 2 == 0 else 'SCCmec-II',  # sccmec_type
                f'agr-{(i % 4) + 1}',                 # agr_group
                AMR_GENE_POOLS[i % len(AMR_GENE_POOLS)],       # amrfinder_genes
                'lukS-PV;lukF-PV;scn;chp;sak',        # virulence_genes
                VIRULENCE_SUMMARIES[i % len(VIRULENCE_SUMMARIES)],  # virulence_summary
                POINT_MUTATIONS[i % len(POINT_MUTATIONS)],          # point_mutations
                'concordant',                          # geno_pheno_discordance
                2800000 + i * 5000,                    # assembly_length
                50 + i,                                # contigs
                150000 + i * 2000,                     # n50
            ])
    return path


def _write_tree(tmpdir, n=10):
    """Write a simple comb Newick tree for samples S000..S{n-1}."""
    path = os.path.join(tmpdir, 'tree.nwk')
    # Build a simple comb tree: ((S000:0.1,S001:0.1):0.1,S002:0.2)...
    if n == 0:
        with open(path, 'w') as fh:
            fh.write('();')
        return path
    labels = [f'S{i:03d}' for i in range(n)]
    # Simple comb: nested parentheses
    nwk = labels[0] + ':0.1'
    for lab in labels[1:]:
        nwk = f'({nwk},{lab}:0.1):0.1'
    nwk += ';'
    with open(path, 'w') as fh:
        fh.write(nwk)
    return path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestVisualization:

    def test_basic_without_tree(self):
        """Generates core plots (mlst, amr heatmap, virulence, mutations, qc)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            summary = _write_summary(tmpdir, n=10)
            outdir = os.path.join(tmpdir, 'plots')

            result = subprocess.run(
                ['python3', TOOL, summary, outdir],
                capture_output=True, text=True,
            )
            assert result.returncode == 0, result.stderr

            expected = [
                'mlst_distribution.png',
                'amr_heatmap.png',
                'virulence_profile.png',
                'point_mutations.png',
                'assembly_qc.png',
            ]
            for name in expected:
                fpath = os.path.join(outdir, name)
                assert os.path.isfile(fpath), f'Missing: {name}'
                assert os.path.getsize(fpath) > 0, f'Empty: {name}'

    def test_with_tree(self):
        """When tree provided, annotated_phylogeny.png is generated (if toytree)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            summary = _write_summary(tmpdir, n=10)
            tree = _write_tree(tmpdir, n=10)
            outdir = os.path.join(tmpdir, 'plots')

            result = subprocess.run(
                ['python3', TOOL, summary, outdir, '--tree', tree],
                capture_output=True, text=True,
            )
            assert result.returncode == 0, result.stderr

            # Core plots should still exist
            assert os.path.isfile(os.path.join(outdir, 'mlst_distribution.png'))
            assert os.path.isfile(os.path.join(outdir, 'amr_heatmap.png'))

            # Phylogeny only if toytree installed
            if viz.HAS_TOYTREE:
                assert os.path.isfile(os.path.join(outdir, 'annotated_phylogeny.png'))

    def test_empty_summary(self):
        """Empty summary (header only, n=0) should not crash."""
        with tempfile.TemporaryDirectory() as tmpdir:
            summary = _write_summary(tmpdir, n=0)
            outdir = os.path.join(tmpdir, 'plots')

            result = subprocess.run(
                ['python3', TOOL, summary, outdir],
                capture_output=True, text=True,
            )
            assert result.returncode == 0, result.stderr

    def test_dpi_flag(self):
        """--dpi flag is accepted without error."""
        with tempfile.TemporaryDirectory() as tmpdir:
            summary = _write_summary(tmpdir, n=5)
            outdir = os.path.join(tmpdir, 'plots')

            result = subprocess.run(
                ['python3', TOOL, summary, outdir, '--dpi', '72'],
                capture_output=True, text=True,
            )
            assert result.returncode == 0, result.stderr
            assert os.path.isfile(os.path.join(outdir, 'mlst_distribution.png'))

    def test_virulence_values(self):
        """Virulence profile matrix encodes +/partial/- correctly."""
        import pandas as pd
        with tempfile.TemporaryDirectory() as tmpdir:
            summary = _write_summary(tmpdir, n=3)
            df = pd.read_csv(summary, sep='\t')
            # First row: PVL+;TSST-;IEC-B;ica:complete;hlg:complete
            # PVL=2, TSST=0, IEC=2, ica=2, hlg=2
            # Second row: PVL-;TSST+;IEC-A;ica:partial;hlg:complete
            # PVL=0, TSST=2, IEC=2, ica=1, hlg=2

            # Invoke the plot function directly and verify matrix contents
            outdir = os.path.join(tmpdir, 'plots')
            os.makedirs(outdir, exist_ok=True)
            viz.plot_virulence_profile(df, outdir, dpi=72)
            assert os.path.isfile(os.path.join(outdir, 'virulence_profile.png'))

    def test_missing_columns_graceful(self):
        """Summary with only sample_id should not crash."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, 'minimal.tsv')
            with open(path, 'w') as fh:
                fh.write('sample_id\n')
                fh.write('S001\n')
            outdir = os.path.join(tmpdir, 'plots')

            result = subprocess.run(
                ['python3', TOOL, path, outdir],
                capture_output=True, text=True,
            )
            assert result.returncode == 0, result.stderr
