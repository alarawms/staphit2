"""Tests for dashboard.data_loader and dashboard.utils."""

import json
import os
import sys
import tempfile

import pandas as pd
import pytest

# Ensure the project root is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dashboard.data_loader import DashboardData, load_all
from dashboard.utils import count_frequencies, parse_virulence_summary


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tmpdir_with_summary(tmp_path):
    """Create a minimal results directory with a 3-row final_summary.tsv."""
    agg = tmp_path / "aggregated"
    agg.mkdir()

    summary = pd.DataFrame(
        {
            "sample_id": ["S001", "S002", "S003"],
            "mlst_st": ["ST5", "ST8", "ST22"],
            "assembly_length": [2800000, 2750000, 2900000],
            "n50": [300000, 250000, 280000],
            "virulence_summary": [
                "PVL+;TSST-;IEC-B;ica:complete;hlg:complete",
                "PVL-;TSST+;IEC-E;ica:partial;hlg:complete",
                "PVL-;TSST-;IEC-B;ica:complete;hlg:partial",
            ],
        }
    )
    summary.to_csv(agg / "final_summary.tsv", sep="\t", index=False)

    return tmp_path


@pytest.fixture
def tmpdir_with_metadata(tmpdir_with_summary):
    """Add a sample_metadata.csv alongside the results directory."""
    meta = pd.DataFrame(
        {
            "sample_id": ["S001", "S002", "S004"],
            "collection_date": ["2024-01-01", "2024-02-15", "2024-03-10"],
            "geo_loc_region": ["Riyadh", "Jeddah", "Dammam"],
        }
    )
    meta_path = tmpdir_with_summary / "sample_metadata.csv"
    meta.to_csv(meta_path, index=False)
    return tmpdir_with_summary, str(meta_path)


@pytest.fixture
def tmpdir_with_antibiogram(tmp_path):
    """Create a minimal antibiogram CSV."""
    abg = pd.DataFrame(
        {
            "sample_id": ["S001", "S001", "S002"],
            "antibiotic": ["Cefoxitin", "Vancomycin", "Cefoxitin"],
            "resistance_phenotype": ["R", "S", "R"],
        }
    )
    abg_path = tmp_path / "antibiogram.csv"
    abg.to_csv(abg_path, index=False)
    return tmp_path, str(abg_path)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestLoadsSummary:
    """Test loading final_summary.tsv."""

    def test_loads_summary(self, tmpdir_with_summary):
        data = load_all(str(tmpdir_with_summary))
        assert len(data.summary) == 3
        assert "sample_id" in data.summary.columns
        assert "mlst_st" in data.summary.columns
        assert len(data.merged) == 3


class TestLoadsAntibiogram:
    """Test loading antibiogram.csv."""

    def test_loads_antibiogram(self, tmpdir_with_antibiogram):
        tmp_path, abg_path = tmpdir_with_antibiogram

        # Need a results dir with at least aggregated/
        agg = tmp_path / "results" / "aggregated"
        agg.mkdir(parents=True)
        summary = pd.DataFrame({"sample_id": ["S001"]})
        summary.to_csv(agg.parent / "aggregated" / "final_summary.tsv", sep="\t", index=False)

        data = load_all(str(tmp_path / "results"), antibiogram_path=abg_path)
        assert len(data.antibiogram) == 3
        assert "antibiotic" in data.antibiogram.columns
        assert "resistance_phenotype" in data.antibiogram.columns


class TestMergesMetadata:
    """Test that summary and metadata are merged on sample_id."""

    def test_merges_metadata(self, tmpdir_with_metadata):
        results_dir, meta_path = tmpdir_with_metadata
        data = load_all(str(results_dir), metadata_path=meta_path)

        # Merged should have summary columns AND metadata columns
        assert "collection_date" in data.merged.columns
        assert "geo_loc_region" in data.merged.columns
        assert "mlst_st" in data.merged.columns

        # Left join: should have 3 rows (from summary), S004 not in summary
        assert len(data.merged) == 3

        # S001 should have metadata
        s001 = data.merged[data.merged["sample_id"] == "S001"].iloc[0]
        assert s001["geo_loc_region"] == "Riyadh"

        # S003 not in metadata, so should be NaN
        s003 = data.merged[data.merged["sample_id"] == "S003"].iloc[0]
        assert pd.isna(s003["collection_date"])


class TestHandlesMissingFiles:
    """Test that load_all handles a completely empty directory gracefully."""

    def test_handles_missing_files(self, tmp_path):
        data = load_all(str(tmp_path))

        assert isinstance(data, DashboardData)
        assert data.summary.empty
        assert data.metadata.empty
        assert data.antibiogram.empty
        assert data.clusters == []
        assert data.cluster_assignments.empty
        assert data.snp_distances.empty
        assert data.tree_newick is None
        assert data.qc.empty
        assert data.plasmids.empty
        assert data.sample_reports == {}
        assert data.merged.empty


class TestParsesVirulenceSummary:
    """Test parse_virulence_summary from utils."""

    def test_parses_virulence_summary(self):
        df = pd.DataFrame(
            {
                "virulence_summary": [
                    "PVL+;TSST-;IEC-B;ica:complete",
                    "PVL-;TSST+;IEC-E;ica:partial;hlg:complete",
                    None,
                ]
            }
        )
        result = parse_virulence_summary(df)

        # Row 0: PVL+, TSST-, IEC-B, ica:complete
        assert result.loc[0, "pvl_positive"] == True
        assert result.loc[0, "tsst_positive"] == False
        assert result.loc[0, "iec_type"] == "B"
        assert result.loc[0, "ica_status"] == "complete"

        # Row 1: PVL-, TSST+, IEC-E, ica:partial, hlg:complete
        assert result.loc[1, "pvl_positive"] == False
        assert result.loc[1, "tsst_positive"] == True
        assert result.loc[1, "iec_type"] == "E"
        assert result.loc[1, "ica_status"] == "partial"
        assert result.loc[1, "hlg_status"] == "complete"

        # Row 2: None -> all defaults
        assert result.loc[2, "pvl_positive"] == False
        assert result.loc[2, "tsst_positive"] == False
        assert pd.isna(result.loc[2, "iec_type"])
        assert pd.isna(result.loc[2, "ica_status"])


class TestCountFrequencies:
    """Test count_frequencies helper."""

    def test_count_frequencies(self):
        series = pd.Series(["geneA;geneB", "geneB;geneC", None, "-", "geneA"])
        result = count_frequencies(series)
        assert result["geneA"] == 2
        assert result["geneB"] == 2
        assert result["geneC"] == 1
        assert "-" not in result
