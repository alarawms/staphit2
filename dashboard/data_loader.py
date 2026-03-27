"""Data loader for the Staphit surveillance dashboard.

Loads all pipeline outputs from a results directory into a single
``DashboardData`` dataclass for use by the dashboard views.
"""

import json
import os
import sys
from dataclasses import dataclass, field
from glob import glob

import pandas as pd

from dashboard.utils import parse_virulence_summary


@dataclass
class DashboardData:
    """Container for all dashboard data sources."""

    summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    metadata: pd.DataFrame = field(default_factory=pd.DataFrame)
    antibiogram: pd.DataFrame = field(default_factory=pd.DataFrame)
    clusters: list = field(default_factory=list)
    cluster_assignments: pd.DataFrame = field(default_factory=pd.DataFrame)
    snp_distances: pd.DataFrame = field(default_factory=pd.DataFrame)
    tree_newick: str = None
    qc: pd.DataFrame = field(default_factory=pd.DataFrame)
    plasmids: pd.DataFrame = field(default_factory=pd.DataFrame)
    sample_reports: dict = field(default_factory=dict)
    merged: pd.DataFrame = field(default_factory=pd.DataFrame)


def load_all(results_dir, metadata_path=None, antibiogram_path=None):
    """Load all pipeline outputs and return a ``DashboardData`` instance.

    Parameters
    ----------
    results_dir : str
        Path to the pipeline results directory (e.g. ``test_hail/``).
    metadata_path : str, optional
        Path to ``sample_metadata.csv``.  Falls back to
        ``<results_dir>/../sample_metadata.csv``.
    antibiogram_path : str, optional
        Path to ``antibiogram.csv``.  Falls back to
        ``<results_dir>/../antibiogram.csv``.

    Returns
    -------
    DashboardData
    """
    data = DashboardData()

    # -- final_summary.tsv --------------------------------------------------
    try:
        path = os.path.join(results_dir, "aggregated", "final_summary.tsv")
        if os.path.isfile(path):
            data.summary = pd.read_csv(path, sep="\t")
    except Exception as exc:
        print(f"Warning: failed to load summary: {exc}", file=sys.stderr)

    # -- sample_metadata.csv ------------------------------------------------
    try:
        if metadata_path and os.path.isfile(metadata_path):
            data.metadata = pd.read_csv(metadata_path)
        else:
            fallback = os.path.join(results_dir, "..", "sample_metadata.csv")
            if os.path.isfile(fallback):
                data.metadata = pd.read_csv(fallback)
    except Exception as exc:
        print(f"Warning: failed to load metadata: {exc}", file=sys.stderr)

    # -- antibiogram.csv ----------------------------------------------------
    try:
        if antibiogram_path and os.path.isfile(antibiogram_path):
            data.antibiogram = pd.read_csv(antibiogram_path)
        else:
            fallback = os.path.join(results_dir, "..", "antibiogram.csv")
            if os.path.isfile(fallback):
                data.antibiogram = pd.read_csv(fallback)
    except Exception as exc:
        print(f"Warning: failed to load antibiogram: {exc}", file=sys.stderr)

    # -- cluster_report.json ------------------------------------------------
    try:
        path = os.path.join(results_dir, "clusters", "cluster_report.json")
        if os.path.isfile(path):
            with open(path) as fh:
                data.clusters = json.load(fh)
    except Exception as exc:
        print(f"Warning: failed to load cluster report: {exc}", file=sys.stderr)

    # -- clusters.tsv -------------------------------------------------------
    try:
        path = os.path.join(results_dir, "clusters", "clusters.tsv")
        if os.path.isfile(path):
            data.cluster_assignments = pd.read_csv(path, sep="\t")
    except Exception as exc:
        print(f"Warning: failed to load cluster assignments: {exc}", file=sys.stderr)

    # -- snp_distances.tsv --------------------------------------------------
    try:
        path = os.path.join(results_dir, "snp_dists", "snp_distances.tsv")
        if os.path.isfile(path):
            data.snp_distances = pd.read_csv(path, sep="\t")
    except Exception as exc:
        print(f"Warning: failed to load SNP distances: {exc}", file=sys.stderr)

    # -- *.treefile (IQ-TREE newick) ----------------------------------------
    try:
        pattern = os.path.join(results_dir, "iqtree", "*.treefile")
        treefiles = glob(pattern)
        if treefiles:
            with open(treefiles[0]) as fh:
                data.tree_newick = fh.read().strip()
    except Exception as exc:
        print(f"Warning: failed to load tree: {exc}", file=sys.stderr)

    # -- qc_report.tsv ------------------------------------------------------
    try:
        path = os.path.join(results_dir, "qc_gate", "qc_report.tsv")
        if os.path.isfile(path):
            data.qc = pd.read_csv(path, sep="\t")
    except Exception as exc:
        print(f"Warning: failed to load QC report: {exc}", file=sys.stderr)

    # -- plasmid_summary.tsv ------------------------------------------------
    try:
        path = os.path.join(results_dir, "mob_recon", "plasmid_summary.tsv")
        if os.path.isfile(path):
            data.plasmids = pd.read_csv(path, sep="\t")
    except Exception as exc:
        print(f"Warning: failed to load plasmid summary: {exc}", file=sys.stderr)

    # -- sample report JSON paths -------------------------------------------
    try:
        pattern = os.path.join(results_dir, "aggregated", "*_report.json")
        for rpath in glob(pattern):
            basename = os.path.basename(rpath)
            sample_id = basename.replace("_report.json", "")
            data.sample_reports[sample_id] = rpath
    except Exception as exc:
        print(f"Warning: failed to index sample reports: {exc}", file=sys.stderr)

    # -- Merge summary + metadata -------------------------------------------
    try:
        if not data.summary.empty:
            if not data.metadata.empty and "sample_id" in data.metadata.columns:
                data.merged = data.summary.merge(
                    data.metadata, on="sample_id", how="left", suffixes=("", "_meta")
                )
            else:
                data.merged = data.summary.copy()
        else:
            data.merged = pd.DataFrame()
    except Exception as exc:
        print(f"Warning: failed to merge summary + metadata: {exc}", file=sys.stderr)
        data.merged = data.summary.copy() if not data.summary.empty else pd.DataFrame()

    # -- Parse virulence_summary into boolean columns -----------------------
    try:
        if not data.merged.empty and "virulence_summary" in data.merged.columns:
            parse_virulence_summary(data.merged)
    except Exception as exc:
        print(f"Warning: failed to parse virulence summary: {exc}", file=sys.stderr)

    # -- Print load summary -------------------------------------------------
    n_samples = len(data.merged)
    n_abg = len(data.antibiogram)
    n_clusters = len(data.clusters)
    print(
        f"Loaded {n_samples} samples, {n_abg} antibiogram rows, {n_clusters} clusters",
        file=sys.stderr,
    )

    return data
