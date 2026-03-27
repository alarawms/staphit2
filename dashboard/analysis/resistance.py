"""Resistance profiles analysis tab for the Staphit surveillance dashboard."""

import json
from collections import Counter

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, dcc, html

from dashboard.utils import SIR_COLORS

# ---------------------------------------------------------------------------
# Gene-prefix to drug-class mapping (for colouring mutation bars)
# ---------------------------------------------------------------------------
GENE_DRUG_CLASS = {
    "gyrA": "Fluoroquinolones",
    "gyrB": "Fluoroquinolones",
    "grlA": "Fluoroquinolones",
    "grlB": "Fluoroquinolones",
    "rpoB": "Rifampicin",
    "rpsJ": "Tetracyclines",
    "dfrB": "Trimethoprim",
    "folA": "Trimethoprim",
    "fusA": "Fusidic acid",
    "ileS": "Mupirocin",
    "murA": "Fosfomycin",
    "pbp2": "Beta-lactams",
    "23S": "Linezolid",
    "rplC": "Linezolid",
    "rplD": "Linezolid",
    "rlmN": "Linezolid",
}

DRUG_CLASS_COLORS = {
    "Fluoroquinolones": "#636efa",
    "Rifampicin": "#ef553b",
    "Tetracyclines": "#00cc96",
    "Trimethoprim": "#ab63fa",
    "Fusidic acid": "#ffa15a",
    "Mupirocin": "#19d3f3",
    "Fosfomycin": "#ff6692",
    "Beta-lactams": "#b6e880",
    "Linezolid": "#ff97ff",
    "Other": "#888888",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _parse_point_mutations(series):
    """Parse semicolon-separated gene:mutation entries and return a Counter."""
    counts = Counter()
    for val in series.dropna():
        val = str(val).strip()
        if not val or val == "-":
            continue
        for item in val.split(";"):
            item = item.strip()
            if item:
                counts[item] += 1
    return counts


def _gene_prefix(mutation_str):
    """Extract gene prefix from a gene:mutation or gene_mutation string."""
    for sep in (":", "_"):
        if sep in mutation_str:
            return mutation_str.split(sep, 1)[0]
    return mutation_str


def _drug_class_for_mutation(mutation_str):
    """Return the drug class for a mutation based on its gene prefix."""
    prefix = _gene_prefix(mutation_str)
    return GENE_DRUG_CLASS.get(prefix, "Other")


def _empty_fig(message="No data available"):
    """Return a placeholder figure with a centred message."""
    fig = go.Figure()
    fig.add_annotation(
        text=message,
        xref="paper", yref="paper",
        x=0.5, y=0.5,
        showarrow=False,
        font=dict(size=16, color="#999"),
    )
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        margin=dict(l=20, r=20, t=30, b=20),
        height=300,
    )
    return fig


_GRAPH_LAYOUT = dict(
    template="plotly_dark",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    margin=dict(l=20, r=20, t=40, b=20),
    font=dict(size=11),
)


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


def create_resistance_layout(data):
    """Build the resistance profiles analysis tab layout.

    Parameters
    ----------
    data : dashboard.data_loader.DashboardData

    Returns
    -------
    dbc.Container
    """
    # Build antibiotic options from antibiogram data
    antibiotics = []
    default_antibiotic = None
    if not data.antibiogram.empty and "antibiotic" in data.antibiogram.columns:
        antibiotics = sorted(data.antibiogram["antibiotic"].dropna().unique().tolist())
        if antibiotics:
            default_antibiotic = antibiotics[0]

    return dbc.Container(
        [
            html.H5("Resistance Profiles", className="mt-3 mb-2"),
            # Row 1: MIC distribution
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.Label("Antibiotic", className="mb-1"),
                            dcc.Dropdown(
                                id="resistance-antibiotic-dropdown",
                                options=[{"label": a, "value": a} for a in antibiotics],
                                value=default_antibiotic,
                                clearable=False,
                                style={"color": "#222"},
                            ),
                        ],
                        width=3,
                    ),
                    dbc.Col(
                        dcc.Graph(
                            id="mic-distribution-histogram",
                            config={"displayModeBar": False},
                            style={"height": "320px"},
                        ),
                        width=9,
                    ),
                ],
                className="mb-3",
            ),
            # Row 2: Resistance rate bar chart & mutation frequency
            dbc.Row(
                [
                    dbc.Col(
                        dcc.Graph(
                            id="resistance-rate-bar",
                            config={"displayModeBar": False},
                            style={"height": "450px"},
                        ),
                        width=6,
                    ),
                    dbc.Col(
                        dcc.Graph(
                            id="mutation-frequency-bar",
                            config={"displayModeBar": False},
                            style={"height": "450px"},
                        ),
                        width=6,
                    ),
                ],
                className="mb-3",
            ),
            # Row 3: Co-resistance heatmap
            dbc.Row(
                [
                    dbc.Col(
                        dcc.Graph(
                            id="coresistance-heatmap",
                            config={"displayModeBar": False},
                            style={"height": "500px"},
                        ),
                        width=12,
                    ),
                ],
                className="mb-3",
            ),
        ],
        fluid=True,
    )


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------


def register_resistance_callbacks(app, data):
    """Register all resistance-tab callbacks.

    Parameters
    ----------
    app : dash.Dash
    data : dashboard.data_loader.DashboardData
    """

    # -- MIC distribution histogram -----------------------------------------
    @app.callback(
        Output("mic-distribution-histogram", "figure"),
        [
            Input("resistance-antibiotic-dropdown", "value"),
            Input("filtered-sample-ids", "data"),
        ],
    )
    def _update_mic_histogram(antibiotic, filtered_ids_json):
        if data.antibiogram.empty or not antibiotic:
            return _empty_fig("No antibiogram data")

        abg = data.antibiogram.copy()

        # Filter to selected antibiotic
        abg = abg[abg["antibiotic"] == antibiotic]

        # Apply sample filter
        if filtered_ids_json:
            try:
                ids = (
                    json.loads(filtered_ids_json)
                    if isinstance(filtered_ids_json, str)
                    else filtered_ids_json
                )
            except (json.JSONDecodeError, TypeError):
                ids = None
            if ids is not None and "sample_id" in abg.columns:
                id_set = set(str(i) for i in ids)
                abg = abg[abg["sample_id"].astype(str).isin(id_set)]

        if abg.empty or "measurement" not in abg.columns:
            return _empty_fig(f"No MIC data for {antibiotic}")

        # Drop rows without numeric measurement
        plot_df = abg.dropna(subset=["measurement"]).copy()
        plot_df["measurement"] = pd.to_numeric(plot_df["measurement"], errors="coerce")
        plot_df = plot_df.dropna(subset=["measurement"])

        if plot_df.empty:
            return _empty_fig(f"No numeric MIC values for {antibiotic}")

        # Ensure phenotype column exists
        if "resistance_phenotype" not in plot_df.columns:
            plot_df["resistance_phenotype"] = "Unknown"

        fig = px.histogram(
            plot_df,
            x="measurement",
            color="resistance_phenotype",
            color_discrete_map={**SIR_COLORS, "Unknown": "#888"},
            category_orders={"resistance_phenotype": ["S", "I", "R"]},
            nbins=30,
            title=f"MIC Distribution — {antibiotic}",
            labels={"measurement": "MIC", "resistance_phenotype": "Phenotype"},
        )
        fig.update_layout(**_GRAPH_LAYOUT, height=300)
        return fig

    # -- Resistance rate bar chart, mutation chart, co-resistance heatmap ---
    @app.callback(
        [
            Output("resistance-rate-bar", "figure"),
            Output("mutation-frequency-bar", "figure"),
            Output("coresistance-heatmap", "figure"),
        ],
        Input("filtered-sample-ids", "data"),
    )
    def _update_resistance_panels(filtered_ids_json):
        # ------ Resolve filtered sample IDs ------
        ids = None
        if filtered_ids_json:
            try:
                ids = (
                    json.loads(filtered_ids_json)
                    if isinstance(filtered_ids_json, str)
                    else filtered_ids_json
                )
            except (json.JSONDecodeError, TypeError):
                ids = None

        id_set = set(str(i) for i in ids) if ids else None

        # ==================================================================
        # 1. Resistance rate bar chart
        # ==================================================================
        rate_fig = _empty_fig("No antibiogram data")
        if not data.antibiogram.empty and "antibiotic" in data.antibiogram.columns:
            abg = data.antibiogram.copy()
            if id_set is not None and "sample_id" in abg.columns:
                abg = abg[abg["sample_id"].astype(str).isin(id_set)]

            if (
                not abg.empty
                and "resistance_phenotype" in abg.columns
            ):
                # Compute per-antibiotic SIR counts
                sir_counts = (
                    abg.groupby(["antibiotic", "resistance_phenotype"])
                    .size()
                    .reset_index(name="count")
                )
                total_per_ab = (
                    abg.groupby("antibiotic")
                    .size()
                    .reset_index(name="total")
                )
                # Filter to antibiotics with >=5 samples
                total_per_ab = total_per_ab[total_per_ab["total"] >= 5]
                valid_abs = set(total_per_ab["antibiotic"])
                sir_counts = sir_counts[sir_counts["antibiotic"].isin(valid_abs)]

                if not sir_counts.empty:
                    sir_counts = sir_counts.merge(total_per_ab, on="antibiotic")
                    sir_counts["pct"] = 100.0 * sir_counts["count"] / sir_counts["total"]

                    # Sort antibiotics by %R descending
                    r_pct = sir_counts[sir_counts["resistance_phenotype"] == "R"].set_index(
                        "antibiotic"
                    )["pct"]
                    ab_order = r_pct.sort_values(ascending=True).index.tolist()
                    # Include antibiotics with 0% R
                    all_abs = sir_counts["antibiotic"].unique().tolist()
                    for ab in all_abs:
                        if ab not in ab_order:
                            ab_order.insert(0, ab)

                    rate_fig = px.bar(
                        sir_counts,
                        y="antibiotic",
                        x="pct",
                        color="resistance_phenotype",
                        color_discrete_map=SIR_COLORS,
                        category_orders={
                            "antibiotic": ab_order,
                            "resistance_phenotype": ["R", "I", "S"],
                        },
                        orientation="h",
                        title="Resistance Rates by Antibiotic",
                        labels={
                            "pct": "% of Isolates",
                            "antibiotic": "Antibiotic",
                            "resistance_phenotype": "Phenotype",
                        },
                    )
                    rate_fig.update_layout(
                        **_GRAPH_LAYOUT,
                        barmode="stack",
                        height=max(350, len(ab_order) * 22 + 80),
                        legend=dict(orientation="h", yanchor="bottom", y=1.02),
                    )

        # ==================================================================
        # 2. Point mutation frequency
        # ==================================================================
        mut_fig = _empty_fig("No point mutation data")
        if not data.merged.empty and "point_mutations" in data.merged.columns:
            merged = data.merged.copy()
            if id_set is not None and "sample_id" in merged.columns:
                merged = merged[merged["sample_id"].astype(str).isin(id_set)]

            if not merged.empty:
                mut_counts = _parse_point_mutations(merged["point_mutations"])
                if mut_counts:
                    top20 = mut_counts.most_common(20)
                    mutations = [m for m, _ in top20]
                    counts = [c for _, c in top20]
                    drug_classes = [_drug_class_for_mutation(m) for m in mutations]

                    mut_df = pd.DataFrame(
                        {
                            "mutation": mutations,
                            "count": counts,
                            "drug_class": drug_classes,
                        }
                    )
                    # Reverse so highest is on top in horizontal bar
                    mut_df = mut_df.iloc[::-1]

                    mut_fig = px.bar(
                        mut_df,
                        y="mutation",
                        x="count",
                        color="drug_class",
                        color_discrete_map=DRUG_CLASS_COLORS,
                        orientation="h",
                        title="Top Point Mutations",
                        labels={
                            "count": "Frequency",
                            "mutation": "Mutation",
                            "drug_class": "Drug Class",
                        },
                    )
                    mut_fig.update_layout(
                        **_GRAPH_LAYOUT,
                        height=max(350, len(mutations) * 22 + 80),
                        legend=dict(orientation="h", yanchor="bottom", y=1.02),
                    )

        # ==================================================================
        # 3. Co-resistance heatmap
        # ==================================================================
        heat_fig = _empty_fig("No co-resistance data")
        if not data.antibiogram.empty and "antibiotic" in data.antibiogram.columns:
            abg = data.antibiogram.copy()
            if id_set is not None and "sample_id" in abg.columns:
                abg = abg[abg["sample_id"].astype(str).isin(id_set)]

            if (
                not abg.empty
                and "resistance_phenotype" in abg.columns
                and "sample_id" in abg.columns
            ):
                # Encode R=1, S=0, I=0.5 for correlation
                encode = {"R": 1, "I": 0.5, "S": 0}
                abg["r_val"] = abg["resistance_phenotype"].map(encode)
                abg = abg.dropna(subset=["r_val"])

                # Pivot: samples x antibiotics
                pivot = abg.pivot_table(
                    index="sample_id",
                    columns="antibiotic",
                    values="r_val",
                    aggfunc="max",
                )

                # Only keep antibiotics with >=5 non-null samples
                col_counts = pivot.notna().sum()
                valid_cols = col_counts[col_counts >= 5].index.tolist()
                pivot = pivot[valid_cols]

                if len(valid_cols) >= 2:
                    corr = pivot.corr()
                    heat_fig = px.imshow(
                        corr,
                        x=corr.columns.tolist(),
                        y=corr.index.tolist(),
                        color_continuous_scale="RdBu_r",
                        zmin=-1,
                        zmax=1,
                        title="Co-resistance Correlation",
                        labels={"color": "Correlation"},
                    )
                    heat_fig.update_layout(
                        **_GRAPH_LAYOUT,
                        height=max(400, len(valid_cols) * 25 + 100),
                    )

        return rate_fig, mut_fig, heat_fig
