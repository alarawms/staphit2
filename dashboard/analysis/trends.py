"""Trends analysis tab for the Staphit surveillance dashboard."""

import json

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, dcc, html

from dashboard.utils import SIR_COLORS, ST_COLORS


def _extract_year(date_series):
    """Parse collection_date to year, handling both 'YYYY' and 'YYYY-MM-DD'."""
    return pd.to_numeric(date_series.astype(str).str[:4], errors="coerce")


def _get_antibiotic_columns(df):
    """Return columns that look like antibiotic SIR results (values are R/I/S)."""
    candidates = []
    for col in df.columns:
        if col.startswith("sir_") or col.endswith("_sir"):
            candidates.append(col)
            continue
        vals = df[col].dropna().astype(str).str.upper().unique()
        if len(vals) > 0 and set(vals).issubset({"R", "I", "S", ""}):
            if len(set(vals) & {"R", "I", "S"}) >= 2:
                candidates.append(col)
    return candidates


def create_trends_layout(data):
    """Build the trends tab layout.

    Parameters
    ----------
    data : dashboard.data_loader.DashboardData

    Returns
    -------
    html.Div
    """
    merged = data.merged
    abg = data.antibiogram

    # Build antibiotic dropdown options from antibiogram or merged columns
    abx_options = []
    if not abg.empty and "antibiotic" in abg.columns:
        abx_list = sorted(abg["antibiotic"].dropna().unique().tolist())
        abx_options = [{"label": a, "value": a} for a in abx_list]
    else:
        sir_cols = _get_antibiotic_columns(merged) if not merged.empty else []
        abx_options = [{"label": c.replace("sir_", "").replace("_sir", ""), "value": c} for c in sir_cols]

    default_abx = abx_options[0]["value"] if abx_options else None

    return html.Div(
        [
            html.H5("Trends Over Time", className="mt-3 mb-3"),
            # Row 1: Resistance over time
            dbc.Card(
                dbc.CardBody(
                    [
                        dbc.Row(
                            [
                                dbc.Col(
                                    html.Label("Antibiotic:", className="small text-muted"),
                                    width="auto",
                                ),
                                dbc.Col(
                                    dcc.Dropdown(
                                        id="trends-abx-dropdown",
                                        options=abx_options,
                                        value=default_abx,
                                        placeholder="Select antibiotic",
                                        className="dash-bootstrap",
                                    ),
                                    md=4,
                                ),
                            ],
                            className="mb-2 align-items-center",
                        ),
                        dcc.Graph(id="trends-resistance-line", config={"displayModeBar": False}),
                    ]
                ),
                className="mb-3",
            ),
            # Row 2: ST frequency + accumulation
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                dcc.Graph(id="trends-st-area", config={"displayModeBar": False})
                            ),
                        ),
                        md=6,
                    ),
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                dcc.Graph(id="trends-accumulation", config={"displayModeBar": False})
                            ),
                        ),
                        md=6,
                    ),
                ],
                className="mb-3",
            ),
        ]
    )


def register_trends_callbacks(app, data):
    """Register callbacks for the trends tab.

    Parameters
    ----------
    app : dash.Dash
    data : dashboard.data_loader.DashboardData
    """
    merged = data.merged
    abg = data.antibiogram

    @app.callback(
        [
            Output("trends-resistance-line", "figure"),
            Output("trends-st-area", "figure"),
            Output("trends-accumulation", "figure"),
        ],
        [
            Input("filtered-sample-ids", "data"),
            Input("trends-abx-dropdown", "value"),
        ],
    )
    def _update_trends(filtered_ids_json, selected_abx):
        dark_template = "plotly_dark"

        # Parse filtered IDs
        try:
            ids = json.loads(filtered_ids_json) if isinstance(filtered_ids_json, str) else filtered_ids_json
        except (json.JSONDecodeError, TypeError):
            ids = None

        if ids is not None and merged is not None and not merged.empty:
            id_set = set(str(i) for i in ids)
            df = merged[merged["sample_id"].astype(str).isin(id_set)].copy()
        else:
            df = merged.copy() if not merged.empty else pd.DataFrame()

        # --- Chart 1: Resistance over time ---
        res_fig = go.Figure()
        res_fig.update_layout(
            template=dark_template,
            title="Resistance Rate Over Time",
            xaxis_title="Year",
            yaxis_title="% Resistant",
            yaxis=dict(range=[0, 100]),
        )

        if selected_abx and not df.empty and "collection_date" in df.columns:
            df_r = df.copy()
            df_r["year"] = _extract_year(df_r["collection_date"])
            df_r = df_r.dropna(subset=["year"])
            df_r["year"] = df_r["year"].astype(int)

            # Get resistance data from antibiogram or merged columns
            if not abg.empty and "antibiotic" in abg.columns and "sample_id" in abg.columns:
                abg_sub = abg[abg["antibiotic"] == selected_abx].copy()
                sir_col = "sir" if "sir" in abg_sub.columns else "resistance_phenotype" if "resistance_phenotype" in abg_sub.columns else None
                if sir_col and not abg_sub.empty:
                    df_r = df_r.merge(
                        abg_sub[["sample_id", sir_col]].rename(columns={sir_col: "_sir"}),
                        on="sample_id",
                        how="inner",
                    )
                    yearly = df_r.groupby("year")["_sir"].apply(
                        lambda x: (x.str.upper() == "R").sum() / len(x) * 100 if len(x) > 0 else 0
                    ).reset_index(name="pct_r")
                    yearly = yearly.sort_values("year")

                    res_fig.add_trace(go.Scatter(
                        x=yearly["year"],
                        y=yearly["pct_r"],
                        mode="lines+markers",
                        name=selected_abx,
                        line=dict(color=SIR_COLORS["R"]),
                        marker=dict(size=8),
                    ))
            elif selected_abx in df_r.columns:
                yearly = df_r.groupby("year")[selected_abx].apply(
                    lambda x: (x.str.upper() == "R").sum() / len(x) * 100 if len(x) > 0 else 0
                ).reset_index(name="pct_r")
                yearly = yearly.sort_values("year")

                res_fig.add_trace(go.Scatter(
                    x=yearly["year"],
                    y=yearly["pct_r"],
                    mode="lines+markers",
                    name=selected_abx,
                    line=dict(color=SIR_COLORS["R"]),
                    marker=dict(size=8),
                ))

        # --- Chart 2: ST frequency over time (stacked area) ---
        st_fig = go.Figure()
        st_fig.update_layout(
            template=dark_template,
            title="ST Frequency Over Time",
            xaxis_title="Year",
            yaxis_title="Count",
        )

        if not df.empty and "collection_date" in df.columns and "mlst_st" in df.columns:
            df_st = df.copy()
            df_st["year"] = _extract_year(df_st["collection_date"])
            df_st = df_st.dropna(subset=["year"])
            df_st["year"] = df_st["year"].astype(int)

            top10_sts = df_st["mlst_st"].value_counts().head(10).index.tolist()
            df_st["st_group"] = df_st["mlst_st"].apply(
                lambda x: str(x) if x in top10_sts else "Other"
            )
            pivot = df_st.groupby(["year", "st_group"]).size().unstack(fill_value=0)
            # Ensure "Other" is last
            cols = [c for c in pivot.columns if c != "Other"]
            if "Other" in pivot.columns:
                cols.append("Other")
            pivot = pivot[cols]

            for col in pivot.columns:
                color = ST_COLORS.get(f"ST{col}", ST_COLORS.get(str(col), None))
                if col == "Other":
                    color = "#888888"
                st_fig.add_trace(go.Scatter(
                    x=pivot.index,
                    y=pivot[col],
                    name=str(col),
                    mode="lines",
                    stackgroup="one",
                    line=dict(color=color) if color else {},
                ))

        # --- Chart 3: Sample accumulation by region ---
        acc_fig = go.Figure()
        acc_fig.update_layout(
            template=dark_template,
            title="Sample Accumulation by Region",
            xaxis_title="Year",
            yaxis_title="Cumulative Count",
        )

        if not df.empty and "collection_date" in df.columns:
            df_acc = df.copy()
            df_acc["year"] = _extract_year(df_acc["collection_date"])
            df_acc = df_acc.dropna(subset=["year"])
            df_acc["year"] = df_acc["year"].astype(int)

            region_col = "geo_loc_region" if "geo_loc_region" in df_acc.columns else None
            if region_col:
                for region in sorted(df_acc[region_col].dropna().unique()):
                    region_data = df_acc[df_acc[region_col] == region]
                    yearly_counts = region_data.groupby("year").size().sort_index().cumsum().reset_index(name="cumulative")
                    acc_fig.add_trace(go.Scatter(
                        x=yearly_counts["year"],
                        y=yearly_counts["cumulative"],
                        mode="lines+markers",
                        name=str(region),
                    ))
            else:
                yearly_counts = df_acc.groupby("year").size().sort_index().cumsum().reset_index(name="cumulative")
                acc_fig.add_trace(go.Scatter(
                    x=yearly_counts["year"],
                    y=yearly_counts["cumulative"],
                    mode="lines+markers",
                    name="All samples",
                ))

        return res_fig, st_fig, acc_fig
