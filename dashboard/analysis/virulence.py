"""Virulence analysis tab for the Staphit surveillance dashboard."""

import json

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, dcc, html


def _parse_enterotoxin_count(row):
    """Count enterotoxin genes from virulence_summary or virulence_genes."""
    count = 0
    for col in ("virulence_genes", "virulence_summary"):
        val = row.get(col)
        if isinstance(val, str):
            tokens = val.split(";")
            for t in tokens:
                t_upper = t.strip().upper()
                if t_upper.startswith("SE") and len(t_upper) <= 5:
                    count += 1
                elif "ENTEROTOXIN" in t_upper:
                    count += 1
    return count


def _parse_exfoliatin(row):
    """Check for exfoliatin genes from virulence_summary or virulence_genes."""
    for col in ("virulence_genes", "virulence_summary"):
        val = row.get(col)
        if isinstance(val, str):
            upper = val.upper()
            if "ETA" in upper or "ETB" in upper or "EXFOLIATIN" in upper:
                return True
    return False


def create_virulence_layout(data):
    """Build the virulence tab layout.

    Parameters
    ----------
    data : dashboard.data_loader.DashboardData

    Returns
    -------
    html.Div
    """
    return html.Div(
        [
            html.H5("Virulence Factors", className="mt-3 mb-3"),
            # Row 1: Toxin prevalence + IEC donut
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                dcc.Graph(id="vir-toxin-bar", config={"displayModeBar": False})
                            ),
                        ),
                        md=4,
                    ),
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                dcc.Graph(id="vir-iec-donut", config={"displayModeBar": False})
                            ),
                        ),
                        md=4,
                    ),
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                dcc.Graph(id="vir-operon-bar", config={"displayModeBar": False})
                            ),
                        ),
                        md=4,
                    ),
                ],
                className="mb-3",
            ),
            # Row 2: Heatmap + enterotoxin histogram
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                dcc.Graph(id="vir-st-heatmap", config={"displayModeBar": False})
                            ),
                        ),
                        md=7,
                    ),
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                dcc.Graph(id="vir-entero-hist", config={"displayModeBar": False})
                            ),
                        ),
                        md=5,
                    ),
                ],
                className="mb-3",
            ),
        ]
    )


def register_virulence_callbacks(app, data):
    """Register callbacks for the virulence tab.

    Parameters
    ----------
    app : dash.Dash
    data : dashboard.data_loader.DashboardData
    """
    merged = data.merged

    @app.callback(
        [
            Output("vir-toxin-bar", "figure"),
            Output("vir-iec-donut", "figure"),
            Output("vir-operon-bar", "figure"),
            Output("vir-st-heatmap", "figure"),
            Output("vir-entero-hist", "figure"),
        ],
        Input("filtered-sample-ids", "data"),
    )
    def _update_virulence(filtered_ids_json):
        dark_template = "plotly_dark"

        # Parse filtered IDs
        try:
            ids = json.loads(filtered_ids_json) if isinstance(filtered_ids_json, str) else filtered_ids_json
        except (json.JSONDecodeError, TypeError):
            ids = None

        if ids is not None and not merged.empty:
            id_set = set(str(i) for i in ids)
            df = merged[merged["sample_id"].astype(str).isin(id_set)].copy()
        else:
            df = merged.copy() if not merged.empty else pd.DataFrame()

        # --- Chart 1: Toxin prevalence bar ---
        toxin_fig = go.Figure()
        toxin_fig.update_layout(
            template=dark_template,
            title="Toxin Prevalence",
            xaxis_title="Toxin",
            yaxis_title="Count",
        )

        if not df.empty:
            pvl_count = int(df["pvl_positive"].sum()) if "pvl_positive" in df.columns else 0
            tsst_count = int(df["tsst_positive"].sum()) if "tsst_positive" in df.columns else 0
            exfol_count = int(df.apply(_parse_exfoliatin, axis=1).sum()) if len(df) > 0 else 0

            toxin_fig.add_trace(go.Bar(
                x=["PVL+", "TSST+", "Exfoliatin+"],
                y=[pvl_count, tsst_count, exfol_count],
                marker_color=["#e74c3c", "#f39c12", "#9b59b6"],
            ))

        # --- Chart 2: IEC type donut ---
        iec_fig = go.Figure()
        iec_fig.update_layout(template=dark_template, title="IEC Type Distribution")

        if not df.empty and "iec_type" in df.columns:
            iec_counts = df["iec_type"].fillna("None").value_counts().reset_index()
            iec_counts.columns = ["iec_type", "count"]
            if len(iec_counts) > 0:
                iec_fig = px.pie(
                    iec_counts,
                    names="iec_type",
                    values="count",
                    hole=0.4,
                    title="IEC Type Distribution",
                    template=dark_template,
                )

        # --- Chart 3: Operon status stacked bar ---
        operon_fig = go.Figure()
        operon_fig.update_layout(
            template=dark_template,
            title="Operon Status",
            xaxis_title="Operon",
            yaxis_title="Count",
            barmode="stack",
        )

        if not df.empty:
            statuses = ["complete", "partial", "absent"]
            status_colors = {"complete": "#2ecc71", "partial": "#f39c12", "absent": "#e74c3c"}

            for operon_col, operon_label in [("ica_status", "ica"), ("hlg_status", "hlg")]:
                if operon_col in df.columns:
                    val_counts = df[operon_col].fillna("absent").value_counts()
                    for status in statuses:
                        count = int(val_counts.get(status, 0))
                        operon_fig.add_trace(go.Bar(
                            x=[operon_label],
                            y=[count],
                            name=status,
                            marker_color=status_colors.get(status, "#888"),
                            showlegend=(operon_label == "ica"),  # legend once
                        ))

        # --- Chart 4: Virulence by ST heatmap ---
        heatmap_fig = go.Figure()
        heatmap_fig.update_layout(template=dark_template, title="Virulence by ST")

        if not df.empty and "mlst_st" in df.columns:
            top10_sts = df["mlst_st"].value_counts().head(10).index.tolist()
            df_top = df[df["mlst_st"].isin(top10_sts)].copy()

            if len(df_top) > 0:
                vir_cols = {
                    "PVL": "pvl_positive",
                    "TSST": "tsst_positive",
                    "IEC+": "iec_type",
                    "ica_complete": "ica_status",
                    "hlg_complete": "hlg_status",
                }

                heatmap_data = []
                for st in top10_sts:
                    st_df = df_top[df_top["mlst_st"] == st]
                    n = len(st_df)
                    row = {}
                    for label, col in vir_cols.items():
                        if col not in st_df.columns:
                            row[label] = 0.0
                        elif col in ("pvl_positive", "tsst_positive"):
                            row[label] = round(st_df[col].sum() / n * 100, 1) if n > 0 else 0.0
                        elif col == "iec_type":
                            row[label] = round(st_df[col].notna().sum() / n * 100, 1) if n > 0 else 0.0
                        elif col in ("ica_status", "hlg_status"):
                            row[label] = round((st_df[col] == "complete").sum() / n * 100, 1) if n > 0 else 0.0
                    heatmap_data.append(row)

                z = [[heatmap_data[i][c] for c in vir_cols.keys()] for i in range(len(top10_sts))]
                st_labels = [f"ST{st}" for st in top10_sts]

                heatmap_fig = px.imshow(
                    z,
                    x=list(vir_cols.keys()),
                    y=st_labels,
                    color_continuous_scale="YlOrRd",
                    labels=dict(color="% Positive"),
                    title="Virulence by ST (% Positive)",
                    template=dark_template,
                    aspect="auto",
                )

        # --- Chart 5: Enterotoxin count distribution ---
        entero_fig = go.Figure()
        entero_fig.update_layout(
            template=dark_template,
            title="Enterotoxin Count Distribution",
            xaxis_title="Enterotoxin Count",
            yaxis_title="Samples",
        )

        if not df.empty:
            entero_counts = df.apply(_parse_enterotoxin_count, axis=1)
            if entero_counts.sum() > 0:
                entero_fig = px.histogram(
                    x=entero_counts,
                    nbins=max(int(entero_counts.max()) + 1, 5),
                    title="Enterotoxin Count Distribution",
                    labels={"x": "Enterotoxin Count", "y": "Samples"},
                    template=dark_template,
                )
            else:
                entero_fig.add_annotation(
                    text="No enterotoxin data detected",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font=dict(size=14, color="#888"),
                )

        return toxin_fig, iec_fig, operon_fig, heatmap_fig, entero_fig
