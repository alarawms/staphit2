"""Plasmids analysis tab for the Staphit surveillance dashboard."""

import json
from collections import Counter

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, dcc, html


def _parse_mobility_entries(mobility_series):
    """Parse plasmid_mobility column into mobility and replicon counts.

    Expected format per cell: ``"AA411:mobilizable;AE018:non-mobilizable"``

    Returns
    -------
    tuple[Counter, Counter]
        (mobility_counts, replicon_counts)
    """
    mobility_counts = Counter()
    replicon_counts = Counter()

    for val in mobility_series.dropna():
        val = str(val).strip()
        if not val or val == "-":
            continue
        for entry in val.split(";"):
            entry = entry.strip()
            if not entry:
                continue
            parts = entry.split(":", 1)
            if len(parts) == 2:
                replicon = parts[0].strip()
                mobility = parts[1].strip().lower()
                mobility_counts[mobility] += 1
                if replicon:
                    replicon_counts[replicon] += 1
            else:
                # Single value -- treat as mobility type
                mobility_counts[entry.strip().lower()] += 1

    return mobility_counts, replicon_counts


def create_plasmids_layout(data):
    """Build the plasmids tab layout.

    Parameters
    ----------
    data : dashboard.data_loader.DashboardData

    Returns
    -------
    html.Div
    """
    plasmids = data.plasmids

    if plasmids.empty:
        return html.Div(
            [
                html.H5("Plasmid Analysis", className="mt-3 mb-3"),
                dbc.Alert(
                    "No plasmid data available. Run the pipeline with mob_recon enabled.",
                    color="info",
                    className="mt-3",
                ),
            ]
        )

    return html.Div(
        [
            html.H5("Plasmid Analysis", className="mt-3 mb-3"),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                dcc.Graph(id="plasmid-count-hist", config={"displayModeBar": False})
                            ),
                        ),
                        md=4,
                    ),
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                dcc.Graph(id="plasmid-mobility-pie", config={"displayModeBar": False})
                            ),
                        ),
                        md=4,
                    ),
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                dcc.Graph(id="plasmid-replicon-bar", config={"displayModeBar": False})
                            ),
                        ),
                        md=4,
                    ),
                ],
                className="mb-3",
            ),
        ]
    )


def register_plasmids_callbacks(app, data):
    """Register callbacks for the plasmids tab.

    Parameters
    ----------
    app : dash.Dash
    data : dashboard.data_loader.DashboardData
    """
    plasmids = data.plasmids

    if plasmids.empty:
        return

    @app.callback(
        [
            Output("plasmid-count-hist", "figure"),
            Output("plasmid-mobility-pie", "figure"),
            Output("plasmid-replicon-bar", "figure"),
        ],
        Input("filtered-sample-ids", "data"),
    )
    def _update_plasmids(filtered_ids_json):
        dark_template = "plotly_dark"

        # Parse filtered IDs
        try:
            ids = json.loads(filtered_ids_json) if isinstance(filtered_ids_json, str) else filtered_ids_json
        except (json.JSONDecodeError, TypeError):
            ids = None

        df = plasmids.copy()
        if ids is not None and "sample_id" in df.columns:
            id_set = set(str(i) for i in ids)
            df = df[df["sample_id"].astype(str).isin(id_set)]

        # --- Chart 1: Plasmid count distribution ---
        if "plasmid_count" in df.columns and not df.empty:
            count_fig = px.histogram(
                df,
                x="plasmid_count",
                nbins=max(int(df["plasmid_count"].max()) + 1, 5) if len(df) > 0 else 10,
                title="Plasmid Count per Sample",
                labels={"plasmid_count": "Number of Plasmids", "count": "Samples"},
                template=dark_template,
            )
        else:
            count_fig = go.Figure()
            count_fig.update_layout(template=dark_template, title="Plasmid Count per Sample")
            count_fig.add_annotation(
                text="No plasmid count data",
                xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="#888"),
            )

        # --- Chart 2: Mobility breakdown pie ---
        mob_fig = go.Figure()
        mob_fig.update_layout(template=dark_template, title="Plasmid Mobility")

        if "plasmid_mobility" in df.columns and not df.empty:
            mobility_counts, replicon_counts = _parse_mobility_entries(df["plasmid_mobility"])
            if mobility_counts:
                mob_df = pd.DataFrame(
                    [{"mobility": k, "count": v} for k, v in mobility_counts.most_common()]
                )
                mob_fig = px.pie(
                    mob_df,
                    names="mobility",
                    values="count",
                    title="Plasmid Mobility Breakdown",
                    template=dark_template,
                    color_discrete_map={
                        "mobilizable": "#f39c12",
                        "non-mobilizable": "#3498db",
                        "conjugative": "#e74c3c",
                    },
                )

        # --- Chart 3: Top replicon types ---
        rep_fig = go.Figure()
        rep_fig.update_layout(
            template=dark_template,
            title="Top Replicon Types",
            xaxis_title="Replicon",
            yaxis_title="Count",
        )

        if "plasmid_mobility" in df.columns and not df.empty:
            _, replicon_counts = _parse_mobility_entries(df["plasmid_mobility"])
            if replicon_counts:
                top_reps = replicon_counts.most_common(15)
                rep_fig.add_trace(go.Bar(
                    x=[r[0] for r in top_reps],
                    y=[r[1] for r in top_reps],
                    marker_color="#3498db",
                ))

        return count_fig, mob_fig, rep_fig
