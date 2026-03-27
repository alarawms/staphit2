"""Clusters analysis tab for the Staphit surveillance dashboard."""

import json
import os

import dash_bootstrap_components as dbc
import pandas as pd
from dash import Input, Output, State, dash_table, dcc, html, no_update


def _build_cluster_records(clusters):
    """Convert cluster report list-of-dicts into table records.

    Parameters
    ----------
    clusters : list[dict]
        Each dict may have keys like cluster_id, tier, members, shared_st,
        locations, date_range, etc.

    Returns
    -------
    list[dict]
        Flattened records for the DataTable.
    """
    records = []
    for c in clusters:
        members = c.get("members", [])
        locations = c.get("locations", [])
        if isinstance(locations, list):
            locations = ", ".join(str(loc) for loc in locations)
        record = {
            "cluster_id": str(c.get("cluster_id", c.get("id", ""))),
            "tier": str(c.get("tier", "")),
            "size": len(members) if isinstance(members, list) else c.get("size", 0),
            "shared_st": str(c.get("shared_st", c.get("st", ""))),
            "locations": str(locations),
            "date_range": str(c.get("date_range", "")),
        }
        records.append(record)
    return records


def _load_transmission_pairs(data):
    """Attempt to load transmission_pairs.tsv from the results directory.

    Looks in data attributes or common relative paths.

    Returns
    -------
    pd.DataFrame or None
    """
    # Try to find the file relative to cluster data paths
    try:
        if hasattr(data, "_results_dir"):
            path = os.path.join(data._results_dir, "clusters", "transmission_pairs.tsv")
            if os.path.isfile(path):
                return pd.read_csv(path, sep="\t")
    except Exception:
        pass
    return None


CLUSTER_TABLE_COLUMNS = [
    {"name": "Cluster ID", "id": "cluster_id"},
    {"name": "Tier", "id": "tier"},
    {"name": "Size", "id": "size", "type": "numeric"},
    {"name": "Shared ST", "id": "shared_st"},
    {"name": "Locations", "id": "locations"},
    {"name": "Date Range", "id": "date_range"},
]

_DARK_TABLE_STYLE_HEADER = {
    "backgroundColor": "#303030",
    "color": "#fff",
    "fontWeight": "bold",
    "border": "1px solid #444",
}

_DARK_TABLE_STYLE_CELL = {
    "backgroundColor": "#222",
    "color": "#ddd",
    "border": "1px solid #444",
    "padding": "6px 10px",
    "fontSize": "0.85rem",
    "whiteSpace": "nowrap",
    "overflow": "hidden",
    "textOverflow": "ellipsis",
    "maxWidth": "220px",
}


def create_clusters_layout(data):
    """Build the clusters tab layout.

    Parameters
    ----------
    data : dashboard.data_loader.DashboardData

    Returns
    -------
    html.Div
    """
    clusters = data.clusters if data.clusters else []
    records = _build_cluster_records(clusters)

    # Tier options for dropdown
    tiers = sorted(set(r["tier"] for r in records if r["tier"]))
    tier_options = [{"label": t, "value": t} for t in tiers]

    # Check for transmission pairs
    trans_df = _load_transmission_pairs(data)
    has_trans = trans_df is not None and not trans_df.empty

    if not records:
        cluster_content = html.Div(
            dbc.Alert(
                "Run the pipeline to generate cluster data.",
                color="info",
                className="mt-3",
            )
        )
    else:
        cluster_content = html.Div(
            [
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                html.Label("Filter by tier:", className="small text-muted"),
                                dcc.Dropdown(
                                    id="cluster-tier-dropdown",
                                    options=tier_options,
                                    placeholder="All tiers",
                                    className="dash-bootstrap",
                                ),
                            ],
                            md=3,
                        ),
                    ],
                    className="mb-2",
                ),
                dash_table.DataTable(
                    id="cluster-table",
                    columns=CLUSTER_TABLE_COLUMNS,
                    data=records,
                    page_size=15,
                    sort_action="native",
                    row_selectable="single",
                    style_table={"overflowX": "auto"},
                    style_header=_DARK_TABLE_STYLE_HEADER,
                    style_cell=_DARK_TABLE_STYLE_CELL,
                    style_data_conditional=[
                        {
                            "if": {"state": "active"},
                            "backgroundColor": "#375a7f",
                            "border": "1px solid #4a90d9",
                        },
                    ],
                ),
            ]
        )

    # Transmission pairs section
    if has_trans:
        trans_cols = [{"name": c, "id": c} for c in trans_df.columns]
        trans_section = html.Div(
            [
                html.H6("Transmission Pairs", className="mt-3 mb-2"),
                dash_table.DataTable(
                    id="transmission-table",
                    columns=trans_cols,
                    data=trans_df.to_dict("records"),
                    page_size=10,
                    sort_action="native",
                    style_table={"overflowX": "auto"},
                    style_header=_DARK_TABLE_STYLE_HEADER,
                    style_cell=_DARK_TABLE_STYLE_CELL,
                ),
            ]
        )
    else:
        trans_section = html.Div(id="transmission-table-container")

    # Cluster detail section
    detail_section = html.Div(id="cluster-detail", className="mt-3")

    return html.Div(
        [
            html.H5("Genomic Clusters", className="mt-3 mb-3"),
            cluster_content,
            trans_section,
            detail_section,
        ]
    )


def register_clusters_callbacks(app, data):
    """Register callbacks for the clusters tab.

    Parameters
    ----------
    app : dash.Dash
    data : dashboard.data_loader.DashboardData
    """
    clusters = data.clusters if data.clusters else []
    all_records = _build_cluster_records(clusters)
    merged = data.merged

    if not clusters:
        # No clusters -- register placeholder callbacks to avoid errors
        # with suppress_callback_exceptions=True, these are optional,
        # but we still register them if the tier dropdown exists.
        return

    @app.callback(
        Output("cluster-table", "data"),
        Input("cluster-tier-dropdown", "value"),
    )
    def _filter_by_tier(tier):
        if not tier:
            return all_records
        return [r for r in all_records if r["tier"] == tier]

    @app.callback(
        Output("cluster-detail", "children"),
        Input("cluster-table", "selected_rows"),
        State("cluster-table", "data"),
    )
    def _show_cluster_detail(selected_rows, table_data):
        if not selected_rows or not table_data:
            return html.Div()

        row = table_data[selected_rows[0]]
        cluster_id = row.get("cluster_id", "")

        # Find the full cluster object
        cluster_obj = None
        for c in clusters:
            cid = str(c.get("cluster_id", c.get("id", "")))
            if cid == cluster_id:
                cluster_obj = c
                break

        if not cluster_obj:
            return html.Div("Cluster not found.", className="text-muted")

        members = cluster_obj.get("members", [])
        member_ids = [str(m) if isinstance(m, str) else str(m.get("sample_id", m)) for m in members]

        # Build member details from merged
        member_details = []
        if not merged.empty and "sample_id" in merged.columns:
            member_df = merged[merged["sample_id"].astype(str).isin(member_ids)]
            for _, mrow in member_df.iterrows():
                member_details.append({
                    "sample_id": str(mrow.get("sample_id", "")),
                    "ST": str(mrow.get("mlst_st", "")),
                    "region": str(mrow.get("geo_loc_region", "")),
                    "date": str(mrow.get("collection_date", "")),
                    "MRSA": str(mrow.get("mrsa_status", "")),
                })
        else:
            member_details = [{"sample_id": mid} for mid in member_ids]

        # Shared resistance summary
        resistance_info = ""
        if not merged.empty and "sample_id" in merged.columns:
            member_df = merged[merged["sample_id"].astype(str).isin(member_ids)]
            if "mrsa_status" in member_df.columns:
                mrsa_counts = member_df["mrsa_status"].value_counts().to_dict()
                resistance_info = ", ".join(f"{k}: {v}" for k, v in mrsa_counts.items())

        detail_cols = [{"name": k, "id": k} for k in (member_details[0].keys() if member_details else [])]

        return dbc.Card(
            dbc.CardBody(
                [
                    html.H6(f"Cluster {cluster_id} Details"),
                    html.P(f"Members: {len(member_ids)}", className="small text-muted"),
                    html.P(f"Resistance: {resistance_info}", className="small text-muted") if resistance_info else html.Div(),
                    dash_table.DataTable(
                        columns=detail_cols,
                        data=member_details,
                        page_size=10,
                        style_table={"overflowX": "auto"},
                        style_header=_DARK_TABLE_STYLE_HEADER,
                        style_cell=_DARK_TABLE_STYLE_CELL,
                    ),
                ]
            ),
            className="mt-2",
        )
