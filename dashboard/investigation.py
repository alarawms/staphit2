"""Sample investigation table for the Staphit surveillance dashboard."""

import json

import dash_bootstrap_components as dbc
from dash import Input, Output, dash_table, html

from dashboard.components.map_viewer import create_map_layout
from dashboard.components.tree_viewer import create_tree_layout

# Columns shown in the sample table
TABLE_COLUMNS = [
    {"name": "Sample", "id": "sample_id"},
    {"name": "ST", "id": "mlst_st"},
    {"name": "spa", "id": "spa_type"},
    {"name": "SCCmec", "id": "sccmec_type"},
    {"name": "agr", "id": "agr_group"},
    {"name": "MRSA", "id": "mrsa_status"},
    {"name": "Date", "id": "collection_date"},
    {"name": "Region", "id": "geo_loc_region"},
    {"name": "PVL", "id": "pvl_positive"},
]


def _safe_col(df, col):
    """Return the column values as list of strings, or empty strings when missing."""
    if col in df.columns:
        return df[col].fillna("").astype(str).tolist()
    return [""] * len(df)


def _build_records(merged):
    """Build the table records from the merged DataFrame."""
    if merged.empty:
        return []
    records = []
    for idx in range(len(merged)):
        row = {}
        for c in TABLE_COLUMNS:
            col = c["id"]
            if col in merged.columns:
                val = merged.iloc[idx][col]
                if col == "pvl_positive":
                    val = "Yes" if val is True or val == "True" else "No"
                row[col] = "" if val is None or (isinstance(val, float) and str(val) == "nan") else str(val)
            else:
                row[col] = ""
        records.append(row)
    return records


def create_investigation_layout(data):
    """Build the investigation layout containing the sample DataTable.

    Parameters
    ----------
    data : dashboard.data_loader.DashboardData

    Returns
    -------
    html.Div
    """
    records = _build_records(data.merged)

    sample_table = dash_table.DataTable(
        id="sample-table",
        columns=TABLE_COLUMNS,
        data=records,
        page_size=25,
        sort_action="native",
        filter_action="native",
        row_selectable=False,
        cell_selectable=True,
        style_table={"overflowX": "auto"},
        style_header={
            "backgroundColor": "#303030",
            "color": "#fff",
            "fontWeight": "bold",
            "border": "1px solid #444",
        },
        style_cell={
            "backgroundColor": "#222",
            "color": "#ddd",
            "border": "1px solid #444",
            "padding": "6px 10px",
            "fontSize": "0.85rem",
            "whiteSpace": "nowrap",
            "overflow": "hidden",
            "textOverflow": "ellipsis",
            "maxWidth": "180px",
        },
        style_filter={
            "backgroundColor": "#2a2a2a",
            "color": "#ddd",
        },
        style_data_conditional=[
            {
                "if": {"state": "active"},
                "backgroundColor": "#375a7f",
                "border": "1px solid #4a90d9",
            },
        ],
    )

    return html.Div(
        [
            html.H5("Sample Investigation", className="mt-2 mb-2"),
            dbc.Row(
                [
                    dbc.Col(create_tree_layout(), width=4),
                    dbc.Col(create_map_layout(), width=4),
                    dbc.Col(sample_table, width=4),
                ],
                className="mb-3",
            ),
        ]
    )


def register_investigation_callbacks(app, data):
    """Register callbacks that update the sample table when filters change.

    Parameters
    ----------
    app : dash.Dash
    data : dashboard.data_loader.DashboardData
    """
    all_records = _build_records(data.merged)

    @app.callback(
        Output("sample-table", "data"),
        Input("filtered-sample-ids", "data"),
    )
    def _update_table(filtered_ids_json):
        if not filtered_ids_json:
            return all_records
        try:
            ids = json.loads(filtered_ids_json) if isinstance(filtered_ids_json, str) else filtered_ids_json
        except (json.JSONDecodeError, TypeError):
            return all_records
        if not ids:
            return []
        id_set = set(str(i) for i in ids)
        return [r for r in all_records if r.get("sample_id") in id_set]
