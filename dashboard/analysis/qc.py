"""QC analysis tab for the Staphit surveillance dashboard."""

import json

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, dash_table, dcc, html

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


def create_qc_layout(data):
    """Build the QC tab layout.

    Parameters
    ----------
    data : dashboard.data_loader.DashboardData

    Returns
    -------
    html.Div
    """
    qc = data.qc
    merged = data.merged
    has_qc = not qc.empty

    children = [html.H5("Quality Control", className="mt-3 mb-3")]

    # Row 1: Assembly scatter + completeness vs contamination
    children.append(
        dbc.Row(
            [
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody(
                            dcc.Graph(id="qc-assembly-scatter", config={"displayModeBar": False})
                        ),
                    ),
                    md=6,
                ),
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody(
                            dcc.Graph(id="qc-completeness-scatter", config={"displayModeBar": False})
                        ),
                    ),
                    md=6,
                ) if has_qc else dbc.Col(md=6),
            ],
            className="mb-3",
        )
    )

    # Row 2: Failed samples table + read survival histogram
    if has_qc or (not merged.empty and "survival_rate" in merged.columns):
        children.append(
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H6("Failed Samples", className="mb-2"),
                                    html.Div(id="qc-failed-table-container"),
                                ]
                            ),
                        ),
                        md=6,
                    ) if has_qc else dbc.Col(md=6),
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                dcc.Graph(id="qc-read-survival", config={"displayModeBar": False})
                            ),
                        ),
                        md=6,
                    ) if not merged.empty and "survival_rate" in merged.columns else dbc.Col(md=6),
                ],
                className="mb-3",
            )
        )

    if not has_qc and (merged.empty or "assembly_length" not in merged.columns):
        children.append(
            dbc.Alert(
                "No QC or assembly data available.",
                color="info",
                className="mt-3",
            )
        )

    return html.Div(children)


def register_qc_callbacks(app, data):
    """Register callbacks for the QC tab.

    Parameters
    ----------
    app : dash.Dash
    data : dashboard.data_loader.DashboardData
    """
    qc = data.qc
    merged = data.merged
    has_qc = not qc.empty

    outputs = [Output("qc-assembly-scatter", "figure")]
    if has_qc:
        outputs.append(Output("qc-completeness-scatter", "figure"))
        outputs.append(Output("qc-failed-table-container", "children"))
    if not merged.empty and "survival_rate" in merged.columns:
        outputs.append(Output("qc-read-survival", "figure"))

    if not outputs:
        return

    @app.callback(
        outputs,
        Input("filtered-sample-ids", "data"),
    )
    def _update_qc(filtered_ids_json):
        dark_template = "plotly_dark"

        # Parse filtered IDs
        try:
            ids = json.loads(filtered_ids_json) if isinstance(filtered_ids_json, str) else filtered_ids_json
        except (json.JSONDecodeError, TypeError):
            ids = None

        # Filter merged data
        df = merged.copy() if not merged.empty else pd.DataFrame()
        if ids is not None and not df.empty and "sample_id" in df.columns:
            id_set = set(str(i) for i in ids)
            df = df[df["sample_id"].astype(str).isin(id_set)]

        # Filter QC data
        qc_df = qc.copy() if has_qc else pd.DataFrame()
        if ids is not None and not qc_df.empty and "sample_id" in qc_df.columns:
            id_set = set(str(i) for i in ids)
            qc_df = qc_df[qc_df["sample_id"].astype(str).isin(id_set)]

        results = []

        # --- Chart 1: Assembly scatter (length vs N50) ---
        asm_fig = go.Figure()
        asm_fig.update_layout(
            template=dark_template,
            title="Assembly Length vs N50",
            xaxis_title="Assembly Length (bp)",
            yaxis_title="N50 (bp)",
        )

        if not df.empty and "assembly_length" in df.columns and "n50" in df.columns:
            # Determine color column
            if has_qc and "sample_id" in df.columns and "sample_id" in qc_df.columns:
                # Merge QC pass/fail into assembly data
                df_asm = df.merge(
                    qc_df[["sample_id", "pass"]].rename(columns={"pass": "qc_pass"}),
                    on="sample_id",
                    how="left",
                )
                df_asm["qc_pass"] = df_asm["qc_pass"].fillna("unknown").astype(str)
                color_col = "qc_pass"
                color_map = {"True": "#2ecc71", "true": "#2ecc71", "False": "#e74c3c", "false": "#e74c3c", "unknown": "#888"}
            elif "mlst_st" in df.columns:
                df_asm = df.copy()
                df_asm["st_label"] = df_asm["mlst_st"].astype(str)
                color_col = "st_label"
                color_map = None
            else:
                df_asm = df.copy()
                color_col = None
                color_map = None

            if color_col and color_map:
                asm_fig = px.scatter(
                    df_asm,
                    x="assembly_length",
                    y="n50",
                    color=color_col,
                    color_discrete_map=color_map,
                    hover_data=["sample_id"] if "sample_id" in df_asm.columns else None,
                    title="Assembly Length vs N50",
                    template=dark_template,
                )
            elif color_col:
                asm_fig = px.scatter(
                    df_asm,
                    x="assembly_length",
                    y="n50",
                    color=color_col,
                    hover_data=["sample_id"] if "sample_id" in df_asm.columns else None,
                    title="Assembly Length vs N50",
                    template=dark_template,
                )
            else:
                asm_fig = px.scatter(
                    df_asm,
                    x="assembly_length",
                    y="n50",
                    hover_data=["sample_id"] if "sample_id" in df_asm.columns else None,
                    title="Assembly Length vs N50",
                    template=dark_template,
                )

        results.append(asm_fig)

        # --- Chart 2: Completeness vs Contamination ---
        if has_qc:
            comp_fig = go.Figure()
            comp_fig.update_layout(
                template=dark_template,
                title="Completeness vs Contamination",
                xaxis_title="Completeness (%)",
                yaxis_title="Contamination (%)",
            )

            if not qc_df.empty and "completeness" in qc_df.columns and "contamination" in qc_df.columns:
                qc_plot = qc_df.copy()
                qc_plot["qc_status"] = qc_plot["pass"].astype(str) if "pass" in qc_plot.columns else "unknown"
                color_map = {"True": "#2ecc71", "true": "#2ecc71", "False": "#e74c3c", "false": "#e74c3c", "unknown": "#888"}

                comp_fig = px.scatter(
                    qc_plot,
                    x="completeness",
                    y="contamination",
                    color="qc_status",
                    color_discrete_map=color_map,
                    hover_data=["sample_id"] if "sample_id" in qc_plot.columns else None,
                    title="Completeness vs Contamination",
                    template=dark_template,
                )

                # Add threshold lines
                comp_fig.add_hline(y=5, line_dash="dash", line_color="#e74c3c",
                                   annotation_text="5% contamination", annotation_position="top right")
                comp_fig.add_vline(x=90, line_dash="dash", line_color="#2ecc71",
                                   annotation_text="90% completeness", annotation_position="top left")

            results.append(comp_fig)

            # --- Failed samples table ---
            if not qc_df.empty and "pass" in qc_df.columns:
                failed = qc_df[qc_df["pass"].astype(str).str.lower() == "false"]
                if not failed.empty:
                    display_cols = ["sample_id", "completeness", "contamination", "reason"]
                    display_cols = [c for c in display_cols if c in failed.columns]
                    table_cols = [{"name": c, "id": c} for c in display_cols]
                    failed_table = dash_table.DataTable(
                        columns=table_cols,
                        data=failed[display_cols].to_dict("records"),
                        page_size=10,
                        sort_action="native",
                        style_table={"overflowX": "auto"},
                        style_header=_DARK_TABLE_STYLE_HEADER,
                        style_cell=_DARK_TABLE_STYLE_CELL,
                    )
                else:
                    failed_table = html.P("All samples passed QC.", className="text-success small")
            else:
                failed_table = html.P("No QC pass/fail data available.", className="text-muted small")

            results.append(failed_table)

        # --- Chart 3: Read survival histogram ---
        if not merged.empty and "survival_rate" in merged.columns:
            surv_df = df if not df.empty and "survival_rate" in df.columns else pd.DataFrame()
            if not surv_df.empty:
                surv_fig = px.histogram(
                    surv_df,
                    x="survival_rate",
                    nbins=20,
                    title="Read Survival Rate Distribution",
                    labels={"survival_rate": "Survival Rate (%)", "count": "Samples"},
                    template=dark_template,
                )
            else:
                surv_fig = go.Figure()
                surv_fig.update_layout(
                    template=dark_template,
                    title="Read Survival Rate Distribution",
                )
                surv_fig.add_annotation(
                    text="No read survival data",
                    xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False,
                    font=dict(size=14, color="#888"),
                )
            results.append(surv_fig)

        return results if len(results) > 1 else results[0]
