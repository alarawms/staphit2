"""Sample detail drawer (Offcanvas) for the Staphit surveillance dashboard."""

import json

import dash_bootstrap_components as dbc
import pandas as pd
from dash import Input, Output, State, html, no_update

from dashboard.utils import SIR_COLORS


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


def create_detail_drawer():
    """Return a hidden dbc.Offcanvas that shows per-sample details.

    Returns
    -------
    dbc.Offcanvas
    """
    return dbc.Offcanvas(
        id="detail-drawer",
        title="Sample Detail",
        is_open=False,
        placement="end",
        style={"width": "520px", "backgroundColor": "#222", "color": "#ddd"},
        children=[html.Div(id="detail-drawer-content")],
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _kv_row(label, value, bold=False):
    """Single key-value row."""
    val_style = {"fontWeight": "bold"} if bold else {}
    return html.Tr(
        [
            html.Td(label, className="text-muted pe-3", style={"whiteSpace": "nowrap"}),
            html.Td(str(value) if value not in (None, "", "nan") else "-", style=val_style),
        ]
    )


def _section_card(title, body, color="secondary"):
    return dbc.Card(
        [
            dbc.CardHeader(title, className="py-1 px-2 fw-bold small"),
            dbc.CardBody(body, className="py-2 px-2"),
        ],
        className="mb-2",
        color=color,
        outline=True,
    )


def _sir_badge(sir):
    """Colored badge for R / I / S."""
    sir = str(sir).strip().upper() if sir else ""
    color = SIR_COLORS.get(sir, "#888")
    return html.Span(
        sir or "-",
        style={
            "backgroundColor": color,
            "color": "#fff",
            "padding": "2px 8px",
            "borderRadius": "4px",
            "fontSize": "0.8rem",
            "fontWeight": "bold",
        },
    )


# ---------------------------------------------------------------------------
# Build drawer contents from report JSON + antibiogram
# ---------------------------------------------------------------------------


def _build_detail(sample_id, report, abg_rows):
    """Return a list of cards for the drawer body.

    Parameters
    ----------
    sample_id : str
    report : dict
        Loaded _report.json for this sample.
    abg_rows : pd.DataFrame
        Antibiogram rows for this sample.
    """
    cards = []

    # 1. Metadata card -------------------------------------------------------
    meta = report.get("metadata", {})
    meta_items = [
        ("Organism", meta.get("organism")),
        ("Date", meta.get("collection_date")),
        ("Region", meta.get("geo_loc_region")),
        ("Country", meta.get("geo_loc_country")),
        ("Host age", meta.get("host_age")),
        ("Host sex", meta.get("host_sex")),
        ("Patient status", meta.get("patient_status")),
        ("Host disease", meta.get("host_disease")),
        ("Body site", meta.get("host_body_site")),
    ]
    cards.append(
        _section_card(
            "Metadata",
            html.Table(
                [_kv_row(k, v) for k, v in meta_items],
                className="table table-sm table-borderless mb-0",
            ),
        )
    )

    # 2. Typing card ---------------------------------------------------------
    typing = report.get("typing", {})
    mlst = typing.get("mlst", {})
    spa = typing.get("spa", {})
    sccmec = typing.get("sccmec", {})
    agr = typing.get("agr", {})
    alleles = ", ".join(mlst.get("alleles", []))
    typing_items = [
        ("ST", mlst.get("st")),
        ("Alleles", alleles or None),
        ("spa type", spa.get("type")),
        ("SCCmec", sccmec.get("type")),
        ("agr group", agr.get("group")),
        ("agr confidence", f"{agr.get('confidence')}%" if agr.get("confidence") else None),
    ]
    cards.append(
        _section_card(
            "Typing",
            html.Table(
                [_kv_row(k, v) for k, v in typing_items],
                className="table table-sm table-borderless mb-0",
            ),
        )
    )

    # 3. Virulence card ------------------------------------------------------
    vp = report.get("virulence_profile", {})
    pvl = vp.get("pvl", {})
    tsst = vp.get("tsst", {})
    iec = vp.get("iec", {})
    entero = vp.get("enterotoxins", {})
    biofilm = vp.get("biofilm", {})
    hemo = vp.get("hemolysins", {})
    vir_items = [
        ("PVL", pvl.get("status", "-")),
        ("TSST", tsst.get("status", "-")),
        ("IEC type", iec.get("type", "-")),
        ("Enterotoxins", f"{entero.get('count', 0)} ({', '.join(entero.get('genes', []))})"),
        ("ica operon", biofilm.get("ica_operon", "-")),
        ("hlg operon", hemo.get("hlg_operon", "-")),
    ]
    cards.append(
        _section_card(
            "Virulence",
            html.Table(
                [_kv_row(k, v) for k, v in vir_items],
                className="table table-sm table-borderless mb-0",
            ),
        )
    )

    # 4. AMR card ------------------------------------------------------------
    amrfinder = report.get("resistance", {}).get("amrfinder", [])
    point_muts = report.get("point_mutations", {}).get("mutations", [])
    if amrfinder or point_muts:
        amr_rows = []
        for hit in amrfinder:
            amr_rows.append(
                html.Tr([
                    html.Td(hit.get("gene") or hit.get("product") or "-", className="small"),
                    html.Td(hit.get("class", "-"), className="small"),
                    html.Td(hit.get("subclass", "-"), className="small"),
                ])
            )
        for mut in point_muts:
            amr_rows.append(
                html.Tr([
                    html.Td(str(mut) if isinstance(mut, str) else mut.get("gene", "-"), className="small text-warning"),
                    html.Td("POINT", className="small"),
                    html.Td(mut.get("class", "-") if isinstance(mut, dict) else "-", className="small"),
                ])
            )
        cards.append(
            _section_card(
                "AMR genes",
                html.Table(
                    [
                        html.Thead(html.Tr([
                            html.Th("Gene", className="small"),
                            html.Th("Class", className="small"),
                            html.Th("Subclass", className="small"),
                        ])),
                        html.Tbody(amr_rows),
                    ],
                    className="table table-sm table-borderless mb-0",
                ),
            )
        )
    else:
        cards.append(_section_card("AMR genes", html.P("No AMR genes detected.", className="text-muted small mb-0")))

    # 5. Antibiogram card ----------------------------------------------------
    if not abg_rows.empty:
        abg_table_rows = []
        for _, row in abg_rows.iterrows():
            abg_table_rows.append(
                html.Tr([
                    html.Td(row.get("antibiotic", "-"), className="small"),
                    html.Td(_sir_badge(row.get("resistance_phenotype", row.get("sir", "")))),
                    html.Td(str(row.get("measurement", row.get("mic", "-"))), className="small"),
                    html.Td(str(row.get("measurement_sign", row.get("mic_sign", ""))), className="small"),
                    html.Td(str(row.get("measurement_units", row.get("units", ""))), className="small"),
                ])
            )
        cards.append(
            _section_card(
                "Antibiogram",
                html.Table(
                    [
                        html.Thead(html.Tr([
                            html.Th("Antibiotic", className="small"),
                            html.Th("SIR", className="small"),
                            html.Th("MIC", className="small"),
                            html.Th("Sign", className="small"),
                            html.Th("Units", className="small"),
                        ])),
                        html.Tbody(abg_table_rows),
                    ],
                    className="table table-sm table-borderless mb-0",
                ),
            )
        )
    else:
        cards.append(_section_card("Antibiogram", html.P("No antibiogram data.", className="text-muted small mb-0")))

    # 6. Plasmid card --------------------------------------------------------
    plasmids = report.get("plasmids", [])
    if plasmids:
        plasm_items = [
            ("Count", len(plasmids)),
        ]
        for i, p in enumerate(plasmids):
            plasm_items.append((f"#{i+1}", p.get("product", p.get("gene", "-"))))
        cards.append(
            _section_card(
                "Plasmids",
                html.Table(
                    [_kv_row(k, v) for k, v in plasm_items],
                    className="table table-sm table-borderless mb-0",
                ),
            )
        )
    else:
        cards.append(_section_card("Plasmids", html.P("No plasmids detected.", className="text-muted small mb-0")))

    # 7. QC card -------------------------------------------------------------
    qc = report.get("qc", {})
    assembly = report.get("assembly", {})
    qc_items = [
        ("Raw reads", qc.get("raw_reads")),
        ("Trimmed reads", qc.get("trimmed_reads")),
        ("Survival rate", f"{qc.get('survival_rate', 0):.1f}%" if qc.get("survival_rate") else None),
        ("Q30 rate", qc.get("q30_rate")),
        ("Assembly length", assembly.get("length")),
        ("Contigs", assembly.get("contigs")),
        ("N50", assembly.get("n50")),
        ("GC%", assembly.get("gc")),
    ]
    cards.append(
        _section_card(
            "QC / Assembly",
            html.Table(
                [_kv_row(k, v) for k, v in qc_items],
                className="table table-sm table-borderless mb-0",
            ),
        )
    )

    return [html.H5(sample_id, className="mb-3")] + cards


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------


def register_detail_callbacks(app, data):
    """Register the callback that opens the detail drawer on row click.

    Parameters
    ----------
    app : dash.Dash
    data : dashboard.data_loader.DashboardData
    """

    @app.callback(
        [
            Output("detail-drawer", "is_open"),
            Output("detail-drawer-content", "children"),
        ],
        Input("sample-table", "active_cell"),
        State("sample-table", "data"),
        prevent_initial_call=True,
    )
    def _open_drawer(active_cell, table_data):
        if not active_cell or not table_data:
            return False, []

        row_idx = active_cell.get("row")
        if row_idx is None or row_idx >= len(table_data):
            return no_update, no_update

        sample_id = table_data[row_idx].get("sample_id")
        if not sample_id:
            return no_update, no_update

        # Load report JSON
        report = {}
        report_path = data.sample_reports.get(sample_id)
        if report_path:
            try:
                with open(report_path) as fh:
                    report = json.load(fh)
            except Exception:
                pass

        # Get antibiogram rows for this sample
        abg_rows = pd.DataFrame()
        if not data.antibiogram.empty and "sample_id" in data.antibiogram.columns:
            abg_rows = data.antibiogram[data.antibiogram["sample_id"] == sample_id]

        content = _build_detail(sample_id, report, abg_rows)
        return True, content
