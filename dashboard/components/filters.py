"""Global filter bar for the Staphit surveillance dashboard."""

import json

import dash_bootstrap_components as dbc
import pandas as pd
from dash import Input, Output, State, callback_context, dcc, html


def _year_range(merged: pd.DataFrame):
    """Return (min_year, max_year) from the collection_date column."""
    if merged.empty or "collection_date" not in merged.columns:
        return 2000, 2026
    years = pd.to_numeric(
        merged["collection_date"].astype(str).str[:4], errors="coerce"
    ).dropna()
    if years.empty:
        return 2000, 2026
    return int(years.min()), int(years.max())


def _top_values(series: pd.Series, n: int = 20):
    """Return the top *n* most frequent non-null values from *series*."""
    if series is None or series.empty:
        return []
    counts = series.dropna().astype(str).value_counts()
    return counts.head(n).index.tolist()


def create_filter_bar(data):
    """Build the global filter bar card.

    Parameters
    ----------
    data : dashboard.data_loader.DashboardData

    Returns
    -------
    dbc.Card
    """
    merged = data.merged

    # Date range
    yr_min, yr_max = _year_range(merged)
    year_marks = {y: str(y) for y in range(yr_min, yr_max + 1)}

    # Location options
    locations = sorted(
        merged["geo_loc_region"].dropna().unique().tolist()
    ) if not merged.empty and "geo_loc_region" in merged.columns else []

    # ST options (top 20 + "Other")
    st_options = []
    if not merged.empty and "mlst_st" in merged.columns:
        top_sts = _top_values(merged["mlst_st"], 20)
        st_options = [{"label": f"ST{v}", "value": v} for v in top_sts]
        st_options.append({"label": "Other", "value": "__other__"})

    # SCCmec options
    sccmec_options = []
    if not merged.empty and "sccmec_type" in merged.columns:
        sccmec_vals = sorted(
            merged["sccmec_type"].dropna().unique().tolist()
        )
        sccmec_options = [{"label": v, "value": v} for v in sccmec_vals if v and str(v).strip() not in ("", "-")]

    total = len(merged)

    return dbc.Card(
        dbc.CardBody(
            [
                dbc.Row(
                    [
                        # Date range slider
                        dbc.Col(
                            [
                                html.Label("Date range", className="small text-muted mb-1"),
                                dcc.RangeSlider(
                                    id="filter-date-range",
                                    min=yr_min,
                                    max=yr_max,
                                    step=1,
                                    value=[yr_min, yr_max],
                                    marks=year_marks,
                                    tooltip={"placement": "bottom", "always_visible": False},
                                ),
                            ],
                            md=3,
                            className="mb-2",
                        ),
                        # Location dropdown
                        dbc.Col(
                            [
                                html.Label("Location", className="small text-muted mb-1"),
                                dcc.Dropdown(
                                    id="filter-location",
                                    options=[{"label": v, "value": v} for v in locations],
                                    multi=True,
                                    placeholder="All regions",
                                    className="dash-bootstrap",
                                ),
                            ],
                            md=2,
                            className="mb-2",
                        ),
                        # ST dropdown
                        dbc.Col(
                            [
                                html.Label("ST", className="small text-muted mb-1"),
                                dcc.Dropdown(
                                    id="filter-st",
                                    options=st_options,
                                    multi=True,
                                    placeholder="All STs",
                                    className="dash-bootstrap",
                                ),
                            ],
                            md=2,
                            className="mb-2",
                        ),
                        # MRSA / MSSA radio
                        dbc.Col(
                            [
                                html.Label("Resistance", className="small text-muted mb-1"),
                                dbc.RadioItems(
                                    id="filter-mrsa",
                                    options=[
                                        {"label": "All", "value": "All"},
                                        {"label": "MRSA", "value": "MRSA"},
                                        {"label": "MSSA", "value": "MSSA"},
                                    ],
                                    value="All",
                                    inline=True,
                                ),
                            ],
                            md=1,
                            className="mb-2",
                        ),
                        # PVL radio
                        dbc.Col(
                            [
                                html.Label("PVL", className="small text-muted mb-1"),
                                dbc.RadioItems(
                                    id="filter-pvl",
                                    options=[
                                        {"label": "All", "value": "All"},
                                        {"label": "PVL+", "value": "PVL+"},
                                        {"label": "PVL-", "value": "PVL-"},
                                    ],
                                    value="All",
                                    inline=True,
                                ),
                            ],
                            md=1,
                            className="mb-2",
                        ),
                        # SCCmec dropdown
                        dbc.Col(
                            [
                                html.Label("SCCmec", className="small text-muted mb-1"),
                                dcc.Dropdown(
                                    id="filter-sccmec",
                                    options=sccmec_options,
                                    multi=True,
                                    placeholder="All types",
                                    className="dash-bootstrap",
                                ),
                            ],
                            md=2,
                            className="mb-2",
                        ),
                        # Sample count badge
                        dbc.Col(
                            [
                                html.Label("\u00a0", className="small mb-1"),
                                html.Div(
                                    dbc.Badge(
                                        f"{total} samples",
                                        id="filter-sample-count",
                                        color="info",
                                        className="fs-6 p-2",
                                    ),
                                    className="text-center",
                                ),
                            ],
                            md=1,
                            className="mb-2 d-flex flex-column align-items-center",
                        ),
                    ],
                    className="align-items-end",
                ),
            ]
        ),
        className="mb-3 mt-3",
    )


def register_filter_callbacks(app, data):
    """Register the filter callback that populates ``filtered-sample-ids`` store.

    Parameters
    ----------
    app : dash.Dash
    data : dashboard.data_loader.DashboardData
    """
    merged = data.merged

    @app.callback(
        [
            Output("filtered-sample-ids", "data"),
            Output("filter-sample-count", "children"),
        ],
        [
            Input("filter-date-range", "value"),
            Input("filter-location", "value"),
            Input("filter-st", "value"),
            Input("filter-mrsa", "value"),
            Input("filter-pvl", "value"),
            Input("filter-sccmec", "value"),
        ],
    )
    def _apply_filters(date_range, locations, sts, mrsa, pvl, sccmec):
        if merged.empty:
            return json.dumps([]), "0 samples"

        mask = pd.Series(True, index=merged.index)

        # Date range
        if date_range and "collection_date" in merged.columns:
            years = pd.to_numeric(
                merged["collection_date"].astype(str).str[:4], errors="coerce"
            )
            mask &= years.between(date_range[0], date_range[1]) | years.isna()

        # Location
        if locations and "geo_loc_region" in merged.columns:
            mask &= merged["geo_loc_region"].isin(locations)

        # ST
        if sts and "mlst_st" in merged.columns:
            if "__other__" in sts:
                top20 = _top_values(merged["mlst_st"], 20)
                explicit = [s for s in sts if s != "__other__"]
                mask &= merged["mlst_st"].isin(explicit) | ~merged["mlst_st"].isin(top20)
            else:
                mask &= merged["mlst_st"].astype(str).isin([str(s) for s in sts])

        # MRSA / MSSA
        if mrsa and mrsa != "All" and "mrsa_status" in merged.columns:
            mask &= merged["mrsa_status"].str.upper() == mrsa.upper()

        # PVL
        if pvl and pvl != "All" and "pvl_positive" in merged.columns:
            if pvl == "PVL+":
                mask &= merged["pvl_positive"] == True  # noqa: E712
            else:
                mask &= merged["pvl_positive"] == False  # noqa: E712

        # SCCmec
        if sccmec and "sccmec_type" in merged.columns:
            mask &= merged["sccmec_type"].isin(sccmec)

        filtered = merged.loc[mask]
        ids = filtered["sample_id"].tolist() if "sample_id" in filtered.columns else []
        return json.dumps(ids), f"{len(ids)} samples"
