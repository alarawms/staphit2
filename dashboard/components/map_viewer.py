"""Geographic map component for the Staphit surveillance dashboard.

Displays sample locations on a Scattermapbox map with markers sized by
sample count and colored by dominant ST.  Clicking a marker filters the
investigation table to that region.
"""

import json
import math

import plotly.graph_objects as go
from dash import Input, Output, State, callback_context, dcc, no_update

from dashboard.utils import ST_COLORS

# Hardcoded coordinates for known Saudi regions
REGION_COORDS = {
    "Riyadh": (24.7136, 46.6753),
    "Jeddah": (21.4858, 39.1925),
    "Hail": (27.5114, 41.7208),
    "Makkah": (21.3891, 39.8579),
    "Madinah": (24.5247, 39.5692),
}

# Fallback color when ST is not in ST_COLORS
_DEFAULT_MARKER_COLOR = "#888888"


def _parse_lat_lon(value):
    """Try to parse a 'lat lon' string into (lat, lon) floats.

    Expected format: ``"24.774265 46.738586"``

    Returns
    -------
    tuple[float, float] or None
    """
    if not isinstance(value, str):
        return None
    parts = value.strip().split()
    if len(parts) != 2:
        return None
    try:
        return float(parts[0]), float(parts[1])
    except (ValueError, TypeError):
        return None


def _aggregate_regions(merged, filtered_ids=None):
    """Aggregate sample data by region.

    Parameters
    ----------
    merged : pd.DataFrame
    filtered_ids : list[str] or None
        If provided, only these sample IDs are included.

    Returns
    -------
    list[dict]
        Each dict has keys: region, lat, lon, count, dominant_st, st_counts,
        mrsa_rate, pvl_rate.
    """
    if merged.empty:
        return []

    df = merged.copy()

    if filtered_ids is not None:
        id_set = set(str(i) for i in filtered_ids)
        if "sample_id" in df.columns:
            df = df[df["sample_id"].astype(str).isin(id_set)]

    if df.empty or "geo_loc_region" not in df.columns:
        return []

    results = []
    for region, group in df.groupby("geo_loc_region"):
        region_str = str(region).strip()
        if not region_str or region_str == "nan":
            continue

        # Resolve coordinates
        lat, lon = None, None
        if region_str in REGION_COORDS:
            lat, lon = REGION_COORDS[region_str]
        elif "lat_lon" in group.columns:
            # Use first valid lat_lon from the group
            for val in group["lat_lon"].dropna():
                parsed = _parse_lat_lon(str(val))
                if parsed:
                    lat, lon = parsed
                    break

        if lat is None or lon is None:
            continue

        count = len(group)

        # ST distribution
        st_counts = {}
        if "mlst_st" in group.columns:
            st_counts = (
                group["mlst_st"]
                .dropna()
                .astype(str)
                .value_counts()
                .to_dict()
            )
        dominant_st = max(st_counts, key=st_counts.get) if st_counts else "Unknown"

        # MRSA rate
        mrsa_rate = 0.0
        if "mrsa_status" in group.columns:
            mrsa_count = (group["mrsa_status"].str.upper() == "MRSA").sum()
            mrsa_rate = (mrsa_count / count) * 100 if count > 0 else 0.0

        # PVL rate
        pvl_rate = 0.0
        if "pvl_positive" in group.columns:
            pvl_count = group["pvl_positive"].sum()
            pvl_rate = (pvl_count / count) * 100 if count > 0 else 0.0

        results.append(
            {
                "region": region_str,
                "lat": lat,
                "lon": lon,
                "count": count,
                "dominant_st": dominant_st,
                "st_counts": st_counts,
                "mrsa_rate": mrsa_rate,
                "pvl_rate": pvl_rate,
            }
        )

    return results


def create_map_figure(data, filtered_ids=None):
    """Build a Scattermapbox figure showing sample locations.

    Parameters
    ----------
    data : dashboard.data_loader.DashboardData
    filtered_ids : list[str] or None

    Returns
    -------
    plotly.graph_objects.Figure
    """
    regions = _aggregate_regions(data.merged, filtered_ids)

    fig = go.Figure()

    if not regions:
        # Empty map centered on Saudi Arabia
        fig.update_layout(
            mapbox=dict(
                style="carto-darkmatter",
                center=dict(lat=24.7, lon=44.0),
                zoom=5,
            ),
            margin=dict(l=0, r=0, t=0, b=0),
            paper_bgcolor="#222",
            height=450,
        )
        return fig

    lats = []
    lons = []
    sizes = []
    colors = []
    hover_texts = []
    custom_data = []

    for r in regions:
        lats.append(r["lat"])
        lons.append(r["lon"])
        sizes.append(math.log(r["count"] + 1) * 15)

        st_key = f"ST{r['dominant_st']}" if not r["dominant_st"].startswith("ST") else r["dominant_st"]
        colors.append(ST_COLORS.get(st_key, _DEFAULT_MARKER_COLOR))

        # Build hover text with top 3 STs
        top3 = sorted(r["st_counts"].items(), key=lambda x: x[1], reverse=True)[:3]
        st_lines = ", ".join(f"ST{k}: {v}" if not str(k).startswith("ST") else f"{k}: {v}" for k, v in top3)

        hover = (
            f"<b>{r['region']}</b><br>"
            f"N = {r['count']}<br>"
            f"Top STs: {st_lines}<br>"
            f"MRSA: {r['mrsa_rate']:.0f}%<br>"
            f"PVL+: {r['pvl_rate']:.0f}%"
        )
        hover_texts.append(hover)
        custom_data.append(r["region"])

    fig.add_trace(
        go.Scattermapbox(
            lat=lats,
            lon=lons,
            mode="markers",
            marker=dict(
                size=sizes,
                color=colors,
                opacity=0.85,
            ),
            hovertext=hover_texts,
            hoverinfo="text",
            customdata=custom_data,
        )
    )

    fig.update_layout(
        mapbox=dict(
            style="carto-darkmatter",
            center=dict(lat=24.7, lon=44.0),
            zoom=5,
        ),
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor="#222",
        height=450,
        showlegend=False,
    )

    return fig


def create_map_layout():
    """Return the map Graph component for embedding in the layout.

    Returns
    -------
    dcc.Graph
    """
    return dcc.Graph(
        id="geo-map",
        config={"scrollZoom": True},
        style={"height": "450px"},
    )


def register_map_callbacks(app, data):
    """Register callbacks that link the map to the filtered-sample-ids store.

    Callbacks
    ---------
    1. filtered-sample-ids -> update geo-map figure
    2. clickData on geo-map -> update filtered-sample-ids (filter to region)

    Parameters
    ----------
    app : dash.Dash
    data : dashboard.data_loader.DashboardData
    """

    @app.callback(
        Output("geo-map", "figure"),
        Input("filtered-sample-ids", "data"),
    )
    def _update_map(filtered_ids_json):
        if not filtered_ids_json:
            return create_map_figure(data)
        try:
            ids = (
                json.loads(filtered_ids_json)
                if isinstance(filtered_ids_json, str)
                else filtered_ids_json
            )
        except (json.JSONDecodeError, TypeError):
            ids = None
        return create_map_figure(data, filtered_ids=ids)

    @app.callback(
        Output("filter-location", "value"),
        Input("geo-map", "clickData"),
        State("filter-location", "value"),
        prevent_initial_call=True,
    )
    def _map_click_to_filter(click_data, current_locations):
        """When a map marker is clicked, add its region to the location filter."""
        if not click_data:
            return no_update

        try:
            point = click_data["points"][0]
            region = point.get("customdata")
            if not region:
                return no_update
        except (KeyError, IndexError, TypeError):
            return no_update

        # Toggle: if region already selected, remove it; otherwise add it
        if current_locations and region in current_locations:
            updated = [loc for loc in current_locations if loc != region]
            return updated if updated else None
        else:
            existing = current_locations or []
            return existing + [region]
