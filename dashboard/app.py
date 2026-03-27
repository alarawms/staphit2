"""Dash application factory for the Staphit surveillance dashboard."""

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, dcc, html

from dashboard.components.filters import create_filter_bar, register_filter_callbacks
from dashboard.components.map_viewer import register_map_callbacks
from dashboard.components.tree_viewer import register_tree_callbacks
from dashboard.components.sample_detail import (
    create_detail_drawer,
    register_detail_callbacks,
)
from dashboard.investigation import (
    create_investigation_layout,
    register_investigation_callbacks,
)
from dashboard.analysis.resistance import (
    create_resistance_layout,
    register_resistance_callbacks,
)
from dashboard.analysis.trends import create_trends_layout, register_trends_callbacks
from dashboard.analysis.virulence import (
    create_virulence_layout,
    register_virulence_callbacks,
)
from dashboard.analysis.clusters import (
    create_clusters_layout,
    register_clusters_callbacks,
)
from dashboard.analysis.plasmids import (
    create_plasmids_layout,
    register_plasmids_callbacks,
)
from dashboard.analysis.qc import create_qc_layout, register_qc_callbacks


def create_app(data):
    """Create and return a configured Dash application.

    Parameters
    ----------
    data : dashboard.data_loader.DashboardData
        Loaded pipeline data.

    Returns
    -------
    dash.Dash
    """
    app = dash.Dash(
        __name__,
        external_stylesheets=[dbc.themes.DARKLY],
        suppress_callback_exceptions=True,
        title="Staphit Surveillance",
    )

    sample_ids = (
        data.merged["sample_id"].tolist()
        if not data.merged.empty and "sample_id" in data.merged.columns
        else []
    )

    app.layout = dbc.Container(
        [
            # Header
            dbc.Navbar(
                dbc.Container(
                    [
                        dbc.NavbarBrand(
                            "Staphit Surveillance Dashboard", className="ms-2"
                        ),
                        dbc.Nav(
                            [
                                dbc.NavItem(
                                    dbc.RadioItems(
                                        id="mode-toggle",
                                        options=[
                                            {
                                                "label": "Investigation",
                                                "value": "investigation",
                                            },
                                            {
                                                "label": "Analysis",
                                                "value": "analysis",
                                            },
                                        ],
                                        value="investigation",
                                        inline=True,
                                        className="ms-4 text-light",
                                    )
                                ),
                                dbc.NavItem(
                                    dbc.Badge(
                                        f"{len(data.merged)} samples",
                                        color="info",
                                        className="ms-2",
                                    )
                                ),
                            ]
                        ),
                    ]
                ),
                dark=True,
                color="dark",
            ),
            # Store for filtered sample IDs (JSON list)
            dcc.Store(id="filtered-sample-ids", data=sample_ids),
            # Filter bar
            create_filter_bar(data),
            # Investigation container (tree + map + table)
            html.Div(
                id="investigation-container",
                children=create_investigation_layout(data),
            ),
            # Analysis container (tabbed views)
            html.Div(
                id="analysis-container",
                children=dcc.Tabs(
                    [
                        dcc.Tab(
                            label="Resistance",
                            children=create_resistance_layout(data),
                            className="custom-tab",
                            selected_className="custom-tab--selected",
                        ),
                        dcc.Tab(
                            label="Trends",
                            children=create_trends_layout(data),
                            className="custom-tab",
                            selected_className="custom-tab--selected",
                        ),
                        dcc.Tab(
                            label="Virulence",
                            children=create_virulence_layout(data),
                            className="custom-tab",
                            selected_className="custom-tab--selected",
                        ),
                        dcc.Tab(
                            label="Clusters",
                            children=create_clusters_layout(data),
                            className="custom-tab",
                            selected_className="custom-tab--selected",
                        ),
                        dcc.Tab(
                            label="Plasmids",
                            children=create_plasmids_layout(data),
                            className="custom-tab",
                            selected_className="custom-tab--selected",
                        ),
                        dcc.Tab(
                            label="QC",
                            children=create_qc_layout(data),
                            className="custom-tab",
                            selected_className="custom-tab--selected",
                        ),
                    ],
                    id="analysis-tabs",
                    className="custom-tabs mt-3",
                ),
                style={"display": "none"},
            ),
            # Detail drawer (hidden until row click)
            create_detail_drawer(),
            # Footer
            dbc.Container(
                html.P(
                    "Staphit Surveillance Dashboard | Powered by Plotly Dash",
                    className="text-muted text-center mt-4 mb-2",
                )
            ),
        ],
        fluid=True,
    )

    # --- Mode toggle callback ---
    @app.callback(
        Output("investigation-container", "style"),
        Output("analysis-container", "style"),
        Input("mode-toggle", "value"),
    )
    def toggle_mode(mode):
        if mode == "investigation":
            return {"display": "block"}, {"display": "none"}
        return {"display": "none"}, {"display": "block"}

    # Register callbacks
    register_filter_callbacks(app, data)
    register_investigation_callbacks(app, data)
    register_map_callbacks(app, data)
    register_tree_callbacks(app, data)
    register_detail_callbacks(app, data)
    register_resistance_callbacks(app, data)
    register_trends_callbacks(app, data)
    register_virulence_callbacks(app, data)
    register_clusters_callbacks(app, data)
    register_plasmids_callbacks(app, data)
    register_qc_callbacks(app, data)

    return app
