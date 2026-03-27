"""Phylogenetic tree viewer for the Staphit surveillance dashboard.

Parses Newick strings, computes cladogram layouts, and renders interactive
Plotly scatter-based tree figures with tip coloring by metadata fields.
"""

import json

import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash import Input, Output, dcc, html

from dashboard.utils import ST_COLORS

# Fields available for tip coloring
COLOR_BY_OPTIONS = [
    {"label": "ST", "value": "mlst_st"},
    {"label": "MRSA / MSSA", "value": "mrsa_status"},
    {"label": "PVL", "value": "pvl_positive"},
    {"label": "SCCmec", "value": "sccmec_type"},
    {"label": "agr group", "value": "agr_group"},
    {"label": "spa type", "value": "spa_type"},
]

# Qualitative palette for non-ST fields
_PALETTE = [
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52",
]


# ---------------------------------------------------------------------------
# Newick parser
# ---------------------------------------------------------------------------


def parse_newick(newick_str):
    """Parse a Newick-format string into a nested dict tree.

    Parameters
    ----------
    newick_str : str
        Newick tree string, e.g.
        ``"(A:0.1,(B:0.2,C:0.3):0.4);"``

    Returns
    -------
    dict
        Tree with keys ``name``, ``branch_length``, ``children``.
    """
    s = newick_str.strip().rstrip(";").strip()
    tree, _ = _parse_subtree(s, 0)
    # Ensure root has a name
    if not tree.get("name"):
        tree["name"] = "root"
    return tree


_internal_counter = 0


def _parse_subtree(s, pos):
    """Recursively parse a subtree starting at *pos*.

    Returns (node_dict, next_pos).
    """
    global _internal_counter

    if pos < len(s) and s[pos] == "(":
        # Internal node: read children
        children = []
        pos += 1  # skip '('
        while True:
            child, pos = _parse_subtree(s, pos)
            children.append(child)
            if pos < len(s) and s[pos] == ",":
                pos += 1  # skip ','
            else:
                break
        if pos < len(s) and s[pos] == ")":
            pos += 1  # skip ')'

        # Read optional label (bootstrap) and branch length after ')'
        name, branch_length, pos = _read_label_and_length(s, pos)
        if not name:
            name = f"internal_{_internal_counter}"
            _internal_counter += 1

        return {
            "name": name,
            "branch_length": branch_length,
            "children": children,
        }, pos
    else:
        # Leaf node
        name, branch_length, pos = _read_label_and_length(s, pos)
        return {
            "name": name,
            "branch_length": branch_length,
            "children": [],
        }, pos


def _read_label_and_length(s, pos):
    """Read a node label and optional ``:length`` starting at *pos*.

    Returns (name, branch_length, next_pos).
    """
    name_chars = []
    while pos < len(s) and s[pos] not in (":", ",", "(", ")", ";"):
        name_chars.append(s[pos])
        pos += 1

    name = "".join(name_chars).strip()

    branch_length = 0.0
    if pos < len(s) and s[pos] == ":":
        pos += 1  # skip ':'
        length_chars = []
        while pos < len(s) and s[pos] not in (",", "(", ")", ";"):
            length_chars.append(s[pos])
            pos += 1
        try:
            branch_length = float("".join(length_chars).strip())
        except ValueError:
            branch_length = 0.0

    return name, branch_length, pos


# ---------------------------------------------------------------------------
# Layout computation
# ---------------------------------------------------------------------------


def compute_layout(tree, leaf_count=None):
    """Compute x, y coordinates for every node in a cladogram layout.

    Parameters
    ----------
    tree : dict
        Parsed Newick tree (from :func:`parse_newick`).
    leaf_count : int or None
        Used internally for recursion; do not supply.

    Returns
    -------
    tuple(list[dict], list[dict])
        ``(nodes, branches)`` where each node is
        ``{'name', 'x', 'y', 'is_leaf', 'branch_length'}`` and each branch
        is ``{'x0', 'y0', 'x1', 'y1'}``.
    """
    # Assign leaf y-positions first via an in-order traversal
    leaf_y = [0]  # mutable counter

    def assign_y(node):
        if not node["children"]:
            y = leaf_y[0]
            leaf_y[0] += 1
            node["_y"] = y
        else:
            for child in node["children"]:
                assign_y(child)
            child_ys = [c["_y"] for c in node["children"]]
            node["_y"] = sum(child_ys) / len(child_ys)

    assign_y(tree)

    # Assign x-positions (cumulative branch length from root)
    nodes = []
    branches = []

    def walk(node, x_parent):
        x = x_parent + node["branch_length"]
        nodes.append({
            "name": node["name"],
            "x": x,
            "y": node["_y"],
            "is_leaf": len(node["children"]) == 0,
            "branch_length": node["branch_length"],
        })
        for child in node["children"]:
            child_x = x + child["branch_length"]
            # Horizontal segment (parent x to child x at child y)
            branches.append({
                "x0": x, "y0": child["_y"],
                "x1": child_x, "y1": child["_y"],
            })
            # Vertical segment (parent x, from parent y-range to child y)
            branches.append({
                "x0": x, "y0": node["_y"],
                "x1": x, "y1": child["_y"],
            })
            walk(child, x)

    walk(tree, 0.0)

    return nodes, branches


# ---------------------------------------------------------------------------
# Figure builder
# ---------------------------------------------------------------------------


def create_tree_figure(nodes, branches, tip_metadata, color_by="mlst_st",
                       filtered_ids=None):
    """Build a Plotly figure showing the phylogenetic tree.

    Parameters
    ----------
    nodes : list[dict]
        Output of :func:`compute_layout` (node list).
    branches : list[dict]
        Output of :func:`compute_layout` (branch list).
    tip_metadata : dict
        ``{sample_id: {mlst_st, spa_type, ...}}``.
    color_by : str
        Metadata field to color tips by.
    filtered_ids : set or None
        Sample IDs currently passing the global filter.  Tips not in this
        set are rendered small and grey.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    fig = go.Figure()

    # Draw branches as individual line segments
    for b in branches:
        fig.add_trace(go.Scatter(
            x=[b["x0"], b["x1"]],
            y=[b["y0"], b["y1"]],
            mode="lines",
            line=dict(color="#666", width=1),
            hoverinfo="skip",
            showlegend=False,
        ))

    # Separate leaves by color group
    leaves = [n for n in nodes if n["is_leaf"]]

    if not leaves:
        _style_figure(fig)
        return fig

    # Build a value -> color mapping
    all_values = set()
    for leaf in leaves:
        meta = tip_metadata.get(leaf["name"], {})
        val = str(meta.get(color_by, "unknown") or "unknown")
        all_values.add(val)

    value_color_map = _build_color_map(color_by, sorted(all_values))

    # Group leaves by (color_value, is_filtered)
    groups = {}
    for leaf in leaves:
        meta = tip_metadata.get(leaf["name"], {})
        val = str(meta.get(color_by, "unknown") or "unknown")

        if filtered_ids is not None:
            in_filter = leaf["name"] in filtered_ids
        else:
            in_filter = True

        key = (val, in_filter)
        groups.setdefault(key, {"x": [], "y": [], "text": [], "ids": []})
        groups[key]["x"].append(leaf["x"])
        groups[key]["y"].append(leaf["y"])
        groups[key]["ids"].append(leaf["name"])

        # Hover text
        st = meta.get("mlst_st", "-")
        spa = meta.get("spa_type", "-")
        mrsa = meta.get("mrsa_status", "-")
        groups[key]["text"].append(
            f"{leaf['name']}<br>ST{st} / {spa}<br>{mrsa}"
        )

    # Plot each group
    for (val, in_filter), g in sorted(groups.items()):
        color = value_color_map.get(val, "#888")
        fig.add_trace(go.Scatter(
            x=g["x"],
            y=g["y"],
            mode="markers",
            marker=dict(
                size=9 if in_filter else 5,
                color=color if in_filter else "#555",
                line=dict(width=1, color="#fff") if in_filter else dict(width=0),
                opacity=1.0 if in_filter else 0.35,
            ),
            text=g["text"],
            hoverinfo="text",
            name=val if in_filter else f"{val} (filtered)",
            legendgroup=val,
            showlegend=in_filter,
            customdata=g["ids"],
        ))

    _style_figure(fig)
    return fig


def _build_color_map(color_by, values):
    """Map categorical values to colors."""
    cmap = {}
    if color_by == "mlst_st":
        for v in values:
            key = f"ST{v}" if not v.startswith("ST") else v
            cmap[v] = ST_COLORS.get(key, "#888")
        # Assign palette colors to any unmapped STs
        idx = 0
        for v in values:
            if cmap.get(v) == "#888":
                cmap[v] = _PALETTE[idx % len(_PALETTE)]
                idx += 1
    else:
        for i, v in enumerate(values):
            cmap[v] = _PALETTE[i % len(_PALETTE)]
    return cmap


def _style_figure(fig):
    """Apply dark-theme styling to the tree figure."""
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=5, r=5, t=5, b=5),
        xaxis=dict(
            showgrid=False, zeroline=False, showticklabels=False,
            title=None,
        ),
        yaxis=dict(
            showgrid=False, zeroline=False, showticklabels=False,
            title=None, autorange="reversed",
        ),
        legend=dict(
            orientation="h", yanchor="top", y=-0.02,
            xanchor="left", x=0, font=dict(size=10),
        ),
        height=max(300, 25 * sum(1 for _ in fig.data if getattr(_, "mode", "") == "markers")),
    )


# ---------------------------------------------------------------------------
# Dash layout
# ---------------------------------------------------------------------------


def create_tree_layout():
    """Return a Div containing the phylogenetic tree graph and color dropdown.

    Returns
    -------
    html.Div
    """
    return html.Div([
        html.H6("Phylogenetic Tree", className="mb-1 mt-1"),
        dbc.Row([
            dbc.Col(
                dcc.Dropdown(
                    id="tree-color-by",
                    options=COLOR_BY_OPTIONS,
                    value="mlst_st",
                    clearable=False,
                    className="dash-bootstrap",
                    style={"fontSize": "0.8rem"},
                ),
                width=6,
            ),
        ], className="mb-1"),
        html.Div(
            id="phylo-tree-container",
            children=[
                dcc.Graph(
                    id="phylo-tree",
                    config={"displayModeBar": False},
                    style={"height": "65vh"},
                ),
            ],
        ),
    ])


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------


def register_tree_callbacks(app, data):
    """Register callbacks that update the phylogenetic tree.

    Parameters
    ----------
    app : dash.Dash
    data : dashboard.data_loader.DashboardData
    """
    # Pre-parse tree and compute layout once
    tree_nodes = []
    tree_branches = []
    tip_meta = {}

    if data.tree_newick:
        global _internal_counter
        _internal_counter = 0
        try:
            tree = parse_newick(data.tree_newick)
            tree_nodes, tree_branches = compute_layout(tree)
        except Exception as exc:
            import sys
            print(f"Warning: failed to parse Newick tree: {exc}", file=sys.stderr)

    # Build tip metadata lookup from merged DataFrame
    if not data.merged.empty and "sample_id" in data.merged.columns:
        for _, row in data.merged.iterrows():
            sid = str(row["sample_id"])
            tip_meta[sid] = {
                "mlst_st": row.get("mlst_st", ""),
                "spa_type": row.get("spa_type", ""),
                "mrsa_status": row.get("mrsa_status", ""),
                "sccmec_type": row.get("sccmec_type", ""),
                "agr_group": row.get("agr_group", ""),
                "pvl_positive": "PVL+" if row.get("pvl_positive") is True else "PVL-",
            }

    @app.callback(
        Output("phylo-tree", "figure"),
        [
            Input("tree-color-by", "value"),
            Input("filtered-sample-ids", "data"),
        ],
    )
    def _update_tree(color_by, filtered_ids_json):
        if not tree_nodes:
            # Return an empty figure with a message
            fig = go.Figure()
            fig.add_annotation(
                text="No tree available",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16, color="#888"),
            )
            fig.update_layout(
                template="plotly_dark",
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                height=300,
            )
            return fig

        # Parse filtered IDs
        filtered_ids = None
        if filtered_ids_json:
            try:
                ids = (
                    json.loads(filtered_ids_json)
                    if isinstance(filtered_ids_json, str)
                    else filtered_ids_json
                )
                filtered_ids = set(str(i) for i in ids)
            except (json.JSONDecodeError, TypeError):
                pass

        color_field = color_by or "mlst_st"

        fig = create_tree_figure(
            tree_nodes, tree_branches, tip_meta,
            color_by=color_field,
            filtered_ids=filtered_ids,
        )

        # Dynamic height based on number of leaves
        n_leaves = sum(1 for n in tree_nodes if n["is_leaf"])
        fig.update_layout(height=max(300, n_leaves * 28))

        return fig
