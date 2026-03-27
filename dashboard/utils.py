"""Utility functions and constants for the Staphit surveillance dashboard."""

from collections import Counter

import plotly.express as px

# ---------------------------------------------------------------------------
# Color palettes
# ---------------------------------------------------------------------------

# Map common STs to Plotly qualitative colors (first 20)
_QUALITATIVE = px.colors.qualitative.Plotly + px.colors.qualitative.D3
_COMMON_STS = [
    "ST5", "ST8", "ST22", "ST30", "ST36", "ST45", "ST59", "ST72",
    "ST80", "ST88", "ST93", "ST97", "ST105", "ST121", "ST149",
    "ST188", "ST239", "ST240", "ST398", "ST1",
]
ST_COLORS = {st: _QUALITATIVE[i % len(_QUALITATIVE)] for i, st in enumerate(_COMMON_STS)}

SIR_COLORS = {"R": "#ff4757", "I": "#ffa502", "S": "#2ed573"}

# Plotly figure template matching the CSS theme
PLOTLY_TEMPLATE = {
    'layout': {
        'paper_bgcolor': '#141920',
        'plot_bgcolor': '#141920',
        'font': {'family': 'JetBrains Mono, monospace', 'color': '#8b95a5', 'size': 11},
        'title': {'font': {'family': 'Outfit, sans-serif', 'color': '#e8ecf0', 'size': 14}},
        'xaxis': {
            'gridcolor': 'rgba(255,255,255,0.04)',
            'linecolor': 'rgba(255,255,255,0.08)',
            'zerolinecolor': 'rgba(255,255,255,0.08)',
        },
        'yaxis': {
            'gridcolor': 'rgba(255,255,255,0.04)',
            'linecolor': 'rgba(255,255,255,0.08)',
            'zerolinecolor': 'rgba(255,255,255,0.08)',
        },
        'colorway': ['#00d2d3', '#f0a500', '#ff4757', '#2ed573', '#a55eea',
                      '#45aaf2', '#fd9644', '#fc5c65', '#26de81', '#4b7bec'],
        'margin': {'l': 50, 'r': 20, 't': 40, 'b': 40},
    }
}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def count_frequencies(series):
    """Count items from a pandas Series of semicolon-separated strings.

    Parameters
    ----------
    series : pd.Series
        Each element is a string like ``"gene1;gene2;gene3"`` or NaN.

    Returns
    -------
    collections.Counter
        Mapping of item -> count across all rows.
    """
    counts = Counter()
    for val in series.dropna():
        val = str(val).strip()
        if not val or val == "-":
            continue
        for item in val.split(";"):
            item = item.strip()
            if item:
                counts[item] += 1
    return counts


def parse_virulence_summary(df):
    """Parse the ``virulence_summary`` column into boolean/status columns.

    Expected format: ``"PVL+;TSST-;IEC-B;ica:complete;hlg:complete"``

    Adds columns to *df* (in-place):
    - pvl_positive  (bool)
    - tsst_positive (bool)
    - iec_type      (str, e.g. "B", "E", or None)
    - ica_status    (str, e.g. "complete", "partial", or None)
    - hlg_status    (str, e.g. "complete", "partial", or None)

    Parameters
    ----------
    df : pd.DataFrame
        Must contain a ``virulence_summary`` column.

    Returns
    -------
    pd.DataFrame
        The same DataFrame with new columns added.
    """
    pvl = []
    tsst = []
    iec = []
    ica = []
    hlg = []

    col = df.get("virulence_summary")
    if col is None:
        df["pvl_positive"] = False
        df["tsst_positive"] = False
        df["iec_type"] = None
        df["ica_status"] = None
        df["hlg_status"] = None
        return df

    for val in col:
        _pvl = False
        _tsst = False
        _iec = None
        _ica = None
        _hlg = None

        if not isinstance(val, str) or not val.strip() or val.strip() == "-":
            pvl.append(_pvl)
            tsst.append(_tsst)
            iec.append(_iec)
            ica.append(_ica)
            hlg.append(_hlg)
            continue

        for token in val.split(";"):
            token = token.strip()
            if not token:
                continue

            upper = token.upper()
            if upper.startswith("PVL"):
                _pvl = upper.endswith("+")
            elif upper.startswith("TSST"):
                _tsst = upper.endswith("+")
            elif upper.startswith("IEC"):
                # IEC-B, IEC-E, etc.
                parts = token.split("-", 1)
                if len(parts) == 2:
                    _iec = parts[1].strip()
            elif token.startswith("ica:"):
                _ica = token.split(":", 1)[1].strip()
            elif token.startswith("hlg:"):
                _hlg = token.split(":", 1)[1].strip()

        pvl.append(_pvl)
        tsst.append(_tsst)
        iec.append(_iec)
        ica.append(_ica)
        hlg.append(_hlg)

    df["pvl_positive"] = pvl
    df["tsst_positive"] = tsst
    df["iec_type"] = iec
    df["ica_status"] = ica
    df["hlg_status"] = hlg

    return df
