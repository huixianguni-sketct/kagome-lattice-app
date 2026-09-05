from __future__ import annotations

import json

import plotly.graph_objects as go

from lattice import KagomeLattice


# Pastel shades chosen to resemble the paper while keeping labels readable.
COLOR_HEX = {
    "Red": "#F3A2A8",
    "Green": "#9FD8BC",
    "Blue": "#A8BCEB",
}


def _build_figure(lattice: KagomeLattice) -> go.Figure:
    edge_x = []
    edge_y = []
    for i, j in lattice.edges:
        x0, y0 = lattice.positions[i]
        x1, y1 = lattice.positions[j]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(width=1.2, color="rgba(70, 70, 70, 0.55)"),
        hoverinfo="skip",
        showlegend=False,
        name="Kagome bonds",
    )

    node_colors = [COLOR_HEX[site.color] for site in lattice.sites]
    customdata = [
        [
            site.site_id,
            site.cell_i,
            site.cell_j,
            site.sublattice,
            site.color,
        ]
        for site in lattice.sites
    ]

    node_trace = go.Scatter(
        x=[site.x for site in lattice.sites],
        y=[site.y for site in lattice.sites],
        mode="markers+text",
        text=["" for _ in lattice.sites],
        textposition="middle center",
        textfont=dict(size=11, color="#111827", family="Arial Black, Arial, sans-serif"),
        customdata=customdata,
        marker=dict(
            size=17,
            color=node_colors,
            line=dict(width=1.2, color="#475569"),
        ),
        hovertemplate=(
            "<b>Qubit %{customdata[0]}</b>"
            "<extra></extra>"
        ),
        showlegend=False,
        name="Qubits",
    )

    # This trace is filled dynamically in JavaScript after a CZ pair is made.
    cz_trace = go.Scatter(
        x=[],
        y=[],
        mode="lines",
        line=dict(width=4, color="#6D28D9"),
        hoverinfo="skip",
        showlegend=False,
        name="CZ pairs",
    )

    # Midpoint labels make it visually obvious which two qubits form a CZ gate.
    cz_label_trace = go.Scatter(
        x=[],
        y=[],
        mode="text",
        text=[],
        textfont=dict(size=11, color="#5B21B6", family="Arial Black, Arial, sans-serif"),
        hoverinfo="skip",
        showlegend=False,
        name="CZ labels",
    )

    fig = go.Figure(data=[edge_trace, node_trace, cz_trace, cz_label_trace])
    fig.update_layout(
        height=720,
        margin=dict(l=5, r=5, t=10, b=5),
        plot_bgcolor="white",
        paper_bgcolor="white",
        dragmode="pan",
        hovermode="closest",
        uirevision="kagome-operator-interface",
        xaxis=dict(
            visible=False,
            showgrid=False,
            zeroline=False,
            fixedrange=False,
        ),
        yaxis=dict(
            visible=False,
            showgrid=False,
            zeroline=False,
            scaleanchor="x",
            scaleratio=1,
            fixedrange=False,
        ),
    )
    return fig


def build_interactive_lattice_html(lattice: KagomeLattice) -> str:
    """
    Build a client-side Plotly interface.

    Vertex clicks are handled entirely in JavaScript, so X/Z/CZ annotations do
    not trigger a Streamlit rerun.  Streamlit only reruns when nx or ny changes.
    """
    fig = _build_figure(lattice)
    figure_json = fig.to_json()

    positions = {
        str(site.site_id): {"x": site.x, "y": site.y}
        for site in lattice.sites
    }
    positions_json = json.dumps(positions)

    # Plotly.js is loaded in the iframe.  All operation state is intentionally
    # client-side for the current UI-only phase of the project.
    return f"""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8" />
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
    * {{ box-sizing: border-box; }}
    body {{
        margin: 0;
        background: white;
        color: #0f172a;
        font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont,
                     "Segoe UI", sans-serif;
    }}
    .toolbar {{
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: 8px;
        padding: 4px 2px 10px 2px;
    }}
    .toolbar-label {{
        font-size: 0.95rem;
        font-weight: 650;
        margin-right: 4px;
    }}
    button {{
        border: 1px solid #cbd5e1;
        background: #ffffff;
        color: #0f172a;
        border-radius: 8px;
        padding: 7px 14px;
        cursor: pointer;
        font-weight: 650;
        font-size: 0.92rem;
    }}
    button:hover {{ background: #f8fafc; }}
    button.operator.active {{
        border-color: #111827;
        background: #111827;
        color: white;
    }}
    button.utility {{
        margin-left: 4px;
        font-weight: 550;
        color: #475569;
    }}
    #status {{
        margin-left: 6px;
        font-size: 0.88rem;
        color: #64748b;
        min-height: 1.2em;
    }}
    #plot {{
        width: 100%;
        height: 720px;
    }}
</style>
</head>
<body>
    <div class="toolbar">
        <span class="toolbar-label">Operator</span>
        <button id="btn-Z" class="operator active" onclick="setOperator('Z')">Z</button>
        <button id="btn-X" class="operator" onclick="setOperator('X')">X</button>
        <button id="btn-CZ" class="operator" onclick="setOperator('CZ')">CZ</button>
        <button class="utility" onclick="undoLast()">Undo</button>
        <button class="utility" onclick="clearAll()">Clear</button>
        <span id="status">Z mode: click a qubit to toggle a Z label.</span>
    </div>
    <div id="plot"></div>

<script>
const fig = {figure_json};
const positions = {positions_json};
const plotDiv = document.getElementById('plot');

const config = {{
    responsive: true,
    displaylogo: false,
    scrollZoom: true,
    modeBarButtonsToRemove: ['select2d', 'lasso2d'],
}};

// UI state.  This lives in the browser, so a qubit click does not rerun Streamlit.
let activeOperator = 'Z';
let localOps = new Map();      // site_id -> 'X' or 'Z'
let czPairs = [];              // array of [site_a, site_b]
let pendingCZ = null;          // first endpoint while constructing a CZ
let history = [];              // {{state, description}} entries for Undo

function snapshot() {{
    return {{
        localOps: Array.from(localOps.entries()),
        czPairs: czPairs.map(pair => [...pair]),
        pendingCZ: pendingCZ,
    }};
}}

function restore(state) {{
    localOps = new Map(state.localOps);
    czPairs = state.czPairs.map(pair => [...pair]);
    pendingCZ = state.pendingCZ;
    renderOperations();
}}

function pushHistory(description, stateOverride = null) {{
    history.push({{
        state: stateOverride || snapshot(),
        description: description,
    }});
    if (history.length > 100) history.shift();
}}

function setOperator(op) {{
    activeOperator = op;
    pendingCZ = null;

    document.querySelectorAll('button.operator').forEach(btn => btn.classList.remove('active'));
    document.getElementById('btn-' + op).classList.add('active');

    if (op === 'CZ') {{
        setStatus('CZ mode: click the first qubit, then click the second qubit.');
    }} else {{
        setStatus(op + ' mode: click a qubit to toggle a ' + op + ' label.');
    }}
    renderOperations();
}}

function setStatus(message) {{
    document.getElementById('status').textContent = message;
}}

function normalizedPair(a, b) {{
    return a < b ? [a, b] : [b, a];
}}

function pairIndex(a, b) {{
    const [u, v] = normalizedPair(a, b);
    return czPairs.findIndex(pair => pair[0] === u && pair[1] === v);
}}

function applyLocalOperator(siteId) {{
    const current = localOps.get(siteId);

    // Clicking the same operator again removes it; choosing the other Pauli
    // replaces the previous local label on that qubit.
    if (current === activeOperator) {{
        pushHistory('removing ' + activeOperator + ' from qubit ' + siteId);
        localOps.delete(siteId);
        setStatus(activeOperator + ' removed from qubit ' + siteId + '.');
    }} else if (current) {{
        pushHistory('replacing ' + current + ' with ' + activeOperator + ' on qubit ' + siteId);
        localOps.set(siteId, activeOperator);
        setStatus(activeOperator + ' applied to qubit ' + siteId + ', replacing ' + current + '.');
    }} else {{
        pushHistory('applying ' + activeOperator + ' to qubit ' + siteId);
        localOps.set(siteId, activeOperator);
        setStatus(activeOperator + ' applied to qubit ' + siteId + '.');
    }}
    renderOperations();
}}

function applyCZClick(siteId) {{
    if (pendingCZ === null) {{
        pendingCZ = siteId;
        setStatus('CZ: first qubit ' + siteId + ' selected. Choose the second qubit.');
        renderOperations();
        return;
    }}

    if (pendingCZ === siteId) {{
        pendingCZ = null;
        setStatus('CZ selection cancelled. Choose the first qubit again.');
        renderOperations();
        return;
    }}

    const first = pendingCZ;
    const existingIndex = pairIndex(first, siteId);

    // Save the state as it was before the first CZ endpoint was selected.
    // This makes one press of Undo reverse the complete two-qubit CZ action.
    const beforeCZ = snapshot();
    beforeCZ.pendingCZ = null;

    if (existingIndex >= 0) {{
        pushHistory('removing CZ between qubits ' + first + ' and ' + siteId, beforeCZ);
        czPairs.splice(existingIndex, 1);
        setStatus('CZ between qubits ' + first + ' and ' + siteId + ' removed.');
    }} else {{
        pushHistory('applying CZ to qubits ' + first + ' and ' + siteId, beforeCZ);
        czPairs.push(normalizedPair(first, siteId));
        setStatus('CZ created between qubits ' + first + ' and ' + siteId + '.');
    }}
    pendingCZ = null;
    renderOperations();
}}

function renderOperations() {{
    const labels = [];
    const nSites = fig.data[1].x.length;

    // Count how many completed CZ gates touch each site.
    const czCount = new Map();
    czPairs.forEach(([a, b]) => {{
        czCount.set(a, (czCount.get(a) || 0) + 1);
        czCount.set(b, (czCount.get(b) || 0) + 1);
    }});

    for (let siteId = 0; siteId < nSites; siteId++) {{
        const local = localOps.get(siteId) || '';
        const hasCZ = (czCount.get(siteId) || 0) > 0;
        const isPending = pendingCZ === siteId;

        if (local && hasCZ) {{
            labels.push('<b>' + local + '</b><br><span style="font-size:8px">CZ</span>');
        }} else if (local) {{
            labels.push('<b>' + local + '</b>');
        }} else if (hasCZ) {{
            labels.push('<b>CZ</b>');
        }} else if (isPending) {{
            labels.push('<b>CZ?</b>');
        }} else {{
            labels.push('');
        }}
    }}

    Plotly.restyle(plotDiv, {{text: [labels]}}, [1]);

    // Rebuild the CZ connection trace and its midpoint labels.
    const lineX = [];
    const lineY = [];
    const labelX = [];
    const labelY = [];
    const labelText = [];

    czPairs.forEach(([a, b]) => {{
        const pa = positions[String(a)];
        const pb = positions[String(b)];
        lineX.push(pa.x, pb.x, null);
        lineY.push(pa.y, pb.y, null);
        labelX.push((pa.x + pb.x) / 2.0);
        labelY.push((pa.y + pb.y) / 2.0);
        labelText.push('CZ');
    }});

    Plotly.restyle(plotDiv, {{x: [lineX], y: [lineY]}}, [2]);
    Plotly.restyle(plotDiv, {{x: [labelX], y: [labelY], text: [labelText]}}, [3]);
}}

function undoLast() {{
    // If the user has only selected the first endpoint of a CZ gate, Undo
    // simply cancels that pending selection.
    if (pendingCZ !== null) {{
        const first = pendingCZ;
        pendingCZ = null;
        renderOperations();
        setStatus('Undone selecting qubit ' + first + ' as the first CZ qubit.');
        return;
    }}

    if (history.length === 0) {{
        setStatus('Nothing to undo.');
        return;
    }}

    const entry = history.pop();
    restore(entry.state);
    setStatus('Undone ' + entry.description + '.');
}}

function clearAll() {{
    if (localOps.size === 0 && czPairs.length === 0 && pendingCZ === null) {{
        setStatus('Nothing to clear.');
        return;
    }}
    pushHistory('clearing all operator annotations');
    localOps.clear();
    czPairs = [];
    pendingCZ = null;
    renderOperations();
    setStatus('All operator annotations cleared.');
}}

Plotly.newPlot(plotDiv, fig.data, fig.layout, config).then(() => {{
    plotDiv.on('plotly_click', eventData => {{
        // Only the qubit trace carries our customdata payload.
        const point = eventData.points.find(point => point.customdata && point.customdata.length > 0);
        if (!point) return;

        const siteId = Number(point.customdata[0]);
        if (activeOperator === 'CZ') {{
            applyCZClick(siteId);
        }} else {{
            applyLocalOperator(siteId);
        }}
    }});
}});

window.addEventListener('resize', () => Plotly.Plots.resize(plotDiv));
</script>
</body>
</html>
"""
