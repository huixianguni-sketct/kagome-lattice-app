from __future__ import annotations

import json
import math


def _coerce_hex_color(value: str) -> str:
    """
    Accept either a hex string or a color name and return a hex color.
    """
    color_map = {
        "red": "#ef4444",
        "green": "#22c55e",
        "blue": "#3b82f6",
        "Red": "#ef4444",
        "Green": "#22c55e",
        "Blue": "#3b82f6",
    }

    if not isinstance(value, str):
        return "#9ca3af"

    if value.startswith("#") or value.startswith("rgb"):
        return value

    return color_map.get(value, "#9ca3af")


def _site_color(site) -> str:
    """
    Try a few possible attribute names for color.
    """
    for attr in ("color_hex", "hex_color", "plot_color", "display_color", "color", "color_name"):
        if hasattr(site, attr):
            return _coerce_hex_color(getattr(site, attr))
    return "#9ca3af"


def _canonical_cycle(cycle):
    """Canonicalize a cyclic list up to rotation and reversal."""
    seq = list(cycle)
    n = len(seq)
    variants = []

    for arr in (seq, list(reversed(seq))):
        for shift in range(n):
            variants.append(tuple(arr[shift:] + arr[:shift]))

    return min(variants)


def _detect_complete_stars(lattice):
    """
    Find complete 12-body D4 star operators in the visible open patch.

    Fig. 2 of the paper shows that a star is NOT simply the six vertices of a
    Kagome hexagon.  Instead it contains:

      * inner_qubits: the six vertices of the central Kagome hexagon.  The
        six CZ gates form a closed ring on these vertices.

      * x_qubits: the six outward-pointing triangle tips surrounding that
        hexagon.  X acts on these six qubits.

    All six x_qubits have the same Red/Green/Blue colour, which is the colour
    of the star.  If even one of the six outward tips is missing at an open
    boundary, that star is excluded.
    """
    sites = sorted(lattice.sites, key=lambda s: int(s.site_id))
    site_by_id = {int(s.site_id): s for s in sites}
    neighbors = {
        int(q): {int(v) for v in values}
        for q, values in lattice.neighbors.items()
    }

    cycles = set()

    # First locate elementary Kagome hexagons as induced 6-cycles.
    for start in sorted(neighbors):

        def dfs(current, path):
            if len(path) == 6:
                if start in neighbors[current]:
                    cyc = _canonical_cycle(path)
                    subset = set(cyc)

                    internal_edges = (
                        sum(
                            1
                            for q in subset
                            for v in neighbors[q]
                            if v in subset
                        )
                        // 2
                    )

                    if internal_edges == 6:
                        cycles.add(cyc)

                return

            for nxt in neighbors[current]:
                if nxt == start:
                    continue
                if nxt in path:
                    continue
                if nxt < start:
                    continue

                dfs(nxt, path + [nxt])

        dfs(start, [start])

    stars = []

    for cyc in cycles:
        cx = sum(float(site_by_id[q].x) for q in cyc) / 6.0
        cy = sum(float(site_by_id[q].y) for q in cyc) / 6.0

        # Geometric cyclic ordering of the six inner hexagon vertices.
        inner = sorted(
            cyc,
            key=lambda q: math.atan2(
                float(site_by_id[q].y) - cy,
                float(site_by_id[q].x) - cx,
            ),
        )

        inner_set = set(inner)
        x_qubits = []
        complete = True

        # Every edge of a complete Kagome hexagon belongs to a small triangle.
        # The third vertex of that triangle is one of the six X-support tips.
        for idx, a in enumerate(inner):
            b = inner[(idx + 1) % 6]

            common_external = (
                neighbors[a]
                .intersection(neighbors[b])
                .difference(inner_set)
            )

            if len(common_external) != 1:
                complete = False
                break

            x_qubits.append(next(iter(common_external)))

        # Open-boundary rule: ignore any star that does not have all 12 qubits.
        if not complete or len(set(x_qubits)) != 6:
            continue

        tip_colors = {str(site_by_id[q].color) for q in x_qubits}

        # In the paper's three-colouring, all six X tips of one star have the
        # same colour.  Reject anything inconsistent with that geometry.
        if len(tip_colors) != 1:
            continue

        star_color = next(iter(tip_colors))

        stars.append(
            {
                "inner_qubits": [int(q) for q in inner],
                "x_qubits": [int(q) for q in x_qubits],
                "cz_edges": [
                    [
                        int(inner[idx]),
                        int(inner[(idx + 1) % 6]),
                    ]
                    for idx in range(6)
                ],
                "color": star_color,
                "center": [cx, cy],
            }
        )

    stars.sort(key=lambda item: (item["center"][1], item["center"][0]))

    for idx, star in enumerate(stars):
        star["id"] = idx

    return stars

def _detect_complete_triangles(lattice):
    """
    Find every complete 3-body Bt triangle in the visible open patch.

    The Bt terms are the nearest equilateral triangles formed by three qubits
    of the same Red/Green/Blue colour.  Their side length is sqrt(3) when the
    Kagome nearest-neighbour spacing is 1.
    """
    sites = sorted(lattice.sites, key=lambda s: int(s.site_id))
    site_by_id = {int(s.site_id): s for s in sites}

    ids_by_color = {"Red": [], "Green": [], "Blue": []}

    for site in sites:
        color = str(site.color)
        if color in ids_by_color:
            ids_by_color[color].append(int(site.site_id))

    triangles = []
    target_distance_squared = 3.0
    tol = 1e-8

    for color, ids in ids_by_color.items():
        same_color_neighbors = {q: set() for q in ids}

        for i, a in enumerate(ids):
            sa = site_by_id[a]

            for b in ids[i + 1 :]:
                sb = site_by_id[b]

                dx = float(sa.x) - float(sb.x)
                dy = float(sa.y) - float(sb.y)
                d2 = dx * dx + dy * dy

                if abs(d2 - target_distance_squared) < tol:
                    same_color_neighbors[a].add(b)
                    same_color_neighbors[b].add(a)

        # Every 3-clique of this nearest-same-colour graph is one Bt triangle.
        for a in ids:
            for b in sorted(q for q in same_color_neighbors[a] if q > a):
                common = same_color_neighbors[a].intersection(
                    same_color_neighbors[b]
                )

                for c in sorted(q for q in common if q > b):
                    cx = (
                        float(site_by_id[a].x)
                        + float(site_by_id[b].x)
                        + float(site_by_id[c].x)
                    ) / 3.0

                    cy = (
                        float(site_by_id[a].y)
                        + float(site_by_id[b].y)
                        + float(site_by_id[c].y)
                    ) / 3.0

                    ordered = sorted(
                        [a, b, c],
                        key=lambda q: math.atan2(
                            float(site_by_id[q].y) - cy,
                            float(site_by_id[q].x) - cx,
                        ),
                    )

                    triangles.append(
                        {
                            "qubits": [int(q) for q in ordered],
                            "color": color,
                            "center": [cx, cy],
                        }
                    )

    triangles.sort(
        key=lambda item: (
            item["center"][1],
            item["center"][0],
            item["color"],
        )
    )

    for idx, triangle in enumerate(triangles):
        triangle["id"] = idx

    return triangles

def build_interactive_lattice_html(lattice) -> str:
    """
    Build the interactive Plotly/HTML component used inside Streamlit.

    Expected lattice interface:
      - lattice.sites : iterable of site objects with attributes
            site_id, x, y, and a color-ish attribute
      - lattice.edges : iterable of (i, j) nearest-neighbour pairs
    """

    sites = sorted(lattice.sites, key=lambda s: int(s.site_id))
    edge_pairs = [(int(i), int(j)) for i, j in lattice.edges]

    # Open-boundary physics objects: only complete stars/triangles are kept.
    # Guided electric strings use complete same-colour star dual edges.
    # Guided magnetic strings use the graph of complete monochromatic Bt
    # triangles: an eligible X qubit joins two of those triangles. Both
    # guides respect the requested open-boundary convention.
    star_data = _detect_complete_stars(lattice)
    triangle_data = _detect_complete_triangles(lattice)

    site_data = [
        {
            "id": int(site.site_id),
            "x": float(site.x),
            "y": float(site.y),
            "color": _site_color(site),
            "colorName": str(site.color),
        }
        for site in sites
    ]

    # Base lattice edges for the permanent grey Kagome bonds
    base_edge_x = []
    base_edge_y = []
    for i, j in edge_pairs:
        si = sites[i]
        sj = sites[j]
        base_edge_x.extend([float(si.x), float(sj.x), None])
        base_edge_y.extend([float(si.y), float(sj.y), None])

    template = r"""
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8"/>
  <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
  <style>
    html, body {
      margin: 0;
      padding: 0;
      background: white;
      font-family: Arial, Helvetica, sans-serif;
    }

    .wrapper {
      width: 100%;
      box-sizing: border-box;
      padding: 8px 8px 0 8px;
    }

    .toolbar {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      align-items: center;
      margin-bottom: 10px;
    }

    .toolbar-label {
      font-weight: 700;
      font-size: 14px;
      margin-right: 4px;
    }

    .anyon-toolbar {
      margin-bottom: 12px;
    }

    .anyon-selector {
      border: 1px solid #cbd5e1;
      border-radius: 8px;
      padding: 7px 10px;
      min-width: 186px;
      background: #ffffff;
      color: #111827;
      font-size: 14px;
      font-weight: 600;
      cursor: pointer;
    }

    .anyon-hint {
      color: #64748b;
      font-size: 12px;
    }

    .anyon-summary {
      color: #0f766e;
      font-size: 12px;
      font-weight: 700;
      white-space: nowrap;
    }

    .visibility-toolbar {
      margin-top: -4px;
      margin-bottom: 10px;
    }

    .violation-toolbar {
      margin-top: -4px;
      margin-bottom: 10px;
    }

    .btn {
      border: 1px solid #cbd5e1;
      background: #f8fafc;
      color: #111827;
      padding: 7px 12px;
      border-radius: 8px;
      cursor: pointer;
      font-size: 14px;
      font-weight: 600;
    }

    .btn:hover {
      background: #eef2ff;
    }

    .btn.active {
      background: #dbeafe;
      border-color: #60a5fa;
      color: #1d4ed8;
    }

    .btn.action {
      background: #f9fafb;
    }

    .status {
      margin: 4px 0 8px 0;
      min-height: 20px;
      font-size: 14px;
      color: #374151;
    }

    .content-layout {
      display: grid;
      grid-template-columns: minmax(0, 1fr) 340px;
      gap: 16px;
      align-items: start;
    }

    .plot-pane {
      min-width: 0;
    }

    .qubit-panel {
      border: 1px solid #e2e8f0;
      border-radius: 10px;
      background: #ffffff;
      padding: 12px;
      box-sizing: border-box;
      max-height: 760px;
      overflow-y: auto;
    }

    .qubit-panel h3 {
      margin: 0 0 10px 0;
      font-size: 16px;
      color: #111827;
    }

    .tracker-controls {
      display: flex;
      gap: 8px;
      margin-bottom: 12px;
    }

    .tracker-controls select {
      flex: 1;
      min-width: 0;
      border: 1px solid #cbd5e1;
      border-radius: 8px;
      padding: 7px 9px;
      background: white;
      color: #111827;
      font-size: 14px;
    }

    .tracker-table {
      width: 100%;
      border-collapse: collapse;
      table-layout: fixed;
      font-size: 13px;
    }

    .tracker-table th,
    .tracker-table td {
      border-bottom: 1px solid #e5e7eb;
      padding: 8px 6px;
      text-align: left;
      vertical-align: top;
      overflow-wrap: anywhere;
    }

    .tracker-table th {
      color: #475569;
      font-weight: 700;
      background: #f8fafc;
    }

    .tracker-table th:first-child,
    .tracker-table td:first-child {
      width: 58px;
    }

    .tracker-table th:last-child,
    .tracker-table td:last-child {
      width: 32px;
      text-align: center;
    }

    .op-sequence {
      display: flex;
      flex-direction: column;
      gap: 4px;
    }

    .op-entry {
      line-height: 1.3;
      color: #111827;
    }

    .empty-ops {
      color: #94a3b8;
    }

    .tracker-empty {
      margin: 8px 0 0 0;
      color: #64748b;
      font-size: 13px;
    }

    .remove-track {
      border: none;
      background: transparent;
      color: #64748b;
      cursor: pointer;
      font-size: 18px;
      line-height: 1;
      padding: 0 2px;
    }

    .remove-track:hover {
      color: #dc2626;
    }

    .expectation-panel {
      margin-top: 20px;
      padding-top: 15px;
      border-top: 1px solid #dbe3ee;
      color: #1f2937;
      font-size: 12px;
      line-height: 1.55;
    }

    .expectation-panel h3 { margin-bottom: 8px; }
    .expectation-intro { color: #64748b; font-size: 12px; }
    .expectation-result { margin: 10px 0; font-size: 16px; font-weight: 800; }
    .expectation-result.good { color: #166534; }
    .expectation-result.bad { color: #991b1b; }
    .expectation-panel code { background: #f1f5f9; padding: 1px 3px; border-radius: 3px; overflow-wrap: anywhere; }
    .expectation-panel details { margin: 8px 0; }
    .expectation-panel summary { cursor: pointer; font-weight: 700; }
    .expectation-steps { list-style: decimal; padding-left: 22px; max-height: 235px; overflow-y: auto; }
    .expectation-steps li { margin: 5px 0; overflow-wrap: anywhere; }
    .expectation-equation { background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 9px; overflow-wrap: anywhere; }
    .expectation-qubits { font-size: 11px; overflow-wrap: anywhere; color: #475569; }
    .expectation-note { color: #64748b; font-size: 11px; }

    #plot {
      width: 100%;
      height: 760px;
    }

    @media (max-width: 900px) {
      .content-layout {
        grid-template-columns: 1fr;
      }

      .qubit-panel {
        max-height: none;
      }
    }
 
    /* Kawagoe-Levin T-junction experiment (independent second canvas). */
    .r-studio[hidden] { display: none !important; }
    .r-studio {margin:20px 0 28px; border-top:2px solid #e2e8f0; padding-top:15px;}
    .r-toolbar {display:flex;align-items:center;flex-wrap:wrap;gap:10px;margin:10px 0;}
    .r-heading {font-size:17px;font-weight:800;color:#0f172a;}
    .r-step-row {display:flex;align-items:center;flex-wrap:wrap;gap:7px;margin:8px 0;}
    .r-process-title {font-size:13px;color:#334155;font-weight:800;min-width:105px;}
    .r-step {font-size:13px;padding:7px 9px;}
    .r-step.r-done:not(.active){border-color:#86efac;background:#f0fdf4;}
    .r-overview,.r-stage-caption {color:#475569;font-size:12px;line-height:1.55;margin:9px 0;}
    .r-layout {display:grid;grid-template-columns:minmax(0,1fr) 350px;gap:14px;align-items:start;}
    #r-plot {height:530px;min-width:0;width:100%;border:1px solid #dbe3ee;border-radius:9px;}
    .r-panel {border:1px solid #dbe3ee;border-radius:9px;padding:12px;max-height:545px;overflow-y:auto;background:#fff;font-size:12px;line-height:1.55;overflow-wrap:anywhere;}
    .r-panel h3 {font-size:16px;color:#0f172a;margin:0 0 10px;}
    .r-panel label{display:block;margin:8px 0 3px;font-weight:700;}
    .r-panel input{box-sizing:border-box;width:100%;padding:6px 8px;border:1px solid #cbd5e1;border-radius:6px;font-size:12px;}
    .r-panel button{margin-top:8px;}
    .r-math{background:#f8fafc;border:1px solid #e2e8f0;padding:9px;border-radius:7px;margin:10px 0;}
    .r-error{color:#9f1239;font-weight:700;}
    .r-caveat{color:#64748b;font-size:11px;}
    .r-fixed-note{font-size:12px;font-weight:700;color:#065f46;background:#ecfdf5;border:1px solid #a7f3d0;padding:7px;border-radius:6px;}
    .r-validation{padding-bottom:10px;}
    .r-validation h3{font-size:14px;margin:0 0 7px;}
    .r-validation select{width:100%;padding:7px;border:1px solid #cbd5e1;border-radius:6px;background:#fff;color:#0f172a;font-size:12px;}
    .r-validation table{width:100%;border-collapse:collapse;font-size:11px;margin:8px 0;}
    .r-validation th,.r-validation td{border-bottom:1px solid #e2e8f0;padding:5px 3px;text-align:left;vertical-align:top;}
    .r-validation th{color:#475569;font-weight:700;}
    .r-v-ok{color:#166534;font-weight:700;}.r-v-bad{color:#9f1239;font-weight:700;}.r-v-neutral{color:#64748b;}
    .r-validation .r-validation-codes{font-family:ui-monospace,SFMono-Regular,Consolas,monospace;overflow-wrap:anywhere;font-size:10px;}
    .r-validation .r-v-row{margin:6px 0;}
    .r-validation .r-v-hint{color:#64748b;font-size:11px;margin:5px 0;}

    .r-symbolic-badge{color:#075985;background:#f0f9ff;border:1px solid #bae6fd;padding:6px 8px;border-radius:6px;font-weight:700;}
    .r-symbolic-table{width:100%;border-collapse:collapse;font-size:11px;margin:8px 0;}
    .r-symbolic-table th,.r-symbolic-table td{padding:5px 3px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}
    .r-symbolic-table small{color:#64748b;}
    #r-calc details{margin:9px 0;}
    #r-calc summary{cursor:pointer;}
    @media(max-width:900px){.r-layout{grid-template-columns:1fr}.r-panel{max-height:none}}

  </style>
</head>
<body>
  <div class="wrapper">
    <div class="toolbar anyon-toolbar">
      <span class="toolbar-label">Anyon</span>
      <select class="anyon-selector" id="anyon-select" aria-label="Choose an anyon to construct">
        <option value="">Select an anyon…</option>
        <option value="e_R">e_R</option>
        <option value="e_B">e_B</option>
        <option value="e_G">e_G</option>
        <option value="m_R">m_R</option>
        <option value="m_B">m_B</option>
        <option value="m_G">m_G</option>
        <option value="m_BG">m_BG (Blue then Green)</option>
        <option value="m_BR">m_BR (Blue then Red)</option>
        <option value="m_RG">m_RG (Red then Green)</option>
        <option value="s1">s1 (straight: horizontal / diagonal)</option>
        <option value="s2">s2 (mirror: horizontal / diagonal)</option>
      </select>
      <span class="anyon-hint" id="anyon-hint">Choose an anyon, then click a highlighted qubit.</span>
      <span class="anyon-summary" id="anyon-summary"></span>
    </div>

    <div class="toolbar">
      <span class="toolbar-label">Operator</span>

      <button class="btn mode-btn active" data-mode="Z">Z</button>
      <button class="btn mode-btn" data-mode="X">X</button>
      <button class="btn mode-btn" data-mode="R" title="R = exp(+iπ Z/4)">R (π/4)</button>
      <button class="btn mode-btn" data-mode="CZ">CZ</button>
      <button class="btn mode-btn" data-mode="ERASE">Erase</button>

      <button class="btn action" id="zoom-btn">Zoom</button>
      <button class="btn action" id="save-btn">Save PNG</button>
      <button class="btn action" id="undo-btn">Undo</button>
      <button class="btn action" id="clear-btn">Clear</button>
    </div>

    <div class="toolbar visibility-toolbar">
      <span class="toolbar-label">Visibility</span>

      <button class="btn visibility-btn active" data-visibility="Z"
              title="Show or hide Z labels">Z</button>

      <button class="btn visibility-btn active" data-visibility="X"
              title="Show or hide X labels">X</button>
      <button class="btn visibility-btn active" data-visibility="R"
              title="Show or hide R labels (also within XR/RX)">R</button>

      <button class="btn visibility-btn active" data-visibility="CZ"
              title="Show or hide CZ connections and labels">CZ</button>
    </div>

    <div class="toolbar violation-toolbar">
      <span class="toolbar-label">Violation</span>

      <button
        class="btn violation-btn"
        data-violation="STAR"
        title="Show or hide violated star operators">
        Star (Aₛ)
      </button>

      <button
        class="btn violation-btn"
        data-violation="TRIANGLE"
        title="Show or hide violated triangle operators">
        Triangle (Bₜ)
      </button>

      <button class="btn expectation-btn" data-expectation="STAR"
        title="Click a complete hexagonal star to inspect its expectation">Star exp.</button>
      <button class="btn expectation-btn" data-expectation="TRIANGLE"
        title="Click a complete monochromatic triangle to inspect its expectation">Triangle exp.</button>
    </div>

    <div class="status" id="status">Mode: Z</div>

    <div class="content-layout">
      <div class="plot-pane">
        <div id="plot"></div>
      </div>

      <aside class="qubit-panel">
        <h3>Tracked qubits</h3>

        <div class="tracker-controls">
          <select id="qubit-select" aria-label="Choose qubit to track"></select>
          <button class="btn action" id="add-qubit-btn">Add</button>
        </div>

        <table class="tracker-table" aria-label="Tracked qubit operations">
          <thead>
            <tr>
              <th>Qubit</th>
              <th>Applied operators</th>
              <th></th>
            </tr>
          </thead>
          <tbody id="tracker-table-body"></tbody>
        </table>

        <p class="tracker-empty" id="tracker-empty">
          No qubits added yet.
        </p>

        <section id="expectation-panel" class="expectation-panel" aria-live="polite">
          <h3>Expectation calculation</h3>
          <p class="expectation-intro" id="expectation-content">
            Press <strong>Star exp.</strong> or <strong>Triangle exp.</strong>, then click inside a complete polygon on the lattice.
          </p>
        </section>
      </aside>
    </div>
    <div class="r-toolbar">
      <span class="toolbar-label">R calculation</span>
      <select class="anyon-selector" id="r-select" aria-label="Choose a candidate string for microscopic R calculation">
        <option value="">Select s1 or s2…</option>
        <option value="s1">s1</option>
        <option value="s2">s2</option>
      </select>
      <span class="anyon-hint">Shows the independent Kawagoe–Levin T-junction protocol below.</span>
    </div>

    <section class="r-studio" id="r-studio" hidden>
      <div class="r-toolbar">
        <span class="r-heading">Microscopic R calculation — symbolic</span>
      </div>
      <p id="r-overview" class="r-overview"></p>
      <div class="r-step-row"><span class="r-process-title">Process (1)</span>
        <button class="btn r-step" data-r-process="A" data-r-step="1">1. Split S(1,2)</button>
        <button class="btn r-step" data-r-process="A" data-r-step="2">2. Move 1 → Y</button>
        <button class="btn r-step" data-r-process="A" data-r-step="3">3. Move 2 → 1</button>
        <button class="btn r-step" data-r-process="A" data-r-step="4">4. Move 1 → 0</button>
      </div>
      <div class="r-step-row"><span class="r-process-title">Process (2)</span>
        <button class="btn r-step" data-r-process="B" data-r-step="1">1. Split S(1,2)</button>
        <button class="btn r-step" data-r-process="B" data-r-step="2">2. Move 1 → 0</button>
        <button class="btn r-step" data-r-process="B" data-r-step="3">3. Move 2 → 1</button>
        <button class="btn r-step" data-r-process="B" data-r-step="4">4. Move 1 → Y</button>
      </div>
      <div id="r-stage-caption" class="r-stage-caption"></div>
      <div class="r-layout">
        <div id="r-plot" aria-label="Selected microscopic T-junction process lattice"></div>
        <aside class="r-panel">
          <section class="r-validation" aria-label="Validate selected microscopic movement">
            <h3>Interface validation</h3>
            <p class="r-v-hint">Inspect each stage of the currently selected process. This only tests the specified finite-patch transport circuit; it does not assume a semion phase.</p>
            <label for="r-checkpoint">Display on lattice</label>
            <select id="r-checkpoint" aria-label="Displayed checkpoint within selected movement">
              <option value="before">Before this movement</option>
              <option value="ribbon">After fresh s ribbon</option>
              <option value="x">After interface X</option>
              <option value="projected" selected>After interface projectors</option>
            </select>
            <div id="r-validation-content" aria-live="polite">Choose a movement button to inspect its interface.</div>
          </section>
          <h3 style="margin-top:14px;border-top:1px solid #e2e8f0;padding-top:10px;">R overlap / consistency checks</h3>
          <div id="r-calc" aria-live="polite">Complete both processes to evaluate their overlap.</div>
          <details open>
            <summary><strong>Interface corrections for displayed move</strong></summary>
            <p class="r-caveat">The suggested X and star-projector IDs are <strong>not</strong> guaranteed by the drawing. Verify them against the physical interface. A +1 star projection is a postselected measurement, not an X/Z/CZ gate.</p>
            <p class="r-fixed-note" id="r-fixed-x-note" aria-live="polite"></p>
            <label for="r-edit-x">Additional interface X qubit IDs, comma-separated</label>
            <input id="r-edit-x" placeholder="e.g. 42, 73" aria-label="Interface X correction qubit IDs"/>
            <label for="r-edit-stars">Complete star IDs to project onto A=+1</label>
            <input id="r-edit-stars" placeholder="e.g. 5, 8, 9" aria-label="Interface projected star IDs"/>
            <button class="btn" id="r-save-adjust" disabled>Apply to this movement (both processes)</button>
          </details>
          <p class="r-caveat">The calculation uses the exact chronological X/R/CZ gates and ordered P=(I+A)/2 projectors. It does not presuppose that s1 or s2 is a semion. The symbolic engine preserves complex phases and normalizes every successful projection without enumerating 2^N amplitudes. The finite-open-patch reference state is unchanged.</p>
        </aside>
      </div>
    </section>
  </div>

  <script>
    const siteData = __SITE_DATA__;
    const edgePairs = __EDGE_PAIRS__;
    const baseEdgeX = __BASE_EDGE_X__;
    const baseEdgeY = __BASE_EDGE_Y__;
    const starData = __STAR_DATA__;
    const triangleData = __TRIANGLE_DATA__;

    const plotDiv = document.getElementById("plot");
    const statusDiv = document.getElementById("status");
    const anyonSelect = document.getElementById("anyon-select");
    const anyonHint = document.getElementById("anyon-hint");
    const anyonSummary = document.getElementById("anyon-summary");
    const modeButtons = Array.from(document.querySelectorAll(".mode-btn"));
    const visibilityButtons = Array.from(document.querySelectorAll(".visibility-btn"));
    const violationButtons = Array.from(document.querySelectorAll(".violation-btn"));
    const expectationButtons = Array.from(document.querySelectorAll(".expectation-btn"));
    const expectationContent = document.getElementById("expectation-content");
    const expectationPanel = document.getElementById("expectation-panel");
    const zoomBtn = document.getElementById("zoom-btn");
    const saveBtn = document.getElementById("save-btn");
    const undoBtn = document.getElementById("undo-btn");
    const clearBtn = document.getElementById("clear-btn");
    const qubitSelect = document.getElementById("qubit-select");
    const addQubitBtn = document.getElementById("add-qubit-btn");
    const trackerTableBody = document.getElementById("tracker-table-body");
    const trackerEmpty = document.getElementById("tracker-empty");

    const siteById = new Map(siteData.map(s => [s.id, s]));

    const nearestNeighborKeys = new Set(
      edgePairs.map(([a, b]) => pairKey(a, b))
    );
    
    const mNeighborIds = new Map(
      siteData.map(s => [s.id, new Set()])
    );
    
    for (const [a, b] of edgePairs) {
      mNeighborIds.get(a).add(b);
      mNeighborIds.get(b).add(a);
    }

    let currentMode = "Z";
    let localOps = new Map();   // site_id -> "X" or "Z"
    let czPairs = [];           // array of [a, b]
    let pendingCZ = null;       // first endpoint for CZ selection
    let history = [];           // array of { state, description }
    let zoomMode = false;

    // Guided Anyon mode: electric (Z), magnetic (X + automatic CZ), or
    // three-color straight s1/s2 diagrams with XR/RX on the spine.
    // A single-colour m circuit is recompiled CZ first, then X.
    // A two-colour m2 ribbon instead keeps TWO ordered factor blocks:
    // first-colour CZ -> X, THEN second-colour CZ -> X.
    let anyonKind = null;      // "e", "m", "m2", "s1" or "s2"; null disables guide
    let anyonColor = null;     // "Red", "Blue", "Green", or null
    let anyonQubits = [];      // selected Z or X physical qubits
    let anyonStarPath = [];    // oriented e-star path (after 2nd Z)
    let mTrianglePath = [];    // oriented monochromatic triangle path
    let mAutoCZ = [];         // automatic CZ pairs in current m ribbon
    let mOptionCache = null;  // computed options for the current m path
    let mGuideId = null;      // unique ID for the current m guide session
    let nextMGuideId = 1;
    // Composite magnetic ribbons: one selected two-colour bond per step.
    // Do not merge their two chronological CZ/X blocks.
    let mPairColors = [];
    let mPairEdges = [];
    let mPairGuideId = null;
    let nextMPairGuideId = 1;
    let mPairBoard = null;
    let mPairOptionCache = null;
    let sBasePath = [];       // clicked straight-line spine qubits, ordered
    let sGuidedX = [];        // all physical X qubits, including triangle tips
    let sAutoCZ = [];         // compiled three-colour magnetic decorations
    let sGuideId = null;
    let nextSGuideId = 1;
    let sCandidateCache = null;

    // A physical qubit carries one colour of e anyon and belongs to zero,
    // one, or two complete stars of that colour in an open-boundary patch.
    const eStarIncidence = new Map();
    for (const star of starData) {
      for (const siteId of star.x_qubits) {
        if (!eStarIncidence.has(siteId)) {
          eStarIncidence.set(siteId, []);
        }
        eStarIncidence.get(siteId).push(star.id);
      }
    }

    // Triangle incidence and the corresponding complete encompassing star.
    // An incomplete star/triangle is not eligible for guided m construction.
    const triangleById = new Map(triangleData.map(t => [t.id, t]));
    const starById = new Map(starData.map(s => [s.id, s]));
    const centerKey = (point) => `${point[0].toFixed(6)}:${point[1].toFixed(6)}`;
    const starByCenter = new Map(starData.map(s => [centerKey(s.center), s.id]));
    const starForTriangle = new Map();
    const mTriangleIncidence = new Map();
    for (const triangle of triangleData) {
      const starId = starByCenter.get(centerKey(triangle.center));
      if (starId === undefined) continue;
      starForTriangle.set(triangle.id, starId);
      for (const id of triangle.qubits) {
        if (!mTriangleIncidence.has(id)) mTriangleIncidence.set(id, []);
        mTriangleIncidence.get(id).push(triangle.id);
      }
    }

    // Display-only visibility state. Turning one of these off does not
    // remove the operator from the lattice state or tracked-qubit table.
    let operatorVisibility = {
      Z: true,
      X: true,
      R: true,
      CZ: true,
    };

    // Violation overlays are display-only.  The physics is always evaluated
    // from the active ordered operationLog.
    let violationVisibility = {
      STAR: false,
      TRIANGLE: false,
    };

    // Inspection is display-only: it never modifies the ordered circuit.
    let expectationMode = null; // "STAR", "TRIANGLE" or null
    let inspectedObject = null; // {kind, id}


    // Ordered list of currently active operations.
    // Each entry is { id, type, qubits }.
    let operationLog = [];
    let nextOperationId = 1;

    // Qubits chosen by the user for the right-hand table.
    let trackedQubits = [];

    function setStatus(text) {
      statusDiv.textContent = text;
    }

    function initializeQubitSelector() {
      qubitSelect.innerHTML = "";

      // Keep the bulk-add action as the first dropdown option.
      const allActiveOption = document.createElement("option");
      allActiveOption.value = "__ALL_ACTIVE__";
      allActiveOption.textContent = "Generate all available option";
      qubitSelect.appendChild(allActiveOption);

      for (const site of siteData) {
        const option = document.createElement("option");
        option.value = String(site.id);
        option.textContent = `Qubit ${site.id}`;
        qubitSelect.appendChild(option);
      }
    }

    function activeQubitIds() {
      const active = new Set();

      // operationLog contains only operations that are still active after
      // erase / undo / clear, so it is the source of truth for this action.
      for (const event of operationLog) {
        for (const siteId of event.qubits) {
          active.add(siteId);
        }
      }

      return Array.from(active).sort((a, b) => a - b);
    }

    function addAllActiveQubits() {
      const activeIds = activeQubitIds();

      if (activeIds.length === 0) {
        setStatus("No qubits currently have operators applied.");
        return;
      }

      let addedCount = 0;

      for (const siteId of activeIds) {
        if (!trackedQubits.includes(siteId)) {
          trackedQubits.push(siteId);
          addedCount += 1;
        }
      }

      updateTrackedTable();

      if (addedCount === 0) {
        setStatus("All qubits with active operators are already being tracked.");
      } else {
        setStatus(`Added ${addedCount} qubit${addedCount === 1 ? "" : "s"} with active operators to the table.`);
      }
    }

    function addTrackedQubit(siteId) {
      if (!siteById.has(siteId)) {
        return;
      }

      if (!trackedQubits.includes(siteId)) {
        trackedQubits.push(siteId);
      }

      updateTrackedTable();
    }

    function removeTrackedQubit(siteId) {
      trackedQubits = trackedQubits.filter(id => id !== siteId);
      updateTrackedTable();
    }

    function operationsForQubit(siteId) {
      return operationLog.filter(event => event.qubits.includes(siteId));
    }

    function operationTextForQubit(event, siteId) {
      if (event.type === "CZ") {
        const partner = event.qubits.find(id => id !== siteId);
        return `CZ with qubit ${partner}`;
      }

      return event.type;
    }

    function updateTrackedTable() {
      trackerTableBody.innerHTML = "";
      trackerEmpty.style.display = trackedQubits.length === 0 ? "block" : "none";

      for (const siteId of trackedQubits) {
        const row = document.createElement("tr");

        const qubitCell = document.createElement("td");
        qubitCell.textContent = String(siteId);

        const opsCell = document.createElement("td");
        const events = operationsForQubit(siteId);

        if (events.length === 0) {
          const empty = document.createElement("span");
          empty.className = "empty-ops";
          empty.textContent = "—";
          opsCell.appendChild(empty);
        } else {
          const sequence = document.createElement("div");
          sequence.className = "op-sequence";

          events.forEach((event, index) => {
            const entry = document.createElement("div");
            entry.className = "op-entry";
            entry.textContent = `${index + 1}. ${operationTextForQubit(event, siteId)}`;
            sequence.appendChild(entry);
          });

          opsCell.appendChild(sequence);
        }

        const removeCell = document.createElement("td");
        const removeButton = document.createElement("button");
        removeButton.className = "remove-track";
        removeButton.type = "button";
        removeButton.title = `Remove qubit ${siteId} from table`;
        removeButton.setAttribute("aria-label", `Remove qubit ${siteId} from table`);
        removeButton.textContent = "×";
        removeButton.addEventListener("click", () => removeTrackedQubit(siteId));
        removeCell.appendChild(removeButton);

        row.appendChild(qubitCell);
        row.appendChild(opsCell);
        row.appendChild(removeCell);
        trackerTableBody.appendChild(row);
      }
    }

    function recordOperation(type, qubits, source = "manual", guideId = null) {
      operationLog.push({
        id: nextOperationId++,
        type: type,
        qubits: [...qubits],
        source: source,
        guideId: guideId,
      });
    }

    function removeLatestLocalOperation(siteId, type) {
      for (let i = operationLog.length - 1; i >= 0; i--) {
        const event = operationLog[i];
        if (event.type === type && event.qubits.length === 1 && event.qubits[0] === siteId) {
          operationLog.splice(i, 1);
          return;
        }
      }
    }

    function removeCZOperation(a, b) {
      const key = pairKey(a, b);

      for (let i = operationLog.length - 1; i >= 0; i--) {
        const event = operationLog[i];
        if (event.type !== "CZ" || event.qubits.length !== 2) {
          continue;
        }

        if (pairKey(event.qubits[0], event.qubits[1]) === key) {
          operationLog.splice(i, 1);
          return;
        }
      }
    }

    function syncVisibleLocalOp(siteId) {
      for (let i = operationLog.length - 1; i >= 0; i--) {
        const event = operationLog[i];
        if ((event.type === "X" || event.type === "Z") && event.qubits.length === 1 && event.qubits[0] === siteId) {
          localOps.set(siteId, event.type);
          return;
        }
      }

      localOps.delete(siteId);
    }

    function pairKey(a, b) {
      return a < b ? `${a}-${b}` : `${b}-${a}`;
    }

    function hasCZPair(a, b) {
      const key = pairKey(a, b);
      return czPairs.some(([u, v]) => pairKey(u, v) === key);
    }

    function isNearestNeighbor(a, b) {
      return nearestNeighborKeys.has(pairKey(a, b));
    }

    function snapshot() {
      return {
        localOps: Array.from(localOps.entries()),
        czPairs: czPairs.map(([a, b]) => [a, b]),
        pendingCZ: pendingCZ,
        currentMode: currentMode,
        anyonKind: anyonKind,
        anyonColor: anyonColor,
        anyonQubits: [...anyonQubits],
        anyonStarPath: [...anyonStarPath],
        mTrianglePath: [...mTrianglePath],
        mAutoCZ: mAutoCZ.map(pair => [...pair]),
        mGuideId: mGuideId,
        nextMGuideId: nextMGuideId,
        mPairColors: [...mPairColors],
        mPairEdges: [...mPairEdges],
        mPairGuideId: mPairGuideId,
        nextMPairGuideId: nextMPairGuideId,
        sBasePath: [...sBasePath],
        sGuidedX: [...sGuidedX],
        sAutoCZ: sAutoCZ.map(pair => [...pair]),
        sGuideId: sGuideId,
        nextSGuideId: nextSGuideId,
        operationLog: operationLog.map(event => ({
          id: event.id,
          type: event.type,
          qubits: [...event.qubits],
          source: event.source || "manual",
          guideId: event.guideId ?? null,
        })),
        nextOperationId: nextOperationId,
      };
    }

    function restoreState(state) {
      localOps = new Map(state.localOps);
      czPairs = state.czPairs.map(([a, b]) => [a, b]);
      pendingCZ = state.pendingCZ;
      currentMode = state.currentMode;
      anyonKind = state.anyonKind || null;
      anyonColor = state.anyonColor || null;
      anyonQubits = [...(state.anyonQubits || [])];
      anyonStarPath = [...(state.anyonStarPath || [])];
      mTrianglePath = [...(state.mTrianglePath || [])];
      mAutoCZ = (state.mAutoCZ || []).map(pair => [...pair]);
      mGuideId = state.mGuideId ?? null;
      nextMGuideId = state.nextMGuideId ?? nextMGuideId;
      mPairColors = [...(state.mPairColors || [])];
      mPairEdges = [...(state.mPairEdges || [])];
      mPairGuideId = state.mPairGuideId ?? null;
      nextMPairGuideId = state.nextMPairGuideId ?? nextMPairGuideId;
      mPairBoard = mPairColors.length === 2 ? buildMPairBoard(mPairColors) : null;
      mPairOptionCache = null;
      sBasePath = [...(state.sBasePath || [])];
      sGuidedX = [...(state.sGuidedX || [])];
      sAutoCZ = (state.sAutoCZ || []).map(pair => [...pair]);
      sGuideId = state.sGuideId ?? null;
      nextSGuideId = state.nextSGuideId ?? nextSGuideId;
      sCandidateCache = null;
      mOptionCache = null;
      anyonSelect.value = anyonKind
        ? (anyonKind === "s1" || anyonKind === "s2"
          ? anyonKind : (anyonKind === "m2"
            ? `m_${mPairColors.map(c => c[0]).join("")}`
            : `${anyonKind}_${anyonColor[0]}`))
        : "";
      operationLog = state.operationLog.map(event => ({
        id: event.id,
        type: event.type,
        qubits: [...event.qubits],
        source: event.source || "manual",
        guideId: event.guideId ?? null,
      }));
      nextOperationId = state.nextOperationId;
      mOptionCache = null;
      updateAnyonSummary();

      updateActiveModeButtons();
      render();
    }

    function pushHistory(description) {
      history.push({
        state: snapshot(),
        description: description,
      });
    }

    // ---------------------------------------------------------------
    // Electric and magnetic string guidance, with open boundaries.
    // ---------------------------------------------------------------
    function anyonName() {
      if (!anyonKind) return "anyon";
      if (anyonKind === "m2") return `m_${mPairColors.map(c => c[0]).join("")}`;
      if (anyonKind === "s1" || anyonKind === "s2") return anyonKind;
      return `${anyonKind}_${anyonColor[0]}`;
    }

    function currentAnyonGate() {
      return (anyonKind === "m" || anyonKind === "m2") ? "X" : (anyonKind === "e" ? "Z" : "XR/RX");
    }

    function sDirectionLabel() {
      if (sBasePath.length < 2) return "choose horizontal or diagonal";
      const [dx,dy]=sStepVector(sBasePath[0],sBasePath[1]);
      if (Math.abs(dy)<1e-5) return "horizontal";
      // Direction is described left-to-right regardless of click order.
      return dx*dy > 0 ? "diagonal rising right" : "diagonal falling right";
    }

    function updateAnyonSummary() {
      if (!anyonKind) {
        anyonHint.textContent = "Choose an anyon, then click a highlighted qubit.";
        anyonSummary.textContent = "";
        return;
      }
      if (anyonKind === "m2") {
        const [first, second] = mPairColors;
        anyonHint.textContent = `Click a gold diamond on a ${first.toLowerCase()}-${second.toLowerCase()} bond. Each step places both Xs on one shared ribbon.`;
        const n = mPairEdges.length;
        const countCZ = operationLog.filter(ev => ev.source === "guided-m2" &&
          ev.guideId === mPairGuideId && ev.type === "CZ").length;
        anyonSummary.textContent = `${n} X per colour · ${countCZ} CZ · ${first}: CZ → X; then ${second}: CZ → X`;
      } else if (anyonKind === "s1" || anyonKind === "s2") {
        anyonHint.textContent = "Select a start, then choose a horizontal or diagonal direction. Existing X/XR/RX strings may overlap; previous gates are retained. Tips and CZs are automatic.";
        anyonSummary.textContent = sBasePath.length < 2
          ? (sBasePath.length ? "Select direction from the highlighted neighbors" : "Up/down triangles: mirrored XR ↔ RX")
          : `${sBasePath.length} spine · ${sDirectionLabel()} · ${sGuidedX.length} X · ${sBasePath.length} R · ${sAutoCZ.length} CZ · gate order: CZ → X/R`;
      } else if (anyonKind === "e") {
        anyonHint.textContent = "Choose highlighted qubits to extend the Z-string.";
        anyonSummary.textContent = anyonQubits.length
          ? `${anyonQubits.length} Z` : "";
      } else {
        anyonHint.textContent = "Choose only highlighted X qubits; CZ links are automatic.";
        anyonSummary.textContent = anyonQubits.length
          ? `${anyonQubits.length} X + ${mAutoCZ.length} auto CZ · gate order: CZ → X`
          : "Gate order: CZ → X";
      }
    }

    function completeEStarEnds(siteId) {
      const site = siteById.get(siteId);
      if (!site || site.colorName !== anyonColor) return [];
      return eStarIncidence.get(siteId) || [];
    }

    function availableEAnyonQubits() {
      const isAvailable = (site) => (
        site.colorName === anyonColor &&
        !localOps.has(site.id) &&
        !anyonQubits.includes(site.id) &&
        completeEStarEnds(site.id).length === 2
      );
      const candidates = siteData.filter(isAvailable);
      if (!anyonQubits.length) return candidates.map(s => s.id);
      if (anyonQubits.length === 1) {
        const starts = new Set(completeEStarEnds(anyonQubits[0]));
        return candidates.filter(site => {
          const ends = completeEStarEnds(site.id);
          return ends.filter(id => starts.has(id)).length === 1;
        }).map(s => s.id);
      }
      const moving = anyonStarPath[anyonStarPath.length - 1];
      const visited = new Set(anyonStarPath);
      return candidates.filter(site => {
        const ends = completeEStarEnds(site.id);
        if (!ends.includes(moving)) return false;
        const other = ends.find(id => id !== moving);
        return other !== undefined && !visited.has(other);
      }).map(s => s.id);
    }

    function completeMTriangleEnds(siteId, color = anyonColor) {
      const site = siteById.get(siteId);
      if (!site || site.colorName !== color) return [];
      return (mTriangleIncidence.get(siteId) || []).filter(
        tId => triangleById.get(tId).color === color
      );
    }

    // Canonical GF(2) reduction of Z strings modulo complete Bt triangles.
    // It preserves the residual vector, unlike the boolean membership test.
    function reduceTriangleMask(mask) {
      let value = mask;
      for (let pivot = maxQubitId; pivot >= 0; pivot--) {
        const bit = qubitBit(pivot);
        if ((value & bit) && triangleZBasis[pivot] !== 0n) {
          value ^= triangleZBasis[pivot];
        }
      }
      return value;
    }

    let reducedZSingle = null;

    // For one single-colour X path, conjugation of an As star produces only
    // diagonal Z residues. All selected X qubits have the same colour, so
    // adjacent Xs never appear on the alternating-colour inner ring of one
    // star and no bare-X sign is generated. Cross-colour CZs do not generate
    // a sign either. We solve for CZs that eliminate all bulk Z residues
    // modulo Bt (not by the erroneous literal requirement Z(mask) === I).
    let mSingleXStarChanges = null;  // filled after triangle basis initialization

    function mStarResiduals(xQubits, czEdges) {
      const res = new Map(starData.map(s => [s.id, 0n]));
      for (const q of xQubits) {
        for (const [id, reduced] of mSingleXStarChanges.get(q) || []) {
          res.set(id, res.get(id) ^ reduced);
        }
      }
      for (const [a,b] of czEdges) {
        for (const id of eStarIncidence.get(a) || []) {
          res.set(id, res.get(id) ^ reducedZSingle.get(b));
        }
        for (const id of eStarIncidence.get(b) || []) {
          res.set(id, res.get(id) ^ reducedZSingle.get(a));
        }
      }
      return res;
    }

    // n=1: two terminal stars; n=2: three terminal stars (1+2).
    // n=3: 1+1, n=4: 1+2, etc. This implements the agreed parity
    // convention: the same CZ decoration is retained for n=2m-1 and 2m.
    function mAllowedStarIds(xQubits, trianglePath) {
      const selected = new Set([trianglePath[0], trianglePath[trianglePath.length-1]]);
      if (xQubits.length > 1 && xQubits.length % 2 === 0) {
        selected.add(trianglePath[trianglePath.length-2]);
      }
      return new Set(Array.from(selected, t => starForTriangle.get(t)));
    }

    function validMResidues(residuals, expectedStarIds) {
      for (const star of starData) {
        const nonzero = residuals.get(star.id) !== 0n;
        if (nonzero !== expectedStarIds.has(star.id)) return false;
      }
      return true;
    }

    function mPathConnectors(xs, color = anyonColor) {
      if (!xs.length || new Set(xs).size !== xs.length) return null;

      if (xs.some(q => siteById.get(q)?.colorName !== color)) {
        return null;
      }

      const connectors = [];

      for (let j = 0; j + 1 < xs.length; j++) {
        const left = mNeighborIds.get(xs[j]);
        const right = mNeighborIds.get(xs[j + 1]);

        if (!left || !right) return null;

        // The connector is directly bonded to BOTH consecutive Xs.
        const shared = [...left].filter(q => right.has(q));

        // Require one unambiguous physical connector.
        if (shared.length !== 1) return null;

        const q = shared[0];

        if (siteById.get(q).colorName === color) return null;

        connectors.push(q);
      }

      // Require consistent alternating connector colors:
      // Green, Red, Green, Red, ... in your example.
      if (connectors.length > 1) {
        const firstColor = siteById.get(connectors[0]).colorName;
        const secondColor = siteById.get(connectors[1]).colorName;

        if (firstColor === secondColor) return null;

        for (let j = 0; j < connectors.length; j++) {
          const expectedColor =
            j % 2 === 0 ? firstColor : secondColor;

          if (siteById.get(connectors[j]).colorName !== expectedColor) {
            return null;
          }
        }
      }

      return connectors;
    }

    function mSupportCZPairs(xs, color = anyonColor) {
      const connectors = mPathConnectors(xs, color);

      if (connectors === null) return null;

      const pairCount = Math.floor(connectors.length / 2);
      const pairs = [];

      // Zero-based indexing:
      // A_i = connectors[2*i]
      // B_i = connectors[2*i + 1]
      //
      // Add CZ(A_i, B_j) for i <= j.
      for (let i = 0; i < pairCount; i++) {
        for (let j = i; j < pairCount; j++) {
          pairs.push([
            connectors[2 * i],
            connectors[2 * j + 1]
          ]);
        }
      }

      return pairs;
    }

    function mNewTrianglePath(newQubit, qs = anyonQubits, path = mTrianglePath, color = anyonColor) {
      const ends = completeMTriangleEnds(newQubit, color);
      if (ends.length !== 2 || qs.includes(newQubit)) return null;
      if (!qs.length) return [...ends];
      if (qs.length === 1) {
        const firstEnds = completeMTriangleEnds(qs[0], color);
        const shared = ends.filter(id => firstEnds.includes(id));
        if (shared.length !== 1) return null;
        const stationary = firstEnds.find(id => id !== shared[0]);
        const moving = ends.find(id => id !== shared[0]);
        return [stationary,shared[0],moving];
      }
      const last = path[path.length-1];
      if (!ends.includes(last)) return null;
      const moving = ends.find(id => id !== last);
      if (path.includes(moving)) return null;
      return [...path,moving];
    }

    function findMDecoration(xs, ts) {
      // Construct the exact CZ list from the selected X path.
      const cz = mSupportCZPairs(xs);

      if (cz === null) return null;

      // Preserve the existing rule against reusing a CZ
      // that belongs to another circuit.
      const own = new Set(
        mAutoCZ.map(pair => pairKey(...pair))
      );

      const occupied = new Set(
        czPairs
          .filter(pair => !own.has(pairKey(...pair)))
          .map(pair => pairKey(...pair))
      );

      if (cz.some(pair => occupied.has(pairKey(...pair)))) {
        return null;
      }

      // Validate the required construction using the existing
      // star-residue check.
      const expected = mAllowedStarIds(xs, ts);

      return validMResidues(
        mStarResiduals(xs, cz),
        expected
      ) ? cz : null;
    }

    function availableMPlans() {
      if (mOptionCache !== null) return mOptionCache;
      const options = new Map();
      for (const site of siteData) {
        if (site.colorName !== anyonColor || localOps.has(site.id)) continue;
        const ts = mNewTrianglePath(site.id);
        if (!ts) continue;
        const qs = [...anyonQubits,site.id];
        const plan = findMDecoration(qs,ts);
        if (plan) options.set(site.id,{triangles:ts,cz:plan});
      }
      mOptionCache=options;
      return options;
    }


    // ===================================================================
    // Two-colour m ribbons. A "rung" is a physical bond joining the selected
    // two colours. Consecutive rungs either cross one remaining-colour
    // hexagon, or the two triangles sharing one remaining-colour site.
    // The former has precisely the user's 111100 / 110110 X support.
    // ===================================================================
    function mPairHexPattern(ring, ids, colors) {
      const chosen = new Set(ids);
      if (ring.length !== 6 || chosen.size !== 4 ||
          [...chosen].some(q => !ring.includes(q))) return null;
      if (colors.some(c => [...chosen].filter(q =>
          siteById.get(q)?.colorName === c).length !== 2)) return null;
      const bits = ring.map(q => chosen.has(q) ? "1" : "0").join("");
      for (let shift = 0; shift < 6; shift++) {
        const rotated = bits.slice(shift) + bits.slice(0, shift);
        if (rotated === "111100") return "four-consecutive";
        if (rotated === "110110") return "opposite-pairs";
      }
      return null;
    }

    function buildMPairBoard(colors) {
      const rungs = new Map(), adjacency = new Map();
      const third = ["Red", "Green", "Blue"].find(c => !colors.includes(c));
      if (colors.length !== 2 || colors[0] === colors[1] || !third) {
        return {rungs, adjacency};
      }
      const hexByEdge = new Map();
      for (const star of starData) {
        if (star.color !== third) continue;
        const ring = star.inner_qubits;
        for (let i = 0; i < 6; i++)
          hexByEdge.set(pairKey(ring[i], ring[(i+1)%6]), star.id);
      }
      const byHex = new Map(), byApex = new Map();
      const group = (map, id, key) => {
        if (!map.has(id)) map.set(id, []);
        map.get(id).push(key);
      };
      for (const edge of edgePairs) {
        const qs = colors.map(c => edge.find(q => siteById.get(q).colorName === c));
        if (qs.some(q => q === undefined)) continue;
        if (qs.some((q, i) => completeMTriangleEnds(q, colors[i]).length !== 2)) continue;
        const key = pairKey(...qs), hex = hexByEdge.get(key);
        if (hex === undefined) continue; // existing complete-open-star rule
        const common = [...mNeighborIds.get(qs[0])].filter(q =>
          mNeighborIds.get(qs[1]).has(q) && siteById.get(q).colorName === third);
        if (common.length !== 1) continue;
        const a = siteById.get(qs[0]), b = siteById.get(qs[1]);
        const rung = {key, qs, hex, apex: common[0], x: (a.x+b.x)/2, y: (a.y+b.y)/2};
        rungs.set(key, rung); adjacency.set(key, new Map());
        group(byHex, hex, key); group(byApex, rung.apex, key);
      }
      function connect(groups, kind) {
        for (const [id, keys] of groups) {
          for (let i = 0; i < keys.length; i++) for (let j = i+1; j < keys.length; j++) {
            const a = rungs.get(keys[i]), b = rungs.get(keys[j]);
            const ids = [...a.qs, ...b.qs];
            if (new Set(ids).size !== 4) continue;
            const pattern = kind === "hex"
              ? mPairHexPattern(starById.get(id).inner_qubits, ids, colors)
              : "shared-apex";
            if (!pattern) continue;
            const step = {kind, id, junction: `${kind}:${id}`, pattern};
            adjacency.get(a.key).set(b.key, step);
            adjacency.get(b.key).set(a.key, step);
          }
        }
      }
      connect(byHex, "hex"); connect(byApex, "apex");
      return {rungs, adjacency};
    }

    function compileMPairRibbon(keys, colors = mPairColors, board = mPairBoard) {
      if (!board || colors.length !== 2 || !keys.length || new Set(keys).size !== keys.length) return null;
      const rungs = keys.map(key => board.rungs.get(key));
      if (rungs.some(r => !r)) return null;
      const allQubits = rungs.flatMap(r => r.qs);
      if (new Set(allQubits).size !== allQubits.length) return null;
      const junctions = new Set(), steps = [];
      for (let i = 1; i < rungs.length; i++) {
        const step = board.adjacency.get(keys[i-1]).get(keys[i]);
        // No revisiting a face/apex: otherwise a hexagon could acquire six Xs.
        if (!step || junctions.has(step.junction)) return null;
        junctions.add(step.junction); steps.push(step);
      }
      const factors = [];
      for (let side = 0; side < 2; side++) {
        const color = colors[side], xs = rungs.map(r => r.qs[side]);
        let triangles = [];
        const prefix = [];
        for (const q of xs) {
          const next = mNewTrianglePath(q, prefix, triangles, color);
          if (!next) return null;
          triangles = next; prefix.push(q);
        }
        const cz = mSupportCZPairs(xs, color);
        if (cz === null || !validMResidues(mStarResiduals(xs, cz),
            mAllowedStarIds(xs, triangles))) return null;
        factors.push({color, xs, triangles, cz});
      }
      
      // This is chronological order, not a commuting/canonicalized gate set.
      // Complete second factor, THEN complete first factor. For two colours c1 c2, U_c1c2 = X_c1 CZ_c1 X_c2 CZ_c2
      const events = [];
      for (const factor of [...factors].reverse()) {
        for (const pair of factor.cz)
          events.push({type: "CZ", qubits: [...pair], factorColor: factor.color});
        for (const q of factor.xs)
          events.push({type: "X", qubits: [q], factorColor: factor.color});
      }
      return {keys: [...keys], factors, events, steps, allQubits};
    }

    function availableMPairPlans() {
      if (mPairOptionCache !== null) return mPairOptionCache;
      const options = new Map();
      if (!mPairBoard) return options;
      const candidateKeys = mPairEdges.length
        ? [...mPairBoard.adjacency.get(mPairEdges.at(-1)).keys()]
        : [...mPairBoard.rungs.keys()];
      for (const key of candidateKeys) {
        const rung = mPairBoard.rungs.get(key);
        if (rung.qs.some(q => localOps.has(q) || anyonQubits.includes(q))) continue;
        const plan = compileMPairRibbon([...mPairEdges, key]);
        if (plan) options.set(key, plan);
      }
      mPairOptionCache = options;
      return options;
    }

    function handleMPairBond(key) {
      const plan = availableMPairPlans().get(key);
      if (!plan) {setStatus(`${anyonName()}: choose a highlighted bond midpoint.`); return false;}
      pushHistory(`extending ${anyonName()} to paired bond ${key}`);
      const affected = new Set([...anyonQubits, ...plan.allQubits]);
      // Replace only THIS composite session. Keep other strings' exact gates,
      // including repeated CZs and their original chronological positions.
      operationLog = operationLog.filter(ev => !(ev.source === "guided-m2" && ev.guideId === mPairGuideId));
      for (const ev of plan.events)
        recordOperation(ev.type, ev.qubits, "guided-m2", mPairGuideId);
      czPairs = operationLog.filter(ev => ev.type === "CZ").map(ev => [...ev.qubits]);
      for (const q of affected) syncVisibleLocalOp(q);
      mPairEdges = [...plan.keys];
      anyonQubits = [...plan.allQubits];
      mPairOptionCache = null;
      render();
      const n = mPairEdges.length, next = availableMPairPlans().size;
      const order = plan.factors.map(f => `${f.color}: ${f.cz.length} CZ → ${f.xs.length} X`).join("; then ");
      setStatus(`${anyonName()}: ${n} paired steps. ${order}. ${next} next bond(s). ` +
        (n % 2 === 0 ? "Even prefix: each factor keeps its preceding odd-prefix CZ decoration." :
          "Odd-length factors complete. Use the violation tools to inspect the ordered product."));
      return true;
    }

    function buildMPairGuideTraces() {
      if (anyonKind !== "m2" || !mPairBoard) return [];
      const candidates = [...availableMPairPlans().keys()].map(key => mPairBoard.rungs.get(key));
      const chosen = mPairEdges.map(key => mPairBoard.rungs.get(key));
      function bonds(rungs, name, color, width) {
        const x = [], y = [];
        for (const rung of rungs) {
          for (const q of rung.qs) {const p = siteById.get(q); x.push(p.x); y.push(p.y);}
          x.push(null); y.push(null);
        }
        return {name, type: "scatter", mode: "lines", x, y,
          line: {color, width}, hoverinfo: "skip", showlegend: false};
      }
      return [
        bonds(chosen, "m2-chosen-bonds", "rgba(13,148,136,0.55)", 3),
        bonds(candidates, "m2-candidate-bonds", "rgba(217,119,6,0.65)", 3),
        {name: "m2-options", type: "scatter", mode: "markers",
          x: candidates.map(r => r.x), y: candidates.map(r => r.y),
          customdata: candidates.map(r => r.key),
          text: candidates.map(r => `${mPairColors[0]} q${r.qs[0]} + ${mPairColors[1]} q${r.qs[1]}`),
          marker: {symbol: "diamond", size: 11, color: "#fbbf24", line: {color: "#92400e", width: 1}},
          hovertemplate: "<b>%{text}</b><br>Click to place both Xs<extra></extra>", showlegend: false}
      ];
    }

    // ------------------------------------------------------------------
    // Three-colour straight ribbons (s1 / s2).
    // Accept horizontal and both 60-degree diagonal directions.
    // IMPORTANT: UP / DOWN is defined relative to the direction of travel
    // along the SELECTED EDGE, not relative to the triangle's horizontal base.
    // A triangle whose third vertex lies LEFT of the directed path is UP;
    // RIGHT of the path is DOWN.  With s1, UP => XR on the first selected
    // qubit, RX on the next for the OLD (reversed) labels.  In the user's
    // ORIGINAL naming, s1 instead uses UP => RX first, XR next and
    // DOWN => XR first, RX next; s2 is the complex-conjugate mirror.
    // Reversing travel reverses both orientation and endpoint order, so the
    // final assignment to each PHYSICAL qubit is unchanged.
    // The label correction changes the former "s2" operator into the
    // user-named "s1" (and vice versa) without changing ribbon topology.
    // Reversing travel reverses endpoint order AND orientation, so physical
    // XR/RX assignments remain unchanged on a reversed directed path.
    // The third triangle vertex receives X only.  All CZ gates remain a
    // jointly optimized three-colour decoration applied first.
    // ------------------------------------------------------------------
    const sStraightNeighbors = new Map();
    const sSiteNeighbors = new Map(siteData.map(s => [s.id, []]));
    for (const [a,b] of edgePairs) {
      sSiteNeighbors.get(a).push(b);
      sSiteNeighbors.get(b).push(a);
    }

    // The 2D cross product is positive when the triangle's third vertex is
    // to the LEFT of the directed edge a -> b, negative on its RIGHT.
    // This is the ribbon-relative UP/DOWN convention, valid for all slopes.
    function sDirectedTriangleUp(a, b, apex) {
      const p=siteById.get(a),q=siteById.get(b),t=siteById.get(apex);
      const cross=(q.x-p.x)*(t.y-p.y)-(q.y-p.y)*(t.x-p.x);
      if (Math.abs(cross)<1e-7) return null;
      return cross>0;
    }
    for (const [a,b] of edgePairs) {
      const sa=siteById.get(a),sb=siteById.get(b);
      // Every Kagome nearest-neighbor bond has a nonzero x displacement,
      // so 'left/right' is well-defined also for the two diagonal slopes.
      if (Math.abs(sa.x-sb.x)<1e-5) continue;
      const common=sSiteNeighbors.get(a).filter(id=>sSiteNeighbors.get(b).includes(id));
      if (common.length!==1) continue;
      const apex=common[0];
      if (![a,b,apex].every(id=>(mTriangleIncidence.get(id)||[]).length===2)) continue;
      if (sDirectedTriangleUp(a,b,apex)===null) continue;
      const l=sa.x<sb.x?a:b,r=sa.x<sb.x?b:a;
      sStraightNeighbors.set(pairKey(a,b),{left:l,right:r,apex});
    }

    function sStepVector(a,b) {
      const u=siteById.get(a),v=siteById.get(b);
      return [v.x-u.x,v.y-u.y];
    }
    function sSameVector(a,b) {
      const u=sStepVector(a[0],a[1]),v=sStepVector(b[0],b[1]);
      return Math.abs(u[0]-v[0])<1e-5 && Math.abs(u[1]-v[1])<1e-5;
    }
    // Guided semion ribbons may overlap physical X, XR, or RX gates left by
    // earlier s1/s2 paths.  Gate history is cumulative: the new guide adds a
    // fresh R/X/CZ circuit, instead of replacing or toggling the old one.
    // Z-occupied and R-only sites are not eligible for a new semion spine.
    function canPlaceSOnQubit(id) {
      if (rCompiling) return true; // independent microscopic-R workspace
      const localGates = operationLog.filter(ev=>
        ev.qubits.length===1 && ev.qubits[0]===id
      ).map(ev=>ev.type);
      if (localGates.includes("Z")) return false;
      if (localGates.includes("R") && !localGates.includes("X")) return false;
      return !localOps.has(id) || localOps.get(id)==="X";
    }

    function availableSAnyonQubits() {
      if (!sBasePath.length) {
        const all=new Set();
        for (const edge of sStraightNeighbors.values()) {all.add(edge.left);all.add(edge.right);}
        return [...all].filter(id=>canPlaceSOnQubit(id));
      }
      const last=sBasePath[sBasePath.length-1];
      return (sSiteNeighbors.get(last)||[]).filter(id=>{
        if (sBasePath.includes(id) || !canPlaceSOnQubit(id)) return false;
        if (!sStraightNeighbors.has(pairKey(last,id))) return false;
        // After the second click, there can be no bending or reversing.
        if (sBasePath.length>=2 && !sSameVector(
          [sBasePath[sBasePath.length-2],last],[last,id]
        )) return false;
        return true;
      });
    }

    function buildSSpine(base) {
      if (base.length<2) return {xIds:[],pairs:[],gateOrder:new Map()};
      const xIds=new Set(base),gateOrder=new Map(),pairs=[];
      for (let i=0;i<base.length-1;i++) {
        if (i>0 && !sSameVector([base[i-1],base[i]],[base[i],base[i+1]])) return null;
        const tri=sStraightNeighbors.get(pairKey(base[i],base[i+1]));
        if (!tri) return null;
        if (!canPlaceSOnQubit(tri.apex)) return null;
        xIds.add(tri.apex);pairs.push(tri);
        // Determine UP/DOWN relative to the user-selected travel direction.
        // This corrects the previous diagonal mapping, which used the
        // physical triangle's horizontal-base orientation instead.
        const orientedUp=sDirectedTriangleUp(base[i],base[i+1],tri.apex);
        if (orientedUp===null) return null;
        // User's original String 1 convention: for a directed edge with
        // triangle apex on its left (UP), the first qubit has RX and the
        // second has XR. String 2 is its conjugate pattern on that SAME edge.
        // DOWN swaps RX and XR. This is the sole generator used by both
        // the main editor and the microscopic R workbench.
        const string1First = orientedUp ? 'RX' : 'XR';
        const firstOrder = anyonKind==='s1' ? string1First
          : (string1First==='RX' ? 'XR' : 'RX');
        const nextOrder=firstOrder==='XR'?'RX':'XR';
        for (const [q,order] of [[base[i],firstOrder],[base[i+1],nextOrder]]) {
          if (gateOrder.has(q)&&gateOrder.get(q)!==order) return null;
          gateOrder.set(q,order);
        }
      }
      return {xIds:[...xIds],pairs,gateOrder};
    }

    // Joint three-colour CZ search for s1/s2.
    // Every CZ joins TWO X sites selected on the physical ribbon, with no
    // distance restriction.  This differs from merely concatenating three
    // separate single-colour m-string decorations.
    //
    // Exact criterion: for each non-endpoint complete star s, the ordered
    // circuit U must obey U† A_s U = + A_s Π_t B_t.  Residual Z strings are
    // reduced modulo the full complete-Bt span, and the overall sign is
    // included in the GF(2) equations.  The original XR/RX phase pattern is
    // held fixed.  A star touching an odd number of R sites is classified
    // as an endpoint-region star and is not required to be +1.
    function compileSDecorations(built) {
      const xs=built.xIds;
      if (!xs.length) return [];

      // With just one triangular segment (two base sites), the string is
      // an unresolved short prefix.  Not all non-endpoint constraints can
      // yet be satisfied by CZs whose endpoints lie on this tiny ribbon.
      if (built.gateOrder.size < 3) return [];

      const localEvents=[];
      for (const q of xs) {
        const order=built.gateOrder.get(q);
        if (order==="XR") {
          localEvents.push({type:"R",q});
          localEvents.push({type:"X",q});
        } else if (order==="RX") {
          localEvents.push({type:"X",q});
          localEvents.push({type:"R",q});
        } else {
          localEvents.push({type:"X",q});
        }
      }
      const phaseSites=new Set(built.gateOrder.keys());
      const bulkStars=[];
      for (const star of starData) {
        const rotations=star.x_qubits.filter(q=>phaseSites.has(q)).length;
        if ((rotations%2)===0) bulkStars.push(star);
      }
      // Each star occupies its own disjoint bit block.  One extra bit in
      // that block records the sign (-1) of the scalar i^phase.
      const width=maxQubitId+2;
      const offsets=new Map(bulkStars.map((s,i)=>[s.id, BigInt(i*width)]));
      const signOffset=BigInt(maxQubitId+1);
      let target=0n;
      for (const star of bulkStars) {
        const ring=star.inner_qubits;
        const ringIndex=new Map(ring.map((q,i)=>[q,i]));
        const tips=new Set(star.x_qubits);
        let zMask=0n, phase=0;
        for (let j=localEvents.length-1;j>=0;j--) {
          const ev=localEvents[j], q=ev.q;
          if (ev.type==="R" && tips.has(q)) {
            zMask ^= qubitBit(q);
            phase=(phase+1)%4;
          } else if (ev.type==="X") {
            if ((zMask & qubitBit(q))!==0n) phase=(phase+2)%4;
            if (ringIndex.has(q)) {
              const i=ringIndex.get(q);
              zMask ^= qubitBit(ring[(i+ring.length-1)%ring.length]);
              zMask ^= qubitBit(ring[(i+1)%ring.length]);
            }
          }
        }
        // For an even number of R rotations the total phase must be real;
        // CZ conjugations can only alter it by two powers of i.
        if (phase%2) return null;
        const shift=offsets.get(star.id);
        target ^= reduceTriangleMask(zMask) << shift;
        if (phase===2) target ^= 1n<<(shift+signOffset);
      }
      if (target===0n) return [];

      const columns=[];
      const ordered=[...xs].sort((a,b)=>a-b);
      for (let i=0;i<ordered.length;i++) {
        const a=ordered[i];
        for (let j=i+1;j<ordered.length;j++) {
          const b=ordered[j];
          let column=0n;
          const intersectA=new Set(eStarIncidence.get(a)||[]);
          const intersectB=new Set(eStarIncidence.get(b)||[]);
          for (const starId of intersectA) {
            const shift=offsets.get(starId);
            if (shift!==undefined) {
              column ^= reducedZSingle.get(b)<<shift;
              if (intersectB.has(starId)) column ^= 1n<<(shift+signOffset);
            }
          }
          for (const starId of intersectB) {
            const shift=offsets.get(starId);
            if (shift!==undefined) column ^= reducedZSingle.get(a)<<shift;
          }
          if (column!==0n) {
            const sa=siteById.get(a),sb=siteById.get(b);
            const dist=Math.hypot(sa.x-sb.x,sa.y-sb.y);
            columns.push({pair:[a,b],column,dist});
          }
        }
      }

      // Solve the GF(2) equation columns * c = target.  Each selected
      // variable indicates one CZ.  All phase and Bt-reduction constraints
      // are enforced simultaneously for every complete bulk star.
      const pivot=(value)=>value.toString(2).length-1;
      const basis=new Map(), nullVectors=[];
      for (let i=0;i<columns.length;i++) {
        let value=columns[i].column,combo=1n<<BigInt(i);
        while (value!==0n) {
          const p=pivot(value),row=basis.get(p);
          if (row) {value ^= row.vec;combo ^= row.combo;}
          else {basis.set(p,{vec:value,combo});break;}
        }
        if (value===0n) nullVectors.push(combo);
      }
      let value=target,selection=0n;
      while (value!==0n) {
        const row=basis.get(pivot(value));
        if (!row) return null;
        value ^= row.vec;
        selection ^= row.combo;
      }

      // Minimize the CZ count within the nullspace.  For short paths we can
      // enumerate every free-variable choice exactly.  For long ribbons use
      // a bounded deterministic greedy reduction to keep the UI responsive.
      function score(mask) {
        let n=0, length=0;
        for (let i=0;i<columns.length;i++) {
          if (mask & (1n<<BigInt(i))) {n++;length+=columns[i].dist;}
        }
        return 10000*n+length;
      }
      let best=selection,bestScore=score(best);
      if (nullVectors.length<=11) {
        for (let pattern=1;pattern<(1<<nullVectors.length);pattern++) {
          let trial=selection;
          for (let j=0;j<nullVectors.length;j++) {
            if (pattern & (1<<j)) trial ^= nullVectors[j];
          }
          const s=score(trial);
          if (s<bestScore) {best=trial;bestScore=s;}
        }
      } else {
        for (let pass=0;pass<3;pass++) {
          let improved=false;
          for (const free of nullVectors) {
            const trial=best^free;
            const s=score(trial);
            if (s<bestScore) {best=trial;bestScore=s;improved=true;}
          }
          if (!improved)break;
        }
      }
      const result=[];
      for (let i=0;i<columns.length;i++) {
        if (best & (1n<<BigInt(i))) result.push(columns[i].pair);
      }
      return result;
    }

    function handleSAnyonQubit(siteId) {
      const nextBase=[...sBasePath,siteId];
      const built=buildSSpine(nextBase);
      if(!built){setStatus(`${anyonName()}: this extension cannot form a straight, consistently oriented alternating triangle ribbon.`);return false;}
      // The first click sets an anchor; operations begin when the first
      // triangular segment and its orientation have been chosen.
      const cz=built.xIds.length?compileSDecorations(built):[];
      if(cz===null){setStatus(`${anyonName()}: this extension cannot satisfy the joint bulk-star constraints with CZs restricted to ribbon qubits. Choose a different start or use Undo.`);return false;}
      pushHistory(nextBase.length===1 ? `choosing ${anyonName()} start ${siteId}` : `extending ${anyonName()} to ${siteId}`);
      // Remove only the *current guide session's* previous compilation.
      // Never remove a CZ just because another ribbon uses the same pair:
      // overlapping ribbons must retain their independent physical gates.
      operationLog=operationLog.filter(ev=>!(ev.source==="guided-s"&&ev.guideId===sGuideId));
      czPairs=operationLog.filter(ev=>ev.type==="CZ").map(ev=>[...ev.qubits]);
      sBasePath=nextBase;
      anyonQubits=[...nextBase];
      sGuidedX=built.xIds;
      sAutoCZ=cz;
      for (const [a,b] of cz) {
        czPairs.push([a,b]);recordOperation("CZ",[a,b],"guided-s",sGuideId);
      }
      // Chronology: existing ribbons, then this guide's CZ first, followed
      // by its R/X operations.  A repeated X, R, or CZ remains in the log
      // (with its own ID) because it can change the physical combined state.
      // For XR apply R then X; for RX apply X then R.
      for (const q of sGuidedX) {
        const order=built.gateOrder.get(q);
        if (order==="XR") {
          recordOperation("R",[q],"guided-s",sGuideId);
          recordOperation("X",[q],"guided-s",sGuideId);
        } else if (order==="RX") {
          recordOperation("X",[q],"guided-s",sGuideId);
          recordOperation("R",[q],"guided-s",sGuideId);
        } else {
          recordOperation("X",[q],"guided-s",sGuideId);
        }
        localOps.set(q,"X");
      }
      sCandidateCache=null;
      return true;
    }

    function availableAnyonQubits() {
      if (!anyonKind) return [];
      if (anyonKind === "m2") return [...new Set([...availableMPairPlans().keys()]
        .flatMap(key => mPairBoard.rungs.get(key).qs))];
      if (anyonKind === "s1" || anyonKind === "s2") return availableSAnyonQubits();
      if (!anyonColor) return [];
      return anyonKind === "m"
        ? [...availableMPlans().keys()]
        : availableEAnyonQubits();
    }

    function clearAnyonGuidance() {
      anyonKind=null;anyonColor=null;anyonQubits=[];
      anyonStarPath=[];mTrianglePath=[];mAutoCZ=[];
      sBasePath=[];sGuidedX=[];sAutoCZ=[];sGuideId=null;sCandidateCache=null;
      mOptionCache=null;mGuideId=null;
      mPairColors=[];mPairEdges=[];mPairGuideId=null;mPairBoard=null;mPairOptionCache=null;
      anyonSelect.value="";
      updateAnyonSummary();
    }

    function setAnyon(value) {
      const match = /^(e|m)_([RGB])$/.exec(value);
      const pairMatch = /^m_([RGB])([RGB])$/.exec(value);
      const validPair = pairMatch && pairMatch[1] !== pairMatch[2];
      const names = {R:"Red",G:"Green",B:"Blue"};
      anyonKind = validPair ? "m2" : ((value === "s1" || value === "s2") ? value : (match ? match[1] : null));
      anyonColor=match ? names[match[2]]:null;
      mPairColors=validPair ? [names[pairMatch[1]],names[pairMatch[2]]] : [];
      mPairEdges=[];mPairOptionCache=null;
      mPairGuideId=validPair ? nextMPairGuideId++ : null;
      mPairBoard=validPair ? buildMPairBoard(mPairColors) : null;
      anyonQubits=[];anyonStarPath=[];
      mTrianglePath=[];mAutoCZ=[];mOptionCache=null;
      sBasePath=[];sGuidedX=[];sAutoCZ=[];sCandidateCache=null;
      sGuideId = (anyonKind === "s1" || anyonKind === "s2") ? nextSGuideId++ : null;
      mGuideId = anyonKind === "m" ? nextMGuideId++ : null;
      pendingCZ=null;
      if (zoomMode) setZoomMode(false);
      currentMode=anyonKind ? "ANYON":"Z";
      updateActiveModeButtons();
      updateAnyonSummary();
      render();
      const options=availableAnyonQubits();
      if (anyonKind==="m2") {
        setStatus(options.length ? `${anyonName()}: click a gold diamond between paired qubits. Choose the shared ribbon, not two independent paths.`
          : `${anyonName()}: no complete interior paired bonds are available. Increase the lattice size or clear conflicting X/Z operations.`);
      } else if (anyonKind==="s1" || anyonKind==="s2") {
        setStatus(options.length ? `${anyonName()}: select a highlighted start, including qubits on earlier X/XR/RX strings. Extend horizontally or diagonally; overlapping gates accumulate. XR/RX, tip X, and joint CZ dressing are automatic.`
          : `${anyonName()}: no interior straight Kagome line is available.`);
      } else if (anyonKind==="m") {
        setStatus(options.length
          ? `${anyonName()}: click a highlighted ${anyonColor.toLowerCase()} qubit to start the X-string. The app automatically builds the CZ-first circuit.`
          : `${anyonName()}: no available start qubits with two complete triangles and stars.`);
      } else if (anyonKind==="e") {
        setStatus(options.length
          ? `${anyonName()}: click any highlighted ${anyonColor.toLowerCase()} qubit to place the starting Z. Only complete interior stars are used.`
          : `${anyonName()}: no available start qubits connecting two complete ${anyonColor.toLowerCase()} stars.`);
      } else {
        setStatus("Anyon guide off. Mode: Z.");
      }
    }

    function handleEAnyonQubit(siteId) {
      const newEnds=completeEStarEnds(siteId);
      const isStart=anyonQubits.length===0;
      const oldTail=anyonQubits[anyonQubits.length-1];
      pushHistory(isStart
        ? `starting ${anyonName()} at qubit ${siteId}`
        : `extending ${anyonName()} to qubit ${siteId}`);
      if (anyonQubits.length===1) {
        const prevEnds=completeEStarEnds(oldTail);
        const shared=prevEnds.find(id=>newEnds.includes(id));
        anyonStarPath=[prevEnds.find(id=>id!==shared),shared,newEnds.find(id=>id!==shared)];
      } else if (anyonQubits.length>=2) {
        const moving=anyonStarPath[anyonStarPath.length-1];
        anyonStarPath.push(newEnds.find(id=>id!==moving));
      }
      anyonQubits.push(siteId);
      localOps.set(siteId,"Z");
      recordOperation("Z",[siteId],"guided-e");
    }

    function handleMAnyonQubit(siteId) {
      const plan=availableMPlans().get(siteId);
      if (!plan) return;
      pushHistory(anyonQubits.length===0
        ? `starting ${anyonName()} at qubit ${siteId}`
        : `extending ${anyonName()} to qubit ${siteId}`);

      // Replace just the current guided-m block; previous manual/e gates
      // are untouched. A snapshot before each click makes Undo atomic.
      const oldPairs = new Set(mAutoCZ.map(pair => pairKey(...pair)));
      czPairs=czPairs.filter(pair=>!oldPairs.has(pairKey(...pair)));
      operationLog=operationLog.filter(ev=>!(ev.source==="guided-m" && ev.guideId===mGuideId));
      for (const [a,b] of plan.cz) {
        czPairs.push([a,b]);
        recordOperation("CZ",[a,b],"guided-m",mGuideId);
      }
      anyonQubits.push(siteId);
      mTrianglePath=[...plan.triangles];
      mAutoCZ=plan.cz.map(pair=>[...pair]);
      // The X operations are logged AFTER the CZs (chronological order).
      // Accordingly U = X_n ... X_1 CZ_k ... CZ_1 as requested.
      for (const q of anyonQubits) {
        localOps.set(q,"X");
        recordOperation("X",[q],"guided-m",mGuideId);
      }
      mOptionCache=null;
    }

    function handleAnyonQubit(siteId) {
      if (anyonKind === "m2") {
        const choices = [...availableMPairPlans().keys()].filter(key => mPairBoard.rungs.get(key).qs.includes(siteId));
        if (choices.length === 1) handleMPairBond(choices[0]);
        else setStatus(`${anyonName()}: click a gold bond-midpoint diamond to choose both qubits unambiguously.`);
        return;
      }
      const options=new Set(availableAnyonQubits());
      if (!options.has(siteId)) {
        setStatus(`${anyonName()}: qubit ${siteId} is not a valid next step. Choose a highlighted qubit.`);
        return;
      }
      if (anyonKind==="m") handleMAnyonQubit(siteId);
      else if (anyonKind==="s1" || anyonKind==="s2") {
        if (!handleSAnyonQubit(siteId)) return;
      } else handleEAnyonQubit(siteId);
      updateAnyonSummary();
      render();
      const next=availableAnyonQubits();
      const pathText=anyonQubits.join(" → ");
      if (anyonKind==="s1" || anyonKind==="s2") {
        setStatus(`${anyonName()} (${sDirectionLabel()}): ${sBasePath.length} straight-line points, ${sGuidedX.length} X + ${(sBasePath.length>1?sBasePath.length:0)} R + ${sAutoCZ.length} automatic CZ. ${next.length} next step(s). Joint three-color CZ decoration checked against all complete bulk stars; use Violation to inspect the endpoint state.`);
      } else if (anyonKind==="m") {
        const n=anyonQubits.length;
        const starCount= n%2===0 && n>1 ? 3:2;
        setStatus(`${anyonName()} X-string: ${pathText}. ${n} X + ${mAutoCZ.length} automatic CZ. Circuit is rebuilt CZ first, then X (written X…X CZ…CZ). Expected ${starCount} endpoint-region star violations. ${next.length} valid next X qubit${next.length===1?"":"s"}.`);
      } else {
        setStatus(`${anyonName()} Z-string: ${pathText}. ${next.length} highlighted next qubit${next.length===1?"":"s"}. ${next.length?"Click one to extend the path.":"No further interior extensions; use Undo, Clear, or choose an anyon to restart."}`);
      }
    }

    function updateActiveModeButtons() {
      for (const btn of modeButtons) {
        const mode = btn.dataset.mode;
        if (mode === currentMode) {
          btn.classList.add("active");
        } else {
          btn.classList.remove("active");
        }
      }
    }
    
    function updateVisibilityButtons() {
      for (const btn of visibilityButtons) {
        const key = btn.dataset.visibility;
        const isVisible = operatorVisibility[key];

        if (isVisible) {
          btn.classList.add("active");
          btn.setAttribute("aria-pressed", "true");
        } else {
          btn.classList.remove("active");
          btn.setAttribute("aria-pressed", "false");
        }
      }
    }


    // ================================================================
    // Local D4 Hamiltonian checks
    //
    // Bt = Z⊗Z⊗Z.
    //
    // As = (Π CZ around the six INNER hexagon vertices)
    //      × (X on the six OUTWARD same-colour tips).
    //
    // The reference ground state obeys As = Bt = +1.  For a conjugated
    // star, U† As U can always be written as
    //
    //      sign × As × Z(mask)
    //
    // for the allowed X/Z/CZ gates.  Its expectation is:
    //      sign, if Z(mask) is generated by the complete Bt constraints;
    //      0,    otherwise.
    //
    // This gives the characteristic As = 0 after a local X excitation,
    // while a local Z on a star qubit gives As = -1.
    // ================================================================

    const maxQubitId = Math.max(...siteData.map(site => site.id));

    function qubitBit(siteId) {
      return 1n << BigInt(siteId);
    }

    function qubitMask(qubits) {
      let mask = 0n;

      for (const siteId of qubits) {
        mask ^= qubitBit(siteId);
      }

      return mask;
    }

    function buildTriangleZBasis() {
      const basis = Array(maxQubitId + 1).fill(0n);

      for (const triangle of triangleData) {
        let value = qubitMask(triangle.qubits);

        for (let pivot = maxQubitId; pivot >= 0; pivot--) {
          const bit = qubitBit(pivot);

          if ((value & bit) === 0n) {
            continue;
          }

          if (basis[pivot] !== 0n) {
            value ^= basis[pivot];
          } else {
            basis[pivot] = value;
            break;
          }
        }
      }

      return basis;
    }

    const triangleZBasis = buildTriangleZBasis();
    reducedZSingle = new Map(
      siteData.map(s => [s.id, reduceTriangleMask(qubitBit(s.id))])
    );
    mSingleXStarChanges = new Map(siteData.map(s => [s.id, []]));
    for (const star of starData) {
      const ring = star.inner_qubits;
      for (let i = 0; i < ring.length; i++) {
        const mask = qubitBit(ring[(i+5)%6]) ^ qubitBit(ring[(i+1)%6]);
        mSingleXStarChanges.get(ring[i]).push([star.id,reduceTriangleMask(mask)]);
      }
    }

    function isTriangleStabilizerMask(mask) {
      let value = mask;

      for (let pivot = maxQubitId; pivot >= 0; pivot--) {
        const bit = qubitBit(pivot);

        if ((value & bit) === 0n) {
          continue;
        }

        if (triangleZBasis[pivot] !== 0n) {
          value ^= triangleZBasis[pivot];
        } else {
          return false;
        }
      }

      return value === 0n;
    }

    // Trace the same conjugation rules used by starExpectation(). The additional
    // certificate tracks which complete Bt triangles were used for the Z reduction.
    function maskSiteIds(mask) {
      return siteData.filter(site => (mask & qubitBit(site.id)) !== 0n).map(site => site.id);
    }

    function zExpression(mask) {
      const ids = maskSiteIds(mask);
      return ids.length ? ids.map(q => `Z<sub>${q}</sub>`).join(" · ") : "I";
    }

    function triangleDecomposition(mask) {
      const rows = Array(maxQubitId + 1).fill(0n);
      const combos = Array(maxQubitId + 1).fill(0n);
      for (let i=0; i<triangleData.length; i++) {
        let v = qubitMask(triangleData[i].qubits);
        let combination = 1n << BigInt(i);
        for (let pivot=maxQubitId; pivot>=0; pivot--) {
          const bit=qubitBit(pivot);
          if ((v & bit) === 0n) continue;
          if (rows[pivot]) {
            v ^= rows[pivot];
            combination ^= combos[pivot];
          } else {
            rows[pivot]=v;
            combos[pivot]=combination;
            break;
          }
        }
      }
      let remainder=mask, combination=0n;
      for (let pivot=maxQubitId; pivot>=0; pivot--) {
        const bit=qubitBit(pivot);
        if ((remainder & bit) && rows[pivot]) {
          remainder ^= rows[pivot];
          combination ^= combos[pivot];
        }
      }
      const indices=[];
      for (let i=0; i<triangleData.length; i++) {
        if (combination & (1n<<BigInt(i))) indices.push(triangleData[i].id);
      }
      return {remainder, triangles:indices};
    }

    function formatPhase(p) {
      return ["+1", "+i", "−1", "−i"][((p%4)+4)%4];
    }

    function traceStarExpectation(star) {
      const inner=star.inner_qubits, xTips=star.x_qubits;
      const innerSet=new Set(inner), tipSet=new Set(xTips);
      const ringIndex=new Map(inner.map((q,i)=>[q,i]));
      const steps=[];
      let phasePower=0, mask=0n, unaffected=0;
      for (let index=operationLog.length-1; index>=0; index--) {
        const ev=operationLog[index];
        const beforeMask=mask, beforePhase=phasePower;
        let explanation="";
        if (ev.qubits.length===1) {
          const q=ev.qubits[0];
          if (ev.type==="Z" && tipSet.has(q)) {
            phasePower+=2;
            explanation=`Z<sub>${q}</sub> anticommutes with the outer X<sub>${q}</sub> → phase × (−1).`;
          } else if (ev.type==="X") {
            if (mask & qubitBit(q)) {
              phasePower+=2;
              explanation=`X<sub>${q}</sub> flips the sign of a residual Z<sub>${q}</sub>. `;
            }
            if (innerSet.has(q)) {
              const i=ringIndex.get(q), a=inner[(i+5)%6], b=inner[(i+1)%6];
              mask ^= qubitBit(a)^qubitBit(b);
              explanation+=`X<sub>${q}</sub> conjugates adjacent ring CZ gates → Z<sub>${a}</sub>Z<sub>${b}</sub>.`;
            }
          } else if (ev.type==="R" && tipSet.has(q)) {
            phasePower+=1;
            mask ^= qubitBit(q);
            explanation=`R<sub>${q}</sub>† X<sub>${q}</sub> R<sub>${q}</sub> = iX<sub>${q}</sub>Z<sub>${q}</sub>.`;
          }
        } else if (ev.type==="CZ" && ev.qubits.length===2) {
          const [a,b]=ev.qubits;
          const aTip=tipSet.has(a),bTip=tipSet.has(b);
          if (aTip && bTip) phasePower+=2;
          if (aTip) mask ^= qubitBit(b);
          if (bTip) mask ^= qubitBit(a);
          if (aTip || bTip) {
            explanation=`CZ<sub>${a},${b}</sub> adds ${aTip?`Z<sub>${b}</sub>`:""}${aTip&&bTip?" · ":""}${bTip?`Z<sub>${a}</sub>`:""}${aTip&&bTip?" and phase (−1)":""}.`;
          }
        }
        if (beforeMask!==mask || beforePhase!==phasePower) {
          steps.push({number:index+1, operator:ev.type, qubits:ev.qubits, explanation,
            phase:formatPhase(phasePower), mask:zExpression(mask)});
        } else { unaffected++; }
      }
      const decomposition=triangleDecomposition(mask);
      const inSpan=decomposition.remainder===0n;
      const p=((phasePower%4)+4)%4;
      const value=!inSpan ? 0 : (p===0?1:(p===2?-1:0));
      // If a residual is not generated by complete Bt's, a star with odd
      // overlap gives an explicit orthogonality/zero-expectation witness.
      let witness=null;
      if (!inSpan) {
        witness=starData.find(other=>
          other.x_qubits.reduce((k,q)=>k+Number((decomposition.remainder&qubitBit(q))!==0n),0)%2===1
        ) || null;
      }
      return {value, mask, phasePower:p, phase:formatPhase(p), steps, unaffected,
        decomposition, inSpan, witness};
    }

    function traceTriangleExpectation(triangle) {
      const support=new Set(triangle.qubits), hits=[];
      for (let i=0; i<operationLog.length; i++) {
        const ev=operationLog[i];
        if (ev.type==="X" && ev.qubits.length===1 && support.has(ev.qubits[0])) {
          hits.push({number:i+1,q:ev.qubits[0]});
        }
      }
      return {hits, value:hits.length%2?-1:1, ignored:operationLog.length-hits.length};
    }

    function expectationResultHTML(value) {
      return `<p class="expectation-result ${value===1?"good":"bad"}">⟨${expectationMode==="STAR"?"Aₛ":"Bₜ"}⟩ = ${value===1?"+1":(value===-1?"−1":"0")}</p>`;
    }

    function updateExpectationPanel() {
      if (!inspectedObject) {
        expectationContent.innerHTML=expectationMode
          ? `<span class="expectation-intro">Click inside a complete ${expectationMode==="STAR"?"hexagonal star":"triangle"} on the plot to see the calculation. You can click another to compare.</span>`
          : `<span class="expectation-intro">Press <b>Star exp.</b> or <b>Triangle exp.</b>, then click a complete star or triangle.</span>`;
        return;
      }
      const isStar=inspectedObject.kind==="STAR";
      const obj=isStar?starById.get(inspectedObject.id):triangleById.get(inspectedObject.id);
      if (!obj) {
        inspectedObject=null;
        expectationContent.textContent="That star/triangle is outside the complete open-boundary set.";
        return;
      }
      if (isStar) {
        const calc=traceStarExpectation(obj);
        const triangleNames=calc.decomposition.triangles.map(id=>`B<sub>${id}</sub>`);
        const steps=calc.steps.map(st=>`<li><strong>#${st.number} ${st.operator}(${st.qubits.join(", ")})</strong>: ${st.explanation}<div class="expectation-qubits">After: phase ${st.phase}; residual ${st.mask}</div></li>`).join("");
        const reduction=calc.inSpan
          ? `Q = ${triangleNames.length?triangleNames.join(" · "):"I"} (product of complete B<sub>t</sub>), hence ⟨Ω|Q|Ω⟩ = +1.`
          : `Q = ${triangleNames.length?triangleNames.join(" · ")+" · ":""}${zExpression(calc.decomposition.remainder)}. The remaining Z-string is not generated by the complete B<sub>t</sub> constraints. ${calc.witness?`It intersects the outer X-support of A<sub>${calc.witness.id}</sub> on an odd number of qubits, so {A<sub>${calc.witness.id}</sub>, Z(remainder)} = 0 and ⟨Ω|Q|Ω⟩ = 0.`:"The current open-boundary checker assigns ⟨Ω|Q|Ω⟩ = 0; no complete-star anticommutation witness was found, so this is a checker assumption rather than a proof."}`;
        expectationContent.innerHTML=`<strong>${obj.color} star A<sub>${obj.id}</sub></strong> · center (${obj.center.map(x=>x.toFixed(2)).join(", ")})
          ${expectationResultHTML(calc.value)}
          <div class="expectation-qubits">Outer X tips: ${obj.x_qubits.join(", ")}<br>Inner CZ ring: ${obj.inner_qubits.join(" → ")} → ${obj.inner_qubits[0]}</div>
          <p>U is the ordered gate sequence (${operationLog.length} gates). Conjugate in <b>reverse chronological order</b>:</p>
          <div class="expectation-equation">U† A<sub>${obj.id}</sub> U = ${calc.phase} · A<sub>${obj.id}</sub> · Q<br>Q = ${zExpression(calc.mask)}</div>
          <details ${calc.steps.length<=12?"open":""}><summary>Gate-by-gate conjugation (${calc.steps.length} contributing; ${calc.unaffected} commuting)</summary>
          <ol class="expectation-steps">${steps||"<li>All applied gates commute with this star.</li>"}</ol></details>
          <p><strong>Triangle reduction:</strong> ${reduction}</p>
          <div class="expectation-equation">⟨Ω|U† A<sub>${obj.id}</sub> U|Ω⟩ = ${calc.phase} × ⟨Ω|Q|Ω⟩ = <b>${calc.value===1?"+1":(calc.value===-1?"−1":"0")}</b></div>
          ${calc.inSpan&&calc.phasePower%2!==0?'<p class="expectation-note">Consistency warning: an odd imaginary phase survived the B reduction. Check the reference state/boundary assumptions.</p>':''}
          <p class="expectation-note">Uses all complete triangle constraints on the open patch, matching the existing violation checker. Star IDs are internal to this lattice size.</p>`;
      } else {
        const calc=traceTriangleExpectation(obj);
        const items=calc.hits.map(hit=>`<li>Gate #${hit.number}: X<sub>${hit.q}</sub> anticommutes with B<sub>${obj.id}</sub> → sign flips.</li>`).join("");
        expectationContent.innerHTML=`<strong>${obj.color} triangle B<sub>${obj.id}</sub></strong> · center (${obj.center.map(x=>x.toFixed(2)).join(", ")})
          ${expectationResultHTML(calc.value)}
          <div class="expectation-qubits">B<sub>${obj.id}</sub> = ${obj.qubits.map(q=>`Z<sub>${q}</sub>`).join(" · ")}</div>
          <div class="expectation-equation">U† B<sub>${obj.id}</sub> U = (−1)<sup>${calc.hits.length}</sup> B<sub>${obj.id}</sub> = ${calc.value===1?"+":"−"}B<sub>${obj.id}</sub></div>
          <details open><summary>Contributing gates (${calc.hits.length})</summary><ol class="expectation-steps">${items||"<li>No X gates on this triangle.</li>"}</ol></details>
          <p>${calc.ignored} other gate${calc.ignored===1?"":"s"} commute with this B<sub>t</sub> (Z, R, CZ, or X on other qubits).</p>
          <div class="expectation-equation">⟨Ω|U† B<sub>${obj.id}</sub> U|Ω⟩ = (−1)<sup>${calc.hits.length}</sup> = <b>${calc.value===1?"+1":"−1"}</b></div>`;
      }
    }

    function polygonContains(x,y,ids) {
      const vertices=ids.map(q=>siteById.get(q));
      let inside=false;
      for(let i=0,j=vertices.length-1;i<vertices.length;j=i++){
        const a=vertices[i],b=vertices[j];
        if ((a.y>y)!==(b.y>y) && x < (b.x-a.x)*(y-a.y)/(b.y-a.y)+a.x) inside=!inside;
      }
      return inside;
    }

    function selectExpectationPoint(x,y) {
      if (!expectationMode) return false;
      const objects=expectationMode==="STAR"?starData:triangleData;
      const eligible=objects.filter(obj=>polygonContains(x,y,
        expectationMode==="STAR"?obj.inner_qubits:obj.qubits));
      if (!eligible.length) return false;
      eligible.sort((a,b)=>Math.hypot(a.center[0]-x,a.center[1]-y)-Math.hypot(b.center[0]-x,b.center[1]-y));
      const obj=eligible[0];
      inspectedObject={kind:expectationMode,id:obj.id};
      updateExpectationPanel();
      render();
      setStatus(`Inspecting ${expectationMode==="STAR"?"star A":"triangle B"}${obj.id}: expectation ${expectationMode==="STAR"?starExpectation(obj):triangleExpectation(obj)}.`);
      // The right panel may already have a tall tracking table.
      const panel=expectationPanel.closest(".qubit-panel");
      if (panel) panel.scrollTop=Math.max(0,expectationPanel.offsetTop-panel.offsetTop-10);
      return true;
    }

    function updateExpectationButtons() {
      for (const btn of expectationButtons) {
        const on=expectationMode===btn.dataset.expectation;
        btn.classList.toggle("active",on);
        btn.setAttribute("aria-pressed",String(on));
      }
    }

    function updateInspectionShapes(shapes) {
      if (!expectationMode) return;
      const objects=expectationMode==="STAR"?starData:triangleData;
      for(const obj of objects){
        const selected=inspectedObject&&inspectedObject.kind===expectationMode&&inspectedObject.id===obj.id;
        shapes.push({type:"path",path:polygonPath(expectationMode==="STAR"?obj.inner_qubits:obj.qubits),
          xref:"x",yref:"y",layer:"above",
          fillcolor: selected?"rgba(37,99,235,0.17)":"rgba(14,165,233,0.025)",
          line:{color:selected?"#1d4ed8":"rgba(14,116,144,0.45)",width:selected?3:1.25,dash:selected?"solid":"dot"}});
      }
    }

    function triangleExpectation(triangle) {
      const triangleSet = new Set(triangle.qubits);
      let sign = 1;

      // Z and CZ commute with Bt.  Each active X on one of the triangle's
      // three qubits anticommutes with Bt and flips its eigenvalue.
      for (const event of operationLog) {
        if (event.type !== "X" || event.qubits.length !== 1) {
          continue;
        }

        if (triangleSet.has(event.qubits[0])) {
          sign *= -1;
        }
      }

      return sign;
    }

    function starExpectation(star) {
      /*
       * Correct 12-body star operator from Fig. 2:
       *
       *   As = [ product of six CZ gates around the INNER hexagon ]
       *        [ product of X on the six OUTWARD tips ].
       *
       * The inner and X-support qubits are disjoint.
       *
       * After conjugating by the user's ordered Clifford circuit U,
       *
       *   U† As U = sign × As × Z(residualMask).
       *
       * Conjugation alone is NOT the final expectation-value test.  Because
       * Bt = +1 in the ground state, a residual Z string that is a product of
       * complete Bt triangle operators is equivalent to +1.  We therefore
       * reduce residualMask modulo the Bt stabilizer span below.
       */

      const inner = star.inner_qubits;
      const xTips = star.x_qubits;

      const innerSet = new Set(inner);
      const xTipSet = new Set(xTips);
      const innerIndex = new Map(
        inner.map((q, index) => [q, index])
      );

      let phasePower = 0;  // i^phasePower, needed for R = exp(+i pi Z/4)
      let zMask = 0n;

      // For U = G_n ... G_2 G_1, evaluate U† As U by conjugating in
      // reverse chronological order through the active operation log.
      for (let opIndex = operationLog.length - 1; opIndex >= 0; opIndex--) {
        const event = operationLog[opIndex];

        // --------------------------------------------------------------
        // User Z_q
        //
        // Z anticommutes only with an X factor of As.  Therefore a Z on
        // one of the six outward tips flips this star's sign.
        // --------------------------------------------------------------
        if (event.type === "Z" && event.qubits.length === 1) {
          const q = event.qubits[0];

          if (xTipSet.has(q)) {
            phasePower += 2;
          }

          continue;
        }

        // --------------------------------------------------------------
        // User X_q
        //
        // If q lies on the INNER CZ ring:
        //
        //   X_q CZ_(p,q) CZ_(q,r) X_q
        //       = CZ_(p,q) CZ_(q,r) Z_p Z_r.
        //
        // Hence X on an inner-hexagon vertex adds Z on its two neighbours
        // around that hexagon.
        //
        // X also flips the sign of any already-existing residual Z_q.
        // --------------------------------------------------------------
        if (event.type === "X" && event.qubits.length === 1) {
          const q = event.qubits[0];

          if ((zMask & qubitBit(q)) !== 0n) {
            phasePower += 2;
          }

          if (innerSet.has(q)) {
            const idx = innerIndex.get(q);
            const previous = inner[
              (idx + inner.length - 1) % inner.length
            ];
            const next = inner[
              (idx + 1) % inner.length
            ];

            zMask ^= qubitBit(previous);
            zMask ^= qubitBit(next);
          }

          continue;
        }

        // --------------------------------------------------------------
        // User CZ_(a,b)
        //
        // CZ commutes with the inner CZ ring.  It only decorates the
        // OUTWARD X support:
        //
        //   CZ_ab X_a CZ_ab = X_a Z_b.
        //
        // If both a and b are X-support tips, both Z factors appear and
        // there is the usual minus sign when put into X...Z normal order.
        // --------------------------------------------------------------
        if (event.type === "R" && event.qubits.length === 1) {
          const q = event.qubits[0];
          if (xTipSet.has(q)) {
            // R† X R = Y = i X Z; R commutes with every diagonal Z/CZ.
            phasePower += 1;
            zMask ^= qubitBit(q);
          }
          continue;
        }
        if (event.type === "CZ" && event.qubits.length === 2) {
          const a = event.qubits[0];
          const b = event.qubits[1];

          const aHasX = xTipSet.has(a);
          const bHasX = xTipSet.has(b);

          if (aHasX && bHasX) {
            phasePower += 2;
          }

          if (aHasX) {
            zMask ^= qubitBit(b);
          }

          if (bHasX) {
            zMask ^= qubitBit(a);
          }
        }
      }

      /*
       * Since As |psi0> = |psi0>,
       *
       *   <psi0| U† As U |psi0>
       *     = sign × <psi0| Z(residualMask) |psi0>.
       *
       * For our open-boundary implementation we keep only complete Bt
       * triangles.  If residualMask is a GF(2) sum/product of those Bt
       * operators, its expectation is +1.  Otherwise the corresponding
       * diagonal Z string is outside the retained local Bt stabilizer group
       * and its expectation is taken to be 0.
       */
      if (!isTriangleStabilizerMask(zMask)) {
        return 0;
      }

      const p=((phasePower % 4)+4)%4;
      // All physical star expectations are real; an imaginary result here
      // indicates an inconsistent stabilizer/boundary reduction.
      return p===0?1:(p===2?-1:0);
    }

    function violatedStars() {
      return starData
        .map(star => ({
          star: star,
          expectation: starExpectation(star),
        }))
        .filter(item => item.expectation !== 1);
    }

    function violatedTriangles() {
      return triangleData
        .map(triangle => ({
          triangle: triangle,
          expectation: triangleExpectation(triangle),
        }))
        .filter(item => item.expectation !== 1);
    }

    function polygonPath(qubits) {
      const points = qubits.map(siteId => siteById.get(siteId));

      if (points.some(point => !point)) {
        return "";
      }

      const first = points[0];
      let path = `M ${first.x},${first.y}`;

      for (let i = 1; i < points.length; i++) {
        path += ` L ${points[i].x},${points[i].y}`;
      }

      path += " Z";
      return path;
    }

    function buildViolationShapes() {
      const shapes = [];

      if (violationVisibility.STAR) {
        for (const item of violatedStars()) {
          const path = polygonPath(item.star.inner_qubits);

          if (!path) {
            continue;
          }

          shapes.push({
            type: "path",
            path: path,
            xref: "x",
            yref: "y",
            layer: "below",
            fillcolor: "rgba(239, 68, 68, 0.11)",
            line: {
              color: "rgba(239, 68, 68, 0.30)",
              width: 2,
            },
          });
        }
      }

      if (violationVisibility.TRIANGLE) {
        for (const item of violatedTriangles()) {
          const path = polygonPath(item.triangle.qubits);

          if (!path) {
            continue;
          }

          shapes.push({
            type: "path",
            path: path,
            xref: "x",
            yref: "y",
            layer: "above",
            fillcolor: "rgba(0, 0, 0, 0)",
            line: {
              color: "rgba(220, 38, 38, 0.72)",
              width: 2,
              dash: "dot",
            },
          });
        }
      }

      updateInspectionShapes(shapes);
      return shapes;
    }

    function updateViolationButtons() {
      for (const btn of violationButtons) {
        const key = btn.dataset.violation;
        const isVisible = violationVisibility[key];

        if (isVisible) {
          btn.classList.add("active");
          btn.setAttribute("aria-pressed", "true");
        } else {
          btn.classList.remove("active");
          btn.setAttribute("aria-pressed", "false");
        }
      }
    }

    function setZoomMode(enabled) {
      zoomMode = enabled;

      if (zoomMode) {
        zoomBtn.classList.add("active");

        Plotly.relayout(plotDiv, {
          dragmode: "zoom"
        });

        setStatus(
          "Zoom mode: click and drag the lattice over area to zoom into."
        );
      } else {
        zoomBtn.classList.remove("active");

        Plotly.relayout(plotDiv, {
          dragmode: "pan"
        });
      }
    }

    function quadraticBezierPoints(x0, y0, x1, y1, cx, cy, n = 32) {
      const xs = [];
      const ys = [];

      for (let k = 0; k < n; k++) {
        const t = k / (n - 1);
        const omt = 1 - t;

        const x =
          omt * omt * x0 +
          2 * omt * t * cx +
          t * t * x1;

        const y =
          omt * omt * y0 +
          2 * omt * t * cy +
          t * t * y1;

        xs.push(x);
        ys.push(y);
      }

      return { xs, ys };
    }

    function buildCZGeometry() {
      const lineX = [];
      const lineY = [];
      const labelX = [];
      const labelY = [];
      const labelText = [];

      for (const [a, b] of czPairs) {
        const sa = siteById.get(a);
        const sb = siteById.get(b);

        if (!sa || !sb) continue;

        const x0 = sa.x;
        const y0 = sa.y;
        const x1 = sb.x;
        const y1 = sb.y;

        if (isNearestNeighbor(a, b)) {
          lineX.push(x0, x1, null);
          lineY.push(y0, y1, null);

          const midX = 0.5 * (x0 + x1);
          const midY = 0.5 * (y0 + y1);

          labelX.push(midX);
          labelY.push(midY + 0.18);
          labelText.push("CZ");
        } else {
          const midX = 0.5 * (x0 + x1);
          const midY = 0.5 * (y0 + y1);
          const dist = Math.hypot(x1 - x0, y1 - y0);

          // "curve up" in +y direction
          const lift = Math.max(0.75, 0.20 * dist);
          const controlX = midX;
          const controlY = midY + lift;

          const curve = quadraticBezierPoints(
            x0, y0,
            x1, y1,
            controlX, controlY,
            36
          );

          lineX.push(...curve.xs, null);
          lineY.push(...curve.ys, null);

          // place label at the top of the arc
          let apexIndex = 0;
          for (let i = 1; i < curve.ys.length; i++) {
            if (curve.ys[i] > curve.ys[apexIndex]) {
              apexIndex = i;
            }
          }

          labelX.push(curve.xs[apexIndex]);
          labelY.push(curve.ys[apexIndex] + 0.18);
          labelText.push("CZ");
        }
      }

      return { lineX, lineY, labelX, labelY, labelText };
    }

    function buildOpLabels() {
      const xs = [],ys = [],texts = [],customdata = [];
      const bySite = new Map();
      for (const ev of operationLog) {
        if (!["X","Z","R"].includes(ev.type) || ev.qubits.length!==1) continue;
        const id=ev.qubits[0];
        if (!bySite.has(id)) bySite.set(id,[]);
        bySite.get(id).push(ev.type);
      }
      for (const [siteId,gates] of bySite) {
        const site=siteById.get(siteId);if(!site)continue;
        const visible=gates.filter(type=>operatorVisibility[type]);
        if(!visible.length)continue;
        xs.push(site.x+0.18);ys.push(site.y+0.25);
        texts.push([...visible].reverse().join(""));
        customdata.push(siteId);
      }
      return {xs,ys,texts,customdata};
    }

    function render() {
      const xs = siteData.map(s => s.x);
      const ys = siteData.map(s => s.y);

      const minX = Math.min(...xs) - 1.2;
      const maxX = Math.max(...xs) + 1.2;
      const minY = Math.min(...ys) - 1.2;
      const maxY = Math.max(...ys) + 1.2;

      const cz = buildCZGeometry();
      const opLabels = buildOpLabels();
      const nextAnyonQubits = availableAnyonQubits();
      updateAnyonSummary();
      const nextAnyonSet = new Set(nextAnyonQubits);
      const usedAnyonSet = new Set((anyonKind==="s1" || anyonKind==="s2") ? [...sGuidedX,...sBasePath] : anyonQubits);

      const edgeTrace = {
        type: "scatter",
        mode: "lines",
        x: baseEdgeX,
        y: baseEdgeY,
        line: {
          color: "rgba(100, 116, 139, 0.65)",
          width: 1.6,
        },
        hoverinfo: "skip",
        showlegend: false,
      };

      const czLineTrace = {
        type: "scatter",
        mode: "lines",
        x: cz.lineX,
        y: cz.lineY,
        line: {
          color: "#7c3aed",
          width: 3.0,
        },
        hoverinfo: "skip",
        showlegend: false,
        visible: operatorVisibility.CZ,
      };

      const czLabelTrace = {
        type: "scatter",
        mode: "text",
        x: cz.labelX,
        y: cz.labelY,
        text: cz.labelText,
        textfont: {
          size: 14,
          color: "#6d28d9",
          family: "Arial, sans-serif",
        },
        hoverinfo: "skip",
        showlegend: false,
        visible: operatorVisibility.CZ,
      };

      const nodeTrace = {
        name: "qubit-nodes",
        type: "scatter",
        mode: "markers",
        x: siteData.map(s => s.x),
        y: siteData.map(s => s.y),
        customdata: siteData.map(s => s.id),
        marker: {
          size: 14,
          color: siteData.map(s => s.color),
          opacity: siteData.map(s => (
            !anyonKind || nextAnyonSet.has(s.id) || usedAnyonSet.has(s.id) ? 1 : 0.3
          )),
          line: {
            width: 1.1,
            color: "#111827",
          },
        },
        hovertemplate: "<b>Qubit %{customdata}</b><extra></extra>",
        showlegend: false,
      };

      const opLabelTrace = {
          name: "operator-labels",
          type: "scatter",
          mode: "text",

          x: opLabels.xs,
          y: opLabels.ys,

          text: opLabels.texts,
          customdata: opLabels.customdata,

          textposition: "middle center",

          textfont: {
            size: 16,
            color: "#111827",
            family: "Arial, sans-serif",
          },

          hoverinfo: "skip",
          showlegend: false,
    };

      const pendingTrace = pendingCZ === null ? null : {
        name: "pending-cz",
        type: "scatter",
        mode: "markers",
        x: [siteById.get(pendingCZ).x],
        y: [siteById.get(pendingCZ).y],
        customdata: [pendingCZ],
        marker: {
          size: 24,
          color: "rgba(0,0,0,0)",
          line: {
            width: 3,
            color: "#7c3aed",
          },
        },
        hovertemplate: "<b>Qubit %{customdata}</b><extra></extra>",
        showlegend: false,
      };

      // Completed guided Z/X string qubits are outlined teal. Valid next
      // qubits receive gold halos; non-options are dimmed.
      const anyonPathTrace = !anyonKind || anyonQubits.length === 0 ? null : {
        name: "anyon-path",
        type: "scatter",
        mode: "markers",
        x: [...usedAnyonSet].map(id => siteById.get(id).x),
        y: [...usedAnyonSet].map(id => siteById.get(id).y),
        customdata: [...usedAnyonSet],
        marker: {
          size: 25,
          color: "rgba(13,148,136,0.10)",
          line: { color: "#0f766e", width: 3 },
        },
        hovertemplate: `<b>Qubit %{customdata}</b><br>${currentAnyonGate()}-string (already selected)<extra></extra>`,
        showlegend: false,
      };

      const anyonNextTrace = !anyonKind ? null : {
        name: "anyon-next",
        type: "scatter",
        mode: "markers",
        x: nextAnyonQubits.map(id => siteById.get(id).x),
        y: nextAnyonQubits.map(id => siteById.get(id).y),
        customdata: [...nextAnyonQubits],
        marker: {
          size: 27,
          color: "rgba(251,191,36,0.18)",
          line: { color: "#d97706", width: 3 },
        },
        hovertemplate: `<b>Qubit %{customdata}</b><br>Available next ${anyonName()} step<extra></extra>`,
        showlegend: false,
      };

      const traces = [
        edgeTrace,
        czLineTrace,
        czLabelTrace,
        nodeTrace,
        opLabelTrace,
      ];

      if (pendingTrace !== null) {
        traces.push(pendingTrace);
      }
      if (anyonPathTrace !== null) {
        traces.push(anyonPathTrace);
      }
      if (anyonNextTrace !== null) {
        traces.push(anyonNextTrace);
      }
      if (expectationMode) {
        const objects=expectationMode==="STAR"?starData:triangleData;
        traces.push({
          name:"inspect-centers", type:"scatter", mode:"markers", showlegend:false,
          x:objects.map(obj=>obj.center[0]), y:objects.map(obj=>obj.center[1]),
          customdata:objects.map(obj=>obj.id),
          marker:{size:12,color:"rgba(14,116,144,0.32)",line:{color:"#0e7490",width:1}},
          hovertemplate:`<b>Click ${expectationMode==="STAR"?"star A":"triangle B"}%{customdata} to calculate expectation</b><extra></extra>`,
        });
      }

      traces.push(...buildMPairGuideTraces());

      const layout = {
          // Keep the user's current zoom/pan whenever Plotly.react()
          // is called after X, Z, CZ, Erase, Undo, etc.
          uirevision: "keep-kagome-view",

          // Set click-and-drag to pan by default
          dragmode: zoomMode ? "zoom" : "pan",

          paper_bgcolor: "white",
          plot_bgcolor: "white",

          // Star glows and dotted Bt triangles are display overlays only.
          shapes: buildViolationShapes(),

          margin: {
            l: 10,
            r: 10,
            t: 10,
            b: 10
          },

          xaxis: {
            visible: false,
            range: [minX, maxX],
            scaleanchor: "y",
            scaleratio: 1,
            fixedrange: false,
          },

          yaxis: {
            visible: false,
            range: [minY, maxY],
            fixedrange: false,
          },
        };

      const config = {
        responsive: true,
        displayModeBar: false,
        scrollZoom: true,
      };

      Plotly.react(plotDiv, traces, layout, config);
      updateTrackedTable();
      updateExpectationPanel();
    }

    function describeConnectedCZs(siteId, pairs) {
      if (pairs.length === 0) return "";

      const parts = pairs.map(([a, b]) => {
        const u = Math.min(a, b);
        const v = Math.max(a, b);
        return `CZ between qubits ${u} and ${v}`;
      });

      if (parts.length === 1) {
        return parts[0];
      }

      return parts.join(", ");
    }

    function handleLocalOp(siteId, newOp) {
      const oldOp = localOps.get(siteId);

      if (oldOp === newOp) {
        pushHistory(`removing ${newOp} from qubit ${siteId}`);
        removeLatestLocalOperation(siteId, newOp);
        syncVisibleLocalOp(siteId);
        render();
        setStatus(`Removed ${newOp} from qubit ${siteId}.`);
        return;
      }

      if (oldOp && oldOp !== newOp) {
        pushHistory(`applying ${newOp} to qubit ${siteId}`);
        recordOperation(newOp, [siteId]);
        localOps.set(siteId, newOp);
        render();
        setStatus(`Applied ${newOp} to qubit ${siteId}.`);
        return;
      }

      pushHistory(`applying ${newOp} to qubit ${siteId}`);
      recordOperation(newOp, [siteId]);
      localOps.set(siteId, newOp);
      render();
      setStatus(`Applied ${newOp} to qubit ${siteId}.`);
    }

    function handleR(siteId) {
      const rEntries=operationLog.filter(ev=>ev.type==="R" && ev.qubits.length===1 && ev.qubits[0]===siteId);
      pushHistory(`${rEntries.length ? "removing" : "applying"} R on qubit ${siteId}`);
      if (rEntries.length) {
        operationLog=operationLog.filter(ev=>!(ev.type==="R" && ev.qubits.length===1 && ev.qubits[0]===siteId));
      } else {
        recordOperation("R",[siteId]);
      }
      render();
      setStatus(`R on qubit ${siteId}: ${rEntries.length ? "removed" : "applied"}.  R = exp(+iπZ/4).`);
    }

    function handleCZ(siteId) {
      if (pendingCZ === null) {
        pushHistory(`selecting qubit ${siteId} as the first CZ qubit`);
        pendingCZ = siteId;
        render();
        setStatus(`CZ: first qubit ${siteId} selected.`);
        return;
      }

      if (pendingCZ === siteId) {
        pushHistory(`cancelling pending CZ selection on qubit ${siteId}`);
        pendingCZ = null;
        render();
        setStatus(`Cancelled pending CZ selection on qubit ${siteId}.`);
        return;
      }

      const a = Math.min(pendingCZ, siteId);
      const b = Math.max(pendingCZ, siteId);

      if (hasCZPair(a, b)) {
        pushHistory(`removing CZ between qubits ${a} and ${b}`);
        czPairs = czPairs.filter(([u, v]) => pairKey(u, v) !== pairKey(a, b));
        removeCZOperation(a, b);
        pendingCZ = null;
        render();
        setStatus(`Removed CZ between qubits ${a} and ${b}.`);
        return;
      }

      pushHistory(`applying CZ to qubits ${a} and ${b}`);
      czPairs.push([a, b]);
      recordOperation("CZ", [a, b]);
      pendingCZ = null;
      render();
      setStatus(`Applied CZ to qubits ${a} and ${b}.`);
    }

    function handleErase(siteId) {
      const existingOp = localOps.get(siteId) || null;
      const hasR = operationLog.some(ev=>ev.type==="R" && ev.qubits.length===1 && ev.qubits[0]===siteId);
      const connectedPairs = czPairs.filter(([a, b]) => a === siteId || b === siteId);
      const wasPending = pendingCZ === siteId;

      if (!existingOp && !hasR && connectedPairs.length === 0 && !wasPending) {
        setStatus(`Qubit ${siteId} has no operations to erase.`);
        return;
      }

      const descriptionParts = [];

      if (existingOp) {
        descriptionParts.push(`${existingOp} on qubit ${siteId}`);
      }

      if (hasR) descriptionParts.push(`R on qubit ${siteId}`);
      if (connectedPairs.length > 0) {
        descriptionParts.push(describeConnectedCZs(siteId, connectedPairs));
      }

      if (wasPending) {
        descriptionParts.push(`pending CZ selection on qubit ${siteId}`);
      }

      pushHistory(`erasing ${descriptionParts.join(" and ")}`);

      if (existingOp) {
        localOps.delete(siteId);
      }

      if (connectedPairs.length > 0) {
        const keysToRemove = new Set(connectedPairs.map(([a, b]) => pairKey(a, b)));
        czPairs = czPairs.filter(([a, b]) => !keysToRemove.has(pairKey(a, b)));
      }

      if (wasPending) {
        pendingCZ = null;
      }

      // Erase invalidates any actively guided path; the old m-specific
      // block must not be rebuilt after some of its gates have been erased.
      if (anyonKind) {
        clearAnyonGuidance();
        currentMode = "ERASE";
        updateActiveModeButtons();
      }

      // Erasing a qubit removes every tracked operation involving it.
      // For CZ this removes the same two-qubit gate from both endpoints' rows.
      operationLog = operationLog.filter(event => !event.qubits.includes(siteId));

      render();
      setStatus(`Erased ${descriptionParts.join(" and ")}.`);
    }

    function handleUndo() {
      if (history.length === 0) {
        setStatus("Nothing to undo.");
        return;
      }

      const last = history.pop();
      restoreState(last.state);
      setStatus(`Undone ${last.description}.`);
    }

    function handleClear() {
      if (operationLog.length === 0 && czPairs.length === 0 && pendingCZ === null && anyonQubits.length === 0) {
        setStatus("Nothing to clear.");
        return;
      }

      pushHistory("clearing all operations");

      localOps.clear();
      czPairs = [];
      pendingCZ = null;
      operationLog = [];
      nextOperationId = 1;
      anyonQubits = [];
      anyonStarPath = [];
      mTrianglePath = [];
      mAutoCZ = [];
      mOptionCache = null;
      mGuideId = anyonKind === "m" ? nextMGuideId++ : null;
      mPairEdges=[];mPairOptionCache=null;
      mPairGuideId = anyonKind === "m2" ? nextMPairGuideId++ : null;
      sBasePath=[];sGuidedX=[];sAutoCZ=[];sCandidateCache=null;
      sGuideId = (anyonKind === "s1" || anyonKind === "s2") ? nextSGuideId++ : null;

      render();
      setStatus("Cleared all operations.");
    }

    function handleSave() {
      const now = new Date();
      const filename =
        "kagome_lattice_" +
        now.getFullYear() +
        String(now.getMonth() + 1).padStart(2, "0") +
        String(now.getDate()).padStart(2, "0") +
        "_" +
        String(now.getHours()).padStart(2, "0") +
        String(now.getMinutes()).padStart(2, "0") +
        String(now.getSeconds()).padStart(2, "0");

      Plotly.downloadImage(plotDiv, {
        format: "png",
        filename: filename,
        scale: 2,
        width: 1800,
        height: 1200,
      });

      setStatus(`Saved current Kagome lattice as ${filename}.png`);
    }

    function handleSiteClick(siteId) {
      if (anyonKind) {
        handleAnyonQubit(siteId);
        return;
      }
      if (currentMode === "R") {
        handleR(siteId);
      } else if (currentMode === "Z") {
        handleLocalOp(siteId, "Z");
      } else if (currentMode === "X") {
        handleLocalOp(siteId, "X");
      } else if (currentMode === "CZ") {
        handleCZ(siteId);
      } else if (currentMode === "ERASE") {
        handleErase(siteId);
      }
    }

    anyonSelect.addEventListener("change", () => {
      setAnyon(anyonSelect.value);
    });

    // Mode button events
    for (const btn of modeButtons) {
      btn.addEventListener("click", () => {

        // Return from navigation mode to operator mode.
        if (zoomMode) {
          setZoomMode(false);
        }

        // Ordinary tools exit Anyon mode, without changing any Z, X, or
        // automatically generated CZ gates already placed by the guide.
        if (anyonKind) {
          clearAnyonGuidance();
        }
        currentMode = btn.dataset.mode;

        updateActiveModeButtons();
        render();

        if (currentMode === "CZ") {
          setStatus(
            "Mode: CZ. Click two qubits to connect them."
          );
        } else if (currentMode === "ERASE") {
          setStatus(
            "Mode: Erase. Click a qubit to remove its local operator and any CZ connected to it."
          );
        } else {
          setStatus(`Mode: ${currentMode}`);
        }
      });
    }
    
    // Visibility controls only affect drawing.
    // Operators remain active and continue to appear in the tracked table.
    for (const btn of visibilityButtons) {
      btn.addEventListener("click", () => {
        const key = btn.dataset.visibility;
        operatorVisibility[key] = !operatorVisibility[key];

        updateVisibilityButtons();
        render();

        const stateText = operatorVisibility[key] ? "shown" : "hidden";

        if (key === "CZ") {
          setStatus(`CZ connections and labels are now ${stateText}.`);
        } else {
          setStatus(`${key} labels are now ${stateText}.`);
        }
      });
    }


    // Violation display controls.  Toggling these buttons does not modify the
    // active gate sequence; it only turns the Hamiltonian diagnostics on/off.
    for (const btn of violationButtons) {
      btn.addEventListener("click", () => {
        const key = btn.dataset.violation;
        violationVisibility[key] = !violationVisibility[key];

        updateViolationButtons();
        render();

        if (violationVisibility[key]) {
          if (key === "STAR") {
            const count = violatedStars().length;
            setStatus(
              `Star violation display on: ${count} violated complete star${count === 1 ? "" : "s"}.`
            );
          } else {
            const count = violatedTriangles().length;
            setStatus(
              `Triangle violation display on: ${count} violated complete triangle${count === 1 ? "" : "s"}.`
            );
          }
        } else {
          setStatus(
            key === "STAR"
              ? "Star violation display off."
              : "Triangle violation display off."
          );
        }
      });
    }

    // Calculation tools do not change any operations. Clicking a second time
    // closes selection mode but leaves the last calculation visible.
    for (const btn of expectationButtons) {
      btn.addEventListener("click", () => {
        const next=btn.dataset.expectation;
        expectationMode=expectationMode===next?null:next;
        if (expectationMode && inspectedObject && inspectedObject.kind!==expectationMode) {
          inspectedObject=null;
        }
        if (zoomMode && expectationMode) setZoomMode(false);
        updateExpectationButtons();
        render();
        setStatus(expectationMode
          ? `Click inside any complete ${expectationMode==="STAR"?"hexagonal star":"triangle"} to calculate its expectation.`
          : "Expectation selection off; the last calculation remains in the right panel.");
      });
    }

    // Unlike Plotly point-click, this also catches empty polygon interiors.
    plotDiv.addEventListener("click", (evt) => {
      if (!expectationMode || zoomMode) return;
      const drag=plotDiv.querySelector(".nsewdrag");
      const fl=plotDiv._fullLayout;
      if (!drag || !fl || !fl.xaxis || !fl.yaxis) return;
      const rect=drag.getBoundingClientRect();
      const px=evt.clientX-rect.left,py=evt.clientY-rect.top;
      if (px<0 || py<0 || px>rect.width || py>rect.height) return;
      const xr=fl.xaxis.range,yr=fl.yaxis.range;
      const x=xr[0]+px/rect.width*(xr[1]-xr[0]);
      const y=yr[1]-py/rect.height*(yr[1]-yr[0]);
      selectExpectationPoint(x,y);
    });

    zoomBtn.addEventListener("click", () => {
      setZoomMode(!zoomMode);
    });
    saveBtn.addEventListener("click", handleSave);
    undoBtn.addEventListener("click", handleUndo);
    clearBtn.addEventListener("click", handleClear);

    addQubitBtn.addEventListener("click", () => {
      if (qubitSelect.value === "__ALL_ACTIVE__") {
        addAllActiveQubits();
        return;
      }

      const siteId = Number(qubitSelect.value);
      if (!Number.isNaN(siteId)) {
        addTrackedQubit(siteId);
      }
    });

    qubitSelect.addEventListener("keydown", (event) => {
      if (event.key === "Enter") {
        event.preventDefault();
        addQubitBtn.click();
      }
    });

    // ====================================================================
    // Microscopic R workbench — independent of the editable lattice.
    // This is a diagnostic protocol, NOT an assumed semion phase.
    // ====================================================================
    const rSelect=document.getElementById('r-select');
    const rStudio=document.getElementById('r-studio');
    const rPlot=document.getElementById('r-plot');
    const rCalc=document.getElementById('r-calc');
    const rCheckpointSelect=document.getElementById('r-checkpoint');
    const rValidationContent=document.getElementById('r-validation-content');
    const rOverview=document.getElementById('r-overview');
    const rEditX=document.getElementById('r-edit-x');
    const rEditStars=document.getElementById('r-edit-stars');
    const rFixedXNote=document.getElementById('r-fixed-x-note');
    const rSaveAdjust=document.getElementById('r-save-adjust');
    const rStepButtons=Array.from(document.querySelectorAll('.r-step'));
    const rSequences={A:['S','Y','M21','M0'],B:['S','M0','M21','Y']};
    const rDescriptions={S:'Vacuum-pair splitting S(1,2)',Y:'Move 1 → Y',M21:'Move 2 → 1',M0:'Move 1 → 0'};
    let rSelectedKind=null,rGeom=null,rSegments=null;
    let rActiveProcess='A',rActiveStep=0;
    let rAdjust={Y:{x:[],stars:[]},M21:{x:[],stars:[]},M0:{x:[],stars:[]}};
    let rCompleted={A:0,B:0};
    let rExactCache=null,rCalcToken=0,rValidationToken=0;
    let rCheckpoint='projected';
    let rCompiling=false;
    let rM21SyndromeSuggestion=null;
    // Fixed DESTINATION-local projectors are separate from the four editable
    // SOURCE-interface projectors.  Both process orders use the same mapping.
    let rEndpointFix={Y:null,M0:null};
    // Fixed splitting-state preparation, calibrated ONLY for the present
    // 12x12 / six-step q223-q232-q241-q133 geometry.  Never silently reuse
    // these numerical star IDs in a lattice with a different numbering.
    let rInitialSplittingStars=[];
    let rExtraM21Star=null;
    function rSplittingCalibration(geom){
      const matches=starData.length===100 && geom &&
        geom.q0===223 && geom.q1===232 && geom.q2===241 && geom.qY===133 &&
        starById.has(53) && starById.has(48);
      return matches?{split:[53,48],extraM21:53}:{split:[],extraM21:null};
    }

    function rResize() {
      // Streamlit HTML components accept the standard component height message.
      if(window.parent && window.parent!==window) {
        const h=Math.max(810,document.documentElement.scrollHeight+24);
        window.parent.postMessage({isStreamlitMessage:true,type:'streamlit:setFrameHeight',height:h},'*');
      }
    }
    function rEscape(t) {return String(t).replace(/[&<>"']/g, c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));}
    function rPoint(id) {return siteById.get(id);}
    function rNearSite(x,y) {return siteData.find(s=>Math.hypot(s.x-x,s.y-y)<1e-5)?.id;}
    function rWalk(start,dx,dy,n) {
      let path=[start];const startSite=rPoint(start);
      for(let k=1;k<n;k++) {
        const next=rNearSite(startSite.x+k*dx,startSite.y+k*dy);
        if(next===undefined||!sStraightNeighbors.has(pairKey(path[path.length-1],next)))return null;
        path.push(next);
      }
      return path;
    }
    function rMakeRibbon(path,kind) {
      const was=anyonKind; const saved= rCompiling;
      anyonKind=kind;rCompiling=true;
      try {
        const built=buildSSpine(path);
        if(!built||!built.xIds.length)return null;
        const cz=compileSDecorations(built);
        if(cz===null)return null;
        const events=cz.map(qs=>({type:'CZ',qubits:[...qs],source:'r-ribbon'}));
        for(const q of built.xIds) {
          const order=built.gateOrder.get(q);
          if(order==='XR') {events.push({type:'R',qubits:[q],source:'r-ribbon'});events.push({type:'X',qubits:[q],source:'r-ribbon'});}
          else if(order==='RX') {events.push({type:'X',qubits:[q],source:'r-ribbon'});events.push({type:'R',qubits:[q],source:'r-ribbon'});}
          else events.push({type:'X',qubits:[q],source:'r-ribbon'});
        }
        return {base:path,events,cz,ordered:[...built.gateOrder.entries()],xSites:built.xIds};
      } finally {anyonKind=was;rCompiling=saved;}
    }
    // Exactly six nearest-neighbor steps per arm.  For these shorter arms,
    // the closest diagonal endpoint pair can contain only ONE complete
    // corridor star with the original stringent 2.5-unit endpoint exclusion.
    // Never silently claim the old two-star criterion still holds.
    const rMinBufferStars=1;
    const rPreferredArmSteps=[6];
    const rRequiredBoundaryStars=1;
    // Geometric boundary of the physical OPEN Kagome patch, not an axis-
    // aligned bounding box (the latter badly misidentifies slanted edges).
    const rBoundaryHull=(()=>{
      const pts=siteData.map(p=>[p.x,p.y]).sort((a,b)=>a[0]-b[0]||a[1]-b[1]);
      const cross=(a,b,c)=>(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
      const lower=[],upper=[];
      for(const p of pts){while(lower.length>1&&cross(lower[lower.length-2],lower.at(-1),p)<=0)lower.pop();lower.push(p);}
      for(let k=pts.length-1;k>=0;k--){const p=pts[k];while(upper.length>1&&cross(upper[upper.length-2],upper.at(-1),p)<=0)upper.pop();upper.push(p);}
      return lower.slice(0,-1).concat(upper.slice(0,-1));
    })();
    function rClosestBoundaryPoint(id){
      const pt=rPoint(id);if(!pt)return null;
      let best=null;
      for(let i=0;i<rBoundaryHull.length;i++){
        const a=rBoundaryHull[i],b=rBoundaryHull[(i+1)%rBoundaryHull.length];
        const dx=b[0]-a[0],dy=b[1]-a[1];
        const frac=Math.max(0,Math.min(1,((pt.x-a[0])*dx+(pt.y-a[1])*dy)/(dx*dx+dy*dy)));
        const x=a[0]+frac*dx,y=a[1]+frac*dy,d=Math.hypot(x-pt.x,y-pt.y);
        if(!best||d<best.distance)best={x,y,distance:d};
      }
      return best;
    }
    function rBoundaryStarsForEndpoint(id){
      const pt=rPoint(id),edge=rClosestBoundaryPoint(id);
      if(!pt||!edge||edge.distance<1e-6)return [];
      const ux=(edge.x-pt.x)/edge.distance,uy=(edge.y-pt.y)/edge.distance;
      // A qualifying star is COMPLETE, sits between the RGB endpoint and
      // the closest physical boundary, and is on the outward radial strip.
      return starData.map(st=>{
        const dx=st.center[0]-pt.x,dy=st.center[1]-pt.y;
        return {star:st,along:dx*ux+dy*uy,across:Math.abs(dx*uy-dy*ux)};
      }).filter(o=>o.along>=1.25-1e-6&&o.along<=edge.distance-.4+1e-6&&o.across<=1.3+1e-6)
        .sort((a,b)=>a.along-b.along||a.across-b.across)
        .map(o=>o.star);
    }
    function rBoundaryReport(events,endpoints,engine=null,state=null){
      const flux=new Set(rFluxIds(rUnitaryOnly(events)).map(t=>t.id));
      const perEndpoint=endpoints.map(q=>{
        const stars=rBoundaryStarsForEndpoint(q);
        const evaluated=stars.map(st=>{
          const tips=new Set(st.x_qubits);
          const fluxFree=triangleData.filter(t=>t.qubits.some(id=>tips.has(id)))
            .every(t=>!flux.has(t.id));
          const value=engine&&state&&!state.dead?engine.starExpectation(state,st.id):null;
          const restored=value===null?null:Number.isFinite(value)&&Math.abs(value-1)<1e-7;
          return {id:st.id,fluxFree,expectation:value,clean:fluxFree&&restored===true};
        });
        return {endpoint:q,starIDs:stars.map(s=>s.id),geometricCount:stars.length,
          fluxFreeCount:evaluated.filter(x=>x.fluxFree).length,
          quantumCleanCount:evaluated.filter(x=>x.clean).length,
          geometryOK:stars.length>=rRequiredBoundaryStars,
          quantumOK:engine&&state&&!state.dead?evaluated.some(x=>x.clean):null,
          evaluations:evaluated};
      });
      return {perEndpoint,geometryOK:perEndpoint.every(x=>x.geometryOK),
        quantumOK:engine&&state&&!state.dead?perEndpoint.every(x=>x.quantumOK):null};
    }
    function rBufferStarsBetween(a,b) {
      const p=rPoint(a),q=rPoint(b);
      if(!p||!q)return [];
      const vx=q.x-p.x,vy=q.y-p.y,L=Math.hypot(vx,vy);
      if(L<5.0)return [];
      const ux=vx/L,uy=vy/L;
      return starData.map(st=>{
        const dx=st.center[0]-p.x,dy=st.center[1]-p.y;
        return {star:st,along:dx*ux+dy*uy,across:Math.abs(dx*uy-dy*ux)};
      }).filter(x=>x.along>=2.5-1e-6&&x.along<=L-2.5+1e-6&&x.across<=1.3+1e-6)
        .sort((x,y)=>x.along-y.along||x.across-y.across)
        .map(x=>x.star);
    }
    function rGapPairs(q0,q1,q2,qY){
      // All pairs that occur after splitting or after any intermediate move.
      return [[q1,q2],[qY,q2],[qY,q1],[q0,q2],[q0,q1],[q0,qY]];
    }
    function rFindGeometry(kind) {
      const xs=siteData.map(s=>s.x),ys=siteData.map(s=>s.y);
      const xmid=(Math.min(...xs)+Math.max(...xs))/2,ymid=(Math.min(...ys)+Math.max(...ys))/2;
      const rt=Math.sqrt(3)/2;
      // Y should lie BELOW the horizontal 0--1--2 axis as in the diagram.
      const dirs=[[.5,-rt],[-.5,-rt]];
      for(const armSteps of rPreferredArmSteps){
        const candidates=[];
        for(const site of siteData){
          const q1=site.id;
          const initial=rWalk(q1,1,0,armSteps+1);
          const left=rWalk(q1,-1,0,armSteps+1);
          if(!initial||!left)continue;
          const q0=left.at(-1),q2=initial.at(-1);
          const move21=rWalk(q2,-1,0,armSteps+1);
          if(!move21)continue;
          for(const diag of dirs){
            const diagPath=rWalk(q1,diag[0],diag[1],armSteps+1);
            if(!diagPath)continue;
            const qY=diagPath.at(-1);
            const gapPairs=rGapPairs(q0,q1,q2,qY);
            const gapCounts=gapPairs.map(([a,b])=>rBufferStarsBetween(a,b).length);
            if(gapCounts.some(n=>n<rMinBufferStars))continue;
            const boundaryStars=[q0,q1,q2,qY].map(q=>rBoundaryStarsForEndpoint(q));
            if(boundaryStars.some(stars=>stars.length<rRequiredBoundaryStars))continue;
            // REQUIRE at least one whole complete star between EACH of the
            // four possible endpoints and its nearest open boundary edge.
            // Prefer maximal clearance, then place near the center.
            // stars, breaking ties by position near the center of the patch.
            const pids=[q0,q1,q2,qY];
            const boundaryPenalty=pids.reduce((sum,id)=>{
              const p=rPoint(id),d=Math.min(p.x-Math.min(...xs),Math.max(...xs)-p.x,
                p.y-Math.min(...ys),Math.max(...ys)-p.y);
              return sum+Math.max(0,3-d)**2;
            },0);
            const centerCost=Math.hypot(site.x-xmid,site.y-ymid);
            const arr={S:initial,Y:diagPath.slice(1),M21:move21.slice(1),M0:left.slice(1)};
            candidates.push({cost:centerCost+2*boundaryPenalty,geom:{q0,q1,q2,qY,
              paths:arr,diagnosticDirection:diag,armSteps,gapCounts,
              minGeometricBuffer:Math.min(...gapCounts),
              boundaryStarIds:boundaryStars.map(xs=>xs.map(st=>st.id))}});
          }
        }
        candidates.sort((a,b)=>a.cost-b.cost||a.geom.q1-b.geom.q1);
        // Only compile physically admissible candidates, rather than
        // running an expensive joint-CZ search at every lattice qubit.
        for(const cand of candidates.slice(0,30)){
          const built={};let valid=true;
          for(const [name,path] of Object.entries(cand.geom.paths)){
            const result=rMakeRibbon(path,kind);
            if(!result){valid=false;break;}
            built[name]=result;
          }
          if(!valid)continue;
          return {geom:cand.geom,segments:built};
        }
      }
      return null; // NEVER silently fall back to a too-short junction.
    }
    function rGapStarReport(events,a,b,engine=null,state=null){
      const stars=rBufferStarsBetween(a,b);
      const flux=new Set(rFluxIds(rUnitaryOnly(events)).map(t=>t.id));
      const evaluated=stars.map(st=>{
        const tips=new Set(st.x_qubits);
        const touching=triangleData.filter(t=>t.qubits.some(q=>tips.has(q)));
        const fluxFree=touching.every(t=>!flux.has(t.id));
        const expectation=(engine&&state&&!state.dead)?engine.starExpectation(state,st.id):null;
        const starClear=expectation===null?null:Number.isFinite(expectation)&&Math.abs(expectation-1)<1e-7;
        return {id:st.id,fluxFree,expectation,starClear,
          clean:fluxFree&&(starClear===true),geometric:fluxFree};
      });
      const geometricClean=evaluated.filter(s=>s.geometric).length;
      const quantumClean=evaluated.filter(s=>s.clean).length;
      return {stars:evaluated,geometricCount:stars.length,geometricClean,
        quantumClean,geometryOK:stars.length>=rMinBufferStars,
        fluxOK:geometricClean>=rMinBufferStars,
        quantumOK:engine&&state&&!state.dead?quantumClean>=rMinBufferStars:null};
    }
    function rFluxIds(events) {
      const xcounts=new Set();for(const ev of events)if(ev.type==='X'&&ev.qubits.length===1){const q=ev.qubits[0];if(xcounts.has(q))xcounts.delete(q);else xcounts.add(q);}
      return triangleData.filter(t=>t.qubits.filter(q=>xcounts.has(q)).length%2===1);
    }
    function rCorrectionAuto(move,geom,segments) {
      // M(2→1) has one mandatory interface X at the original position 2.
      // It is inserted separately by rStepEvents, not guessed as an optional correction.
      // Never place q2 in rAdjust.M21.x: that would apply X twice and cancel it.
      if(move==='M21')return [];
      const srcId=geom.q1, src=rPoint(srcId);
      const flux=new Set(rFluxIds([...segments.S.events,...segments[move].events]).map(t=>t.id));
      const matching=[];
      for(const q of siteData) {
        if(Math.hypot(q.x-src.x,q.y-src.y)>2.3)continue;
        const touched=triangleData.filter(t=>t.qubits.includes(q.id));
        const gain=touched.filter(t=>flux.has(t.id)).length*2-touched.length;
        if(gain<=0)continue;
        matching.push({q:q.id,gain,distance:Math.hypot(q.x-src.x,q.y-src.y)});
      }
      matching.sort((a,b)=>b.gain-a.gain||a.distance-b.distance||a.q-b.q);
      return matching.length?[matching[0].q]:[];
    }
    function rStarsAuto(move,geom) {
      const near=rPoint(move==='M21'?geom.q2:geom.q1),count=4; // Four local source-star projectors (parallel and slanted).
      return starData.slice().sort((a,b)=>Math.hypot(a.center[0]-near.x,a.center[1]-near.y)-Math.hypot(b.center[0]-near.x,b.center[1]-near.y)).slice(0,count).map(s=>s.id);
    }
    // An outward destination-local star must have the new endpoint on its
    // X-support, and its center must lie beyond that endpoint in the directed
    // movement. This gives A50 (1->0) and A18 (1->Y) in the default 12x12
    // junction, without hardcoding star IDs on other lattice sizes.
    // Require uniqueness rather than guessing if the endpoint is near a boundary.
    function rOutwardDestinationStar(move,geom) {
      if(move!=='M0'&&move!=='Y')return null;
      const src=rPoint(geom.q1),q=move==='M0'?geom.q0:geom.qY,dest=rPoint(q);
      if(!src||!dest)return null;
      const dx=dest.x-src.x,dy=dest.y-src.y;
      const candidates=starData.filter(st=>st.x_qubits.includes(q)&&
        (st.center[0]-dest.x)*dx+(st.center[1]-dest.y)*dy>1e-7);
      return candidates.length===1?candidates[0].id:null;
    }
    function rUserXForSixCellJunction(geom) {
      // Qubit IDs are lattice-specific.  These two experimentally identified
      // additional corrections apply ONLY to the matching 6 x 6 junction.
      // In this geometry the existing mandatory 2→1 X is at q67, so q68
      // is a SECOND, distinct interface X and must not replace q67.
      const matches=starData.length===16 &&
        geom.q0===58 && geom.q1===63 && geom.q2===67 && geom.qY===29 &&
        siteById.has(61) && siteById.has(68);
      return matches?{Y:[61],M21:[68],M0:[]}:null;
    }
    // Calibrate only the source-interface X correction for M(2 -> 1), using
    // the triangle syndrome.  This never inspects the exchange phase.
    // The old q2 cluster is replaced by the *actual* q1 cluster obtained from
    // the splitting S; the spectator's triangle IDs are kept exactly as they
    // were before the move.  Both Kawagoe-Levin spectator contexts must give
    // the same, uniquely localized additional X operator.
    function rM21FindOneLocalX() {
      const fail=(reason,more={})=>({ok:false,reason,...more});
      const split=rEndpointFluxReport(rStepEvents('A',1),[rGeom.q1,rGeom.q2]);
      if(!split.ok)return fail('Splitting has no unambiguous two-RGB endpoint syndrome.');
      const targetOne=split.groups[0].items.map(t=>t.id);
      const tests=[];
      for(const which of ['A','B']){
        const spectator=which==='A'?rGeom.qY:rGeom.q0;
        const before=rEndpointFluxReport(rStepEvents(which,2),[spectator,rGeom.q2]);
        if(!before.ok)return fail(`Process (${which==='A'?'1':'2'}) lacks a clean spectator/source syndrome before 2→1.`);
        const wanted=new Set([...targetOne,...before.groups[0].items.map(t=>t.id)]);
        const measured=new Set(rFluxIds(rUnitaryOnly(rStepEvents(which,3))).map(t=>t.id));
        const mismatch=triangleData.filter(t=>wanted.has(t.id)!==measured.has(t.id)).map(t=>t.id).sort((a,b)=>a-b);
        tests.push({which,spectator,expected:[...wanted].sort((a,b)=>a-b),measured:[...measured].sort((a,b)=>a-b),mismatch});
      }
      const delta=tests[0].mismatch;
      if(!rSameIds(delta,tests[1].mismatch))return fail('Spectator contexts require different triangle flips; no shared correction inferred.',{tests});
      if(!delta.length)return {ok:true,qubits:[],reason:'No residual interface triangle mismatch.',tests};
      const p=rPoint(rGeom.q2);
      // The local radius is a strict locality constraint, not a search for
      // an overlap phase.  Candidate support must match the complete delta:
      // no unplanned flips on any other triangle are permitted.
      const matches=siteData.filter(q=>q.id!==rGeom.q2 &&
          Math.hypot(q.x-p.x,q.y-p.y)<=2.1+1e-6 &&
          rSameIds(triangleData.filter(t=>t.qubits.includes(q.id)).map(t=>t.id).sort((a,b)=>a-b),delta))
        .map(q=>q.id);
      if(matches.length!==1)return fail(`Found ${matches.length} local single-X solutions for the required syndrome; no automatic correction applied.`,{tests,candidates:matches});
      return {ok:true,qubits:matches,reason:`One local X(q${matches[0]}) flips precisely B${delta.join(' and B')} in both contexts.`,tests,candidates:matches};
    }
    function rSetup(kind) {
      rSelectedKind=kind;rExactCache=null;rLastSymbolicReport=null;rM21SyndromeSuggestion=null;rCalcToken++;rValidationToken++;rCheckpoint='projected';rCheckpointSelect.value=rCheckpoint;
      rValidationContent.textContent='Choose one of the movement buttons to inspect its checkpoints.';
      const model=rFindGeometry(kind);
      rGeom=model?.geom||null;rSegments=model?.segments||null;
      rEndpointFix=rGeom?{Y:rOutwardDestinationStar('Y',rGeom),M0:rOutwardDestinationStar('M0',rGeom)}:{Y:null,M0:null};
      const splitCalibration=rSplittingCalibration(rGeom);
      rInitialSplittingStars=splitCalibration.split;
      rExtraM21Star=splitCalibration.extraM21;
      rAdjust={Y:{x:[],stars:[]},M21:{x:[],stars:[]},M0:{x:[],stars:[]}};
      rCompleted={A:0,B:0};rActiveProcess='A';rActiveStep=0;
      rStudio.hidden=!kind;setTimeout(rResize,120);
      if(!kind)return;
      if(!model){rOverview.textContent='No six-step T-junction has one complete star in each endpoint-to-boundary strip and one between each endpoint pair. Increase the lattice size; the workbench will not silently shorten the arms or waive the clearance test.';rCalc.textContent='No R result: a complete junction could not be constructed.';rUpdateButtons();return;}
      const userX=rUserXForSixCellJunction(rGeom);
      // Initialize both possible first movements before checking M(2→1) in
      // its TWO spectator contexts. This ensures we don't fit only one order.
      for(const move of ['Y','M0','M21']){
        const automatic=rCorrectionAuto(move,rGeom,rSegments);
        const userAdditions=userX?userX[move]:[];
        rAdjust[move]={
          x:[...new Set([...automatic,...userAdditions])].filter(q=>!(move==='M21'&&q===rGeom.q2)),
          stars:rStarsAuto(move,rGeom)
        };
      }
      // Preserve the manually calibrated 6x6 correction unchanged.  On
      // other admissible geometries, add a uniquely syndrome-verified local X.
      if(!userX){
        rM21SyndromeSuggestion=rM21FindOneLocalX();
        if(rM21SyndromeSuggestion.ok)
          rAdjust.M21.x=[...new Set([...rAdjust.M21.x,...rM21SyndromeSuggestion.qubits])];
      }
      const userInfo=userX
        ? '<strong>6×6 corrections from user inspection:</strong> 1→Y also applies X(q61); 2→1 also applies X(q68), in addition to mandatory X(q67). These are applied in BOTH braiding processes for s1 and s2. '
        : '6×6-specific extra X locations have not been established for this geometry; no q61/q68 overrides applied. ';
      const sourceNote=rM21SyndromeSuggestion?
        `<br><strong>2→1 triangle-only calibration:</strong> ${rEscape(rM21SyndromeSuggestion.reason)} ${rM21SyndromeSuggestion.ok?'Verified against both spectator contexts; this is not a fusion-state validation.':'No correction automatically inserted.'}`:'';
      const splitProjInfo=rInitialSplittingStars.length
        ? `<strong>Splitting preparation:</strong> after the fresh S(1,2) ribbon, apply ${rInitialSplittingStars.map(id=>'P'+id).join(' then ')} (each successful branch normalized). `
        : '<strong>Splitting preparation:</strong> the P53/P48 star IDs are calibrated only for the matching default 12×12 T-junction; none are inserted for this geometry. ';
      const extraMoveInfo=rExtraM21Star===null?'':`<strong>2→1 additional local projection:</strong> fixed P${rExtraM21Star} after its four source-interface projectors, in both processes. `;
      const endProjInfo=`${splitProjInfo}${extraMoveInfo}<strong>Fixed destination-local +1 projections:</strong> M(1→0): ${rEndpointFix.M0===null?'unavailable':'P'+rEndpointFix.M0}; M(1→Y): ${rEndpointFix.Y===null?'unavailable':'P'+rEndpointFix.Y}. These occur AFTER the four source-interface stars and are identical in both process orders. They define chosen endpoint representatives, not a proven anyon transport isometry. `;
      rOverview.innerHTML=`${endProjInfo}T-junction: exactly ${rGeom.armSteps} nearest-neighbor steps per arm; at least ${rGeom.minGeometricBuffer} complete stars between each possible pair of endpoints, plus at least one COMPLETE outward star between each endpoint and the physical open boundary (q0/q1/q2/qY: ${rGeom.boundaryStarIds.map(xs=>xs.map(id=>'A'+id).join(',')).join(' / ')}). Geometric clearance alone does not guarantee Bₜ=+1 or Aₛ=+1; all finished stages are audited below. Fixed junction: 0=q${rGeom.q0}, 1=q${rGeom.q1}, 2=q${rGeom.q2}, Y=q${rGeom.qY}. Each fresh segment begins at the next qubit. <strong>M(2→1) always applies exactly one X at position 2 (q${rGeom.q2}), after the fresh ribbon and before the projectors.</strong> ${userInfo}${sourceNote} Each parallel movement now projects four complete stars around its old endpoint. Other suggested X corrections and the selected star projectors remain editable hypotheses, not independently certified.`;
      rCalc.innerHTML='<p>Complete both four-step processes to compare their final quantum states. The interface corrections must be validated first.</p>';
      rUpdateButtons();rRenderStep();setTimeout(rResize,350);
    }
    function rStepEvents(which,step) {
      const events=[];
      for(let j=0;j<step;j++){
        const move=rSequences[which][j],seg=rSegments[move];
        for(const ev of seg.events)events.push({...ev,step:j+1,move});
        if(move==='S'){
          // Chronological splitting preparation: S -> P53 -> P48 for the
          // default calibrated 12x12 junction.  The same S is shared by both.
          for(const starId of rInitialSplittingStars)
            events.push({type:'P',starId,source:'splitting-local',step:j+1,move});
        }
        if(move!=='S') {
          // User-specified microscopic rule for M(2→1): the interface is at
          // the OLD position 2. Exactly one X there, AFTER the fresh s string
          // and BEFORE any interface star projection (no duplicate X).
          if(move==='M21')events.push({type:'X',qubits:[rGeom.q2],source:'interface-fixed',step:j+1,move});
          for(const q of rAdjust[move].x)events.push({type:'X',qubits:[q],source:'interface',step:j+1,move});
          for(const starId of rAdjust[move].stars)events.push({type:'P',starId,source:'interface',step:j+1,move});
          // The separately requested P53 at M(2→1) follows the four source
          // projectors; it does not replace P48 and it must not be duplicated.
          if(move==='M21'&&rExtraM21Star!==null)
            events.push({type:'P',starId:rExtraM21Star,source:'m21-local',step:j+1,move});
          // Select the same local endpoint representative regardless of
          // whether this movement happens first or last in the T-junction.
          const targetStar=rEndpointFix[move];
          if(targetStar!==null&&targetStar!==undefined)
            events.push({type:'P',starId:targetStar,source:'endpoint-local',step:j+1,move});
        }
      }
      return events;
    }
    function rUnitaryOnly(events){return events.filter(e=>e.type!=='P');}
    function rStageEvents(which,step,checkpoint='projected'){
      if(step<1)return [];
      if(checkpoint==='projected')return rStepEvents(which,step);
      const previous=rStepEvents(which,step-1);
      const current=rStepEvents(which,step).slice(previous.length);
      if(checkpoint==='before')return previous;
      if(checkpoint==='ribbon')return [...previous,...current.filter(e=>e.source==='r-ribbon')];
      if(checkpoint==='x')return [...previous,...current.filter(e=>e.type!=='P')];
      return [...previous,...current];
    }
    function rExpectedPair(which,step){
      const endpoints=[rGeom.q1,rGeom.q2];
      for(let k=1;k<step-1;k++){
        const move=rSequences[which][k];
        const old=move==='M21'?rGeom.q2:rGeom.q1;
        const dest=move==='M21'?rGeom.q1:move==='Y'?rGeom.qY:rGeom.q0;
        const index=endpoints.indexOf(old);
        if(index>=0)endpoints[index]=dest;
      }
      return endpoints;
    }
    function rMoveEndpoints(which,step){
      const before=rExpectedPair(which,step),move=rSequences[which][step-1];
      const source=move==='M21'?rGeom.q2:rGeom.q1;
      const destination=move==='M21'?rGeom.q1:move==='Y'?rGeom.qY:rGeom.q0;
      const after=before.slice(),i=after.indexOf(source);
      if(i>=0)after[i]=destination;
      const spectator=before.find(q=>q!==source);
      return {source,destination,spectator,before,after};
    }
    // The RGB endpoint-locality radius is a GEOMETRIC diagnostic, not a
    // microscopic projector or a proof that the anyon charge is semionic.
    function rEndpointFluxReport(events,endpoints){
      const flux=rFluxIds(rUnitaryOnly(events));
      const groups=endpoints.map(q=>({q,items:[],maxDist:0}));
      for(const t of flux){
        const dd=endpoints.map(q=>Math.hypot(t.center[0]-rPoint(q).x,t.center[1]-rPoint(q).y));
        const j=dd[0]<=dd[1]?0:1;groups[j].items.push(t);groups[j].maxDist=Math.max(groups[j].maxDist,dd[j]);
      }
      const countRGB=groups.every(g=>g.items.length===3&&new Set(g.items.map(t=>t.color)).size===3);
      const withinRadius=groups.every(g=>g.maxDist<=1.85);
      return {flux,groups,countRGB,withinRadius,ok:flux.length===6&&countRGB&&withinRadius};
    }
    function rNearInterfaceTriangles(events,site,radius=1.15){
      const p=rPoint(site);return rFluxIds(rUnitaryOnly(events)).filter(t=>Math.hypot(t.center[0]-p.x,t.center[1]-p.y)<radius);
    }
    function rSpectatorFlux(events,site,radius=1.85){
      const p=rPoint(site);return rFluxIds(rUnitaryOnly(events)).filter(t=>Math.hypot(t.center[0]-p.x,t.center[1]-p.y)<=radius).map(t=>t.id).sort((a,b)=>a-b);
    }
    function rSameIds(a,b){return a.length===b.length&&a.every((v,i)=>v===b[i]);}

    function rTriangleShapes(events) {
      const triangles=rFluxIds(events),shapes=[];
      for(const t of triangles){const path=polygonPath(t.qubits);if(!path)continue;shapes.push({type:'path',path,fillcolor:'rgba(248,113,113,.14)',line:{color:'#ef4444',width:1.7,dash:'dot'},layer:'below'});}
      return shapes;
    }
    function rRenderStep() {
      if(!rSelectedKind||!rGeom)return;
      const active=rSequences[rActiveProcess][Math.max(0,rActiveStep-1)];
      const events=rStageEvents(rActiveProcess,rActiveStep,rCheckpoint),visibleEvents=rUnitaryOnly(events);
      const outlines=rTriangleShapes(visibleEvents);
      const pp=rGeom;
      const points=[['0',pp.q0],['1',pp.q1],['2',pp.q2],['Y',pp.qY]];
      const traceEdges={type:'scatter',mode:'lines',x:baseEdgeX,y:baseEdgeY,line:{width:1,color:'#94a3b8'},hoverinfo:'skip',name:'lattice',showlegend:false};
      const traceSites={type:'scatter',mode:'markers',x:siteData.map(s=>s.x),y:siteData.map(s=>s.y),marker:{size:7,color:siteData.map(s=>s.color),line:{width:.5,color:'#334155'}},customdata:siteData.map(s=>s.id),text:siteData.map(s=>'Qubit '+s.id),hovertemplate:'%{text}<extra></extra>',name:'qubits',showlegend:false};
      const xTrace={type:'scatter',mode:'markers+text',x:[],y:[],text:[],marker:{size:12,symbol:'circle-open',line:{width:2,color:'#047857'}},textposition:'top center',textfont:{size:12,color:'#0f172a'},hovertext:[],hoverinfo:'text',name:'applied single-qubit gates',showlegend:false};
      const flags=new Map();
      for(const ev of visibleEvents)if(ev.qubits.length===1){const q=ev.qubits[0];if(!flags.has(q))flags.set(q,[]);flags.get(q).push(ev.type);}
      for(const [q,ops] of flags){const s=rPoint(q);xTrace.x.push(s.x);xTrace.y.push(s.y);xTrace.text.push(ops.length>3?ops.at(-2)+ops.at(-1)+'×'+Math.floor(ops.length/2):ops.join(''));xTrace.hovertext.push('Q'+q+': '+ops.join(' → '));}
      const czTraces=[];let ci=0;
      for(const ev of visibleEvents)if(ev.type==='CZ'){
        const [a,b]=ev.qubits,s1=rPoint(a),s2=rPoint(b),cx=(s1.x+s2.x)/2,cy=(s1.y+s2.y)/2+.22+.07*(ci%5);ci++;
        const x=[],y=[];for(let k=0;k<=15;k++){let t=k/15;x.push((1-t)**2*s1.x+2*(1-t)*t*cx+t*t*s2.x);y.push((1-t)**2*s1.y+2*(1-t)*t*cy+t*t*s2.y);}
        czTraces.push({type:'scatter',mode:'lines',x,y,line:{color:'#7c3aed',width:1.5},hoverinfo:'text',text:`CZ(Q${a},Q${b})`,showlegend:false});
      }
      const projectorShapes=[],projectorTrace={type:'scatter',mode:'markers+text',x:[],y:[],text:[],textposition:'middle center',textfont:{size:10,color:'#0369a1'},marker:{size:26,color:'rgba(56,189,248,.15)',line:{width:1,color:'#0284c7'}},hovertext:[],hoverinfo:'text',showlegend:false};
      for(const ev of events)if(ev.type==='P'){
        const st=starById.get(ev.starId);if(!st)continue;
        const path=polygonPath(st.inner_qubits);if(path)projectorShapes.push({type:'path',path,fillcolor:'rgba(14,165,233,.10)',line:{width:1.5,color:'#0284c7',dash:'dot'},layer:'below'});
        projectorTrace.x.push(st.center[0]);projectorTrace.y.push(st.center[1]);
        projectorTrace.text.push('P'+st.id+(['endpoint-local','splitting-local','m21-local'].includes(ev.source)?'*':''));
        const meaning=ev.source==='endpoint-local'?'Fixed destination-local ':ev.source==='splitting-local'?'Fixed splitting-state ':ev.source==='m21-local'?'Fixed added M(2→1) ': 'Source-interface ';
        projectorTrace.hovertext.push(meaning+'projector (I+A_'+st.id+')/2');
      }
      const allStarCenters={type:'scatter',mode:'markers',x:starData.map(s=>s.center[0]),y:starData.map(s=>s.center[1]),marker:{size:10,color:'#0284c7',opacity:.28,symbol:'diamond'},customdata:starData.map(s=>s.id),hovertemplate:'Star A%{customdata}<extra></extra>',name:'r-star-centers',showlegend:false};
      const shownEndpoints=rActiveStep===1?[rGeom.q1,rGeom.q2]:rMoveEndpoints(rActiveProcess,rActiveStep).after;
      const gapStars=rBufferStarsBetween(shownEndpoints[0],shownEndpoints[1]);
      const gapTrace={type:'scatter',mode:'markers',x:gapStars.map(s=>s.center[0]),y:gapStars.map(s=>s.center[1]),
        marker:{size:17,color:'#16a34a',opacity:.6,symbol:'diamond-open',line:{width:2,color:'#15803d'}},
        customdata:gapStars.map(s=>s.id),hovertemplate:'Buffer star A%{customdata}: inspect excitation validation<extra></extra>',
        name:'geometric-buffer-stars',showlegend:false};
      const boundaryShown=shownEndpoints.flatMap(q=>rBoundaryStarsForEndpoint(q));
      const boundaryUnique=[...new Map(boundaryShown.map(st=>[st.id,st])).values()];
      const boundaryTrace={type:'scatter',mode:'markers',x:boundaryUnique.map(s=>s.center[0]),y:boundaryUnique.map(s=>s.center[1]),
        marker:{size:19,color:'#f59e0b',symbol:'square-open',line:{width:2,color:'#b45309'}},
        customdata:boundaryUnique.map(s=>s.id),hovertemplate:'Outward boundary-buffer star A%{customdata}<extra></extra>',
        name:'geometric-boundary-stars',showlegend:false};
      let violatedNow=[];
      let violatedStarShapes=[];
      let violatedStarTrace={type:'scatter',mode:'markers+text',x:[],y:[],text:[],textposition:'top center',
        textfont:{size:12,color:'#991b1b'},marker:{size:18,symbol:'diamond-open',color:'rgba(239,68,68,.08)',line:{width:2.3,color:'#dc2626'}},
        hovertext:[],hoverinfo:'text',name:'violated-stars',showlegend:false};
      try{
        if(!rExactCache)rExactCache=new RSymbolicEngine(siteData,starData,triangleData);
        const stageState=rExactCache.run(events);
        if(stageState && !stageState.dead){
          violatedNow=starData
            .map(st=>({star:st,value:rExactCache.starExpectation(stageState,st.id)}))
            .filter(item=>Number.isFinite(item.value)&&Math.abs(item.value-1)>1e-7);
          for(const item of violatedNow){
            const path=polygonPath(item.star.inner_qubits);
            if(path)violatedStarShapes.push({type:'path',path,fillcolor:'rgba(248,113,113,.08)',line:{width:2,color:'#dc2626',dash:'dot'},layer:'below'});
            violatedStarTrace.x.push(item.star.center[0]);
            violatedStarTrace.y.push(item.star.center[1]);
            violatedStarTrace.text.push('A'+item.star.id);
            violatedStarTrace.hovertext.push(`Violated star A${item.star.id}: ⟨A⟩ = ${item.value.toFixed(6)}`);
          }
        }
      }catch(_err){}
      const labels={type:'scatter',mode:'markers+text',x:points.map(p=>rPoint(p[1]).x),y:points.map(p=>rPoint(p[1]).y),text:points.map(p=>p[0]),textposition:'bottom center',textfont:{size:15,color:'#0f172a'},marker:{size:17,symbol:'circle-open',line:{width:2,color:'#0f172a'}},hoverinfo:'skip',showlegend:false};
      const xvals=siteData.map(s=>s.x),yvals=siteData.map(s=>s.y);
      const layout={margin:{l:26,r:18,t:26,b:24},height:530,showlegend:false,hovermode:'closest',uirevision:'r-study',paper_bgcolor:'#fff',plot_bgcolor:'#fff',xaxis:{visible:false,range:[Math.min(...xvals)-.8,Math.max(...xvals)+.8]},yaxis:{visible:false,scaleanchor:'x',scaleratio:1,range:[Math.min(...yvals)-.8,Math.max(...yvals)+.8]},shapes:[...outlines,...projectorShapes,...violatedStarShapes]};
      Promise.resolve(Plotly.react(rPlot,[traceEdges,traceSites,...czTraces,xTrace,allStarCenters,gapTrace,boundaryTrace,projectorTrace,violatedStarTrace,labels],layout,{displaylogo:false,responsive:true,scrollZoom:true})).then(()=>{
        if(!rPlot._rClickBound && typeof rPlot.on==='function') {
          rPlot.on('plotly_click',rHandlePlotClick);rPlot._rClickBound=true;
        }
      }).catch(console.error);
      if(rActiveStep>1){
        const move=rSequences[rActiveProcess][rActiveStep-1];
        const targetStar=rEndpointFix[move];
        const targetText=targetStar===null||targetStar===undefined?'':` Fixed destination-local P${targetStar} is applied AFTER the four editable source-interface projectors in both process orders. It is shown with * on the canvas and is not part of the editable source-star list.`;
        const extraText=move==='M21'&&rExtraM21Star!==null?` Fixed P${rExtraM21Star} is applied once after the four editable interface stars.`:'';
        rFixedXNote.textContent=(move==='M21'
          ? `Mandatory interface operation: X at position 2 = qubit ${rGeom.q2} (applied exactly once, after the fresh 2→1 segment). Do NOT enter this qubit again below.`
          : 'There is no mandatory source X for this movement. Additional interface X corrections are editable below.')+extraText+targetText;
        rEditX.value=rAdjust[move].x.join(', ');
        rEditStars.value=rAdjust[move].stars.join(', ');
        rSaveAdjust.disabled=false;
      }else{
        rFixedXNote.textContent='';rEditX.value='';rEditStars.value='';rSaveAdjust.disabled=true;
      }
      rUpdateStatus(events,active);rValidateSelected();setTimeout(rResize,220);
    }
    function rUpdateButtons(){for(const b of rStepButtons){const seq=b.dataset.rProcess,st=Number(b.dataset.rStep);b.classList.toggle('active',rActiveProcess===seq&&rActiveStep===st);b.classList.toggle('r-done',rCompleted[seq]>=st);}}
    function rUpdateStatus(events,active){
      const counts={X:0,R:0,CZ:0,P:0};events.forEach(e=>counts[e.type]++);
      const fluxes=rFluxIds(rUnitaryOnly(events));
      const label=rActiveStep?`${rActiveProcess} · Step ${rActiveStep}: ${rDescriptions[active]} · checkpoint: ${rCheckpoint}`:'Choose one of the process buttons.';
      const value=`${counts.X} X · ${counts.R} R · ${counts.CZ} CZ · ${counts.P} star projectors; ${fluxes.length} excited triangles before projectors`;
      let violatedStage='';
      try{
        if(!rExactCache)rExactCache=new RSymbolicEngine(siteData,starData,triangleData);
        const stageState=rExactCache.run(events);
        if(stageState && !stageState.dead){
          const ids=starData.filter(st=>{const v=rExactCache.starExpectation(stageState,st.id);return Number.isFinite(v)&&Math.abs(v-1)>1e-7;}).map(st=>'A'+st.id);
          violatedStage=ids.length?`Violated stars shown on the canvas: ${ids.join(', ')}`:'No violated complete stars at this checkpoint.';
        }else if(stageState && stageState.dead){
          violatedStage='This checkpoint has zero postselection weight, so violated-star readout is not shown.';
        }
      }catch(_err){}
      const diff=rActiveStep>1?
        `${active==='M21'?`Fixed X at position 2: q${rGeom.q2}; `:''}Additional interface X: ${rAdjust[active].x.length?rAdjust[active].x.join(', '):'none'}; source-interface stars: ${rAdjust[active].stars.map(id=>'A'+id).join(', ')||'none'}; ${active==='M21'&&rExtraM21Star!==null?'fixed extra P'+rExtraM21Star+'; ':''}${rEndpointFix[active]!=null?'fixed endpoint-local P'+rEndpointFix[active]+' (after source projections); ':''}${violatedStage?`; ${violatedStage}`:''}`:(`Splitting preparation: ${rInitialSplittingStars.length?rInitialSplittingStars.map(id=>'P'+id).join(' → '):'no fixed splitting projections on this geometry'}; ${violatedStage}`);
      document.getElementById('r-stage-caption').innerHTML=`<strong>${rEscape(label)}</strong><div>${rEscape(value)}</div><div>${rEscape(diff)}</div>`;
      if(rCompleted.A===4 && rCompleted.B===4)rComputeR();
      else rCalc.innerHTML=`<p><strong>Not evaluated yet.</strong> Complete both four-stage processes. For each move, examine the suggested interface X and projector star IDs and correct them as needed.</p><div class="r-math">R<sub>candidate</sub> = &lt;Ψ<sub>2</sub>|Ψ<sub>1</sub>&gt; / (‖Ψ<sub>1</sub>‖‖Ψ<sub>2</sub>‖)</div>`;
    }
    function rFormatComplex(z){return `${z.re.toFixed(6)} ${z.im<0?'−':'+'} ${Math.abs(z.im).toFixed(6)} i`;}
    function rWholeStateChecks(events,expected0,expectedY){
      const tris=rFluxIds(rUnitaryOnly(events));
      const colored={L:[],Y:[]};
      const p0=rPoint(expected0),py=rPoint(expectedY);
      for(const t of tris){let d0=Math.hypot(t.center[0]-p0.x,t.center[1]-p0.y),dy=Math.hypot(t.center[0]-py.x,t.center[1]-py.y);colored[d0<=dy?'L':'Y'].push(t);}
      const ok=Object.values(colored).every(arr=>arr.length===3 && new Set(arr.map(t=>t.color)).size===3);
      return {ok,fluxIds:tris.map(t=>t.id),colored};
    }
    let rLastSymbolicReport=null;
    function rWeightGroups(state,which) {
      return rSequences[which].map((move,k)=>{
        const terms=state.projectionLog.filter(x=>x.step===k+1);
        const zero=terms.some(t=>t.p===0);
        const logp=zero?-Infinity:terms.reduce((n,t)=>n+(t.log2p??0),0);
        return {step:k+1,move,stars:terms.map(t=>t.starId),p:zero?0:2**logp,w:zero?null:2**(-logp/2),log2p:zero?null:logp};
      });
    }
    function rNumber(v) {
      if(v===null||!Number.isFinite(v))return 'undefined';
      if(v===0)return '0';if(v===1)return '1';
      return (Math.abs(v)<1e-5||Math.abs(v)>1e6)?v.toExponential(6):Number(v.toPrecision(8)).toString();
    }
    function rExactForm(result) {
      if(result.zero)return '0';
      return `2<sup>${result.power}</sup> exp(iπ × ${result.phase8}/4)`;
    }
    function rWeightTable(state,which) {
      const groups=rWeightGroups(state,which);
      let html=`<details open><summary><strong>Process (${which==='A'?1:2}): projection normalization</strong></summary>`;
      html+='<table class="r-symbolic-table"><thead><tr><th>Step</th><th>p = ⟨P†P⟩</th><th>w = 1/√p</th></tr></thead><tbody>';
      for(const g of groups)html+=`<tr><td>${g.step}. ${rEscape(rDescriptions[g.move])}<br><small>${g.stars.map(id=>'A'+id).join(', ')||'No projections'}</small></td><td>${rNumber(g.p)}</td><td>${rNumber(g.w)}</td></tr>`;
      html+='</tbody></table><p class="r-caveat">p is conditional on all earlier successful, normalized steps. Each movement uses the product of its per-projector weights; no independence of outcomes is assumed.</p></details>';
      html+='<details><summary>Individual star probabilities and weights</summary><table class="r-symbolic-table"><thead><tr><th>Star / step</th><th>p</th><th>w</th></tr></thead><tbody>';
      for(const t of state.projectionLog)html+=`<tr><td>A${t.starId} / ${t.step??'–'}</td><td>${rNumber(t.p)}</td><td>${rNumber(t.w)}</td></tr>`;
      html+='</tbody></table></details>';return html;
    }
    function rSymbolicMethodHTML(engine,contraction) {
      const meta=engine.metadata();
      return `<details><summary><strong>How the symbolic calculation works</strong></summary>
        <p>Evaluates ⟨Ω|K₂†K₁|Ω⟩ and projector norms by exact operator/phase-polynomial contraction. Pauli parities and CZ phases are combined algebraically. Projection branches remain binary variables in a compact quadratic polynomial; they are never expanded into an amplitude vector.</p>
        <div class="r-math">0 amplitude arrays · ${meta.physicalQubits} physical qubits · ${meta.completeStars} independent complete stars<br>Last overlap: ${contraction?.variables??0} binary variables; ${contraction?.rank??0} matching constraints; ${contraction?.freeVariables??0} free variables eliminated algebraically.</div>
        <p>All complex phases, including global R-gate phases and the signs from reordering X/CZ/star operations, are retained. A zero postselection branch has no normalization and no R value.</p>
        <p><strong>Open-boundary reference:</strong> Ω is the same normalized star-generated state from |0…0⟩ as in the earlier amplitude engine. Products of Bₜ have expectation +1, and an anticommuting star proves zero. A residual outside the Bₜ span is <em>not automatically zero</em>: this reference fixes ${meta.extraDiagonalConstraints} additional diagonal degrees of freedom beyond the ${meta.triangleRank} independent triangle constraints.</p>
        <label for="r-z-query">Check a residual Z string on Ω (qubit IDs)</label>
        <input id="r-z-query" placeholder="e.g. 99, 100, 101" aria-label="Residual Z qubit IDs"/>
        <button class="btn" id="r-z-test">Check Z residual</button><div id="r-z-answer" aria-live="polite"></div>
        <p class="r-caveat">This compresses the user's residual-Z summation rather than listing exponentially many terms. It supports the displayed X/Z/R/CZ/star-projector gate set and independent open-patch star basis; it is not a generic polynomial-time simulator for arbitrary quantum gates.</p>
      </details>`;
    }
    function rAllStageBufferAudit(engine){
      const rows=[];
      for(const which of ['A','B'])for(let step=1;step<=4;step++){
        const events=rStepEvents(which,step);
        const endpoints=step===1?[rGeom.q1,rGeom.q2]:rMoveEndpoints(which,step).after;
        const state=engine.run(events);
        const gap=rGapStarReport(events,endpoints[0],endpoints[1],engine,state);
        const boundary=rBoundaryReport(events,endpoints,engine,state);
        const flux=rEndpointFluxReport(events,endpoints);
        rows.push({which,step,endpoints,geometricStars:gap.geometricCount,
          cleanStars:gap.quantumClean,cleanIDs:gap.stars.filter(x=>x.clean).map(x=>x.id),
          separationOK:gap.quantumOK===true,boundary, boundaryOK:boundary.quantumOK===true,
          fluxCount:flux.flux.length,endpointsOK:flux.ok,postselectionAlive:!state.dead});
      }
      return rows;
    }
    function rAllStageAuditHTML(rows){
      let html='<details open><summary><strong>Complete endpoint/boundary buffer audit (all eight stages)</strong></summary>';
      html+='<table class="r-symbolic-table"><thead><tr><th>Process / step</th><th>Clean inter-endpoint stars</th><th>Clean boundary stars</th><th>Endpoint flux</th></tr></thead><tbody>';
      for(const x of rows){
        const move=rSequences[x.which][x.step-1];
        html+=`<tr><td>(${x.which==='A'?'1':'2'}) · ${x.step}: ${rEscape(rDescriptions[move])}</td>`+
          `<td>${x.cleanStars}/${x.geometricStars} ${x.separationOK?'✓':'✗'}<br><small>${x.cleanIDs.map(id=>'A'+id).join(', ')||'none'}</small></td>`+
          `<td>${x.boundary.perEndpoint.map(o=>'q'+o.endpoint+': '+o.quantumCleanCount+'/'+o.geometricCount+(o.quantumOK?' ✓':' ✗')).join('<br>')}</td>`+
          `<td>${x.fluxCount} Bₜ=−1 ${x.endpointsOK?'✓':'✗'}${x.postselectionAlive?'':' (zero branch)'}</td></tr>`;
      }
      html+='</tbody></table><p class="r-caveat">A clean buffer star is a complete Aₛ with ⟨Aₛ⟩=+1 and no excited Bₜ sharing its six X-tip qubits. At least one clean corridor star is required between each actual pair of endpoints, AND one complete clean star is required on the outward boundary side of each occupied endpoint, at every completed stage. This six-step geometry cannot meet the previous two-star diagonal corridor requirement with the same 2.5-unit endpoint exclusion. Passing the gap check alone does not guarantee that all endpoint fluxes or fusion data are correct.</p></details>';
      return html;
    }
    function rComputeR(){
      const A=rStepEvents('A',4),B=rStepEvents('B',4),cA=rWholeStateChecks(A,rGeom.q0,rGeom.qY),cB=rWholeStateChecks(B,rGeom.q0,rGeom.qY);
      let section=`<p><strong>Two protocol outputs:</strong> (1): ${cA.fluxIds.length} excited triangles; (2): ${cB.fluxIds.length}. ${cA.ok&&cB.ok?'Both have three RGB triangles near each endpoint.':'The final triangle pattern is not two clean RGB endpoints.'}</p>`;
      if(!cA.ok||!cB.ok)section+='<p class="r-error">Not a calibrated exchange: the final flux syndrome does not match two RGB anyons at 0 and Y.</p>';
      rLastSymbolicReport=null;
      rCalc.innerHTML=section+'<p>Contracting operators and normalizing projector branches…</p>';
      const token=++rCalcToken;
      setTimeout(()=>{
        if(token!==rCalcToken)return;
        try {
          const start=performance.now();
          if(!rExactCache)rExactCache=new RSymbolicEngine(siteData,starData,triangleData);
          const first=rExactCache.run(A),second=rExactCache.run(B);
          section+='<p class="r-symbolic-badge">Phase-sensitive symbolic calculation · no amplitude enumeration</p>';
          if(first.dead||second.dead){
            section+='<p class="r-error">R is undefined: '+[first.dead?'process (1)':null,second.dead?'process (2)':null].filter(Boolean).join(' and ')+' has a zero-probability projection. A factor w cannot normalize the zero vector.</p>';
            section+=rWeightTable(first,'A')+rWeightTable(second,'B');
            rLastSymbolicReport={status:'zero-postselection',engine:rExactCache.metadata(),kind:rSelectedKind,geometry:rGeom,adjustments:rAdjust,splittingProjectors:rInitialSplittingStars,extraM21Projector:rExtraM21Star,destinationLocalProjectors:rEndpointFix,eventsA:A,eventsB:B,projectorsA:first.projectionLog,projectorsB:second.projectionLog};
            section+='<button class="btn" id="r-export-symbolic">Export calculation JSON</button>';
            rCalc.innerHTML=section;setTimeout(rResize,120);return;
          }
          // Both states are already normalized after every successful P.
          const exact=rExactCache.innerExact(first,second,true),w=rSymComplex(exact);
          const n1=rExactCache.norm(first),n2=rExactCache.norm(second);
          const rawN1=rExactCache.rawNorm(first),rawN2=rExactCache.rawNorm(second);
          const rawExact={...exact,power:exact.power+(first.logSuccess+second.logSuccess)/2};
          const z=rSymComplex(rawExact),mag=Math.hypot(w.re,w.im),phase=mag>1e-10?Math.atan2(w.im,w.re)/Math.PI:null;
          const checked=[...new Set([...A,...B].filter(ev=>ev.type==='P').map(ev=>ev.starId))];
          const starA=checked.map(id=>({id,value:rExactCache.starExpectation(first,id)}));
          const starB=checked.map(id=>({id,value:rExactCache.starExpectation(second,id)}));
          const chosenEndpointStars=[...new Set([rEndpointFix.M0,rEndpointFix.Y].filter(id=>id!==null&&id!==undefined))];
          const endpointStarChecks=chosenEndpointStars.map(id=>({id,
            process1:rExactCache.starExpectation(first,id),
            process2:rExactCache.starExpectation(second,id)}));
          const allSingleStarMismatch=starData.map(st=>({id:st.id,
            process1:rExactCache.starExpectation(first,st.id),process2:rExactCache.starExpectation(second,st.id)}))
            .filter(v=>!Number.isFinite(v.process1)||!Number.isFinite(v.process2)||Math.abs(v.process1-v.process2)>1e-8);
          const brokenA=starA.filter(s=>!Number.isFinite(s.value)||Math.abs(s.value-1)>1e-7).map(s=>s.id);
          const brokenB=starB.filter(s=>!Number.isFinite(s.value)||Math.abs(s.value-1)>1e-7).map(s=>s.id);
          // Extra local-star audit at the *vacated* source q2.  This catches
          // stars that are not among the three heuristic M(2→1) projectors.
          // An endpoint gap can be clean even if the old endpoint is not.
          const srcPoint=rPoint(rGeom.q2);
          const vacatedStars=starData.filter(st=>Math.hypot(
            st.center[0]-srcPoint.x,st.center[1]-srcPoint.y)<=2.1+1e-6);
          const vacatedA=vacatedStars.map(st=>({id:st.id,value:rExactCache.starExpectation(first,st.id)}));
          const vacatedB=vacatedStars.map(st=>({id:st.id,value:rExactCache.starExpectation(second,st.id)}));
          const unrestoredA=vacatedA.filter(st=>!Number.isFinite(st.value)||Math.abs(st.value-1)>1e-7);
          const unrestoredB=vacatedB.filter(st=>!Number.isFinite(st.value)||Math.abs(st.value-1)>1e-7);
          const gapA=rGapStarReport(A,rGeom.q0,rGeom.qY,rExactCache,first);
          const gapB=rGapStarReport(B,rGeom.q0,rGeom.qY,rExactCache,second);
          const stageAudit=rAllStageBufferAudit(rExactCache);
          const allGapsOK=stageAudit.every(x=>x.separationOK&&x.boundaryOK&&x.postselectionAlive);
          const allEndpointsOK=stageAudit.every(x=>x.endpointsOK);
          const phaseOnly=cA.ok&&cB.ok&&rSameIds(cA.fluxIds.slice().sort((a,b)=>a-b),cB.fluxIds.slice().sort((a,b)=>a-b))&&Math.abs(mag-1)<1e-9;
          const passes=phaseOnly&&allGapsOK&&allEndpointsOK&&gapA.quantumOK&&gapB.quantumOK&&brokenA.length===0&&brokenB.length===0&&unrestoredA.length===0&&unrestoredB.length===0;
          const phaseName=phaseOnly&&Math.abs(w.re)<1e-8&&Math.abs(Math.abs(w.im)-1)<1e-8?(w.im>0?'+i':'−i'):'unit-magnitude complex phase';
          const status=phaseOnly
            ?`Phase-only final states: ${phaseName} in this finite-patch calculation. ${passes?'Local audit passes;':'Some strict boundary/intermediate-stage checks remain unresolved;'} Fusion-sector identification and transport isometry are not yet independently certified.`
            :'INCONCLUSIVE — final states are not related by a unit-magnitude phase';
          section+=`<div class="r-math">Normalized ⟨ψ₂|ψ₁⟩ = <strong>${rFormatComplex(w)}</strong><br>Exact form: ${rExactForm(exact)}<br>|overlap| = ${mag.toFixed(8)}${phase===null?'':`; arg / π = ${phase.toFixed(6)}`}<br>Final normalized norms: ${rNumber(n1)}, ${rNumber(n2)}</div>`;
          section+=`<p><strong>${status}</strong></p>`;
          section+=`<p><strong>Chosen local endpoint sectors:</strong> ${endpointStarChecks.map(s=>'A'+s.id+': ⟨A⟩ in process (1) = '+rNumber(s.process1)+', process (2) = '+rNumber(s.process2)).join('; ')||'no unique outward destination-local star available on this geometry'}. The fixed destination projectors follow the four source-interface projectors.</p>`;
          section+=`<p><strong>Single-star comparison:</strong> ${allSingleStarMismatch.length===0?'All '+starData.length+' complete-star expectations agree between the two final states.':'Different expectations at '+allSingleStarMismatch.length+' complete star(s): '+allSingleStarMismatch.map(x=>'A'+x.id).join(', ')+'.'} ${allSingleStarMismatch.length===0&&mag<1-1e-9?'Equal individual Aₛ expectations and identical triangle syndromes do not imply identical endpoint quantum states.':''}</p>`;
          section+=`<p><strong>Required one-star inter-endpoint buffer:</strong> process (1): ${gapA.quantumClean} clean out of ${gapA.geometricCount} complete corridor stars; process (2): ${gapB.quantumClean} clean out of ${gapB.geometricCount}. ${gapA.quantumOK&&gapB.quantumOK?'Both final endpoint gaps have at least one star with flux-free adjacent triangles and ⟨Aₛ⟩=+1.':'Insufficient clean-star separation in one or both final states. Correct the microscopic movement before interpreting a phase.'}</p>`;
          if(!passes&&mag<1-1e-9)section+='<p class="r-error">The output states are not proportional. Normalization is already included; it cannot turn a magnitude below one into a valid exchange phase.</p>';
          if(stageAudit.some(x=>!x.boundaryOK))section+='<p class="r-error">Boundary one-star physical-cleanliness test fails at one or more completed stages. This is independent of the purely geometric one-star clearance.</p>';
          section+=rAllStageAuditHTML(stageAudit);
          section+=rWeightTable(first,'A')+rWeightTable(second,'B');
          section+=`<details><summary><strong>Ground-state sandwich and total normalization</strong></summary><div class="r-math">Kⱼ = ordered circuit with unscaled P=(I+A)/2<br>⟨Ω|K₂†K₁|Ω⟩ = ${rFormatComplex(z)}<br>Exact raw form: ${rExactForm(rawExact)}<br>⟨Ω|K₁†K₁|Ω⟩ = ${rNumber(rawN1**2)}<br>⟨Ω|K₂†K₂|Ω⟩ = ${rNumber(rawN2**2)}<br>Total w₁ = ${rNumber(1/rawN1)}; total w₂ = ${rNumber(1/rawN2)}<br>⟨ψ₂|ψ₁⟩ = (∏ w₂)(∏ w₁) ⟨Ω|K₂†K₁|Ω⟩</div></details>`;
          section+=`<p>Final selected-star check: process (1) ${brokenA.length?'A'+brokenA.join(', A')+' not restored':'all requested stars +1'}; process (2) ${brokenB.length?'A'+brokenB.join(', A')+' not restored':'all requested stars +1'}.</p>`;
          section+=`<p><strong>Vacated position 2 (all ${vacatedStars.length} neighboring complete stars):</strong> process (1): ${unrestoredA.length?unrestoredA.map(x=>'A'+x.id+' = '+rNumber(x.value)).join(', '):'all +1'}; process (2): ${unrestoredB.length?unrestoredB.map(x=>'A'+x.id+' = '+rNumber(x.value)).join(', '):'all +1'}. ${unrestoredA.length||unrestoredB.length?'Not fully restored: a correct triangle count is insufficient for transport certification.':'Locally restored in both outputs.'}</p>`;
          section+=rSymbolicMethodHTML(rExactCache,exact);
          section+='<p class="r-caveat">The positive w factors normalize these specific successful branches. This alone does not establish state-independent success probability or preservation of a multidimensional fusion space. This longer-arm geometry changes the generated ribbons; heuristic interface operations may require recalibration, and no R phase is accepted if the clean-star buffer fails.</p>';
          const elapsed=performance.now()-start;
          section+=`<p class="r-caveat">Calculation and selected-star checks: ${elapsed.toFixed(1)} ms on this browser.</p><button class="btn" id="r-export-symbolic">Export calculation JSON</button>`;
          rLastSymbolicReport={status,kind:rSelectedKind,geometry:rGeom,adjustments:rAdjust,splittingProjectors:rInitialSplittingStars,extraM21Projector:rExtraM21Star,destinationLocalProjectors:rEndpointFix,m21SyndromeCalibration:rM21SyndromeSuggestion,engine:rExactCache.metadata(),normalization:'positive normalization after every successful star projector',eventsA:A,eventsB:B,projectionGroupsA:rWeightGroups(first,'A'),projectionGroupsB:rWeightGroups(second,'B'),projectorsA:first.projectionLog,projectorsB:second.projectionLog,rawNormA:rawN1,rawNormB:rawN2,normalizedNormA:n1,normalizedNormB:n2,rawOverlap:z,overlap:w,exactOverlap:exact,exactRawOverlap:rawExact,magnitude:mag,phasePi:phase,selectedStarExpectationsA:starA,selectedStarExpectationsB:starB,endpointStarChecks,allSingleStarMismatch,vacatedSourceStarsA:vacatedA,vacatedSourceStarsB:vacatedB,bufferA:gapA,bufferB:gapB,allStageBufferAudit:stageAudit,boundaryClearanceAllStages:stageAudit.every(x=>x.boundaryOK),triangleIdsA:cA.fluxIds,triangleIdsB:cB.fluxIds,timeMs:elapsed};
          if(token!==rCalcToken)return;
          rCalc.innerHTML=section;setTimeout(rResize,120);
        }catch(e){rCalc.innerHTML=section+'<p class="r-error">Symbolic evaluation failed: '+rEscape(e.message)+'</p>';console.error(e);}
      },15);
    }

    function rValidationContext(which,step,engine){
      const move=rSequences[which][step-1],info=rMoveEndpoints(which,step);
      const beforeEv=rStageEvents(which,step,'before');
      const ribbonEv=rStageEvents(which,step,'ribbon');
      const xEv=rStageEvents(which,step,'x');
      const projEv=rStageEvents(which,step,'projected');
      const eventSnapshots=[beforeEv,ribbonEv,xEv,projEv];
      const stageEvents=[ribbonEv.slice(beforeEv.length),xEv.slice(ribbonEv.length),projEv.slice(xEv.length)];
      let state=null,stageStates=null,results=[];
      if(engine){
        state=engine.run(beforeEv);stageStates=[engine.copy(state)];
        for(const events of stageEvents){
          for(const ev of events)engine.gate(state,ev);
          stageStates.push(engine.copy(state));
        }
      }
      const projectors=projEv.slice(xEv.length).filter(e=>e.type==='P');
      const checked=[...new Set(projectors.map(e=>e.starId))];
      const nearby=starData.slice().sort((a,b)=>{
        const p=rPoint(info.source);
        return Math.hypot(a.center[0]-p.x,a.center[1]-p.y)-Math.hypot(b.center[0]-p.x,b.center[1]-p.y);
      }).slice(0,Math.max(6,checked.length)).map(s=>s.id);
      const starIds=[...new Set([...checked,...nearby])];
      for(let k=0;k<eventSnapshots.length;k++){
        const report=rEndpointFluxReport(eventSnapshots[k],k===0?info.before:info.after);
        const fluxNearSource=rNearInterfaceTriangles(eventSnapshots[k],info.source);
        const spectatorIds=rSpectatorFlux(eventSnapshots[k],info.spectator);
        const expectations=state?Object.fromEntries(starIds.map(id=>[id,engine.starExpectation(stageStates[k],id)])):null;
        results.push({events:eventSnapshots[k],report,fluxNearSource,spectatorIds,expectations,norm:state?engine.norm(stageStates[k]):null});
      }
      const beforeSpectator=rSpectatorFlux(beforeEv,info.spectator);
      const spectatorUnchanged=rSameIds(beforeSpectator,results[3].spectatorIds);
      const sourceClean=results[2].fluxNearSource.length===0;
      // Check EVERY complete star intersecting the vacated interface region,
      // not only the selected projectors. On the 6-step junction there are
      // four complete source stars inside this 2.1-unit neighborhood.
      const sourceCenter=rPoint(info.source);
      const interfaceStarIds=starData.filter(st=>Math.hypot(st.center[0]-sourceCenter.x,st.center[1]-sourceCenter.y)<=2.1+1e-7).map(st=>st.id);
      const unrestoredInterface=state?interfaceStarIds.filter(id=>{
        const v=results[3].expectations[id];
        return !Number.isFinite(v)||Math.abs(v-1)>1e-7;
      }):null;
      const restored=state?checked.filter(id=>!Number.isFinite(results[3].expectations[id])||Math.abs(results[3].expectations[id]-1)>1e-7):null;
      const nPrior=stageStates?stageStates[0].outcomes.length:0;
      const newRatios=stageStates?stageStates[3].outcomes.slice(nPrior):null;
      const nonzero=newRatios?newRatios.every(x=>x>1e-12):null;
      const pCombined=newRatios?newRatios.reduce((a,b)=>a*b,1):null;
      const gap=rGapStarReport(projEv,info.after[0],info.after[1],engine,stageStates?stageStates[3]:null);
      return {which,step,move,info,results,starIds,checked,projectors,restored,interfaceStarIds,unrestoredInterface,newRatios,nonzero,pCombined,spectatorUnchanged,sourceClean,inputValid:results[0].report.ok,gap};
    }
    function rVFlag(value,yes,no,unknown='not evaluated'){
      if(value===null||value===undefined)return `<span class="r-v-neutral">${unknown}</span>`;
      return `<span class="${value?'r-v-ok':'r-v-bad'}">${value?yes:no}</span>`;
    }
    function rVContextSummary(ctx){
      const fin=ctx.results[3], report=fin.report;
      return `<tr><td>${ctx.which==='A'?'(1): Y spectator':'(2): 0 spectator'}</td><td>${rVFlag(ctx.sourceClean,'clear','residual flux')}</td><td>${rVFlag(report.ok,'RGB × 2','mismatch')}</td><td>${rVFlag(ctx.spectatorUnchanged,'unchanged','changed')}</td><td>${rVFlag(ctx.restored===null?null:ctx.restored.length===0,'all +1','not all +1')}</td></tr>`;
    }
    function rValidationHTML(ctx,other){
      const v=ctx.results, labels=['Before move','After new ribbon','After interface X','After projections'];
      const info=ctx.info,validFlux=v[3].report;
      let html=`<p class="r-v-row"><strong>${ctx.which==='A'?'Process (1)':'Process (2)'}, ${rEscape(rDescriptions[ctx.move])}</strong><br>Source: q${info.source} · Target: q${info.destination} · Spectator: q${info.spectator}</p>`;
      html+=`<div class="r-v-row">${rVFlag(ctx.inputValid,'✓ Input has two expected RGB clusters','⚠ Input already has extra/misplaced flux')}</div>`;
      if(!ctx.inputValid)html+='<p class="r-v-hint">The preceding move has already failed the flux check. A failure here cannot be attributed to this movement alone; fix the earlier move first.</p>';
      html+='<table><thead><tr><th>Checkpoint</th><th>B<sub>t</sub>=−1</th><th>RGB at endpoints</th><th>Source flux*</th></tr></thead><tbody>';
      for(let i=0;i<v.length;i++){
        const out=v[i].report;
        const expectation=i===0?'before move':i===1?'intermediate':(out.ok?'yes':'no');
        html+=`<tr><td>${labels[i]}</td><td>${out.flux.length}</td><td>${expectation}</td><td>${v[i].fluxNearSource.length}</td></tr>`;
      }
      html+='</tbody></table>';
      html+=`<div class="r-v-row">${rVFlag(ctx.sourceClean,'✓ Old interface flux cleared after X','✗ Residual interface flux after X')} · ${rVFlag(validFlux.ok,'✓ Final RGB clusters at targets','✗ Final RGB clusters not clean')} · ${rVFlag(ctx.spectatorUnchanged,'✓ Spectator triangle IDs unchanged','✗ Spectator triangle IDs changed')}</div>`;
      html+=`<p class="r-v-hint">*Source flux means excited triangle centers within 1.15 lattice units of the source. RGB locality uses a 1.85-unit radius and checks three distinct colors at each expected endpoint. These geometric thresholds are diagnostics, not topological proofs.</p>`;
      html+=`<div class="r-v-row">Final triangle IDs: <span class="r-validation-codes">${validFlux.flux.map(t=>'B'+t.id).join(', ')||'none'}</span></div>`;
      if(ctx.move==='M21'&&rM21SyndromeSuggestion){
        const info=rM21SyndromeSuggestion;
        html+=`<div class="r-v-row"><strong>Both-context source syndrome check:</strong> ${rEscape(info.reason)}<br><small>Extra interface X from triangle calibration: ${info.ok?(info.qubits.map(q=>'q'+q).join(', ')||'none'):'unresolved'}. Correction settings remain editable. Matching triangles does not certify fusion-channel transport.</small></div>`;
      }
      // Projected-star success alone can conceal another unresolved complete
      // star near the source. Report it; do NOT silently add projectors.
      if(ctx.restored!==null){
        const unprojected=ctx.starIds.filter(id=>!ctx.checked.includes(id)&&
          Number.isFinite(v[3].expectations[id])&&Math.abs(v[3].expectations[id]-1)>1e-7);
        if(unprojected.length)html+=`<div class="r-v-row r-v-bad"><strong>Warning: unprojected nearby stars are not +1:</strong> ${unprojected.map(id=>'A'+id+' = '+v[3].expectations[id].toFixed(3)).join('; ')}. Flux correction alone does not restore these stars.</div>`;
      }
      const gap=ctx.gap;
      const edge=rBoundaryReport(ctx.results[3].events,ctx.info.after,null,null);
      html+=`<div class="r-v-row"><strong>Open-boundary one-star geometry:</strong> ${edge.perEndpoint.map(e=>'q'+e.endpoint+': '+e.geometricCount+' complete star(s) to nearest edge').join('; ')}. Full Aₛ/Bₜ cleanliness is checked in the eight-stage audit once both processes are completed.</div>`;
      html+=`<div class="r-v-row"><strong>One-star separation (after this move):</strong> ${gap.geometricCount} complete star centers in the endpoint gap; ${gap.geometricClean} have no adjacent Bₜ violations; ${gap.quantumOK===null?'star expectations unavailable':gap.quantumClean+' have both Bₜ=+1 locally and ⟨Aₛ⟩=+1'}. ${rVFlag(gap.quantumOK,'✓ At least one physically clean buffer star','✗ No fully restored buffer star','not yet validated')}.<br><small>Gap-star IDs: ${gap.stars.map(x=>'A'+x.id+(x.clean?' ✓':x.fluxFree?' (star unchecked)':' (flux nearby)')).join(', ')||'none'}.</small></div>`;
      if(ctx.restored!==null){
        const fullInterfaceOK=ctx.unrestoredInterface.length===0;
        html+=`<div class="r-v-row"><strong>Complete source-interface star check:</strong> ${rVFlag(fullInterfaceOK,'✓ Every nearby complete source star is +1','✗ Unrestored source stars: '+ctx.unrestoredInterface.map(id=>'A'+id).join(', '))}<br><small>Source neighborhood: ${ctx.interfaceStarIds.map(id=>'A'+id+' = '+rNumber(v[3].expectations[id])).join('; ')||'none'}.</small> This is a local diagnostic, not a fusion-space isometry proof.</div>`;
        html+=`<div class="r-v-row">${rVFlag(ctx.restored.length===0,'✓ All chosen projected stars have ⟨A⟩=+1','✗ Stars not restored: '+ctx.restored.map(id=>'A'+id).join(', '))}</div>`;
        html+=`<div class="r-v-row">${rVFlag(ctx.nonzero,'✓ Every projector branch has nonzero probability','✗ A projector branch has zero probability')}</div>`;
        html+=`<div class="r-v-row">Conditional projector probabilities: <span class="r-validation-codes">${ctx.checked.length?ctx.projectors.map((e,i)=>'A'+e.starId+': '+ctx.newRatios[i].toFixed(4)).join('; '):'no projectors'}</span><br>Combined conditional probability: <strong>${rNumber(ctx.pCombined)}</strong><br>Movement normalization w = 1/√p = <strong>${ctx.pCombined>0?rNumber(1/Math.sqrt(ctx.pCombined)):'undefined (zero branch)'}</strong><br>Each successful projector is positively normalized; the checkpoint state has norm 1.</div>`;
        html+='<details><summary>Star expectations at each checkpoint</summary><table><thead><tr><th>Star</th><th>Before</th><th>Ribbon</th><th>X</th><th>P</th></tr></thead><tbody>';
        for(const id of ctx.starIds){
          const fmt=(k)=>Number.isFinite(v[k].expectations[id])?v[k].expectations[id].toFixed(3):'undefined';
          html+=`<tr><td>A${id}${ctx.checked.includes(id)?' *':''}</td>${[0,1,2,3].map(k=>`<td>${fmt(k)}</td>`).join('')}</tr>`;
        }
        html+='</tbody></table><p class="r-v-hint">* Selected projectors. The remaining rows are nearby complete stars (diagnostic only). Expectations are conditional on the successful projector branches.</p></details>';
      }else{
        html+=`<p class="r-v-hint"><strong>Star expectations and projector probabilities unavailable:</strong> The symbolic engine could not evaluate the current reference sector. The triangle syndromes above are still evaluated exactly from X parity.</p>`;
      }
      if(other){
        html+='<details open><summary><strong>Same M(2→1) in both spectator contexts</strong></summary>';
        html+=`<p class="r-v-hint">Process (1) input: ${rVFlag(ctx.which==='A'?ctx.inputValid:other.inputValid,'clean RGB','already violated')}; Process (2) input: ${rVFlag(ctx.which==='B'?ctx.inputValid:other.inputValid,'clean RGB','already violated')}.</p>`;
        html+='<table><thead><tr><th>Context</th><th>Source</th><th>Endpoints</th><th>Spectator</th><th>Stars</th></tr></thead><tbody>';
        html+=rVContextSummary(ctx)+rVContextSummary(other)+'</tbody></table>';
        html+='<p class="r-v-hint">The exact same ribbon, one mandatory X(q2), and projector list are used in both contexts; the input quantum states differ. Matching syndromes alone does not prove preservation of the fusion space.</p></details>';
      }
      html+='<p class="r-v-hint"><strong>Not yet certified:</strong> preservation of the non-Abelian fusion space and inverse-path consistency. Passing these local tests is necessary but not sufficient to identify an R-symbol.</p>';
      return html;
    }
    function rValidateSelected(){
      const token=++rValidationToken;
      if(!rGeom||!rActiveStep){rValidationContent.innerHTML='<p>Choose one of the movement buttons to inspect its checkpoints.</p>';return;}
      if(rActiveStep===1){
        const ev=rStepEvents(rActiveProcess,1),report=rEndpointFluxReport(ev,[rGeom.q1,rGeom.q2]);
        let splitGap;
        try{if(!rExactCache)rExactCache=new RSymbolicEngine(siteData,starData,triangleData);
          splitGap=rGapStarReport(ev,rGeom.q1,rGeom.q2,rExactCache,rExactCache.run(ev));}
        catch(e){splitGap=rGapStarReport(ev,rGeom.q1,rGeom.q2);}
        let splitWeights='',splitExpectations='';
        try{
          if(!rExactCache)rExactCache=new RSymbolicEngine(siteData,starData,triangleData);
          const prepared=rExactCache.run(ev);
          const terms=prepared.projectionLog.filter(x=>x.step===1);
          splitWeights=terms.length?terms.map(x=>`P${x.starId}: p=${rNumber(x.p)}, w=${rNumber(x.w)}`).join('; '):'none for this geometry';
          splitExpectations=rInitialSplittingStars.length?rInitialSplittingStars.map(id=>`⟨A${id}⟩ = ${rNumber(rExactCache.starExpectation(prepared,id))}`).join('; '):'';
        }catch(e){splitWeights='symbolic evaluation unavailable: '+rEscape(e.message);}
        rValidationContent.innerHTML=`<p>Initial splitting S(1,2): ${report.flux.length} excited triangles. Expected two RGB endpoints: ${rVFlag(report.ok,'geometrically consistent','not yet consistent')}.</p><p><strong>Splitting-state projections (applied after S, in order):</strong> ${splitWeights}. ${splitExpectations}</p><p><strong>Buffer between 1 and 2:</strong> ${splitGap.geometricCount} complete stars, ${splitGap.geometricClean} flux-free, ${splitGap.quantumOK===null?'star expectation unavailable':splitGap.quantumClean+' fully restored'}. ${rVFlag(splitGap.quantumOK,'At least one clean star','Insufficient clean buffer star','Not fully checked')}.</p><p class="r-v-hint">Choose movement step 2, 3 or 4 to inspect the ribbon, X and projector checkpoints.</p>`;
        return;
      }
      rValidationContent.innerHTML='<p>Evaluating interface checkpoints…</p>';
      setTimeout(()=>{
        if(token!==rValidationToken)return;
        try{
          let engine=null;
          if(!rExactCache)rExactCache=new RSymbolicEngine(siteData,starData,triangleData);engine=rExactCache;
          const ctx=rValidationContext(rActiveProcess,rActiveStep,engine);
          const other=ctx.move==='M21'?rValidationContext(rActiveProcess==='A'?'B':'A',3,engine):null;
          if(token!==rValidationToken)return;
          rValidationContent.innerHTML=rValidationHTML(ctx,other);
          setTimeout(rResize,100);
        }catch(e){rValidationContent.innerHTML=`<p class="r-error">Interface validation failed: ${rEscape(e.message)}</p>`;console.error(e);}
      },15);
    }

    /* Phase-sensitive symbolic contraction for the Kagome complete-star orbit.
     * No 2^N amplitude array. Projection branches are affine variables, not vectors.
     * Gate set: X, Z, R=exp(+i*pi*Z/4), Rdag, CZ, A_s and (I +/- A_s)/2.
     * Reference: normalized product_s (I+A_s)|0...0>, independent star X supports.
     * This is NOT a claim that the physical D4 ground state is a stabilizer state.
     */
    const rSymBitsCache = new Map();
    function rSymIndex(bit) {
      let v=rSymBitsCache.get(bit);
      if(v===undefined){v=bit.toString(2).length-1;rSymBitsCache.set(bit,v);}
      return v;
    }
    function rSymIndices(mask) {
      const out=[];
      while(mask){const bit=mask&-mask;out.push(rSymIndex(bit));mask^=bit;}
      return out;
    }
    function rSymParity(mask) {let p=0;while(mask){mask&=mask-1n;p^=1;}return p;}
    function rSymMod4(x) {return ((x%4)+4)%4;}
    function rSymMod8(x) {return ((x%8)+8)%8;}
    function rSymRank(rows) {
      const basis=new Map();
      for(let v of rows){while(v){const p=v.toString(2).length-1;if(basis.has(p))v^=basis.get(p);else{basis.set(p,v);break;}}}
      return basis;
    }
    function rSymContains(basis,vec) {
      while(vec){const p=vec.toString(2).length-1;if(!basis.has(p))return false;vec^=basis.get(p);}return true;
    }
    class RPhasePolynomial {
      constructor(n=0){this.n=n;this.c=0;this.lin=Array(n).fill(0);this.adj=Array(n).fill(0n);}
      clone(){const p=new RPhasePolynomial(this.n);p.c=this.c;p.lin=this.lin.slice();p.adj=this.adj.slice();return p;}
      grow(){this.lin.push(0);this.adj.push(0n);return this.n++;}
      edge(i,j){if(i===j){this.lin[i]=(this.lin[i]+2)&3;return;}this.adj[i]^=1n<<BigInt(j);this.adj[j]^=1n<<BigInt(i);}
      // Add coeff*(constant XOR parity(mask)) mod 4. Odd parity has pair terms.
      parity(coeff,mask,constant=0){
        coeff=rSymMod4(coeff);if(!coeff)return;
        if(constant){this.c=(this.c+coeff)&3;coeff=(-coeff)&3;}
        const ids=rSymIndices(mask);
        for(const i of ids)this.lin[i]=(this.lin[i]+coeff)&3;
        if(coeff&1)for(let a=0;a<ids.length;a++)for(let b=a+1;b<ids.length;b++)this.edge(ids[a],ids[b]);
      }
      // Add 2*(ca XOR parity(a))*(cb XOR parity(b)) mod 4.
      product(a,ca,b,cb){
        if(ca&&cb)this.c^=2;
        if(ca)this.parity(2,b);if(cb)this.parity(2,a);
        const ai=rSymIndices(a),bi=rSymIndices(b);
        for(const i of ai)for(const j of bi)this.edge(i,j);
      }
      append(other,offset,sign=1){
        this.c=rSymMod4(this.c+sign*other.c);
        for(let i=0;i<other.n;i++){
          this.lin[i+offset]=rSymMod4(this.lin[i+offset]+sign*other.lin[i]);
          // A quadratic coefficient is 2; its negative is also 2 mod 4.
          for(const j of rSymIndices(other.adj[i]))if(j>i)this.edge(i+offset,j+offset);
        }
      }
      substitute(forms,nfree){
        const out=new RPhasePolynomial(nfree);out.c=this.c;
        for(let i=0;i<this.n;i++){
          if(this.lin[i])out.parity(this.lin[i],forms[i].mask,forms[i].c);
          for(const j of rSymIndices(this.adj[i]))if(j>i)out.product(forms[i].mask,forms[i].c,forms[j].mask,forms[j].c);
        }
        return out;
      }
      toJSON(){
        const quadratic=[];
        for(let i=0;i<this.n;i++)for(const j of rSymIndices(this.adj[i]))if(j>i)quadratic.push([i,j]);
        return {variables:this.n,constantMod4:this.c,linearMod4:this.lin.slice(),quadraticCoefficient2:quadratic,meaning:'q(x)=constant+sum linear_i*x_i+2*sum x_i*x_j (mod 4), with x_i in {0,1}'};
      }
      // Exact binary quadratic Gauss sum. Result is 0 or 2^power * exp(i*pi*phase8/4).
      gauss(){
        const q=this.clone();let active=(1n<<BigInt(q.n))-1n,power=0,phase8=0,oddEliminations=0,pairEliminations=0;
        const remove=(i)=>{const bit=1n<<BigInt(i);for(const j of rSymIndices(q.adj[i]))q.adj[j]&=~bit;q.adj[i]=0n;q.lin[i]=0;active&=~bit;};
        while(active){
          const ids=rSymIndices(active),odd=ids.find(i=>q.lin[i]&1);
          if(odd!==undefined){
            const coeff=q.lin[odd],neighbors=q.adj[odd]&active;
            remove(odd);q.parity(-coeff,neighbors);power+=.5;phase8+=coeff===1?1:-1;oddEliminations++;continue;
          }
          const i=ids.find(j=>(q.adj[j]&active)!==0n);
          if(i!==undefined){
            const j=rSymIndex((q.adj[i]&active)&-(q.adj[i]&active));
            const a=q.adj[i]&active&~(1n<<BigInt(j)),b=q.adj[j]&active&~(1n<<BigInt(i));
            const ca=q.lin[i]>>1,cb=q.lin[j]>>1;
            remove(i);remove(j);q.product(a,ca,b,cb);power++;pairEliminations++;continue;
          }
          if(ids.some(j=>q.lin[j]===2))return {zero:true,power:-Infinity,phase8:0,oddEliminations,pairEliminations};
          power+=ids.length;active=0n;
        }
        return {zero:false,power,phase8:rSymMod8(phase8+2*q.c),oddEliminations,pairEliminations};
      }
    }
    function rSymAffineSolve(equations,n){
      const basis=new Map();
      for(const row of equations){
        let mask=row.mask,c=row.c&1;
        while(mask){const p=mask.toString(2).length-1,old=basis.get(p);if(old){mask^=old.mask;c^=old.c;}else{basis.set(p,{mask,c});break;}}
        if(mask===0n&&c)return null;
      }
      const forms=Array(n),free=[];
      for(let i=0;i<n;i++)if(!basis.has(i)){forms[i]={mask:1n<<BigInt(free.length),c:0};free.push(i);}
      for(const i of [...basis.keys()].sort((a,b)=>a-b)){
        const row=basis.get(i);let mask=0n,c=row.c;
        for(const j of rSymIndices(row.mask^(1n<<BigInt(i)))){mask^=forms[j].mask;c^=forms[j].c;}
        forms[i]={mask,c};
      }
      return {forms,nfree:free.length,rank:basis.size};
    }
    function rSymComplex(result){
      if(result.zero)return {re:0,im:0};
      const rt=Math.SQRT1_2,root=[[1,0],[rt,rt],[0,1],[-rt,rt],[-1,0],[-rt,-rt],[0,-1],[rt,-rt]][rSymMod8(result.phase8)];
      const a=2**result.power;return {re:a*root[0],im:a*root[1]};
    }
    class RSymbolicEngine {
      constructor(sites,stars,triangles=[]){
        this.sites=sites;this.stars=stars;this.triangles=triangles;this.n=stars.length;
        if(!this.n)throw Error('The reference patch needs at least one complete star.');
        this.qIndex=new Map(sites.map((q,i)=>[Number(q.id??q.site_id),i]));
        this.sIndex=new Map(stars.map((s,i)=>[Number(s.id),i]));
        if(this.qIndex.size!==sites.length||this.sIndex.size!==stars.length)throw Error('Duplicate site/star IDs.');
        this.siteMask=Array(sites.length).fill(0n);this.starX=[];this.starEdges=[];this.basis=new Map();
        for(let j=0;j<this.n;j++){
          let bits=0n;
          for(const q of stars[j].x_qubits){const i=this.index(q);this.siteMask[i]^=1n<<BigInt(j);bits^=1n<<BigInt(i);}
          this.starX.push(bits);this.starEdges.push(stars[j].cz_edges.map(([a,b])=>[this.index(a),this.index(b)]));
          let v=bits,comb=1n<<BigInt(j);
          while(v){const p=v.toString(2).length-1,prev=this.basis.get(p);if(prev){v^=prev.vec;comb^=prev.comb;}else{this.basis.set(p,{vec:v,comb});break;}}
          if(!v)throw Error('Dependent complete-star X supports: this reference engine requires an independent open-patch basis.');
        }
        // Verify that the assumed orbit remains in the retained triangle +1 sector.
        const triBits=[];
        for(const t of triangles){let p=0n,z=0n;for(const q of t.qubits){p^=this.siteMask[this.index(q)];z^=1n<<BigInt(this.index(q));}if(p)throw Error('A triangle anticommutes with a complete-star X support.');triBits.push(z);}
        this.triangleBasis=rSymRank(triBits);
        this.stats={contractions:0,maxVariables:0,maxFreeVariables:0,zeroContractions:0,oddEliminations:0,pairEliminations:0};
        this.lastContraction=null;this.prefixCache=new Map();
      }
      index(q){const i=this.qIndex.get(Number(q));if(i===undefined)throw Error('Unknown physical qubit '+q);return i;}
      sid(id){const i=this.sIndex.get(Number(id));if(i===undefined)throw Error('Unknown complete star '+id);return i;}
      solveOffset(delta){let shift=0n;while(delta){const p=delta.toString(2).length-1,old=this.basis.get(p);if(!old)return null;delta^=old.vec;shift^=old.comb;}return shift;}
      state(){return {poly:new RPhasePolynomial(this.n),t:Array.from({length:this.n},(_,j)=>1n<<BigInt(j)),tc:0n,off:0n,phase8:0,logScale:-this.n/2,logSuccess:0,outcomes:[],projectionLog:[],dead:false,history:[]};}
      copy(s){return {poly:s.poly.clone(),t:s.t.slice(),tc:s.tc,off:s.off,phase8:s.phase8,logScale:s.logScale,logSuccess:s.logSuccess,outcomes:s.outcomes.slice(),projectionLog:s.projectionLog.map(v=>({...v})),dead:s.dead,history:s.history.map(e=>({...e,qubits:e.qubits?.slice()}))};}
      virtual(s,mask){let m=0n;for(const i of rSymIndices(mask))m^=s.t[i];return {mask:m,c:rSymParity(s.tc&mask)};}
      physical(s,i,extraVirtual=0n,includeOff=true){const f=this.virtual(s,this.siteMask[i]);f.c^=rSymParity(extraVirtual&this.siteMask[i]);if(includeOff)f.c^=Number((s.off>>BigInt(i))&1n);return f;}
      residual(s,sid){
        let m=0n,c=0;
        for(const [a,b] of this.starEdges[sid]){const oa=Number((s.off>>BigInt(a))&1n),ob=Number((s.off>>BigInt(b))&1n);if(oa)m^=this.siteMask[b];if(ob)m^=this.siteMask[a];c^=oa&ob;}
        if(m&(1n<<BigInt(sid)))throw Error('Non-Hermitian star-frame Pauli: inconsistent geometry.');
        const f=this.virtual(s,m);f.c^=c;return f;
      }
      applyStar(s,id){const sid=this.sid(id),f=this.residual(s,sid);s.poly.parity(2,f.mask,f.c);s.tc^=1n<<BigInt(sid);}
      innerExact(A,B,explain=false){
        const stat=this.stats;stat.contractions++;
        if(A.dead||B.dead)return {zero:true,power:-Infinity,phase8:0};
        const shift=this.solveOffset(A.off^B.off);
        if(shift===null){stat.zeroContractions++;return {zero:true,power:-Infinity,phase8:0,reason:'Disjoint computational-basis affine supports'};}
        const na=A.poly.n,nb=B.poly.n,n=na+nb,off=BigInt(na);
        stat.maxVariables=Math.max(stat.maxVariables,n);
        const p=new RPhasePolynomial(n);p.append(A.poly,0,1);p.append(B.poly,na,-1);
        // Align the eta ground-state signs, not just the bit-string supports.
        let previous=0n;
        for(const sid of rSymIndices(shift)){
          for(const [a,b] of this.starEdges[sid]){const fa=this.physical(A,a,previous,false),fb=this.physical(A,b,previous,false);p.product(fa.mask,fa.c,fb.mask,fb.c);}
          previous^=1n<<BigInt(sid);
        }
        const eq=[];
        for(let j=0;j<this.n;j++)eq.push({mask:A.t[j]^(B.t[j]<<off),c:Number(((A.tc^B.tc^shift)>>BigInt(j))&1n)});
        const solved=rSymAffineSolve(eq,n);
        if(!solved){stat.zeroContractions++;return {zero:true,power:-Infinity,phase8:0,reason:'Incompatible affine constraints'};}
        stat.maxFreeVariables=Math.max(stat.maxFreeVariables,solved.nfree);
        const reduced=p.substitute(solved.forms,solved.nfree),g=reduced.gauss();
        stat.oddEliminations+=g.oddEliminations||0;stat.pairEliminations+=g.pairEliminations||0;
        const result={...g,power:g.zero?-Infinity:g.power+A.logScale+B.logScale,phase8:rSymMod8(g.phase8+A.phase8-B.phase8),variables:n,freeVariables:solved.nfree,rank:solved.rank};
        if(result.zero)stat.zeroContractions++;
        if(explain)result.reduction={
          prefactorPower:A.logScale+B.logScale,prefactorPhase8:rSymMod8(A.phase8-B.phase8),
          groundSignAlignmentStars:rSymIndices(shift).map(i=>this.stars[i].id),
          affineConstraintRank:solved.rank,originalVariables:n,freeVariables:solved.nfree,
          reducedPolynomial:reduced.toJSON(),gaussSum:g,
          interpretation:'Sum all binary assignments of the reduced polynomial symbolically, then multiply by the power-of-two and eighth-root-of-unity prefactor. No assignments or state amplitudes are enumerated.'
        };
        this.lastContraction={...result};return result;
      }
      overlap(A,B){return rSymComplex(this.innerExact(A,B));} // <B|A>, full phase retained
      norm(s){if(s.dead)return 0;const n=this.innerExact(s,s);if(n.zero)return 0;if(n.phase8!==0)throw Error('Norm contraction is not positive real.');return 2**(n.power/2);}
      rawNorm(s){return s.dead?0:2**(s.logSuccess/2);}
      project(s,id,meta={}){
        if(s.dead){s.outcomes.push(0);s.projectionLog.push({starId:id,p:0,w:null,step:meta.step??null,move:meta.move??null,unreachable:true});return;}
        const before=this.innerExact(s,s);if(before.zero)throw Error('Cannot project a zero branch.');
        const sid=this.sid(id),f=this.residual(s,sid),h=s.poly.grow(),bit=1n<<BigInt(h);
        // Explicit (I +/- A)/2 branch, represented by one binary variable h.
        s.poly.product(bit,0,f.mask,f.c^(meta.sign===-1?1:0));
        s.t[sid]^=bit;s.logScale-=1;
        const after=this.innerExact(s,s);
        if(after.zero){s.dead=true;s.logSuccess=-Infinity;s.outcomes.push(0);s.projectionLog.push({starId:id,p:0,w:null,step:meta.step??null,move:meta.move??null});return;}
        if(before.phase8!==0||after.phase8!==0)throw Error('Projector norm is not positive real.');
        const logp=after.power-before.power,p=2**logp;
        if(p>1+1e-9)throw Error('Projector probability exceeds one.');
        // A positive real rescaling retains the phase of the successful branch.
        s.logScale-=after.power/2;s.logSuccess+=logp;s.outcomes.push(p);
        s.projectionLog.push({starId:id,p,w:2**(-logp/2),log2p:logp,step:meta.step??null,move:meta.move??null,variables:after.variables,freeVariables:after.freeVariables});
      }
      gate(s,ev){
        s.history.push({...ev,qubits:ev.qubits?.slice()});
        if(ev.type==='P'){this.project(s,ev.starId,ev);return;}
        if(s.dead)return;
        if(ev.type==='A'){this.applyStar(s,ev.starId);return;}
        const ids=(ev.qubits||[]).map(q=>this.index(q));const a=ids[0],b=ids[1];
        if(ev.type==='X'){s.off^=1n<<BigInt(a);return;}
        if(ev.type==='Z'){const f=this.physical(s,a);s.poly.parity(2,f.mask,f.c);return;}
        if(ev.type==='R'||ev.type==='Rdag'){
          const sign=ev.type==='R'?1:-1,f=this.physical(s,a);
          s.phase8=rSymMod8(s.phase8+sign);s.poly.parity(-sign,f.mask,f.c);return;
        }
        if(ev.type==='CZ'){if(a===b)throw Error('CZ needs two distinct qubits.');const fa=this.physical(s,a),fb=this.physical(s,b);s.poly.product(fa.mask,fa.c,fb.mask,fb.c);return;}
        throw Error('Unsupported symbolic gate '+ev.type);
      }
      run(events){
        const key=JSON.stringify(events.map(e=>[e.type,e.qubits??null,e.starId??null,e.sign??1,e.step??null,e.move??null]));
        const old=this.prefixCache.get(key);if(old)return this.copy(old);
        const s=this.state();for(const ev of events)this.gate(s,ev);
        if(this.prefixCache.size>=24)this.prefixCache.delete(this.prefixCache.keys().next().value);
        this.prefixCache.set(key,this.copy(s));return s;
      }
      starExpectation(s,id){if(s.dead)return NaN;const acted=this.copy(s);this.applyStar(acted,id);const z=this.overlap(acted,s);if(Math.abs(z.im)>1e-8)throw Error('Star expectation is not real.');return Math.abs(z.re)<1e-12?0:z.re;}
      zStringReference(qubits){
        let qbits=0n,parity=0n;
        for(const q of qubits){const i=this.index(q);qbits^=1n<<BigInt(i);parity^=this.siteMask[i];}
        const bt=rSymContains(this.triangleBasis,qbits);
        return {value:parity?0:1,triangleProduct:bt,witnessStar:parity?this.stars[rSymIndex(parity&-parity)].id:null,boundarySectorNeeded:!parity&&!bt};
      }
      metadata(){return {method:'Exact phase-polynomial / residual-character contraction',reference:'normalized product_s (I+A_s)|all physical qubits 0>, with retained independent complete stars',completeStars:this.n,physicalQubits:this.sites.length,triangleRank:this.triangleBasis.size,diagonalReferenceRank:this.sites.length-this.n,extraDiagonalConstraints:this.sites.length-this.n-this.triangleBasis.size,amplitudeArraysAllocated:0,...this.stats};}
    }
    if(typeof module!=='undefined'&&module.exports)module.exports={RPhasePolynomial,RSymbolicEngine,rSymAffineSolve,rSymComplex};


    rCalc.addEventListener('click',(ev)=>{
      if(ev.target.id==='r-export-symbolic'&&rLastSymbolicReport){
        const text=JSON.stringify(rLastSymbolicReport,(_,v)=>typeof v==='bigint'?v.toString():v,2);
        const blob=new Blob([text],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');
        a.href=url;a.download='kagome_R_symbolic_'+rSelectedKind+'.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),500);
      }
      if(ev.target.id==='r-z-test'){
        const input=document.getElementById('r-z-query'),answer=document.getElementById('r-z-answer');
        try{
          const qs=input.value.trim()?input.value.trim().split(/[\s,;]+/).map(Number):[];
          if(qs.some(q=>!Number.isInteger(q)||!siteById.has(q)))throw Error('Enter valid qubit IDs; an empty list means identity.');
          const c=rExactCache.zStringReference(qs);
          answer.innerHTML='<p><strong>⟨Ω|Z(residual)|Ω⟩ = '+c.value+'</strong><br>'+(c.triangleProduct?'Residual is in the span of Bₜ (expectation +1).':c.witnessStar!==null?'It anticommutes with complete star A'+c.witnessStar+', proving expectation zero.':'Outside the Bₜ span, but +1 in the chosen open-boundary reference sector. Setting this term to zero would change the reference state.')+'</p>';
        }catch(e){answer.textContent=e.message;}
      }
    });

    rSelect.addEventListener('change',()=>rSetup(rSelect.value));
    // Click the second canvas to edit interface X and star projector choices.
    function rHandlePlotClick(event){
      if(!rGeom||rActiveStep<=1||!event?.points?.length)return;
      const p=event.points[0],name=p.data?.name||p.fullData?.name;
      const move=rSequences[rActiveProcess][rActiveStep-1];
      if(name==='qubits' && Number.isInteger(Number(p.customdata))) {
        const id=Number(p.customdata),arr=rAdjust[move].x;
        // The fixed X(q2) already exists as a mandatory interface gate.
        if(move==='M21'&&id===rGeom.q2)return;
        rAdjust[move].x=arr.includes(id)?arr.filter(q=>q!==id):[...arr,id];
      } else if((name==='r-star-centers'||name==='geometric-buffer-stars') && Number.isInteger(Number(p.customdata))) {
        const id=Number(p.customdata),arr=rAdjust[move].stars;
        if(id===rEndpointFix[move]||(move==='M21'&&id===rExtraM21Star))return; // fixed projector already supplied
        rAdjust[move].stars=arr.includes(id)?arr.filter(q=>q!==id):[...arr,id];
      } else return;
      rCalcToken++;rValidationToken++;rRenderStep();
    }
    for(const btn of rStepButtons){
      btn.addEventListener('click',()=>{
        if(!rGeom)return;
        rActiveProcess=btn.dataset.rProcess;rActiveStep=Number(btn.dataset.rStep);rCheckpoint='projected';rCheckpointSelect.value=rCheckpoint;rValidationToken++;
        rCompleted[rActiveProcess]=Math.max(rCompleted[rActiveProcess],rActiveStep);
        rUpdateButtons();rRenderStep();
      });
    }
    rSaveAdjust.addEventListener('click',()=>{
      if(rActiveStep<=1||!rGeom)return;
      const move=rSequences[rActiveProcess][rActiveStep-1];
      const parse=(s)=>s.split(/[\s,;]+/).map(x=>x.trim()).filter(x=>x!=='').map(Number);
      const qs=parse(rEditX.value),stars=parse(rEditStars.value);
      if(qs.some(q=>!Number.isInteger(q)||!siteById.has(q))||stars.some(q=>!Number.isInteger(q)||!starById.has(q))){
        rCalc.innerHTML='<p class="r-error">Invalid qubit or complete-star ID.</p>';return;
      }
      if(qs.length>8||stars.length>8){rCalc.innerHTML='<p class="r-error">Interface has too many proposed corrections; review the geometry.</p>';return;}
      if((move==='M21'||move==='M0')&&stars.length!==4){rCalc.innerHTML='<p class="r-error">Parallel movement requires exactly four editable SOURCE-interface stars. Any fixed destination-local projector is separate and applied afterward.</p>';return;}
      if(stars.includes(rEndpointFix[move])){rCalc.innerHTML=`<p class="r-error">P${rEndpointFix[move]} is already included as the fixed destination-local projector for this movement; do not duplicate it in the source list.</p>`;return;}
      if(move==='M21'&&rExtraM21Star!==null&&stars.includes(rExtraM21Star)){
        rCalc.innerHTML=`<p class="r-error">P${rExtraM21Star} is already included as an additional fixed M(2→1) projector. Remove it from the four-star source list.</p>`;return;
      }
      if(move==='M21'&&qs.includes(rGeom.q2)){
        rCalc.innerHTML=`<p class="r-error">Position 2 (qubit ${rGeom.q2}) already receives one mandatory interface X. Remove it from the additional X list to avoid X²=I.</p>`;
        return;
      }
      // Duplicate entries would generate unintended pairwise X cancellations.
      if(new Set(qs).size!==qs.length||new Set(stars).size!==stars.length){
        rCalc.innerHTML='<p class="r-error">Please remove duplicate qubit or star IDs.</p>';return;
      }
      rAdjust[move]={x:qs,stars};rCalcToken++;rValidationToken++;rRenderStep();
    });

    rCheckpointSelect.addEventListener('change',()=>{
      rCheckpoint=rCheckpointSelect.value;rValidationToken++;rRenderStep();
    });

    initializeQubitSelector();
    updateTrackedTable();
    updateVisibilityButtons();
    updateViolationButtons();
    updateExpectationButtons();
    render();

    plotDiv.on("plotly_click", (eventData) => {
      if (zoomMode) {
        return;
      }
      if (!eventData || !eventData.points || eventData.points.length === 0) {
        return;
      }

      const point = eventData.points[0];
      if (expectationMode) {
        if (point.data && point.data.name==="inspect-centers") {
          const id=Number(point.customdata);
          const obj=(expectationMode==="STAR"?starById:triangleById).get(id);
          if (obj) selectExpectationPoint(obj.center[0],obj.center[1]);
        }
        // In inspection mode, never apply Z/X/R/CZ by clicking a qubit.
        return;
      }

      // Named traces avoid fragile curve-index dependencies as overlays
      // (violation markers, anyon highlights) are added or removed.
      const traceName = (point.data && point.data.name)
        || (point.fullData && point.fullData.name)
        || "";
      if (traceName === "m2-options") {
        if (anyonKind === "m2") handleMPairBond(String(point.customdata));
        return;
      }
      const clickableTraceNames = new Set([
        "qubit-nodes", "operator-labels", "pending-cz",
        "anyon-path", "anyon-next",
      ]);

      if (!clickableTraceNames.has(traceName)) {
        return;
      }

      const siteId = Number(point.customdata);
      if (Number.isNaN(siteId)) {
        return;
      }

      handleSiteClick(siteId);
    });
  </script>
</body>
</html>
"""

    html = (
        template
        .replace("__SITE_DATA__", json.dumps(site_data))
        .replace("__EDGE_PAIRS__", json.dumps(edge_pairs))
        .replace("__BASE_EDGE_X__", json.dumps(base_edge_x))
        .replace("__BASE_EDGE_Y__", json.dumps(base_edge_y))
        .replace("__STAR_DATA__", json.dumps(star_data))
        .replace("__TRIANGLE_DATA__", json.dumps(triangle_data))
    )

    return html
