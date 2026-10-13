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

      <button class="btn action" id="pan-btn">Pan</button>
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
    const panBtn = document.getElementById("pan-btn");
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

    let currentMode = "Z";
    let localOps = new Map();   // site_id -> "X" or "Z"
    let czPairs = [];           // array of [a, b]
    let pendingCZ = null;       // first endpoint for CZ selection
    let history = [];           // array of { state, description }
    let panMode = false;

    // Guided Anyon mode: electric (Z), magnetic (X + automatic CZ), or
    // three-color straight s1/s2 diagrams with XR/RX on the spine.
    // The magnetic circuit is *recompiled* after every chosen X, with all its
    // CZ gates chronologically before all its X gates. The mathematical
    // operator is therefore (product X) (product CZ), as requested.
    let anyonKind = null;      // "e", "m", "s1" or "s2"; null disables guide
    let anyonColor = null;     // "Red", "Blue", "Green", or null
    let anyonQubits = [];      // selected Z or X physical qubits
    let anyonStarPath = [];    // oriented e-star path (after 2nd Z)
    let mTrianglePath = [];    // oriented monochromatic triangle path
    let mAutoCZ = [];         // automatic CZ pairs in current m ribbon
    let mOptionCache = null;  // computed options for the current m path
    let mGuideId = null;      // unique ID for the current m guide session
    let nextMGuideId = 1;
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
      sBasePath = [...(state.sBasePath || [])];
      sGuidedX = [...(state.sGuidedX || [])];
      sAutoCZ = (state.sAutoCZ || []).map(pair => [...pair]);
      sGuideId = state.sGuideId ?? null;
      nextSGuideId = state.nextSGuideId ?? nextSGuideId;
      sCandidateCache = null;
      mOptionCache = null;
      anyonSelect.value = anyonKind
        ? (anyonKind === "s1" || anyonKind === "s2"
          ? anyonKind : `${anyonKind}_${anyonColor[0]}`)
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
      if (anyonKind === "s1" || anyonKind === "s2") return anyonKind;
      return `${anyonKind}_${anyonColor[0]}`;
    }

    function currentAnyonGate() {
      return anyonKind === "m" ? "X" : (anyonKind === "e" ? "Z" : "XR/RX");
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
      if (anyonKind === "s1" || anyonKind === "s2") {
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

    function completeMTriangleEnds(siteId) {
      const site = siteById.get(siteId);
      if (!site || site.colorName !== anyonColor) return [];
      return (mTriangleIncidence.get(siteId) || []).filter(
        tId => triangleById.get(tId).color === anyonColor
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

    function mCandidateCZEdges(xQubits, trianglePath, currentCZ) {
      // Candidate CZs connect the *other two colours*; both endpoints must
      // be tips of stars along the selected triangle path. They may be
      // non-nearest neighbours (the same upward arcs used in the paper).
      const tips = new Set();
      for (const tId of trianglePath.slice(1, -1)) {
        const sId = starForTriangle.get(tId);
        const star = starById.get(sId);
        if (star) star.x_qubits.forEach(id => tips.add(id));
      }
      const sites = Array.from(tips).filter(id => {
        const s = siteById.get(id);
        if (s.colorName === anyonColor) return false;
        return xQubits.some(x => {
          const p = siteById.get(x);
          const dx = p.x - s.x, dy = p.y - s.y;
          return dx*dx + dy*dy <= 2.7*2.7;
        });
      });
      const colours = ["Red", "Green", "Blue"].filter(c => c !== anyonColor);
      const left = sites.filter(id => siteById.get(id).colorName === colours[0]);
      const right = sites.filter(id => siteById.get(id).colorName === colours[1]);
      const occupied = new Set([
        ...currentCZ.map(pair => pairKey(...pair)),
        ...czPairs.filter(pair => !mAutoCZ.some(m => pairKey(...m) === pairKey(...pair)))
                 .map(pair => pairKey(...pair)),
      ]);
      const edges = [];
      for (const a of left) {
        for (const b of right) {
          if (occupied.has(pairKey(a,b))) continue;
          const sa = siteById.get(a), sb = siteById.get(b);
          const dist = Math.hypot(sa.x - sb.x, sa.y - sb.y);
          edges.push({pair:[a,b],dist});
        }
      }
      edges.sort((a,b) => a.dist - b.dist || a.pair[0]-b.pair[0] || a.pair[1]-b.pair[1]);
      return edges;
    }

    function mNewTrianglePath(newQubit) {
      const ends = completeMTriangleEnds(newQubit);
      if (ends.length !== 2 || anyonQubits.includes(newQubit)) return null;
      if (!anyonQubits.length) return [...ends];
      if (anyonQubits.length === 1) {
        const firstEnds = completeMTriangleEnds(anyonQubits[0]);
        const shared = ends.filter(id => firstEnds.includes(id));
        if (shared.length !== 1) return null;
        const stationary = firstEnds.find(id => id !== shared[0]);
        const moving = ends.find(id => id !== shared[0]);
        return [stationary,shared[0],moving];
      }
      const last = mTrianglePath[mTrianglePath.length-1];
      if (!ends.includes(last)) return null;
      const moving = ends.find(id => id !== last);
      if (mTrianglePath.includes(moving)) return null;
      return [...mTrianglePath,moving];
    }

    function findMDecoration(xs, ts) {
      const previous = mAutoCZ.map(pair => [...pair]);
      const expected = mAllowedStarIds(xs,ts);
      const current = mStarResiduals(xs,previous);
      if (validMResidues(current,expected)) return previous;
      // Every even X-count extends by X alone. No CZ rearrangements.
      if (xs.length % 2 === 0) return null;

      // A successful odd-number extension repairs *two* intermediate stars.
      // The normal straight ribbon needs 1,2,3,... new CZs at odd steps
      // n=3,5,7,...; turns may admit fewer. Search only valid solutions.
      const budget = Math.floor((xs.length-1)/2);
      const otherStars = starData.filter(st => !expected.has(st.id));
      const shift = new Map(otherStars.map((s,i) => [s.id, BigInt(i * (maxQubitId + 1))]));
      let target = 0n;
      for (const s of otherStars) target ^= current.get(s.id) << shift.get(s.id);
      if (target === 0n) return null;

      const candidates = mCandidateCZEdges(xs,ts,previous);
      const columns = [];
      const dict = new Map();
      for (const candidate of candidates) {
        const [a,b] = candidate.pair;
        let column = 0n;
        for (const sId of eStarIncidence.get(a) || []) {
          const pos = shift.get(sId);
          if (pos !== undefined) column ^= reducedZSingle.get(b) << pos;
        }
        for (const sId of eStarIncidence.get(b) || []) {
          const pos = shift.get(sId);
          if (pos !== undefined) column ^= reducedZSingle.get(a) << pos;
        }
        if (column === 0n) continue;
        columns.push({...candidate, column});
        if (!dict.has(column)) dict.set(column, []);
        dict.get(column).push(columns.length - 1);
      }

      const verify = (indexes) => {
        const add = indexes.map(i => columns[i].pair);
        const cz = [...previous,...add];
        return validMResidues(mStarResiduals(xs,cz),expected) ? cz : null;
      };

      if (budget >= 1) {
        for (const i of dict.get(target) || []) {
          const found = verify([i]);
          if (found) return found;
        }
      }
      if (budget >= 2) {
        for (let i=0;i<columns.length;i++) {
          for (const j of dict.get(target ^ columns[i].column) || []) {
            if (i === j) continue;
            const found = verify([i,j]);
            if (found) return found;
          }
        }
      }
      if (budget === 3) {
        for (let i=0;i<columns.length;i++) {
          for (let j=i+1;j<columns.length;j++) {
            const missing = target ^ columns[i].column ^ columns[j].column;
            for (const k of dict.get(missing) || []) {
              if (k===i || k===j) continue;
              const found = verify([i,j,k]);
              if (found) return found;
            }
          }
        }
      }
      // A longer ribbon may need more than four new edges. Solve the same
      // GF(2) equations without a gate-count constraint as a fallback.
      if (budget >= 4) {
        const basis = new Map();
        function pivot(x) { return x.toString(2).length-1; }
        for (let i=0;i<columns.length;i++) {
          let v=columns[i].column, combo=1n<<BigInt(i);
          while (v !== 0n) {
            const p=pivot(v), row=basis.get(p);
            if (row) {v ^=row.vec;combo ^=row.combo;}
            else {basis.set(p,{vec:v,combo});break;}
          }
        }
        let residue=target,combo=0n;
        while(residue!==0n) {
          const row=basis.get(pivot(residue));
          if (!row) break;
          residue ^=row.vec;combo ^=row.combo;
        }
        if (residue===0n) {
          const chosen=[];
          for(let i=0;i<columns.length;i++) {
            if (combo & (1n<<BigInt(i))) chosen.push(i);
          }
          const found=verify(chosen);
          if (found)return found;
        }
      }
      return null;
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

    // ------------------------------------------------------------------
    // Three-colour straight ribbons (s1 / s2).
    // Accept a horizontal spine, or either 60-degree diagonal line.
    // Every selected Kagome edge belongs to one small triangle. Its
    // UP/DOWN orientation is determined geometrically from its horizontal
    // side (not from the apex relative to the selected DIAGONAL edge).
    // The left and right base sites are sorted by global x coordinate,
    // independent of click direction. s1 uses XR(left),RX(right) on UP
    // triangles and RX(left),XR(right) on DOWN triangles. s2 swaps them.
    // The third triangle vertex always carries X only. All CZ gates
    // remain a jointly optimized three-colour decoration applied first.
    // ------------------------------------------------------------------
    const sStraightNeighbors = new Map();
    const sSiteNeighbors = new Map(siteData.map(s => [s.id, []]));
    for (const [a,b] of edgePairs) {
      sSiteNeighbors.get(a).push(b);
      sSiteNeighbors.get(b).push(a);
    }

    // Locate the orientation of a small Kagome triangle by its horizontal
    // side. This works for horizontal AND diagonal selected edges.
    function sTriangleUp(vertices) {
      let horizontal=null;
      for (let i=0;i<3;i++) {
        for (let j=i+1;j<3;j++) {
          const a=siteById.get(vertices[i]), b=siteById.get(vertices[j]);
          if (Math.abs(a.y-b.y)<1e-5) horizontal=[vertices[i],vertices[j]];
        }
      }
      if (!horizontal) return null;
      const tip=vertices.find(id=>!horizontal.includes(id));
      return siteById.get(tip).y>siteById.get(horizontal[0]).y;
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
      const up=sTriangleUp([a,b,apex]);
      if (up===null) continue;
      const l=sa.x<sb.x?a:b,r=sa.x<sb.x?b:a;
      sStraightNeighbors.set(pairKey(a,b),{left:l,right:r,apex,up});
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
        const flip=(anyonKind==='s2');
        const leftOrder=(tri.up!==flip)?'XR':'RX';
        const rightOrder=leftOrder==='XR'?'RX':'XR';
        for (const [q,order] of [[tri.left,leftOrder],[tri.right,rightOrder]]) {
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
      anyonSelect.value="";
      updateAnyonSummary();
    }

    function setAnyon(value) {
      const match = /^(e|m)_([RGB])$/.exec(value);
      anyonKind = (value === "s1" || value === "s2") ? value : (match ? match[1] : null);
      anyonColor=match ? ({R:"Red",G:"Green",B:"Blue"}[match[2]]):null;
      anyonQubits=[];anyonStarPath=[];
      mTrianglePath=[];mAutoCZ=[];mOptionCache=null;
      sBasePath=[];sGuidedX=[];sAutoCZ=[];sCandidateCache=null;
      sGuideId = (anyonKind === "s1" || anyonKind === "s2") ? nextSGuideId++ : null;
      mGuideId = anyonKind === "m" ? nextMGuideId++ : null;
      pendingCZ=null;
      if (panMode) setPanMode(false);
      currentMode=anyonKind ? "ANYON":"Z";
      updateActiveModeButtons();
      updateAnyonSummary();
      render();
      const options=availableAnyonQubits();
      if (anyonKind==="s1" || anyonKind==="s2") {
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

    function setPanMode(enabled) {
      panMode = enabled;

      if (panMode) {
        panBtn.classList.add("active");

        Plotly.relayout(plotDiv, {
          dragmode: "pan"
        });

        setStatus(
          "Pan mode: click and drag the lattice to move around."
        );
      } else {
        panBtn.classList.remove("active");

        Plotly.relayout(plotDiv, {
          dragmode: "zoom"
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

      const layout = {
          // Keep the user's current zoom/pan whenever Plotly.react()
          // is called after X, Z, CZ, Erase, Undo, etc.
          uirevision: "keep-kagome-view",

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
        if (panMode) {
          setPanMode(false);
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
        if (panMode && expectationMode) setPanMode(false);
        updateExpectationButtons();
        render();
        setStatus(expectationMode
          ? `Click inside any complete ${expectationMode==="STAR"?"hexagonal star":"triangle"} to calculate its expectation.`
          : "Expectation selection off; the last calculation remains in the right panel.");
      });
    }

    // Unlike Plotly point-click, this also catches empty polygon interiors.
    plotDiv.addEventListener("click", (evt) => {
      if (!expectationMode || panMode) return;
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

    panBtn.addEventListener("click", () => {
      setPanMode(!panMode);
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

    initializeQubitSelector();
    updateTrackedTable();
    updateVisibilityButtons();
    updateViolationButtons();
    updateExpectationButtons();
    render();

    plotDiv.on("plotly_click", (eventData) => {
      if (panMode) {
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
