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
    Find complete Kagome hexagons (open-boundary stars).

    A valid star is an induced 6-cycle of nearest-neighbour Kagome bonds.
    Any hexagon cut by the finite boundary is therefore automatically ignored.
    """
    sites = sorted(lattice.sites, key=lambda s: int(s.site_id))
    site_by_id = {int(s.site_id): s for s in sites}
    neighbors = {
        int(q): [int(v) for v in values]
        for q, values in lattice.neighbors.items()
    }

    cycles = set()

    for start in sorted(neighbors):

        def dfs(current, path):
            if len(path) == 6:
                if start in neighbors[current]:
                    cyc = _canonical_cycle(path)
                    subset = set(cyc)

                    # An elementary Kagome hexagon is an induced 6-cycle:
                    # exactly six graph edges occur inside the six vertices.
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

                # The smallest vertex ID is used as the DFS start.
                if nxt < start:
                    continue

                dfs(nxt, path + [nxt])

        dfs(start, [start])

    all_colors = {"Red", "Green", "Blue"}
    stars = []

    for cyc in cycles:
        cx = sum(float(site_by_id[q].x) for q in cyc) / 6.0
        cy = sum(float(site_by_id[q].y) for q in cyc) / 6.0

        # Stable cyclic order around the hexagon.
        ordered = sorted(
            cyc,
            key=lambda q: math.atan2(
                float(site_by_id[q].y) - cy,
                float(site_by_id[q].x) - cx,
            ),
        )

        # A star of one colour is surrounded by the other two colours.
        present = {str(site_by_id[q].color) for q in ordered}
        missing = all_colors - present
        star_color = next(iter(missing)) if len(missing) == 1 else "Unknown"

        stars.append(
            {
                "qubits": [int(q) for q in ordered],
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
    star_data = _detect_complete_stars(lattice)
    triangle_data = _detect_complete_triangles(lattice)

    site_data = [
        {
            "id": int(site.site_id),
            "x": float(site.x),
            "y": float(site.y),
            "color": _site_color(site),
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
    <div class="toolbar">
      <span class="toolbar-label">Operator</span>

      <button class="btn mode-btn active" data-mode="Z">Z</button>
      <button class="btn mode-btn" data-mode="X">X</button>
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
    const modeButtons = Array.from(document.querySelectorAll(".mode-btn"));
    const visibilityButtons = Array.from(document.querySelectorAll(".visibility-btn"));
    const violationButtons = Array.from(document.querySelectorAll(".violation-btn"));
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

    // Display-only visibility state. Turning one of these off does not
    // remove the operator from the lattice state or tracked-qubit table.
    let operatorVisibility = {
      Z: true,
      X: true,
      CZ: true,
    };

    // Violation overlays are display-only.  The physics is always evaluated
    // from the active ordered operationLog.
    let violationVisibility = {
      STAR: false,
      TRIANGLE: false,
    };

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

    function recordOperation(type, qubits) {
      operationLog.push({
        id: nextOperationId++,
        type: type,
        qubits: [...qubits],
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
        operationLog: operationLog.map(event => ({
          id: event.id,
          type: event.type,
          qubits: [...event.qubits],
        })),
        nextOperationId: nextOperationId,
      };
    }

    function restoreState(state) {
      localOps = new Map(state.localOps);
      czPairs = state.czPairs.map(([a, b]) => [a, b]);
      pendingCZ = state.pendingCZ;
      currentMode = state.currentMode;
      operationLog = state.operationLog.map(event => ({
        id: event.id,
        type: event.type,
        qubits: [...event.qubits],
      }));
      nextOperationId = state.nextOperationId;

      updateActiveModeButtons();
      render();
    }

    function pushHistory(description) {
      history.push({
        state: snapshot(),
        description: description,
      });
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
    // As = (Π CZ around a complete hexagon) X⊗6.
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
      const cycle = star.qubits;
      const starSet = new Set(cycle);
      const cycleIndex = new Map(cycle.map((q, index) => [q, index]));

      let sign = 1;
      let zMask = 0n;

      // If the user applied gates G1, G2, ..., Gn, then
      //
      // U† As U = G1†(...Gn† As Gn...)G1,
      //
      // so conjugation is evaluated in reverse chronological order.
      for (let opIndex = operationLog.length - 1; opIndex >= 0; opIndex--) {
        const event = operationLog[opIndex];

        if (event.type === "Z" && event.qubits.length === 1) {
          const q = event.qubits[0];

          // Zq anticommutes with the Xq factor of As.
          if (starSet.has(q)) {
            sign *= -1;
          }

          continue;
        }

        if (event.type === "X" && event.qubits.length === 1) {
          const q = event.qubits[0];

          // Conjugating an already-present Zq by Xq contributes a minus sign.
          if ((zMask & qubitBit(q)) !== 0n) {
            sign *= -1;
          }

          if (starSet.has(q)) {
            const idx = cycleIndex.get(q);
            const previous = cycle[(idx + cycle.length - 1) % cycle.length];
            const next = cycle[(idx + 1) % cycle.length];

            // Xq conjugates the two CZs of As touching q, producing Z on
            // the two neighbouring star vertices.
            zMask ^= qubitBit(previous);
            zMask ^= qubitBit(next);
          }

          continue;
        }

        if (event.type === "CZ" && event.qubits.length === 2) {
          const a = event.qubits[0];
          const b = event.qubits[1];

          const aInStar = starSet.has(a);
          const bInStar = starSet.has(b);

          // CZ_ab X_a X_b CZ_ab = - X_a X_b Z_a Z_b.
          if (aInStar && bInStar) {
            sign *= -1;
          }

          // If one endpoint carries an X factor from As, conjugation by CZ
          // decorates it with Z on the opposite endpoint.
          if (aInStar) {
            zMask ^= qubitBit(b);
          }

          if (bInStar) {
            zMask ^= qubitBit(a);
          }
        }
      }

      // The complete open-boundary Bt triangles define the diagonal +1
      // constraints used here.  If the residual Z string is not generated by
      // them, its ground-state expectation vanishes.
      if (!isTriangleStabilizerMask(zMask)) {
        return 0;
      }

      return sign;
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
          const path = polygonPath(item.star.qubits);

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
      const xs = [];
      const ys = [];
      const texts = [];
      const customdata = [];

      // Offset from the center of the qubit.
      // Nearest-neighbour spacing is ~1, so these keep
      // the labels close while clearly outside the vertex.
      const offsetX = 0.18;
      const offsetY = 0.25;

      for (const [siteId, op] of localOps.entries()) {
        const s = siteById.get(siteId);
        if (!s) continue;

        // Visibility only changes what is drawn; the operator remains active.
        if (!operatorVisibility[op]) continue;

        xs.push(s.x + offsetX);
        ys.push(s.y + offsetY);

        texts.push(op);
        customdata.push(siteId);
      }

      return {
        xs,
        ys,
        texts,
        customdata
      };
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
        type: "scatter",
        mode: "markers",
        x: siteData.map(s => s.x),
        y: siteData.map(s => s.y),
        customdata: siteData.map(s => s.id),
        marker: {
          size: 14,
          color: siteData.map(s => s.color),
          line: {
            width: 1.1,
            color: "#111827",
          },
        },
        hovertemplate: "<b>Qubit %{customdata}</b><extra></extra>",
        showlegend: false,
      };

      const opLabelTrace = {
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
      const connectedPairs = czPairs.filter(([a, b]) => a === siteId || b === siteId);
      const wasPending = pendingCZ === siteId;

      if (!existingOp && connectedPairs.length === 0 && !wasPending) {
        setStatus(`Qubit ${siteId} has no operations to erase.`);
        return;
      }

      const descriptionParts = [];

      if (existingOp) {
        descriptionParts.push(`${existingOp} on qubit ${siteId}`);
      }

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
      if (localOps.size === 0 && czPairs.length === 0 && pendingCZ === null) {
        setStatus("Nothing to clear.");
        return;
      }

      pushHistory("clearing all operations");

      localOps.clear();
      czPairs = [];
      pendingCZ = null;
      operationLog = [];
      nextOperationId = 1;

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
      if (currentMode === "Z") {
        handleLocalOp(siteId, "Z");
      } else if (currentMode === "X") {
        handleLocalOp(siteId, "X");
      } else if (currentMode === "CZ") {
        handleCZ(siteId);
      } else if (currentMode === "ERASE") {
        handleErase(siteId);
      }
    }

    // Mode button events
    for (const btn of modeButtons) {
      btn.addEventListener("click", () => {

        // Return from navigation mode to operator mode.
        if (panMode) {
          setPanMode(false);
        }

        currentMode = btn.dataset.mode;

        updateActiveModeButtons();

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
    render();

    plotDiv.on("plotly_click", (eventData) => {
      if (panMode) {
        return;
      }
      if (!eventData || !eventData.points || eventData.points.length === 0) {
        return;
      }

      const point = eventData.points[0];
      const curveNumber = point.curveNumber;

      // Clicks are accepted on:
      // 3 = node markers
      // 4 = local operator text
      // 5 = pending marker (if present)
      const clickableCurveNumbers = new Set([3, 4, 5]);

      if (!clickableCurveNumbers.has(curveNumber)) {
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
