from __future__ import annotations

import json


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

    #plot {
      width: 100%;
      height: 760px;
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

      <button class="btn action" id="save-btn">Save PNG</button>
      <button class="btn action" id="undo-btn">Undo</button>
      <button class="btn action" id="clear-btn">Clear</button>
    </div>

    <div class="status" id="status">Mode: Z</div>
    <div id="plot"></div>
  </div>

  <script>
    const siteData = __SITE_DATA__;
    const edgePairs = __EDGE_PAIRS__;
    const baseEdgeX = __BASE_EDGE_X__;
    const baseEdgeY = __BASE_EDGE_Y__;

    const plotDiv = document.getElementById("plot");
    const statusDiv = document.getElementById("status");
    const modeButtons = Array.from(document.querySelectorAll(".mode-btn"));
    const saveBtn = document.getElementById("save-btn");
    const undoBtn = document.getElementById("undo-btn");
    const clearBtn = document.getElementById("clear-btn");

    const siteById = new Map(siteData.map(s => [s.id, s]));

    const nearestNeighborKeys = new Set(
      edgePairs.map(([a, b]) => pairKey(a, b))
    );

    let currentMode = "Z";
    let localOps = new Map();   // site_id -> "X" or "Z"
    let czPairs = [];           // array of [a, b]
    let pendingCZ = null;       // first endpoint for CZ selection
    let history = [];           // array of { state, description }

    function setStatus(text) {
      statusDiv.textContent = text;
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
      };
    }

    function restoreState(state) {
      localOps = new Map(state.localOps);
      czPairs = state.czPairs.map(([a, b]) => [a, b]);
      pendingCZ = state.pendingCZ;
      currentMode = state.currentMode;

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
        paper_bgcolor: "white",
        plot_bgcolor: "white",
        margin: { l: 10, r: 10, t: 10, b: 10 },
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
        localOps.delete(siteId);
        render();
        setStatus(`Removed ${newOp} from qubit ${siteId}.`);
        return;
      }

      if (oldOp && oldOp !== newOp) {
        pushHistory(`replacing ${oldOp} with ${newOp} on qubit ${siteId}`);
        localOps.set(siteId, newOp);
        render();
        setStatus(`Replaced ${oldOp} with ${newOp} on qubit ${siteId}.`);
        return;
      }

      pushHistory(`applying ${newOp} to qubit ${siteId}`);
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
        pendingCZ = null;
        render();
        setStatus(`Removed CZ between qubits ${a} and ${b}.`);
        return;
      }

      pushHistory(`applying CZ to qubits ${a} and ${b}`);
      czPairs.push([a, b]);
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
        currentMode = btn.dataset.mode;
        updateActiveModeButtons();

        if (currentMode === "CZ") {
          setStatus("Mode: CZ. Click two qubits to connect them.");
        } else if (currentMode === "ERASE") {
          setStatus("Mode: Erase. Click a qubit to remove its local operator and any CZ connected to it.");
        } else {
          setStatus(`Mode: ${currentMode}`);
        }
      });
    }

    saveBtn.addEventListener("click", handleSave);
    undoBtn.addEventListener("click", handleUndo);
    clearBtn.addEventListener("click", handleClear);

    render();

    plotDiv.on("plotly_click", (eventData) => {
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
    )

    return html
