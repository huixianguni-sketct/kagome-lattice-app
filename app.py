from __future__ import annotations

from typing import Iterable, Set

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from lattice import KagomeLattice, generate_kagome
from physics import (
    adjacency_matrix,
    count_selected_bonds,
    tight_binding_hamiltonian,
)


# -----------------------------------------------------------------------------
# Streamlit configuration
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Interactive Kagome Lattice",
    page_icon="🔺",
    layout="wide",
)


# -----------------------------------------------------------------------------
# Session-state helpers
# -----------------------------------------------------------------------------
def initialize_state() -> None:
    defaults = {
        "selected_sites": set(),
        "lattice_signature": None,
        "chart_version": 0,
        "calculation": None,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def reset_selection() -> None:
    st.session_state.selected_sites = set()
    st.session_state.chart_version += 1
    st.session_state.calculation = None


def sync_lattice_signature(nx: int, ny: int) -> None:
    """Clear site-dependent state when lattice dimensions change."""
    signature = (nx, ny)
    if st.session_state.lattice_signature != signature:
        st.session_state.lattice_signature = signature
        reset_selection()


# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------
def build_lattice_figure(
    lattice: KagomeLattice,
    selected_sites: Iterable[int],
    show_labels: bool,
    marker_size: int,
) -> go.Figure:
    selected_sites = sorted(set(int(i) for i in selected_sites))

    # One Plotly line trace for all bonds.
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
        line=dict(width=1.5, color="rgba(80,80,80,0.65)"),
        hoverinfo="skip",
        showlegend=False,
        name="Bonds",
    )

    # customdata travels back with Plotly selection events.
    degrees = lattice.degrees
    customdata = np.array(
        [
            [
                site.site_id,
                site.cell_i,
                site.cell_j,
                site.sublattice,
                int(degrees[site.site_id]),
            ]
            for site in lattice.sites
        ],
        dtype=object,
    )

    node_trace = go.Scatter(
        x=lattice.positions[:, 0],
        y=lattice.positions[:, 1],

        mode="markers+text",

        text=[
            str(site.site_id)
            for site in lattice.sites
        ],

        textposition="top center",

        customdata=customdata,

        selectedpoints=list(selected_sites),

        marker=dict(
            size=14,
            color="#3b82f6",

            # line IS allowed here
            line=dict(
                width=1.3,
                color="#111827"
            ),
        ),

        # But selected.marker can only have
        # color, size and opacity
        selected=dict(
            marker=dict(
                color="#ef4444",
                size=18,
                opacity=1.0,
            )
        ),

        unselected=dict(
            marker=dict(
                opacity=0.75,
            )
        ),

        hovertemplate=(
            "<b>Site %{customdata[0]}</b><br>"
            "cell = (%{customdata[1]}, %{customdata[2]})<br>"
            "sublattice = %{customdata[3]}<br>"
            "degree = %{customdata[4]}<br>"
            "x = %{x:.3f}<br>"
            "y = %{y:.3f}"
            "<extra></extra>"
        ),

        showlegend=False,
    )

    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        height=700,
        margin=dict(l=10, r=10, t=20, b=10),
        plot_bgcolor="white",
        paper_bgcolor="white",
        clickmode="event+select",
        dragmode="pan",
        uirevision="kagome-lattice",
        xaxis=dict(
            visible=False,
            showgrid=False,
            zeroline=False,
        ),
        yaxis=dict(
            visible=False,
            showgrid=False,
            zeroline=False,
            scaleanchor="x",
            scaleratio=1,
        ),
    )

    return fig


def selected_ids_from_plotly_event(event) -> Set[int]:
    """Extract our site IDs from Streamlit's Plotly selection state."""
    try:
        points = event.selection.points
    except (AttributeError, TypeError):
        try:
            points = event.get("selection", {}).get("points", [])
        except AttributeError:
            points = []

    site_ids: Set[int] = set()
    for point in points:
        try:
            customdata = point["customdata"]
        except (KeyError, TypeError):
            continue

        # Bond trace has no customdata; only the site trace reaches here.
        if customdata is not None and len(customdata) > 0:
            site_ids.add(int(customdata[0]))

    return site_ids


def site_table(lattice: KagomeLattice, selected_sites: Iterable[int]) -> pd.DataFrame:
    rows = []
    for site_id in sorted(set(int(i) for i in selected_sites)):
        site = lattice.sites[site_id]
        rows.append(
            {
                "site_id": site.site_id,
                "cell_i": site.cell_i,
                "cell_j": site.cell_j,
                "sublattice": site.sublattice,
                "x": round(site.x, 6),
                "y": round(site.y, 6),
                "degree": len(lattice.neighbors[site_id]),
                "neighbors": ", ".join(map(str, lattice.neighbors[site_id])),
            }
        )
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# App
# -----------------------------------------------------------------------------
initialize_state()

st.title("Interactive Kagome Lattice")
st.caption(
    "Click a site to select it. Use Shift-click to build a multi-site selection. "
    "You can also use Plotly's box/lasso tools from the toolbar."
)

with st.sidebar:
    st.header("Lattice")
    nx = st.slider("Unit cells along $a_1$", min_value=1, max_value=12, value=5)
    ny = st.slider("Unit cells along $a_2$", min_value=1, max_value=12, value=4)
    marker_size = st.slider("Site marker size", min_value=8, max_value=28, value=14)
    show_labels = st.checkbox("Show site IDs", value=True)

    st.divider()
    st.header("Example Hamiltonian")
    hopping = st.number_input("Hopping $t$", value=1.0, step=0.1, format="%.3f")
    onsite = st.number_input("Uniform onsite $\\epsilon$", value=0.0, step=0.1, format="%.3f")
    selected_potential = st.number_input(
        "Extra potential $V$ on selected sites",
        value=2.0,
        step=0.1,
        format="%.3f",
    )

    st.divider()
    if st.button("Clear selection", use_container_width=True):
        reset_selection()
        st.rerun()

sync_lattice_signature(nx, ny)
lattice = generate_kagome(nx, ny)

left, right = st.columns([4.2, 1.8], gap="large")

with left:
    st.subheader("Lattice")
    fig = build_lattice_figure(
        lattice=lattice,
        selected_sites=st.session_state.selected_sites,
        show_labels=show_labels,
        marker_size=marker_size,
    )

    event = st.plotly_chart(
        fig,
        use_container_width=True,
        key=f"kagome_plot_{st.session_state.chart_version}",
        on_select="rerun",
        selection_mode=("points", "box", "lasso"),
        config={
            "displaylogo": False,
            "scrollZoom": True,
        },
    )

    # Treat Plotly's current non-empty selection as authoritative.
    # Empty selections are handled explicitly by the Clear button so that
    # unrelated Streamlit reruns do not accidentally erase user state.
    plot_selection = selected_ids_from_plotly_event(event)
    if plot_selection:
        st.session_state.selected_sites = plot_selection

with right:
    st.subheader("Current state")
    c1, c2 = st.columns(2)
    c1.metric("Sites", lattice.n_sites)
    c2.metric("Bonds", len(lattice.edges))

    selected = sorted(st.session_state.selected_sites)
    c3, c4 = st.columns(2)
    c3.metric("Selected", len(selected))
    c4.metric(
        "Selected bonds",
        count_selected_bonds(lattice.edges, selected),
    )

    st.markdown("**Selected site IDs**")
    if selected:
        st.code(", ".join(map(str, selected)), language=None)
    else:
        st.info("No sites selected yet.")

    st.markdown("**Inspect one site**")
    inspect_id = st.number_input(
        "Site ID",
        min_value=0,
        max_value=lattice.n_sites - 1,
        value=0,
        step=1,
        label_visibility="collapsed",
    )
    inspect_id = int(inspect_id)
    site = lattice.sites[inspect_id]
    st.write(
        {
            "site_id": site.site_id,
            "cell": (site.cell_i, site.cell_j),
            "sublattice": site.sublattice,
            "position": (round(site.x, 4), round(site.y, 4)),
            "neighbors": lattice.neighbors[inspect_id],
        }
    )

st.divider()

selected = sorted(st.session_state.selected_sites)
st.subheader("Selected-site data")
if selected:
    st.dataframe(
        site_table(lattice, selected),
        use_container_width=True,
        hide_index=True,
    )
else:
    st.caption("Select one or more lattice sites to populate this table.")

with st.expander("Graph / adjacency data"):
    A = adjacency_matrix(lattice.n_sites, lattice.edges)
    st.write(
        "The graph representation is what you will normally feed into later "
        "physics calculations. `A[i, j] = 1` means sites `i` and `j` are nearest neighbors."
    )

    if lattice.n_sites <= 90:
        st.dataframe(pd.DataFrame(A), use_container_width=True)
    else:
        st.info(
            f"The full adjacency matrix is {lattice.n_sites}×{lattice.n_sites}. "
            "Showing the top-left 30×30 block to keep the browser responsive."
        )
        st.dataframe(pd.DataFrame(A[:30, :30]), use_container_width=True)

st.subheader("Example calculation")
st.write(
    "This demonstrates how the clicked site IDs can feed directly into a calculation. "
    "Selected sites receive an additional onsite potential $V$."
)

if st.button("Run tight-binding calculation", type="primary"):
    H = tight_binding_hamiltonian(
        n_sites=lattice.n_sites,
        edges=lattice.edges,
        hopping=hopping,
        onsite_energy=onsite,
        selected_sites=selected,
        selected_site_potential=selected_potential,
    )
    eigenvalues = np.linalg.eigvalsh(H)

    st.session_state.calculation = {
        "signature": (
            nx,
            ny,
            float(hopping),
            float(onsite),
            float(selected_potential),
            tuple(selected),
        ),
        "H": H,
        "eigenvalues": eigenvalues,
    }

calc = st.session_state.calculation
current_calc_signature = (
    nx,
    ny,
    float(hopping),
    float(onsite),
    float(selected_potential),
    tuple(selected),
)

if calc is not None:
    if calc["signature"] != current_calc_signature:
        st.warning("The lattice/parameters changed. Run the calculation again for fresh results.")
    else:
        eigenvalues = calc["eigenvalues"]
        H = calc["H"]

        r1, r2, r3 = st.columns(3)
        r1.metric("Lowest eigenvalue", f"{eigenvalues[0]:.6f}")
        r2.metric("Highest eigenvalue", f"{eigenvalues[-1]:.6f}")
        r3.metric("Bandwidth", f"{eigenvalues[-1] - eigenvalues[0]:.6f}")

        spectrum = go.Figure(
            go.Scatter(
                x=np.arange(len(eigenvalues)),
                y=eigenvalues,
                mode="markers",
                marker=dict(size=6),
                hovertemplate="state %{x}<br>E = %{y:.6f}<extra></extra>",
            )
        )
        spectrum.update_layout(
            title="Single-particle spectrum",
            xaxis_title="Eigenstate index",
            yaxis_title="Energy",
            height=420,
            margin=dict(l=30, r=20, t=55, b=35),
        )
        st.plotly_chart(spectrum, use_container_width=True)

        with st.expander("Hamiltonian matrix"):
            if lattice.n_sites <= 90:
                st.dataframe(pd.DataFrame(H), use_container_width=True)
            else:
                st.info(
                    f"The Hamiltonian is {lattice.n_sites}×{lattice.n_sites}; "
                    "showing its top-left 30×30 block."
                )
                st.dataframe(pd.DataFrame(H[:30, :30]), use_container_width=True)

st.divider()
st.caption(
    "Architecture: lattice.py owns geometry/graph data; physics.py owns calculations; "
    "app.py owns only the UI and interaction state."
)
