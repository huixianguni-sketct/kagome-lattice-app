from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np


# Global three-colouring of the Kagome vertices.
#
# IMPORTANT: colour is NOT tied to the fixed A/B/C basis label.  If it were,
# straight Kagome chains would alternate between only two colours.  Instead we
# colour the Kagome sites using integer coordinates (u, v) of the underlying
# triangular grid:
#
#     e1 = (1, 0),  e2 = (1/2, sqrt(3)/2)
#
# and assign colour from (u - v) mod 3.  Therefore moving along any of the
# three Kagome line directions cycles through all three colours.  With the
# mapping below, a horizontal chain reads
#
#     Red -> Blue -> Green -> Red -> ...
#
# matching the visible horizontal colour sequence in Fig. 2 of the paper.
COLOR_BY_INDEX = {
    0: "Red",
    1: "Blue",
    2: "Green",
}


@dataclass(frozen=True)
class Site:
    """Metadata for one qubit on the Kagome lattice."""

    site_id: int
    x: float
    y: float
    cell_i: int
    cell_j: int
    sublattice: str
    color: str


@dataclass
class KagomeLattice:
    """Geometry and graph representation of a finite Kagome-lattice patch."""

    sites: List[Site]
    positions: np.ndarray  # shape (N, 2)
    edges: List[Tuple[int, int]]
    neighbors: Dict[int, List[int]]

    @property
    def n_sites(self) -> int:
        return len(self.sites)


def _vertex_color(cell_i: int, cell_j: int, sublattice: str) -> str:
    """
    Return the global Red/Blue/Green colour of a Kagome vertex.

    In fine triangular-grid coordinates (u, v):
        A(i,j) = (2i,   2j)
        B(i,j) = (2i+1, 2j)
        C(i,j) = (2i,   2j+1)

    The colour index c = (u - v) mod 3 changes nontrivially along every one of
    the three straight lattice directions, so every straight chain cycles
    through all three colours instead of alternating between only two.
    """
    if sublattice == "A":
        u, v = 2 * cell_i, 2 * cell_j
    elif sublattice == "B":
        u, v = 2 * cell_i + 1, 2 * cell_j
    elif sublattice == "C":
        u, v = 2 * cell_i, 2 * cell_j + 1
    else:
        raise ValueError(f"Unknown Kagome sublattice: {sublattice}")

    return COLOR_BY_INDEX[(u - v) % 3]


def generate_kagome(nx: int, ny: int) -> KagomeLattice:
    """
    Generate an open-boundary Kagome lattice.

    Primitive vectors:
        a1 = (2, 0)
        a2 = (1, sqrt(3))

    Three-site geometric basis:
        A = (0, 0)
        B = (1, 0)
        C = (1/2, sqrt(3)/2)

    The A/B/C labels describe geometry only.  The Red/Blue/Green qubit colour
    is assigned globally by _vertex_color(), so the colour pattern shifts from
    unit cell to unit cell.

    The nearest-neighbour distance is 1 in these units.
    """
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must both be >= 1")

    sqrt3 = np.sqrt(3.0)
    a1 = np.array([2.0, 0.0])
    a2 = np.array([1.0, sqrt3])

    basis = {
        "A": np.array([0.0, 0.0]),
        "B": np.array([1.0, 0.0]),
        "C": np.array([0.5, sqrt3 / 2.0]),
    }

    sites: List[Site] = []
    positions: List[np.ndarray] = []

    site_id = 0
    for j in range(ny):
        for i in range(nx):
            cell_origin = i * a1 + j * a2
            for sublattice, offset in basis.items():
                r = cell_origin + offset
                sites.append(
                    Site(
                        site_id=site_id,
                        x=float(r[0]),
                        y=float(r[1]),
                        cell_i=i,
                        cell_j=j,
                        sublattice=sublattice,
                        color=_vertex_color(i, j, sublattice),
                    )
                )
                positions.append(r)
                site_id += 1

    pos = np.asarray(positions, dtype=float)

    # Build the nearest-neighbour graph geometrically.
    delta = pos[:, None, :] - pos[None, :, :]
    dist2 = np.einsum("ijk,ijk->ij", delta, delta)
    nearest = np.isclose(dist2, 1.0, atol=1e-9, rtol=1e-9)
    upper = np.triu(nearest, k=1)
    rows, cols = np.where(upper)
    edges = [(int(i), int(j)) for i, j in zip(rows, cols)]

    neighbors: Dict[int, List[int]] = {i: [] for i in range(len(sites))}
    for i, j in edges:
        neighbors[i].append(j)
        neighbors[j].append(i)

    for site_neighbors in neighbors.values():
        site_neighbors.sort()

    return KagomeLattice(
        sites=sites,
        positions=pos,
        edges=edges,
        neighbors=neighbors,
    )
