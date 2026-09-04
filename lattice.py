from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np


@dataclass(frozen=True)
class Site:
    """Metadata for one Kagome-lattice site."""

    site_id: int
    x: float
    y: float
    cell_i: int
    cell_j: int
    sublattice: str


@dataclass
class KagomeLattice:
    """Geometry + graph representation of a finite open-boundary Kagome lattice."""

    sites: List[Site]
    positions: np.ndarray  # shape (N, 2)
    edges: List[Tuple[int, int]]
    neighbors: Dict[int, List[int]]

    @property
    def n_sites(self) -> int:
        return len(self.sites)

    @property
    def degrees(self) -> np.ndarray:
        return np.array([len(self.neighbors[i]) for i in range(self.n_sites)], dtype=int)


def generate_kagome(nx: int, ny: int) -> KagomeLattice:
    """
    Generate a finite Kagome lattice with open boundaries.

    We use a triangular Bravais lattice with primitive vectors

        a1 = (2, 0)
        a2 = (1, sqrt(3))

    and a three-site basis

        A = (0, 0)
        B = (1, 0)
        C = (1/2, sqrt(3)/2).

    With this convention, nearest-neighbor distance = 1.

    Parameters
    ----------
    nx, ny:
        Number of Bravais unit cells along a1 and a2.

    Returns
    -------
    KagomeLattice
        Site metadata, coordinates, nearest-neighbor edges, and neighbor lists.
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
                    )
                )
                positions.append(r)
                site_id += 1

    pos = np.asarray(positions, dtype=float)

    # Build nearest-neighbor graph from geometry.
    # For the chosen basis, NN distance is exactly 1 up to floating-point error.
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
