from __future__ import annotations

from typing import Iterable, Sequence, Tuple

import numpy as np


def adjacency_matrix(n_sites: int, edges: Sequence[Tuple[int, int]]) -> np.ndarray:
    """Return the symmetric 0/1 nearest-neighbor adjacency matrix."""
    A = np.zeros((n_sites, n_sites), dtype=int)
    for i, j in edges:
        A[i, j] = 1
        A[j, i] = 1
    return A


def tight_binding_hamiltonian(
    n_sites: int,
    edges: Sequence[Tuple[int, int]],
    hopping: float = 1.0,
    onsite_energy: float = 0.0,
    selected_sites: Iterable[int] = (),
    selected_site_potential: float = 0.0,
) -> np.ndarray:
    r"""
    Build a simple single-particle nearest-neighbor tight-binding Hamiltonian.

        H = -t sum_<ij> (|i><j| + |j><i|)
            + eps sum_i |i><i|
            + V sum_{i in selected} |i><i|.

    This is only an example calculation layer. Replace/extend it later with
    your spin, stabilizer, anyon, or many-body Hamiltonian code.
    """
    H = np.eye(n_sites, dtype=float) * onsite_energy

    for i, j in edges:
        H[i, j] = -hopping
        H[j, i] = -hopping

    for i in selected_sites:
        if not 0 <= int(i) < n_sites:
            raise IndexError(f"Selected site {i} is outside 0..{n_sites - 1}")
        H[int(i), int(i)] += selected_site_potential

    return H


def count_selected_bonds(
    edges: Sequence[Tuple[int, int]], selected_sites: Iterable[int]
) -> int:
    """Count bonds whose two endpoints are both selected."""
    selected = set(int(i) for i in selected_sites)
    return sum(1 for i, j in edges if i in selected and j in selected)
