# Interactive Kagome Lattice — Streamlit + Plotly

## Setup

```bash
python -m venv .venv
```

Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
streamlit run app.py
```

macOS/Linux:

```bash
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

The app uses Streamlit's built-in Plotly selection events. Click a site to select it; Shift-click can be used for multi-selection. The toolbar also exposes box/lasso selection.

## Files

- `app.py` — Streamlit UI and interaction state
- `lattice.py` — Kagome geometry, edges, neighbors, site metadata
- `physics.py` — adjacency matrix and example tight-binding Hamiltonian
