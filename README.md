# Interactive Kagome Qubit Lattice

A Streamlit + Plotly interface for annotating qubits on a three-coloured Kagome lattice.

## Current interaction

- Qubits are coloured Red / Green / Blue.
- Default lattice size: 8 unit cells along `a1` and 8 along `a2` (192 qubits).
- `Z`: click a qubit to toggle a Z label.
- `X`: click a qubit to toggle an X label.
- `CZ`: click the first qubit and then the second; both endpoints are labelled and a CZ link is drawn.
- `Undo` and `Clear` are available above the lattice.
- Qubit clicks are handled client-side, so they do not rerun Streamlit.
- Changing `nx` or `ny` rebuilds the lattice and resets annotations.

This version intentionally does **not** implement Hamiltonians, state vectors, stabilizer values, anyons, or tight-binding physics.

## Run locally

```bash
python -m venv .venv
```

Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
streamlit run app.py
```

## Project structure

```text
app.py              Streamlit page and nx/ny sidebar
lattice.py          Kagome geometry and three-colouring
ui_component.py     Plotly figure + client-side X/Z/CZ interaction
requirements.txt
.gitignore
README.md
```
