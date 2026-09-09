from __future__ import annotations

import streamlit as st
import streamlit.components.v1 as components

from lattice import generate_kagome
from ui_component import build_interactive_lattice_html


st.set_page_config(
    page_title="Interactive Kagome Qubits",
    page_icon="🔺",
    layout="wide",
)


with st.sidebar:
    st.header("Lattice")
    nx = st.slider("Unit cells along $a_1$", min_value=1, max_value=20, value=12)
    ny = st.slider("Unit cells along $a_2$", min_value=1, max_value=20, value=12)


st.title("Interactive Kagome Qubit Lattice")
st.caption(
    "Choose Z, X, or CZ above the lattice. X/Z act on one clicked qubit; "
    "CZ is completed after choosing two qubits."
)

lattice = generate_kagome(nx, ny)
component_html = build_interactive_lattice_html(lattice)

# Vertex interaction happens inside this client-side Plotly component.  Qubit
# clicks therefore do not trigger a Streamlit Python rerun.  Changing nx or ny
# does rerun the app and currently resets the operator annotations.
components.html(component_html, height=790, scrolling=False)
