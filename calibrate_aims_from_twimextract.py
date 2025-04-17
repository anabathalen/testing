import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_ciu_heatmap(df):
    st.subheader("CIU Heatmap")

    # === User controls ===
    colormap = st.selectbox("Select Colormap", plt.colormaps(), index=plt.colormaps().index("viridis"))
    font_size = st.slider("Font size", 8, 24, 12)
    fig_width = st.slider("Figure width", 4, 20, 10)
    fig_height = st.slider("Figure height", 4, 20, 6)
    dpi = st.slider("Figure resolution (dpi)", 100, 1000, 300)
    
    # Optional cropping
    crop_x = st.checkbox("Crop X-axis (Collision Voltage)")
    x_min, x_max = None, None
    if crop_x:
        x_min = st.number_input("X min", value=float(df['Collision Voltage'].min()))
        x_max = st.number_input("X max", value=float(df['Collision Voltage'].max()))
    
    crop_y = st.checkbox("Crop Y-axis (CCS)")
    y_min, y_max = None, None
    if crop_y:
        y_min = st.number_input("Y min", value=float(df['CCS'].min()))
        y_max = st.number_input("Y max", value=float(df['CCS'].max()))

    # Labels
    label_cvs = st.multiselect("Label specific Collision Voltages", sorted(df['Collision Voltage'].unique().tolist()))
    label_ccs = st.multiselect("Label specific CCS values", sorted(df['CCS'].unique().tolist()))

    # === Prepare data ===
    pivot = df.pivot(index="CCS", columns="Collision Voltage", values="Intensity")
    pivot = pivot.sort_index(ascending=False)  # Flip CCS axis

    # Optional cropping
    if crop_x:
        pivot = pivot.loc[:, (pivot.columns >= x_min) & (pivot.columns <= x_max)]
    if crop_y:
        pivot = pivot.loc[(pivot.index >= y_min) & (pivot.index <= y_max), :]

    # === Plot ===
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)
    im = ax.imshow(pivot, aspect='auto', cmap=colormap,
                   extent=[pivot.columns.min(), pivot.columns.max(), pivot.index.min(), pivot.index.max()],
                   origin='lower')

    # Axis settings
    ax.set_xlabel("Collision Voltage (V)", fontsize=font_size)
    ax.set_ylabel("CCS (Å²)", fontsize=font_size)
    ax.tick_params(labelsize=font_size - 2)

    # Gridline labels
    for val in label_cvs:
        ax.axvline(x=val, color='white', linestyle='--', linewidth=0.8)
    for val in label_ccs:
        ax.axhline(y=val, color='white', linestyle='--', linewidth=0.8)

    # Smarter ticks
    ax.set_xticks(np.linspace(pivot.columns.min(), pivot.columns.max(), 6))
    ax.set_yticks(np.linspace(pivot.index.min(), pivot.index.max(), 6))

    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=font_size - 2)
    cbar.set_label("Intensity", fontsize=font_size)

    st.pyplot(fig)
