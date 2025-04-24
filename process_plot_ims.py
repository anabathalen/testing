import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO

def plot_and_scale_page():
    st.title("Plot and Scale Calibrated Data")

    # Upload calibrated CSV for a protein
    cal_file = st.file_uploader("Upload a calibrated CSV file for a protein", type="csv")

    # Upload mass spectrum TXT file
    ms_file = st.file_uploader("Upload the mass spectrum TXT file (no headers)", type="txt")

    # Input protein mass
    protein_mass = st.number_input("Enter the protein mass (Da)", min_value=0.0, step=1.0)

    if cal_file and ms_file and protein_mass > 0:
        # Read calibrated drift data
        cal_df = pd.read_csv(cal_file)
        if not {'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'}.issubset(cal_df.columns):
            st.error("The CSV file must contain 'Charge', 'CCS', 'CCS Std.Dev.', and 'Intensity' columns.")
            return

        # Filter where standard deviation < CCS value
        filtered_df = cal_df[cal_df['CCS Std.Dev.'] < cal_df['CCS']].copy()

        # Read mass spectrum data
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df.dropna(inplace=True)

        # Constants
        PROTON_MASS = 1.007276

        # Get unique charge states
        unique_charges = sorted(filtered_df["Charge"].unique())

        # Compute scale factors for each charge state
        scale_factors = {}
        for z in unique_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[z] = intensity_sum

        # Add scale factor and scaled intensity to dataframe
        filtered_df["Scale Factor"] = filtered_df["Charge"].map(scale_factors)
        filtered_df["Scaled Intensity"] = filtered_df["Intensity"] * filtered_df["Scale Factor"]

        st.subheader("Scaled Calibrated Data")
        st.dataframe(filtered_df)

        csv_output = filtered_df.to_csv(index=False).encode("utf-8")
        st.download_button("Download Scaled CSV", data=csv_output, file_name="scaled_calibrated_data.csv", mime="text/csv")

        # Plot Options
        st.subheader("Plot Options")
        palette_choice = st.selectbox("Choose a color palette", list(sns.palettes.SEABORN_PALETTES.keys()))
        fig_width = st.slider("Figure width", min_value=4, max_value=20, value=10)
        fig_height = st.slider("Figure height", min_value=3, max_value=15, value=6)
        fig_dpi = st.slider("Figure DPI", min_value=72, max_value=300, value=100)

        palette = sns.color_palette(palette_choice, n_colors=len(unique_charges))

        # ==== 1. MASS SPECTRUM WITH INTEGRATION WINDOWS ====
        st.subheader("Mass Spectrum with Charge State Integration Regions")

        fig1, ax1 = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)
        ax1.plot(ms_df["m/z"], ms_df["Intensity"], color="gray", label="Mass Spectrum")

        for i, z in enumerate(unique_charges):
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            region = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]
            ax1.fill_between(region["m/z"], region["Intensity"], color=palette[i], alpha=0.5, label=f"Charge {z} window")

        ax1.set_xlabel("m/z")
        ax1.set_ylabel("Intensity")
        ax1.set_title("Mass Spectrum with Integration Windows")
        ax1.legend()

        for spine in ax1.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)

        st.pyplot(fig1)

        # ==== 2. INTERPOLATED CCS PLOT ====
        st.subheader("Scaled Intensity vs CCS")

        # Interpolate to common CCS grid
        ccs_min = np.floor(filtered_df["CCS"].min())
        ccs_max = np.ceil(filtered_df["CCS"].max())
        ccs_grid = np.arange(ccs_min, ccs_max + 1, 1.0)

        interpolated_traces = []

        fig2, ax2 = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)

        for i, (charge, group) in enumerate(filtered_df.groupby("Charge")):
            group_sorted = group.sort_values("CCS")
            interp_intensity = np.interp(ccs_grid, group_sorted["CCS"], group_sorted["Scaled Intensity"], left=0, right=0)
            interpolated_traces.append(interp_intensity)
            ax2.plot(group_sorted["CCS"], group_sorted["Scaled Intensity"], label=f"Charge {charge}", color=palette[i])

        # Plot summed total trace
        total_trace = np.sum(interpolated_traces, axis=0)
        ax2.plot(ccs_grid, total_trace, color="black", linewidth=2.0, label="Total (Interpolated)")

        ax2.set_xlabel("CCS (Å²)")
        ax2.set_title("Scaled Intensity vs CCS")
        ax2.set_yticks([])
        ax2.set_ylabel("")
        ax2.yaxis.set_ticklabels([])
        ax2.grid(False)

        for spine in ax2.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)

        st.pyplot(fig2)

        # Download button for figure
        fig_buffer = BytesIO()
        fig2.savefig(fig_buffer, format='png', dpi=fig_dpi, bbox_inches='tight')
        fig_buffer.seek(0)
        st.download_button("Download CCS Plot as PNG", data=fig_buffer, file_name="ccs_plot.png", mime="image/png")
