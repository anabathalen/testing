import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO

def plot_and_scale_page():
    st.title("Plot and Scale Calibrated Data")

    cal_file = st.file_uploader("Upload a calibrated CSV file for a protein", type="csv")
    ms_file = st.file_uploader("Upload the mass spectrum TXT file (no headers)", type="txt")
    protein_mass = st.number_input("Enter the protein mass (Da)", min_value=0.0, step=1.0)

    if cal_file and ms_file and protein_mass > 0:
        cal_df = pd.read_csv(cal_file)
        if not {'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'}.issubset(cal_df.columns):
            st.error("The CSV must include 'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'.")
            return

        cal_df = cal_df[cal_df['CCS Std.Dev.'] < 0.5 * cal_df['CCS']]
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df.dropna(inplace=True)

        PROTON_MASS = 1.007276
        unique_charges = sorted(cal_df["Charge"].unique())

        scale_factors = {}
        for z in unique_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[z] = intensity_sum

        cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        st.subheader("Scaled Calibrated Data")
        st.dataframe(cal_df)

        csv_output = cal_df.to_csv(index=False).encode("utf-8")
        st.download_button("Download Scaled CSV", data=csv_output, file_name="scaled_calibrated_data.csv", mime="text/csv", key="scaled_csv")

        st.subheader("Plot Options")
        palette_choice = st.selectbox("Choose a color palette", list(sns.palettes.SEABORN_PALETTES.keys()))
        fig_width = st.slider("Figure width", min_value=4, max_value=20, value=10)
        fig_height = st.slider("Figure height", min_value=3, max_value=15, value=6)
        fig_dpi = st.slider("Figure DPI", min_value=72, max_value=300, value=100)
        font_size = st.slider("Font size", min_value=6, max_value=24, value=12)

        palette = sns.color_palette(palette_choice, n_colors=len(unique_charges))

        # === 1. MASS SPECTRUM WITH INTEGRATION WINDOWS ===
        st.subheader("Mass Spectrum with Charge State Integration Regions")
        fig1, ax1 = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)
        ax1.plot(ms_df["m/z"], ms_df["Intensity"], color="gray", label="Mass Spectrum")

        for i, z in enumerate(unique_charges):
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            ax1.axvspan(mz_min, mz_max, color=palette[i], alpha=0.2, label=f"Charge {z} window")
            region = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]
            ax1.fill_between(region["m/z"], region["Intensity"], color=palette[i], alpha=0.5)

        ax1.set_xlabel("m/z", fontsize=font_size)
        ax1.set_ylabel("Intensity", fontsize=font_size)
        ax1.set_title("Mass Spectrum with Integration Windows", fontsize=font_size + 2)
        ax1.tick_params(axis='both', labelsize=font_size)
        ax1.legend(fontsize=font_size)

        for spine in ax1.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)

        st.pyplot(fig1)

        # === 2. INTERPOLATED CCS PLOT ===
        st.subheader("Scaled Intensity vs CCS")
        ccs_min = np.floor(cal_df["CCS"].min())
        ccs_max = np.ceil(cal_df["CCS"].max())
        ccs_grid = np.arange(ccs_min, ccs_max + 1, 1.0)

        interpolated_traces = []
        fig2, ax2 = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)

        for i, (charge, group) in enumerate(cal_df.groupby("Charge")):
            group_sorted = group.sort_values("CCS")
            interp_intensity = np.interp(ccs_grid, group_sorted["CCS"], group_sorted["Scaled Intensity"], left=0, right=0)
            interpolated_traces.append(interp_intensity)
            ax2.plot(group_sorted["CCS"], group_sorted["Scaled Intensity"], label=f"Charge {charge}", color=palette[i])

        total_trace = np.sum(interpolated_traces, axis=0)
        ax2.plot(ccs_grid, total_trace, color="black", linewidth=2.0, label="Total (Interpolated)")

        ax2.set_xlabel("CCS (Å²)", fontsize=font_size)
        ax2.set_title("Scaled Intensity vs CCS", fontsize=font_size + 2)
        ax2.tick_params(axis='both', labelsize=font_size)
        ax2.set_yticks([])
        ax2.set_ylabel("")
        ax2.yaxis.set_ticklabels([])
        ax2.legend(fontsize=font_size)
        ax2.grid(False)

        for spine in ax2.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)

        st.pyplot(fig2)

        # PNG download with unique key
        fig_buffer = BytesIO()
        fig2.savefig(fig_buffer, format='png', dpi=fig_dpi, bbox_inches='tight')
        fig_buffer.seek(0)
        st.download_button("Download CCS Plot as PNG", data=fig_buffer, file_name="ccs_plot.png", mime="image/png", key="ccs_plot_png")
