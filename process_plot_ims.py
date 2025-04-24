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
        if not {"Charge", "Intensity", "CCS", "CCS Std.Dev."}.issubset(cal_df.columns):
            st.error("CSV must contain 'Charge', 'Intensity', 'CCS', 'CCS Std.Dev.' columns.")
            return

        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"]).dropna()

        PROTON_MASS = 1.007276
        unique_charges = sorted(cal_df["Charge"].unique())

        charge_scale_factors = {}
        for z in unique_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            charge_scale_factors[z] = intensity_sum

        cal_df["Scale Factor"] = cal_df["Charge"].map(charge_scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        st.subheader("Scaled Calibrated Data")
        st.dataframe(cal_df)

        csv_output = cal_df.to_csv(index=False).encode("utf-8")
        st.download_button("Download Scaled CSV", data=csv_output, file_name="scaled_calibrated_data.csv", mime="text/csv")

        st.subheader("CCS Plot Options")
        palette_choice = st.selectbox("Choose a color palette", sns.palettes.SEABORN_PALETTES.keys(), index=0)
        fig_width = st.slider("Figure width (CCS plot)", min_value=4, max_value=20, value=10)
        fig_height = st.slider("Figure height (CCS plot)", min_value=3, max_value=15, value=6)
        fig_dpi = st.slider("Figure DPI (CCS plot)", min_value=50, max_value=300, value=150)
        font_size = st.slider("Font size for CCS plot", min_value=8, max_value=20, value=12)

        palette = sns.color_palette(palette_choice, n_colors=len(unique_charges))

        st.subheader("Mass Spectrum with Charge State Integration Regions")
        fig1, ax1 = plt.subplots(figsize=(fig_width, fig_height))
        ax1.plot(ms_df["m/z"], ms_df["Intensity"], color="gray", label="Mass Spectrum")
        for i, z in enumerate(unique_charges):
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            ax1.vlines(x=[mz_min, mz_max], ymin=0, ymax=max(ms_df["Intensity"]), color=palette[i], alpha=0.5, label=f"Charge {z} region")
        ax1.set_xlabel("m/z", fontsize=font_size)
        ax1.set_ylabel("Intensity", fontsize=font_size)
        ax1.set_title("Mass Spectrum with Integration Regions", fontsize=font_size)
        ax1.legend(fontsize=font_size)
        ax1.grid(False)
        for spine in ax1.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)
        st.pyplot(fig1)

        # === CCS Plot ===
        st.subheader("Scaled Intensity vs CCS by Charge State")

        # Filter: keep only good CCS Std.Dev. (< 50%) and non-zero scaled intensities
        filtered_cal_df = cal_df[
            (cal_df["CCS Std.Dev."] < cal_df["CCS"] * 0.5) &
            (cal_df["Scaled Intensity"] > 0)
        ]

        fig2, ax2 = plt.subplots(figsize=(fig_width, fig_height))

        min_ccs = int(filtered_cal_df["CCS"].min())
        max_ccs = int(filtered_cal_df["CCS"].max())
        common_ccs = np.arange(min_ccs, max_ccs + 1, 1)  # 1 Å² resolution

        total_intensity = np.zeros_like(common_ccs, dtype=float)

        for i, (charge, group) in enumerate(filtered_cal_df.groupby("Charge")):
            group = group.sort_values("CCS")
            if len(group) >= 2:
                interp_intensity = np.interp(common_ccs, group["CCS"], group["Scaled Intensity"])
                ax2.plot(common_ccs, interp_intensity, label=f"Charge {charge}", color=palette[i])
                ax2.fill_between(common_ccs, 0, interp_intensity, color=palette[i], alpha=0.3)
                total_intensity += interp_intensity

        ax2.plot(common_ccs, total_intensity, color="black", linewidth=2.0, label="Total")

        ax2.set_xlabel("CCS (Å²)", fontsize=font_size)
        ax2.set_ylabel("", fontsize=font_size)  # No y-axis label
        ax2.set_title("Scaled Intensity vs CCS by Charge State", fontsize=font_size)
        ax2.legend(fontsize=font_size)
        ax2.grid(False)

        # Remove y-axis ticks and labels
        ax2.set_yticks([])
        ax2.set_yticklabels([])

        for spine in ax2.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)

        fig_buffer = BytesIO()
        fig2.savefig(fig_buffer, format='png', dpi=fig_dpi, bbox_inches="tight")
        fig_buffer.seek(0)
        st.download_button("Download CCS Plot as PNG", data=fig_buffer, file_name="ccs_plot.png", mime="image/png")

        st.pyplot(fig2)
