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
        if "Charge" not in cal_df.columns or "Intensity" not in cal_df.columns or "CCS" not in cal_df.columns or "CCS Std. Dev." not in cal_df.columns:
            st.error("The CSV file must contain 'Charge', 'Intensity', 'CCS', and 'CCS Std. Dev.' columns.")
            return

        # Read mass spectrum data
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df = ms_df.dropna()

        # Constants
        PROTON_MASS = 1.007276

        # Get unique charge states from calibrated data
        unique_charges = sorted(cal_df["Charge"].unique())

        # Calculate m/z for each charge state
        charge_scale_factors = {}
        for z in unique_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            charge_scale_factors[z] = intensity_sum

        # Add Scale Factor and Scaled Intensity to the DataFrame
        cal_df["Scale Factor"] = cal_df["Charge"].map(charge_scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        # Show and allow download
        st.subheader("Scaled Calibrated Data")
        st.dataframe(cal_df)

        csv_output = cal_df.to_csv(index=False).encode("utf-8")
        st.download_button("Download Scaled CSV", data=csv_output, file_name="scaled_calibrated_data.csv", mime="text/csv")

        # Plot config controls for CCS plot only
        st.subheader("CCS Plot Options")
        palette_choice = st.selectbox("Choose a color palette", sns.palettes.SEABORN_PALETTES.keys(), index=0)
        fig_width = st.slider("Figure width (CCS plot)", min_value=4, max_value=20, value=10)
        fig_height = st.slider("Figure height (CCS plot)", min_value=3, max_value=15, value=6)
        fig_dpi = st.slider("Figure DPI (CCS plot)", min_value=50, max_value=300, value=150)
        font_size = st.slider("Font size for CCS plot", min_value=8, max_value=20, value=12)

        palette = sns.color_palette(palette_choice, n_colors=len(unique_charges))

        # ==== 1. MASS SPECTRUM WITH INTEGRATION BARS ====
        st.subheader("Mass Spectrum with Charge State Integration Regions")
        
        fig1, ax1 = plt.subplots(figsize=(fig_width, fig_height))
        
        # Plot raw spectrum
        ax1.plot(ms_df["m/z"], ms_df["Intensity"], color="gray", label="Mass Spectrum")
        
        # Plot vertical bars for integration regions for each charge state
        for i, z in enumerate(unique_charges):
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            color = palette[i]
        
            ax1.vlines(x=[mz_min, mz_max], ymin=0, ymax=max(ms_df["Intensity"]), color=color, alpha=0.5, label=f"Charge {z} region")

        ax1.set_xlabel("m/z")
        ax1.set_ylabel("Intensity")
        ax1.set_title("Mass Spectrum with Integration Regions")
        ax1.legend()
        ax1.grid(False)

        # Black border
        for spine in ax1.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)
        
        st.pyplot(fig1)

        # ==== 2. CCS PLOT WITH SHADING UNDER EACH CHARGE STATE ====
        st.subheader("Scaled Intensity vs CCS by Charge State")
        
        # Filter out rows where CCS Std. Dev. >= CCS
        filtered_cal_df = cal_df[cal_df["CCS Std. Dev."] < cal_df["CCS"]]

        # Plotting the CCS plot
        fig2, ax2 = plt.subplots(figsize=(fig_width, fig_height))
        
        # Plot each charge state with semi-transparent shading
        for i, (charge, group) in enumerate(filtered_cal_df.groupby("Charge")):
            # Interpolating each charge state to 1 A2 resolution
            interpolated_ccs = np.arange(group["CCS"].min(), group["CCS"].max(), 1)
            interpolated_intensity = np.interp(interpolated_ccs, group["CCS"], group["Scaled Intensity"])
            ax2.plot(interpolated_ccs, interpolated_intensity, label=f"Charge {charge}", color=palette[i])
            ax2.fill_between(interpolated_ccs, 0, interpolated_intensity, color=palette[i], alpha=0.3)

        # Plot total intensity across all charges
        total_df = filtered_cal_df.groupby("CCS")["Scaled Intensity"].sum().reset_index()
        ax2.plot(total_df["CCS"], total_df["Scaled Intensity"], color="black", linewidth=2.0, label="Total")

        # Set axes labels and title
        ax2.set_xlabel("CCS (s)", fontsize=font_size)
        ax2.set_ylabel("", fontsize=font_size)  # No label for y-axis
        ax2.set_title("Scaled Intensity vs CCS by Charge State", fontsize=font_size)
        ax2.legend(fontsize=font_size)
        ax2.grid(False)

        # Add black border around plot area
        for spine in ax2.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)

        # Allow downloading the figure
        fig_buffer = BytesIO()
        fig2.savefig(fig_buffer, format='png', dpi=fig_dpi, bbox_inches="tight")
        fig_buffer.seek(0)
        
        st.download_button("Download CCS Plot as PNG", data=fig_buffer, file_name="ccs_plot.png", mime="image/png")

        st.pyplot(fig2)

