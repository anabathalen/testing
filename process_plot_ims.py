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
        if "Charge" not in cal_df.columns or "Intensity" not in cal_df.columns:
            st.error("The CSV file must contain 'Charge' and 'Intensity' columns.")
            return

        # Read mass spectrum data
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df = ms_df.dropna()

        # Constants
        PROTON_MASS = 1.007276

        # Get unique charge states from calibrated data
        unique_charges = sorted(cal_df["Charge"].unique())

        # Calculate m/z for each charge state and integrate spectrum
        charge_scale_factors = {}
        for z in unique_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            charge_scale_factors[z] = intensity_sum

        # Add scale factor columns
        cal_df["Scale Factor"] = cal_df["Charge"].map(charge_scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        # Show and download scaled data
        st.subheader("Scaled Calibrated Data")
        st.dataframe(cal_df)

        csv_output = cal_df.to_csv(index=False).encode("utf-8")
        st.download_button("Download Scaled CSV", data=csv_output, file_name="scaled_calibrated_data.csv", mime="text/csv")

        # === Plot configuration ===
        st.subheader("Plot Options")
        palette_choice = st.selectbox("Choose a color palette", list(sns.palettes.SEABORN_PALETTES.keys()), index=0)
        fig_width = st.slider("Figure width", min_value=4, max_value=20, value=10)
        fig_height = st.slider("Figure height", min_value=3, max_value=15, value=6)
        fig_dpi = st.slider("Figure resolution (DPI)", min_value=72, max_value=600, value=150)

        palette = sns.color_palette(palette_choice, n_colors=len(unique_charges))

        # === 1. Mass Spectrum with Integration Windows ===
        st.subheader("Mass Spectrum with Charge State Integration Regions")
        
        fig1, ax1 = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)
        ax1.plot(ms_df["m/z"], ms_df["Intensity"], color="gray", linewidth=1.5, label="Mass Spectrum")
        
        for i, z in enumerate(unique_charges):
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            region = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]
            
            # Make region visible even if empty
            if not region.empty:
                ax1.fill_between(region["m/z"], region["Intensity"], color=palette[i], alpha=0.7, label=f"Charge {z} window", linewidth=0)
                ax1.axvline(mz_min, color=palette[i], linestyle="--", linewidth=1)
                ax1.axvline(mz_max, color=palette[i], linestyle="--", linewidth=1)
        
        ax1.set_xlabel("m/z")
        ax1.set_title("Mass Spectrum with Integration Regions")
        ax1.legend()
        ax1.set_yticks([])
        ax1.set_ylabel("")  # remove y-axis label
        ax1.yaxis.set_ticklabels([])  # ensure tick labels are gone
        ax1.grid(False)
        
        for spine in ax1.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)
        
        st.pyplot(fig1)
        
        # === 2. CCS Plot by Charge State ===
        st.subheader("Scaled Intensity vs CCS by Charge State")
        # Filter data: Only include rows where CCS Std.Dev. < CCS
        filtered_df = cal_df[cal_df["CCS Std.Dev."] < cal_df["CCS"]]

        fig2, ax2 = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)

        for i, (charge, group) in enumerate(filtered_df.groupby("Charge")):
            ax2.plot(group["CCS"], group["Scaled Intensity"], label=f"Charge {charge}", color=palette[i])
        
        total_df = filtered_df.groupby("CCS")["Scaled Intensity"].sum().reset_index()
        ax2.plot(total_df["CCS"], total_df["Scaled Intensity"], color="black", linewidth=2.0, label="Total")
        
        ax2.set_xlabel("CCS (Å²)")
        ax2.set_title("Scaled Intensity vs CCS")
        ax2.set_yticks([])
        ax2.set_ylabel("")  # Remove y-axis label
        ax2.yaxis.set_ticklabels([])
        ax2.grid(False)
        
        for spine in ax2.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)
        
        st.pyplot(fig2)
        
        # === Download Plot as PNG ===
        buf = BytesIO()
        fig2.savefig(buf, format="png", dpi=fig_dpi, bbox_inches="tight")
        st.download_button(
            label="Download CCS Plot as PNG",
            data=buf.getvalue(),
            file_name="ccs_plot.png",
            mime="image/png"
        )
