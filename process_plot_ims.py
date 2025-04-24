import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

        # Plot config controls
        st.subheader("Plot Options")
        palette_choice = st.selectbox("Choose a color palette", sns.palettes.SEABORN_PALETTES.keys(), index=0)
        fig_width = st.slider("Figure width", min_value=4, max_value=20, value=10)
        fig_height = st.slider("Figure height", min_value=3, max_value=15, value=6)
        
        palette = sns.color_palette(palette_choice, n_colors=len(unique_charges))
        
        # ==== 1. MASS SPECTRUM WITH INTEGRATION WINDOWS ====
        st.subheader("Mass Spectrum with Charge State Integration Regions")
        
        fig1, ax1 = plt.subplots(figsize=(fig_width, fig_height))
        
        # Plot raw spectrum
        ax1.plot(ms_df["m/z"], ms_df["Intensity"], color="gray", label="Mass Spectrum")
        
        # Highlight integration regions for each charge state
        for i, z in enumerate(unique_charges):
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            color = palette[i]
        
            region = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]
            ax1.fill_between(region["m/z"], region["Intensity"], color=color, alpha=0.5, label=f"Charge {z} window")
        
        ax1.set_xlabel("m/z")
        ax1.set_ylabel("Intensity")
        ax1.set_title("Mass Spectrum with Integration Windows")
        ax1.legend()
        ax1.grid(True)
        
        # Black border
        for spine in ax1.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)
        
st.pyplot(fig1)

        
        # Add color palette selection
        st.subheader("Plot Options")
        palette_choice = st.selectbox("Choose a color palette", sns.palettes.SEABORN_PALETTES.keys(), index=0)
        fig_width = st.slider("Figure width", min_value=4, max_value=20, value=10)
        fig_height = st.slider("Figure height", min_value=3, max_value=15, value=6)
        
        # Plotting
        palette = sns.color_palette(palette_choice, n_colors=cal_df["Charge"].nunique())
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        # Plot each charge state
        for i, (charge, group) in enumerate(cal_df.groupby("Charge")):
            ax.plot(group["Drift"], group["Scaled Intensity"], label=f"Charge {charge}", color=palette[i])
        
        # Plot total intensity across all charges
        total_df = cal_df.groupby("Drift")["Scaled Intensity"].sum().reset_index()
        ax.plot(total_df["Drift"], total_df["Scaled Intensity"], color="black", linewidth=2.0, label="Total")
        
        # Styling
        ax.set_xlabel("Drift Time (s)")
        ax.set_ylabel("Scaled Intensity")
        ax.set_title("Scaled Intensity vs Drift Time by Charge State")
        ax.legend()
        ax.grid(True)
        
        # Add black border around plot area
        for spine in ax.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)
        
        st.pyplot(fig)

