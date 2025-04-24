import streamlit as st
import pandas as pd
import numpy as np

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
