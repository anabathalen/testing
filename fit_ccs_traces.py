import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from scipy.stats import norm
from scipy.optimize import curve_fit
import io

# Constants
PROTON_MASS = 1.007276

# Function to load the mass spectrum data
def load_mass_spectrum():
    uploaded_file = st.file_uploader("Upload your mass spectrum file", type=["txt"])
    if uploaded_file is not None:
        # Read the mass spectrum data into a pandas DataFrame
        ms_df = pd.read_csv(uploaded_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df.dropna(inplace=True)
        return ms_df
    return None

# Function to calculate scale factors for the selected charge states
def calculate_scale_factors(ms_df, cal_df, selected_charges, protein_mass):
    scale_factors = {}
    for z in selected_charges:
        mz = (protein_mass + z * PROTON_MASS) / z  # Calculate m/z for each charge state
        mz_min = mz * 0.98  # Define the range for finding peaks
        mz_max = mz * 1.02
        # Sum the intensities in the range of the m/z
        intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
        scale_factors[z] = intensity_sum  # Store the scale factor for each charge state

    return scale_factors

# Function to scale the intensities for each charge state
def scale_intensities(cal_df, scale_factors):
    cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)  # Map scale factors to the charge states
    cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]  # Scale intensities by the scale factor
    return cal_df

# Gaussian fitting function
def gaussian(x, mean, stddev, amplitude):
    return amplitude * np.exp(-0.5 * ((x - mean) / stddev) ** 2)

# Function to fit Gaussian to a charge state
def fit_gaussian_to_trace(mz_values, intensities):
    popt, _ = curve_fit(gaussian, mz_values, intensities, p0=[mz_values[np.argmax(intensities)], 0.1, np.max(intensities)])
    return popt  # Return the optimized parameters (mean, stddev, amplitude)

# Function to save data to a CSV file for download
def to_csv(df):
    csv = df.to_csv(index=False)
    return io.StringIO(csv)

# Function to fit and sum the CCS traces and perform Gaussian fitting
def fit_ccs_traces_page():
    st.title("Fit CCS Traces")

    # Load mass spectrum file
    ms_df = load_mass_spectrum()

    if ms_df is not None:
        # Ask for protein mass and selected charge states
        protein_mass = st.number_input("Enter protein mass (Da)", min_value=0.0, step=0.1)
        all_charges = sorted(ms_df['Charge'].unique())  # Assuming charges are already in the mass spectrum
        selected_charges = st.multiselect("Select charge states to include", all_charges, default=all_charges)

        # Perform calculations
        scale_factors = calculate_scale_factors(ms_df, ms_df, selected_charges, protein_mass)
        scaled_df = scale_intensities(ms_df, scale_factors)

        # Fit Gaussians for each charge state
        gaussian_params = {}
        summed_trace = np.zeros_like(ms_df['m/z'].values)

        for z in selected_charges:
            charge_df = scaled_df[scaled_df['Charge'] == z]
            mz_values = charge_df['m/z'].values
            intensities = charge_df['Scaled Intensity'].values

            # Fit Gaussian to the charge state
            popt = fit_gaussian_to_trace(mz_values, intensities)
            gaussian_params[z] = popt

            # Create a fitted Gaussian trace
            fitted_trace = gaussian(mz_values, *popt)
            summed_trace += fitted_trace

            # Plot individual Gaussian fit for the charge state
            plt.plot(mz_values, intensities, label=f"Charge {z} Data", alpha=0.7)
            plt.plot(mz_values, fitted_trace, label=f"Charge {z} Fit", linestyle="--")

        # Plot summed Gaussian fit
        plt.plot(ms_df['m/z'], summed_trace, label="Summed Fit", color='black', linestyle='-')
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title('CCS Trace Fitting')
        plt.legend()
        st.pyplot()

        # Prepare the Gaussian fit parameters for download
        gaussian_params_df = pd.DataFrame.from_dict(gaussian_params, orient='index', columns=["Mean", "Std.Dev.", "Amplitude"])

        # Download button for Gaussian fits
        st.download_button("Download Gaussian Fit Parameters", to_csv(gaussian_params_df), file_name="gaussian_fits.csv")

        # Prepare the summed trace for download
        summed_df = pd.DataFrame({
            "m/z": ms_df['m/z'],
            "Summed Intensity": summed_trace
        })

        # Download button for summed CCS vs fitted intensity
        st.download_button("Download Summed CCS vs Fitted Intensity", to_csv(summed_df), file_name="summed_ccs_fitted.csv")

