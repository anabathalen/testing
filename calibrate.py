import streamlit as st
import zipfile
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import io

# Function to handle .zip file upload and extract folder names
def handle_zip_upload(unique_key):
    uploaded_file = st.file_uploader("Upload a ZIP file", type="zip", key=unique_key)
    
    if uploaded_file is not None:
        try:
            with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
                zip_contents = zip_ref.namelist()
                folders = set()
                for item in zip_contents:
                    if item.endswith('/'):  # Folders end with '/'
                        folders.add(item)
                
                if folders:
                    st.write("Folders in the ZIP file:")
                    for folder in sorted(folders):
                        st.write(folder)
                else:
                    st.write("No folders found in the ZIP file.")
        
        except zipfile.BadZipFile:
            st.error("The file uploaded is not a valid ZIP file.")
        except Exception as e:
            st.error(f"An error occurred: {e}")
        return folders
    return []

# Define the Gaussian function
def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

# R² Calculation
def r_squared(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1 - (ss_res / ss_tot)

# Function for fitting Gaussian and retrying with different initial guesses
def fit_gaussian_with_retries(drift_time, intensity, n_attempts=10):
    best_r2 = -np.inf
    best_params = None
    best_fitted_values = None

    for _ in range(n_attempts):
        initial_guess = [
            np.random.uniform(0.8 * max(intensity), 1.2 * max(intensity)),
            np.random.uniform(np.min(drift_time), np.max(drift_time)),
            np.random.uniform(0.1 * np.std(drift_time), 2 * np.std(drift_time))
        ]

        try:
            params, _ = curve_fit(gaussian, drift_time, intensity, p0=initial_guess)
            fitted_values = gaussian(drift_time, *params)
            r2 = r_squared(intensity, fitted_values)

            if r2 > best_r2:
                best_r2 = r2
                best_params = params
                best_fitted_values = fitted_values
        except RuntimeError:
            continue

    return best_params, best_r2, best_fitted_values

# Process the protein files and return results as DataFrame
def process_protein_files(protein_name, zip_file):
    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
        folder_path = f"/tmp/{protein_name}/"
        zip_ref.extractall(folder_path)
    
    results = []
    plots = []

    protein_folder = os.path.join(folder_path, protein_name)
    
    for filename in os.listdir(protein_folder):
        if filename.endswith('.txt') and re.match(r'^\d', filename):
            file_path = os.path.join(protein_folder, filename)

            data = np.loadtxt(file_path)
            drift_time = data[:, 0]
            intensity = data[:, 1]

            params, r2, fitted_values = fit_gaussian_with_retries(drift_time, intensity)
            if params is not None:
                amp, apex, stddev = params
                file_number = filename.split('.')[0]
                results.append([file_number, apex, r2, amp, stddev])

                plots.append((drift_time, intensity, fitted_values, filename, apex, r2))

    results_df = pd.DataFrame(results, columns=['File Number', 'Apex Drift Time', 'R²', 'Amplitude', 'Standard Deviation'])

    st.write(f"Gaussian Fitting Results for {protein_name}:")
    st.dataframe(results_df)

    n_plots = len(plots)
    n_cols = 3
    n_rows = (n_plots + n_cols - 1) // n_cols

    plt.figure(figsize=(12, 4 * n_rows))
    for i, (drift_time, intensity, fitted_values, filename, apex, r2) in enumerate(plots):
        plt.subplot(n_rows, n_cols, i + 1)
        plt.plot(drift_time, intensity, 'b.', label='Raw Data', markersize=3)
        plt.plot(drift_time, fitted_values, 'r-', label='Gaussian Fit', linewidth=1)
        plt.title(f'{filename}\nApex: {apex:.2f}, R²: {r2:.3f}')
        plt.xlabel('Drift Time')
        plt.ylabel('Intensity')
        plt.legend()
        plt.grid()

    plt.tight_layout()
    st.pyplot(plt)

    return results_df

# Main app code to handle ZIP upload and processing
def calibrate_page():
    st.title("ZIP File Folder Extractor and Gaussian Fitting")

    uploaded_zip_file = st.file_uploader("Upload a ZIP file", type="zip", key="zip_file_uploader_1")
    if uploaded_zip_file is not None:
        folders = handle_zip_upload(unique_key="zip_file_uploader_1")  # Extract folders from the ZIP file
        if folders:
            protein_name = st.selectbox("Select Protein Folder", options=list(folders))
            if protein_name:
                st.write(f"Processing {protein_name}...")
                results_df = process_protein_files(protein_name, uploaded_zip_file)

                # Provide an option to download the results as CSV
                csv_buffer = io.StringIO()
                results_df.to_csv(csv_buffer, index=False)
                st.download_button(
                    label="Download Gaussian Fit Results (CSV)",
                    data=csv_buffer.getvalue(),
                    file_name=f"{protein_name}_gaussian_fit_results.csv",
                    mime="text/csv"
                )


