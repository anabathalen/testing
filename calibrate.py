import streamlit as st
import zipfile
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import io

# Function to handle ZIP file upload and extract folder names
def handle_zip_upload(uploaded_file):
    # Create a temporary directory to extract ZIP contents
    temp_dir = '/tmp/extracted_zip/'
    os.makedirs(temp_dir, exist_ok=True)

    # Extract the ZIP file contents
    with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    # List folders in the extracted ZIP content
    folders = [f for f in os.listdir(temp_dir) if os.path.isdir(os.path.join(temp_dir, f))]
    
    if not folders:
        st.error("No folders found in the ZIP file.")
    return folders, temp_dir

# Gaussian fit function
def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

# R² Calculation
def r_squared(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1 - (ss_res / ss_tot)

# Fit Gaussian and retry with different initial guesses
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

# Function to process each folder and extract data from .txt files
def process_folder_data(folder_name, base_path):
    folder_path = os.path.join(base_path, folder_name)
    folder_data = []

    # Iterate through each file in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt') and filename[0].isdigit():  # Only process .txt files
            file_path = os.path.join(folder_path, filename)

            data = np.loadtxt(file_path)
            drift_time = data[:, 0]
            intensity = data[:, 1]

            # Perform Gaussian fit
            params, r2, fitted_values = fit_gaussian_with_retries(drift_time, intensity)
            if params is not None:
                file_number = int(filename.split('.')[0])  # Extract charge state number from filename
                amp, apex, stddev = params
                # Store results as a row to append to the final dataframe
                folder_data.append({
                    'protein': folder_name,
                    'charge state': file_number,
                    'drift time': apex,  # Storing Apex Drift Time
                    'R²': r2
                })

    return folder_data

# Function to display the final results and plots
def display_results(all_folders_data):
    st.write("Final Gaussian Fit Results for All Folders:")

    # Combine all the results into a single DataFrame for easy viewing
    all_results = []
    for folder, data in all_folders_data.items():
        all_results.extend(data)

    results_df = pd.DataFrame(all_results)
    st.dataframe(results_df)

    # Plot all the fits for each charge state in each folder
    n_plots = len(all_results)
    n_cols = 3
    n_rows = (n_plots + n_cols - 1) // n_cols

    plt.figure(figsize=(12, 4 * n_rows))
    for i, (folder, data) in enumerate(all_folders_data.items()):
        for charge_state, params in data.items():
            # For illustration purposes, we just use arbitrary values for plotting
            st.write(f"Plotting {folder} - Charge State {charge_state}")
            # Actual plotting code can go here if needed

    plt.tight_layout()
    st.pyplot(plt)

# Main function for the Streamlit page
def calibrate_page():
    st.title("ZIP File Folder Extractor and Gaussian Fitting")

    uploaded_zip_file = st.file_uploader("Upload a ZIP file", type="zip")
    if uploaded_zip_file is not None:
        # Step 1: Extract the folders from the ZIP file
        folders, temp_dir = handle_zip_upload(uploaded_zip_file)

        # Step 2: Automatically process all folders
        all_folders_data = {}
        for folder in folders:
            st.write(f"Processing folder: {folder}")
            folder_data = process_folder_data(folder, temp_dir)
            all_folders_data[folder] = folder_data

        # Step 3: Display the results and allow download
        display_results(all_folders_data)

        # Option to download all results as CSV
        csv_buffer = io.StringIO()
        all_results = []
        for folder, data in all_folders_data.items():
            for charge_state, params in data.items():
                row = {'Folder': folder, 'Charge State': charge_state, **params}
                all_results.append(row)

        results_df = pd.DataFrame(all_results)
        results_df.to_csv(csv_buffer, index=False)
        st.download_button(
            label="Download All Gaussian Fit Results (CSV)",
            data=csv_buffer.getvalue(),
            file_name="all_gaussian_fit_results.csv",
            mime="text/csv"
        )
