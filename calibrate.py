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
    results = []
    plots = []

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
                amp, apex, stddev = params
                file_number = filename.split('.')[0]
                results.append([file_number, apex, r2, amp, stddev])
                plots.append((drift_time, intensity, fitted_values, filename, apex, r2))

    # Convert results to DataFrame
    results_df = pd.DataFrame(results, columns=['File Number', 'Apex Drift Time', 'R²', 'Amplitude', 'Standard Deviation'])

    return results_df, plots

# Function to display the data and plots
def display_results(results_df, plots):
    st.write("Gaussian Fit Results:")
    st.dataframe(results_df)

    # Plot all the fits
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

# Main function for the Streamlit page
def calibrate_page():
    st.title("ZIP File Folder Extractor and Gaussian Fitting")

    uploaded_zip_file = st.file_uploader("Upload a ZIP file", type="zip")
    if uploaded_zip_file is not None:
        # Step 1: Extract the folders from the ZIP file
        folders, temp_dir = handle_zip_upload(uploaded_zip_file)

        # Step 2: Let the user choose a folder to process
        if folders:
            selected_folder = st.selectbox("Select Folder to Process", options=folders)
            if selected_folder:
                st.write(f"Processing folder: {selected_folder}")

                # Process the data in the selected folder
                results_df, plots = process_folder_data(selected_folder, temp_dir)

                # Display the results and plots
                display_results(results_df, plots)

                # Option to download the results as CSV
                csv_buffer = io.StringIO()
                results_df.to_csv(csv_buffer, index=False)
                st.download_button(
                    label="Download Gaussian Fit Results (CSV)",
                    data=csv_buffer.getvalue(),
                    file_name=f"{selected_folder}_gaussian_fit_results.csv",
                    mime="text/csv"
                )
