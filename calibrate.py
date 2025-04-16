import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import io
import zipfile
import streamlit as st

# Function to handle ZIP file upload and extract folder names
def handle_zip_upload(uploaded_file):
    temp_dir = '/tmp/extracted_zip/'
    os.makedirs(temp_dir, exist_ok=True)

    with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    folders = [f for f in os.listdir(temp_dir) if os.path.isdir(os.path.join(temp_dir, f))]
    if not folders:
        st.error("No folders found in the ZIP file.")
    return folders, temp_dir

# Read the bush.csv file from the same folder as calibrant.py
def read_bush_csv():
    calibrant_file_path = os.path.join(os.path.dirname(__file__), 'bush.csv')
    if os.path.exists(calibrant_file_path):
        bush_df = pd.read_csv(calibrant_file_path)
    else:
        st.error(f"'{calibrant_file_path}' not found. Make sure 'bush.csv' is in the same directory as the script.")
        bush_df = pd.DataFrame()  # Empty DataFrame if the file isn't found
    return bush_df

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
def process_folder_data(folder_name, base_path, bush_df, calibrant_type):
    folder_path = os.path.join(base_path, folder_name)
    results = []
    plots = []

    # Determine the column for the selected calibrant type
    calibrant_column = 'CCS_he' if calibrant_type == 'He' else 'CCS_n2'

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
                charge_state = filename.split('.')[0]
                # Look up the calibrant data from bush.csv based on protein and charge state
                calibrant_row = bush_df[(bush_df['protein'] == folder_name) & (bush_df['charge'] == int(charge_state)) & (bush_df['mass'] == float(mass))]
                calibrant_value = calibrant_row[calibrant_column].values[0] if not calibrant_row.empty else None
                results.append([folder_name, mass, charge_state, apex, r2, calibrant_value])
                plots.append((drift_time, intensity, fitted_values, filename, apex, r2))

    # Convert results to DataFrame
    results_df = pd.DataFrame(results, columns=['protein', 'mass', 'charge state', 'drift time', 'r2', 'calibrant_value'])

    return results_df, plots

# Function to display the data and plots
def display_results(results_df, plots):
    st.write("Combined Gaussian Fit Results:")
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

# Function to create .dat file from results
def generate_dat_file(results_df, velocity, voltage, pressure, length):
    dat_content = f"# length {length}\n# velocity {velocity}\n# voltage {voltage}\n# pressure {pressure}\n"

    # Create .dat content
    for _, row in results_df.iterrows():
        protein = row['protein']
        charge_state = row['charge state']
        mass = 
        mass = row['calibrant_value']*100 if not pd.isna(row['calibrant_value']) else 0  # Handle missing calibrant value
        drift_time = row['drift time']
        dat_content += f"{protein}_{charge_state} {mass} {charge_state} {drift_time} {drift_time}\n"
    
    return dat_content

# Main function for the Streamlit page
def calibrate_page():
    st.title("Calibrate with IMSCal")

    # Step 1: Upload ZIP file
    uploaded_zip_file = st.file_uploader("Upload a ZIP file", type="zip")
    if uploaded_zip_file is not None:
        # Extract the folders from the ZIP file
        folders, temp_dir = handle_zip_upload(uploaded_zip_file)

        # Step 2: Read bush.csv for calibrant data
        bush_df = read_bush_csv()

        # Step 3: Dropdown for selecting calibrant type (He or N2)
        calibrant_type = st.selectbox("Select Calibrant Type", options=["He", "N2"])

        # Step 4: Get user inputs for parameters
        velocity = st.number_input("Enter velocity", min_value=0.0, value=281.0)
        voltage = st.number_input("Enter voltage", min_value=0.0, value=20.0)
        pressure = st.number_input("Enter pressure", min_value=0.0, value=1.63)
        length = st.number_input("Enter length", min_value=0.0, value=0.980)

        # Step 5: Process all folders and files
        all_results_df = pd.DataFrame(columns=['protein', 'charge state', 'drift time', 'r2', 'calibrant_value'])
        all_plots = []

        # Process each folder
        for folder in folders:
            st.write(f"Processing folder: {folder}")

            # Process the data in the folder
            results_df, plots = process_folder_data(folder, temp_dir, bush_df, calibrant_type)

            # Append to the combined results DataFrame
            all_results_df = pd.concat([all_results_df, results_df], ignore_index=True)

            # Append plots for visualization
            all_plots.extend(plots)

        # Display the combined results and plots
        display_results(all_results_df, all_plots)

        # Option to download the combined results as CSV
        csv_buffer = io.StringIO()
        all_results_df.to_csv(csv_buffer, index=False)
        st.download_button(
            label="Download Combined Gaussian Fit Results (CSV)",
            data=csv_buffer.getvalue(),
            file_name="combined_gaussian_fit_results.csv",
            mime="text/csv"
        )

        # Generate .dat file and allow download
        dat_file_content = generate_dat_file(all_results_df, velocity, voltage, pressure, length)
        st.download_button(
            label="Download .dat File",
            data=dat_file_content,
            file_name="calibration_data.dat",
            mime="text/plain"
        )
