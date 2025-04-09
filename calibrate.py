import os
import pandas as pd
import zipfile
import streamlit as st
import numpy as np

# Function to read bush.csv and get protein data
def read_bush_csv():
    calibrant_file_path = os.path.join(os.path.dirname(__file__), 'bush.csv')
    if os.path.exists(calibrant_file_path):
        bush_df = pd.read_csv(calibrant_file_path)
    else:
        st.error(f"'{calibrant_file_path}' not found. Make sure 'bush.csv' is in the same directory as the script.")
        bush_df = pd.DataFrame()  # Empty DataFrame if the file isn't found
    return bush_df

# Function to handle ZIP file upload and extract folders
def handle_zip_upload(uploaded_file):
    temp_dir = '/tmp/extracted_zip/'
    os.makedirs(temp_dir, exist_ok=True)
    
    with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)
        
    # List all folders in the extracted ZIP
    folders = [f for f in os.listdir(temp_dir) if os.path.isdir(os.path.join(temp_dir, f))]
    return folders, temp_dir

# Function to fit the data (for now we are just simulating)
def fit_data_for_folder(folder_path):
    # Placeholder fitting function, replace with your actual fitting process
    # For now, return simulated drift time and charge state
    return np.random.rand(10), np.random.randint(10, 30, size=10)  # Simulated drift time and charge states

# Function to generate the .dat file
def generate_dat_file(bush_df, velocity, voltage, pressure, length, selected_data):
    # Choose the correct column based on user input for calibrant type
    calibrant_type = st.selectbox("Select Calibrant Type", options=["He", "N2"])
    calibrant_column = 'CCS_he' if calibrant_type == 'He' else 'CCS_n2'

    # Open the file for writing
    dat_file_path = os.path.join(os.path.dirname(__file__), 'calibration_data.dat')
    with open(dat_file_path, 'w') as f:
        # Write metadata information
        f.write(f"# length {length}\n")
        f.write(f"# velocity {velocity}\n")
        f.write(f"# voltage {voltage}\n")
        f.write(f"# pressure {pressure}\n")

        # Iterate through the selected data and write to the .dat file
        for protein, charge_state, mass, reference_value, drift_time in selected_data:
            f.write(f"{protein} {mass} {charge_state} {reference_value} {drift_time}\n")

    st.success(f"Calibration data .dat file generated successfully: {dat_file_path}")

# Streamlit UI for calibration page
def calibrate_page():
    st.title("Calibration Data Generator")

    # Step 1: Upload ZIP file with data
    uploaded_zip_file = st.file_uploader("Upload a ZIP file containing data", type="zip")
    
    if uploaded_zip_file:
        folders, temp_dir = handle_zip_upload(uploaded_zip_file)
        st.write(f"Found the following folders in the ZIP file: {folders}")
        
        # Step 2: Fit the data for each folder
        selected_data = []
        for folder in folders:
            folder_path = os.path.join(temp_dir, folder)
            st.write(f"Processing folder: {folder}")
            
            # Simulate fitting the data for this folder
            drift_times, charge_states = fit_data_for_folder(folder_path)
            for charge_state, drift_time in zip(charge_states, drift_times):
                protein_name = f"{folder}_{charge_state}"  # protein name for the charge state
                # Use mass and CCS from bush.csv (e.g., simulation with 'myoglobin')
                protein_row = bush_df[bush_df['protein'] == 'myoglobin']
                if not protein_row.empty:
                    mass = protein_row['mass'].values[0]
                    reference_value = protein_row['CCS_he'].values[0] * 100  # Use He or N2 based on selection
                    selected_data.append((protein_name, charge_state, mass, reference_value, drift_time))
        
        # Step 3: User inputs for calibration parameters
        velocity = st.number_input("Enter velocity (e.g., 281)", min_value=1.0)
        voltage = st.number_input("Enter voltage (e.g., 20)", min_value=1.0)
        pressure = st.number_input("Enter pressure (e.g., 1.63)", min_value=0.1)
        length = st.number_input("Enter length (e.g., 0.980)", min_value=0.1)

        # Step 4: Generate the .dat file upon button press
        if st.button("Generate .dat File"):
            if selected_data:
                generate_dat_file(bush_df, velocity, voltage, pressure, length, selected_data)
            else:
                st.warning("No data to generate the .dat file. Please ensure the fitting step is complete.")

# Run the Streamlit app
if __name__ == "__main__":
    calibrate_page()


