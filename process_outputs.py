import os
import io
import zipfile
import numpy as np
import pandas as pd
import streamlit as st
from pathlib import Path

# Function to handle ZIP file upload and extract folder names
def handle_zip_upload(uploaded_file):
    temp_dir = '/tmp/protein_data/'
    os.makedirs(temp_dir, exist_ok=True)

    with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    # Filter out system folders like __MACOSX and only keep real folders
    protein_folders = [
        f for f in os.listdir(temp_dir)
        if os.path.isdir(os.path.join(temp_dir, f)) and not f.startswith("__")
    ]
    if not protein_folders:
        st.error("No valid protein folders found in the ZIP file.")
    return protein_folders, temp_dir

# Function to process each protein folder and extract data from input and output files
def process_protein_folder(protein_folder, base_path):
    protein_path = os.path.join(base_path, protein_folder)
    
    input_files = {}
    output_files = {}
    
    # Search for input and output .dat files
    for filename in os.listdir(protein_path):
        if filename.startswith('input_') and filename.endswith('.dat'):
            input_files[filename] = os.path.join(protein_path, filename)
        elif filename.startswith('output_') and filename.endswith('.dat'):
            output_files[filename] = os.path.join(protein_path, filename)
    
    # Skip if no input or output files are found for this protein
    if not input_files or not output_files:
        return pd.DataFrame()

    # Process the data from both input and output files
    protein_data = []

    for input_filename, input_file_path in input_files.items():
        output_filename = f"output_{input_filename[6:]}"  # Match corresponding output file
        if output_filename not in output_files:
            continue  # Skip if no matching output file for this input
        
        # Load input data
        input_data = np.loadtxt(input_file_path)
        drift_time = input_data[:, 3]  # Drift time column (4th column)
        intensity = input_data[:, 2]  # Intensity column (3rd column)

        # Load output data, skipping non-data rows and using the correct columns
        try:
            # Skip rows until the data starts (skip headers, empty rows, or comments)
            output_data = np.loadtxt(output_files[output_filename], dtype=str, skiprows=4, usecols=(3, 4, 5))
            
            # Check if the output file has data
            if output_data.ndim == 1:  # Single column case
                output_data = np.reshape(output_data, (-1, 3))  # Reshape to avoid dimension errors

            # Sort the output data by index (the first column)
            output_data = output_data[output_data[:, 0].argsort()]

            # Convert the CCS and CCS stddev columns to numeric
            ccs = pd.to_numeric(output_data[:, 1], errors='coerce')
            ccs_stddev = pd.to_numeric(output_data[:, 2], errors='coerce')
        except Exception as e:
            st.warning(f"Error processing output file {output_filename}: {str(e)}")
            continue

        # Ensure the lengths match
        if len(drift_time) != len(ccs):
            st.warning(f"Warning: Drift time and CCS lengths mismatch in {protein_folder}/{input_filename}")
            continue
        
        # Combine input and output data
        for i in range(len(drift_time)):
            protein_data.append([protein_folder, int(input_filename[6:-4]), drift_time[i], intensity[i], ccs[i], ccs_stddev[i]])

    # Convert to DataFrame
    protein_df = pd.DataFrame(protein_data, columns=['protein', 'charge_state', 'drift_time', 'intensity', 'ccs', 'ccs_stddev'])
    
    return protein_df

# Function to combine results into a single CSV
def combine_results(protein_folders_data):
    combined_df = pd.concat(protein_folders_data, ignore_index=True)
    return combined_df

# Main Streamlit app function
def analyze_output_dat_files_app():
    st.title("Analyze Protein Data - Combine Input and Output")

    # Step 1: Upload ZIP file
    uploaded_zip_file = st.file_uploader("Upload ZIP containing protein folders", type="zip")

    if uploaded_zip_file is not None:
        protein_folders, base_path = handle_zip_upload(uploaded_zip_file)

        all_protein_data = []
        for protein_folder in protein_folders:
            st.write(f"Processing protein folder: {protein_folder}")
            protein_data = process_protein_folder(protein_folder, base_path)
            
            # Skip empty dataframes (i.e., no valid data)
            if not protein_data.empty:
                all_protein_data.append(protein_data)

        # Combine all results into a single dataframe
        if all_protein_data:
            combined_df = combine_results(all_protein_data)
            st.write("Combined Data:")
            st.dataframe(combined_df)

            # Allow user to download the combined data as a CSV file
            csv_buffer = io.StringIO()
            combined_df.to_csv(csv_buffer, index=False)
            st.download_button(
                label="Download Combined CSV",
                data=csv_buffer.getvalue(),
                file_name="combined_protein_data.csv",
                mime="text/csv"
            )
        else:
            st.error("No valid data extracted from uploaded files.")

# Run the Streamlit app
if __name__ == "__main__":
    analyze_output_dat_files_app()
