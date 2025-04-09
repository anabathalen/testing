import os
import zipfile
import numpy as np
import pandas as pd
import streamlit as st
from pathlib import Path
from tempfile import TemporaryDirectory

def handle_zip_upload(uploaded_file):
    temp_dir = '/tmp/samples_extracted/'
    os.makedirs(temp_dir, exist_ok=True)

    with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    folders = [f for f in os.listdir(temp_dir) if os.path.isdir(os.path.join(temp_dir, f))]
    if not folders:
        st.error("No folders found in the ZIP file.")
    return folders, temp_dir

def process_protein_folder(protein_folder, base_path):
    output_files = {}
    protein_output_folder = os.path.join(base_path, protein_folder)

    # Skip if folder doesn't have any relevant files
    if not any(f.endswith('.dat') for f in os.listdir(protein_output_folder)):
        st.warning(f"No valid output .dat files found in folder: {protein_folder}")
        return None

    # Read all the relevant .dat files in the folder
    for filename in os.listdir(protein_output_folder):
        if filename.endswith('.dat') and 'output' in filename:
            file_path = os.path.join(protein_output_folder, filename)
            try:
                # Try reading file, ignoring first row (PARAMETERS, DIAGNOSTICS sections)
                data = pd.read_csv(file_path, delimiter=',', skiprows=3, header=None)

                # Check if the number of columns matches expectations
                if data.shape[1] != 6:  # Expecting 6 columns (ID, Mass, Z, Drift, CCS, CCS Std.Dev.)
                    st.error(f"Error: unexpected number of columns in file: {filename}")
                    continue

                # Store the dataframe for further processing
                output_files[filename] = data
            except Exception as e:
                st.error(f"Error processing output file {filename}: {e}")

    return output_files

def merge_data(drift_file, output_data):
    # Read drift time and intensity from input file
    drift_data = pd.read_csv(drift_file, sep=' ', header=None, names=["index", "mass", "intensity", "drift_time"])
    
    # Merge the drift data with the output data on the drift time column (using index as a join key)
    merged_data = pd.merge(drift_data, output_data, how="left", on="index")

    # Extract relevant columns to return
    merged_data = merged_data[["protein", "charge_state", "drift_time", "intensity"]]
    return merged_data

def generate_output_zip(final_dataframes):
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w") as zipf:
        for sample_name, final_df in final_dataframes.items():
            # Save each dataframe as a .csv file in the zip
            sample_file_path = f"{sample_name}_merged.csv"
            final_df.to_csv(sample_file_path, index=False)
            zipf.writestr(sample_file_path, final_df.to_csv(index=False))
    
    return zip_buffer

# Main Streamlit app
def analyze_output_dat_files_app():
    st.title("Analyze Output .dat Files")

    uploaded_zip_file = st.file_uploader("Upload ZIP containing protein output folders", type="zip")
    
    if uploaded_zip_file is not None:
        protein_folders, base_path = handle_zip_upload(uploaded_zip_file)

        all_final_dataframes = {}

        for protein_folder in protein_folders:
            st.subheader(f"Processing folder: {protein_folder}")
            
            output_files = process_protein_folder(protein_folder, base_path)
            if output_files is not None:
                for output_filename, output_data in output_files.items():
                    # Here we process each file, assuming the drift files are available
                    drift_file = os.path.join(base_path, protein_folder, "input_charge.dat")

                    # Merge the data
                    final_data = merge_data(drift_file, output_data)
                    all_final_dataframes[protein_folder] = final_data

        # If we have processed data, allow for downloading
        if all_final_dataframes:
            zip_buffer = generate_output_zip(all_final_dataframes)

            st.success("Data merged successfully!")

            st.download_button(
                label="Download Merged Data (ZIP)",
                data=zip_buffer.getvalue(),
                file_name="merged_data.zip",
                mime="application/zip"
            )

# Run the Streamlit app
if __name__ == "__main__":
    analyze_output_dat_files_app()

