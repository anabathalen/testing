import os
import io
import zipfile
import pandas as pd
import streamlit as st

def extract_zip(uploaded_file):
    extract_dir = "/tmp/extracted_outputs"
    if os.path.exists(extract_dir):
        for root, dirs, files in os.walk(extract_dir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
    os.makedirs(extract_dir, exist_ok=True)

    with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
        zip_ref.extractall(extract_dir)

    folders = [f for f in os.listdir(extract_dir) if os.path.isdir(os.path.join(extract_dir, f))]
    return extract_dir, folders

def parse_output_dat(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    try:
        cal_data_index = lines.index('[CALIBRATED DATA]\n') + 1
    except ValueError:
        return pd.DataFrame()

    # Read header
    headers = lines[cal_data_index].strip().split(',')
    data_lines = lines[cal_data_index + 1:]
    data = [line.strip().split(',') for line in data_lines if line.strip()]
    df = pd.DataFrame(data, columns=headers)
    df = df.apply(pd.to_numeric, errors='ignore')
    return df

def read_input_dat(filepath):
    return pd.read_csv(filepath, sep=' ', header=None, names=['index', 'mass', 'intensity', 'drift_time'])

def merge_data(input_df, output_df):
    return pd.merge(input_df, output_df[['ID', 'CCS', 'CCS Std.Dev.']], left_on='index', right_on='ID')

def process_protein_folder(folder_path, folder_name):
    results = []

    dat_files = [f for f in os.listdir(folder_path) if f.startswith("input_") and f.endswith(".dat")]
    if not dat_files:
        return pd.DataFrame()  # Skip folders with no input files

    for file in dat_files:
        charge_state = file.split("_")[1].split(".")[0]
        input_path = os.path.join(folder_path, file)
        output_path = os.path.join(folder_path, f"output_{charge_state}.dat")

        if not os.path.exists(output_path):
            continue  # Skip if there's no matching output file

        input_df = read_input_dat(input_path)
        output_df = parse_output_dat(output_path)

        if output_df.empty:
            continue  # Skip malformed or empty output files

        merged = pd.merge(input_df, output_df, left_on='index', right_on='ID', how='inner')
        if merged.empty:
            continue  # No matching rows, skip this pair

        merged['protein'] = folder_name
        merged['charge state'] = charge_state
        results.append(merged[['protein', 'charge state', 'drift_time', 'intensity', 'CCS', 'CCS Std.Dev.']])

    if results:
        return pd.concat(results, ignore_index=True)
    return pd.DataFrame()  # Return empty if nothing useful was processed


    if results:
        return pd.concat(results, ignore_index=True)
    return pd.DataFrame()

def analyze_output_dat_files_app():
    st.title("Combine Input/Output .dat Data into Summary CSV")

    uploaded_zip = st.file_uploader("Upload ZIP file with input/output .dat files", type="zip")
    if uploaded_zip:
        base_path, protein_folders = extract_zip(uploaded_zip)
        all_results = []

        for folder in protein_folders:
            st.write(f"Processing {folder}...")
            folder_path = os.path.join(base_path, folder)
            result = process_protein_folder(folder_path, folder)
            if not result.empty:
                all_results.append(result)

        if all_results:
            final_df = pd.concat(all_results, ignore_index=True)
            st.dataframe(final_df)

            csv = final_df.to_csv(index=False)
            st.download_button(
                label="Download Combined CSV",
                data=csv,
                file_name="combined_protein_data.csv",
                mime="text/csv"
            )
        else:
            st.warning("No valid data extracted from uploaded files.")

# To run this function in your app
if __name__ == "__main__":
    analyze_output_dat_files_app()
