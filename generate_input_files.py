import os
import io
import zipfile
import numpy as np
import pandas as pd
import streamlit as st
from pathlib import Path
from tempfile import TemporaryDirectory

def handle_zip_upload(uploaded_file):
    temp_dir = '/tmp/samples_extracted/'
    
    # 🚨 Clear old files
    if os.path.exists(temp_dir):
        for root, dirs, files in os.walk(temp_dir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))

    os.makedirs(temp_dir, exist_ok=True)

    with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    folders = [f for f in os.listdir(temp_dir) if os.path.isdir(os.path.join(temp_dir, f))]
    if not folders:
        st.error("No folders found in the ZIP file.")
    return folders, temp_dir


def process_sample_folder(folder_name, folder_path, mass, drift_mode, inject_time=None):
    sample_output = io.BytesIO()
    dat_files = {}

    output_folder = os.path.join(folder_path, folder_name)
    os.makedirs(output_folder, exist_ok=True)

    for filename in os.listdir(os.path.join(folder_path, folder_name)):
        if filename.endswith('.txt') and filename[0].isdigit():
            file_path = os.path.join(folder_path, folder_name, filename)
            data = np.loadtxt(file_path)

            drift_time = data[:, 0]
            intensity = data[:, 1]

            if drift_mode == "Cyclic" and inject_time is not None:
                drift_time = drift_time - inject_time

            drift_time = np.maximum(drift_time, 0)  # avoid negative times

            index = np.arange(len(drift_time))
            df = pd.DataFrame({
                "index": index,
                "mass": mass,
                "intensity": intensity,
                "drift_time": drift_time
            })

            dat_filename = f"input_{os.path.splitext(filename)[0]}.dat"
            dat_path = os.path.join(output_folder, dat_filename)
            df.to_csv(dat_path, sep=' ', index=False, header=False)
    
    return output_folder

def generate_output_zip(sample_folders_paths):
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w") as zipf:
        for sample_folder in sample_folders_paths:
            for root, _, files in os.walk(sample_folder):
                for file in files:
                    full_path = os.path.join(root, file)
                    relative_path = os.path.relpath(full_path, os.path.dirname(sample_folder))
                    zipf.write(full_path, arcname=relative_path)
    return zip_buffer

# Main Streamlit app
def generate_input_dat_files_app():
    st.title("Generate Sample Input .dat Files")

    uploaded_zip_file = st.file_uploader("Upload ZIP containing sample folders", type="zip")

    if uploaded_zip_file is not None:
        sample_folders, base_path = handle_zip_upload(uploaded_zip_file)

        drift_mode = st.radio("What type of data is this?", options=["Cyclic", "Synaptic"])
        inject_time = None
        if drift_mode == "Cyclic":
            inject_time = st.number_input("Enter inject time to subtract (ms)", min_value=0.0, value=0.0)

        all_sample_paths = []

        st.subheader("Input Mass for Each Sample")

        with st.form("sample_mass_form"):
            sample_mass_map = {}
            for sample in sample_folders:
                mass = st.number_input(f"Mass for sample '{sample}'", min_value=0.0, key=sample)
                sample_mass_map[sample] = mass
            submitted = st.form_submit_button("Generate .dat Files")

        if submitted:
            with TemporaryDirectory() as tmp_output_dir:
                for sample in sample_folders:
                    sample_path = process_sample_folder(
                        folder_name=sample,
                        folder_path=base_path,
                        mass=sample_mass_map[sample],
                        drift_mode=drift_mode,
                        inject_time=inject_time
                    )
                    all_sample_paths.append(sample_path)

                zip_buffer = generate_output_zip(all_sample_paths)

                st.success("All .dat files generated and zipped!")

                st.download_button(
                    label="Download All .dat Files (ZIP)",
                    data=zip_buffer.getvalue(),
                    file_name="sample_dat_files.zip",
                    mime="application/zip"
                )
