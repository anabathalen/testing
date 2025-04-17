import streamlit as st
import zipfile
import tempfile
import os
import pandas as pd
from io import BytesIO

def calibrate_drift_files_page():
    # Upload zipped drift files (X.txt files for proteins)
    drift_zip = st.file_uploader("Upload zipped folder of raw drift files", type="zip")

    # Ask if the data is from Synapt or Cyclic
    data_type = st.radio("Is your data from a Synapt or Cyclic instrument?", ["Synapt", "Cyclic"])
    
    # If Cyclic, ask for the injection time to subtract from the drift time
    inject_time = None
    if data_type == "Cyclic":
        inject_time = st.number_input("Enter the injection time (ms)", min_value=0.0, value=0.0, step=0.1)
    
    # Upload the calibration CSV files from the previous step
    cal_csvs = st.file_uploader("Upload the CSV files from the 'Process Output Files' page", type="csv", accept_multiple_files=True)

    # Check if both drift zip and calibration CSVs are uploaded
    if drift_zip and cal_csvs:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Save the uploaded zipped drift files temporarily
            drift_zip_path = os.path.join(tmpdir, "drift.zip")
            with open(drift_zip_path, "wb") as f:
                f.write(drift_zip.getvalue())

            # Extract the drift files
            with zipfile.ZipFile(drift_zip_path, 'r') as zip_ref:
                zip_ref.extractall(tmpdir)

            # Load the calibration data into a dictionary
            calibration_lookup = {}
            for file in cal_csvs:
                protein_name = file.name.replace(".csv", "")
                df = pd.read_csv(file)
                for _, row in df.iterrows():
                    key = (protein_name, int(row["Z"]))
                    if key not in calibration_lookup:
                        calibration_lookup[key] = []
                    calibration_lookup[key].append({
                        "Drift": row["Drift"],
                        "CCS": row["CCS"],
                        "CCS Std.Dev.": row["CCS Std.Dev."]
                    })

            # Prepare to save the output dataframes
            output_buffers = {}

            # Process each drift file
            for root, _, files in os.walk(tmpdir):
                for file in files:
                    if file.endswith(".txt") and file.split(".")[0].isdigit():
                        charge_state = int(file.split(".")[0])
                        protein_name = os.path.basename(root)
                        key = (protein_name, charge_state)
                        cal_data = calibration_lookup.get(key)

                        # Skip if no calibration data is found
                        if not cal_data:
                            continue

                        # Read the raw drift data
                        file_path = os.path.join(root, file)
                        try:
                            raw_df = pd.read_csv(file_path, sep="\t", header=None, names=["Drift", "Intensity"])
                        except Exception as e:
                            st.error(f"Failed to read file {file}: {e}")
                            continue

                        # Adjust drift time for Cyclic data
                        if data_type == "Cyclic":
                            raw_df["Drift"] = raw_df["Drift"] - inject_time

                        # Match calibration drift times to intensities
                        out_rows = []
                        for entry in cal_data:
                            drift_val = entry["Drift"]
                            # Find the row in raw_df with the closest drift time
                            closest_idx = (raw_df["Drift"] - drift_val).abs().idxmin()
                            matched_intensity = raw_df.loc[closest_idx, "Intensity"]

                            out_rows.append({
                                "Charge": charge_state,
                                "Drift": drift_val,
                                "CCS": entry["CCS"],
                                "CCS Std.Dev.": entry["CCS Std.Dev."],
                                "Intensity": matched_intensity
                            })

                        # Create the dataframe for the matched data
                        out_df = pd.DataFrame(out_rows)

                        # Save each protein's data
                        out_key = f"{protein_name}.csv"
                        if out_key not in output_buffers:
                            output_buffers[out_key] = []
                        output_buffers[out_key].append(out_df)

            # If any output data was created, zip it and offer download
            if output_buffers:
                zip_buffer = BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_out:
                    for filename, dfs in output_buffers.items():
                        combined = pd.concat(dfs, ignore_index=True)
                        csv_bytes = combined.to_csv(index=False).encode("utf-8")
                        zip_out.writestr(filename, csv_bytes)

                zip_buffer.seek(0)
                st.download_button(
                    label="Download Calibrated Drift Data (ZIP)",
                    data=zip_buffer,
                    file_name="calibrated_drift_data.zip",
                    mime="application/zip"
                )
            else:
                st.warning("No matching calibration or intensity data found.")
