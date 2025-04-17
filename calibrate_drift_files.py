import streamlit as st
import zipfile
import tempfile
import os
import pandas as pd
from io import BytesIO, StringIO

def calibrate_drift_files_page():
    drift_zip = st.file_uploader("Upload zipped folder of raw drift files", type="zip")

    data_type = st.radio("Is your data from a Synapt or Cyclic instrument?", ["Synapt", "Cyclic"])
    inject_time = None
    if data_type == "Cyclic":
        inject_time = st.number_input("Enter the injection time (ms)", min_value=0.0, value=0.0, step=0.1)

    cal_csvs = st.file_uploader("Upload the CSV files from the 'Process Output Files' page", type="csv", accept_multiple_files=True)

    if drift_zip and cal_csvs:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Extract drift files
            drift_zip_path = os.path.join(tmpdir, "drift.zip")
            with open(drift_zip_path, "wb") as f:
                f.write(drift_zip.getvalue())

            with zipfile.ZipFile(drift_zip_path, 'r') as zip_ref:
                zip_ref.extractall(tmpdir)

            # Load calibration data into dictionary
            calibration_lookup = {}
            for file in cal_csvs:
                name = file.name.replace(".csv", "")
                df = pd.read_csv(file)
                for _, row in df.iterrows():
                    key = (name, int(row["Z"]))
                    calibration_lookup[key] = {
                        "CCS": row["CCS"],
                        "CCS Std.Dev.": row["CCS Std.Dev."]
                    }

            # Process drift files
            output_buffers = {}
            for root, _, files in os.walk(tmpdir):
                for file in files:
                    if file.endswith(".txt") and file.split(".")[0].isdigit():
                        charge_state = int(file.split(".")[0])
                        protein_name = os.path.basename(root)
                        key = (protein_name, charge_state)
                        cal = calibration_lookup.get(key)

                        if not cal:
                            continue  # Skip if no calibration for this protein/charge

                        file_path = os.path.join(root, file)
                        df = pd.read_csv(file_path, sep="\t", header=None, names=["Drift", "Intensity"])
                        if data_type == "Cyclic" and inject_time is not None:
                            df["Drift"] = df["Drift"] - inject_time

                        df["Charge"] = charge_state
                        df["CCS"] = cal["CCS"]
                        df["CCS Std.Dev."] = cal["CCS Std.Dev."]

                        df = df[["Charge", "Drift", "CCS", "CCS Std.Dev.", "Intensity"]]

                        # Save per-protein
                        output_key = f"{protein_name}.csv"
                        if output_key not in output_buffers:
                            output_buffers[output_key] = []
                        output_buffers[output_key].append(df)

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
                st.warning("No matching calibration data found for any files.")
