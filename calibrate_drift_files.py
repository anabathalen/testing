import streamlit as st
import zipfile
import tempfile
import os
import pandas as pd
from io import BytesIO

def calibrate_drift_files_page():
    drift_zip = st.file_uploader("Upload zipped folder of raw drift files", type="zip")
    data_type = st.radio("Is your data from a Synapt or Cyclic instrument?", ["Synapt", "Cyclic"])
    inject_time = None
    if data_type == "Cyclic":
        inject_time = st.number_input("Enter the injection time (ms)", min_value=0.0, value=0.0, step=0.1)
    
    cal_csvs = st.file_uploader("Upload the CSV files from the 'Process Output Files' page", type="csv", accept_multiple_files=True)

    if drift_zip and cal_csvs:
        with tempfile.TemporaryDirectory() as tmpdir:
            drift_zip_path = os.path.join(tmpdir, "drift.zip")
            with open(drift_zip_path, "wb") as f:
                f.write(drift_zip.getvalue())

            with zipfile.ZipFile(drift_zip_path, 'r') as zip_ref:
                zip_ref.extractall(tmpdir)

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

            output_buffers = {}

            for root, _, files in os.walk(tmpdir):
                for file in files:
                    if file.endswith(".txt") and file.split(".")[0].isdigit():
                        charge_state = int(file.split(".")[0])
                        protein_name = os.path.basename(root)
                        key = (protein_name, charge_state)
                        cal_data = calibration_lookup.get(key)

                        if not cal_data:
                            continue

                        file_path = os.path.join(root, file)
                        try:
                            raw_df = pd.read_csv(file_path, sep="\t", header=None, names=["Drift", "Intensity"])
                        except Exception as e:
                            st.error(f"Failed to read file {file}: {e}")
                            continue

                        if data_type == "Cyclic":
                            raw_df["Drift"] = raw_df["Drift"] - inject_time

                        out_rows = []
                        for entry in cal_data:
                            drift_val = entry["Drift"]
                            closest_idx = (raw_df["Drift"] - drift_val).abs().idxmin()
                            matched_intensity = raw_df.loc[closest_idx, "Intensity"]

                            out_rows.append({
                                "Charge": charge_state,
                                "Drift": drift_val,
                                "CCS": entry["CCS"],
                                "CCS Std.Dev.": entry["CCS Std.Dev."],
                                "Intensity": matched_intensity
                            })

                        out_df = pd.DataFrame(out_rows)
                        out_key = f"{protein_name}.csv"
                        if out_key not in output_buffers:
                            output_buffers[out_key] = []
                        output_buffers[out_key].append(out_df)

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
