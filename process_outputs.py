# process_outputs.py
import streamlit as st
import zipfile
import tempfile
import os
import pandas as pd
from io import BytesIO, StringIO

def process_outputs_page():
    
    uploaded_zip = st.file_uploader("Upload a zipped folder containing output files", type="zip")

    if uploaded_zip:
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "uploaded.zip")

            with open(zip_path, "wb") as f:
                f.write(uploaded_zip.getvalue())

            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(tmpdir)

            st.success("Zip file extracted!")

            csv_buffers = {}

            for root, dirs, files in os.walk(tmpdir):
                for file in files:
                    if file.startswith("output_") and file.endswith(".dat"):
                        protein_folder = os.path.basename(root)
                        file_path = os.path.join(root, file)

                        with open(file_path, 'r') as f:
                            lines = f.readlines()

                        try:
                            cal_index = next(i for i, line in enumerate(lines) if line.strip() == "[CALIBRATED DATA]")
                            data_lines = lines[cal_index + 1:]
                        except StopIteration:
                            st.warning(f"No [CALIBRATED DATA] section found in {file}")
                            continue

                        try:
                            df = pd.read_csv(StringIO(''.join(data_lines)))
                            df = df[['Z', 'Drift', 'CCS', 'CCS Std.Dev.']]
                        except Exception as e:
                            st.error(f"Failed to parse calibrated data from {file}: {e}")
                            continue

                        csv_buffer = BytesIO()
                        df.to_csv(csv_buffer, index=False)
                        csv_buffer.seek(0)

                        label = f"{protein_folder}_{file.replace('.dat', '.csv')}"
                        csv_buffers[label] = csv_buffer

            if csv_buffers:
                st.header("Download Processed CSVs")
                for label, buffer in csv_buffers.items():
                    st.download_button(
                        label=f"Download {label}",
                        data=buffer,
                        file_name=label,
                        mime='text/csv'
                    )
            else:
                st.warning("No output_X.dat files with valid calibrated data found.")

