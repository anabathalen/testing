# process_outputs.py
import streamlit as st
import zipfile
import tempfile
import os
import pandas as pd
from io import BytesIO, StringIO

def process_outputs_page():
    st.title("Process Outputs")
    st.subheader("Upload zipped output folders for each protein sample (with output_X.dat files inside).")

    uploaded_zip = st.file_uploader("Upload a zipped folder", type="zip")

    if uploaded_zip:
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "uploaded.zip")

            with open(zip_path, "wb") as f:
                f.write(uploaded_zip.getvalue())

            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(tmpdir)

            st.success("Zip file extracted!")

            protein_data = {}

            # Traverse protein subfolders
            for root, dirs, files in os.walk(tmpdir):
                # Skip root folder — we want only subfolders
                if root == tmpdir:
                    continue

                protein_name = os.path.basename(root)
                protein_df_list = []

                for file in files:
                    if file.startswith("output_") and file.endswith(".dat"):
                        file_path = os.path.join(root, file)
                        with open(file_path, 'r') as f:
                            lines = f.readlines()

                        try:
                            cal_index = next(i for i, line in enumerate(lines) if line.strip() == "[CALIBRATED DATA]")
                            data_lines = lines[cal_index + 1:]
                        except StopIteration:
                            st.warning(f"No [CALIBRATED DATA] section in {file}")
                            continue

                        try:
                            df = pd.read_csv(StringIO(''.join(data_lines)))
                            df = df[['Z', 'Drift', 'CCS', 'CCS Std.Dev.']]
                            protein_df_list.append(df)
                        except Exception as e:
                            st.error(f"Failed to parse data from {file_path}: {e}")
                            continue

                if protein_df_list:
                    combined_df = pd.concat(protein_df_list, ignore_index=True)
                    protein_data[protein_name] = combined_df

            if protein_data:
                st.header("Download CSVs per Protein")
                for protein_name, df in protein_data.items():
                    buffer = BytesIO()
                    df.to_csv(buffer, index=False)
                    buffer.seek(0)
                    st.download_button(
                        label=f"Download {protein_name}.csv",
                        data=buffer,
                        file_name=f"{protein_name}.csv",
                        mime="text/csv"
                    )
            else:
                st.warning("No valid output_X.dat files with calibrated data found.")
