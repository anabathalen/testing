import streamlit as st
import zipfile
import tempfile
import os
import pandas as pd
from io import BytesIO, StringIO

def process_outputs_page():
    uploaded_zip = st.file_uploader("Upload a zipped folder", type="zip")

    if uploaded_zip:
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "uploaded.zip")

            with open(zip_path, "wb") as f:
                f.write(uploaded_zip.getvalue())

            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(tmpdir)

            protein_data = {}

            for root, dirs, files in os.walk(tmpdir):
                for file in files:
                    if file.startswith("output_") and file.endswith(".dat"):
                        file_path = os.path.join(root, file)

                        # Infer protein name from the first folder in the relative path
                        rel_path = os.path.relpath(file_path, tmpdir)
                        parts = rel_path.split(os.sep)
                        if len(parts) < 2:
                            continue
                        protein_name = parts[0]

                        with open(file_path, 'r') as f:
                            lines = f.readlines()

                        try:
                            cal_index = next(i for i, line in enumerate(lines) if line.strip() == "[CALIBRATED DATA]")
                            data_lines = lines[cal_index + 1:]
                        except StopIteration:
                            continue

                        try:
                            df = pd.read_csv(StringIO(''.join(data_lines)))
                            df = df[['Z', 'Drift', 'CCS', 'CCS Std.Dev.']]
                            if protein_name not in protein_data:
                                protein_data[protein_name] = []
                            protein_data[protein_name].append(df)
                        except Exception:
                            continue

            if protein_data:
                for protein_name, dfs in protein_data.items():
                    combined_df = pd.concat(dfs, ignore_index=True)
                    buffer = BytesIO()
                    combined_df.to_csv(buffer, index=False)
                    buffer.seek(0)
                    st.download_button(
                        label=f"Download {protein_name}.csv",
                        data=buffer,
                        file_name=f"{protein_name}.csv",
                        mime="text/csv"
                    )
            else:
                st.warning("No valid output_X.dat files found.")


