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

            st.write("📁 File tree inside ZIP:")
            for root, dirs, files in os.walk(tmpdir):
                for file in files:
                    full_path = os.path.join(root, file)
                    rel_path = os.path.relpath(full_path, tmpdir)
                    st.text(f"  - {rel_path}")


            protein_data = {}

            st.write("🔍 Scanning extracted files...")
            for root, dirs, files in os.walk(tmpdir):
                for file in files:
                    if file.startswith("output_") and file.endswith(".dat"):
                        file_path = os.path.join(root, file)
                        st.write(f"Found: `{file_path}`")

                        # Infer protein name as first folder after root
                        rel_path = os.path.relpath(file_path, tmpdir)
                        parts = rel_path.split(os.sep)
                        if len(parts) < 2:
                            st.warning(f"⚠️ Could not determine protein name from: {rel_path}")
                            continue
                        protein_name = parts[0]

                        with open(file_path, 'r') as f:
                            lines = f.readlines()

                        try:
                            cal_index = next(i for i, line in enumerate(lines) if line.strip() == "[CALIBRATED DATA]")
                            data_lines = lines[cal_index + 1:]
                        except StopIteration:
                            st.warning(f"No [CALIBRATED DATA] in {file}")
                            continue

                        try:
                            df = pd.read_csv(StringIO(''.join(data_lines)))
                            df = df[['Z', 'Drift', 'CCS', 'CCS Std.Dev.']]
                            if protein_name not in protein_data:
                                protein_data[protein_name] = []
                            protein_data[protein_name].append(df)
                        except Exception as e:
                            st.error(f"Failed to parse {file_path}: {e}")
                            continue

            if protein_data:
                st.header("Download CSVs per Protein")
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

