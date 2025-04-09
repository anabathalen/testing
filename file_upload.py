# file_upload.py
import streamlit as st
import pandas as pd

def handle_file_upload():
    uploaded_file = st.file_uploader("Upload a CSV file", type="csv")
    
    if uploaded_file is not None:
        st.session_state["uploaded_file"] = uploaded_file
        st.write("File uploaded successfully!")

    if st.session_state.get("uploaded_file"):
        try:
            df = pd.read_csv(st.session_state["uploaded_file"])
            st.write("Your data:")
            st.dataframe(df)
        except Exception as e:
            st.error(f"Error reading the file: {e}")
    else:
        st.write("No file uploaded yet.")
