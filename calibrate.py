# calibrate.py
import streamlit as st
from file_upload import handle_file_upload

def calibrate_page():
    st.title("Calibrate")
    handle_file_upload()
