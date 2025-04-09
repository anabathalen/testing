import streamlit as st
import zipfile
import os

# Function to handle .zip file upload and extract folder names
def handle_zip_upload():
    # Allow the user to upload a .zip file
    uploaded_file = st.file_uploader("Upload a ZIP file", type="zip")
    
    if uploaded_file is not None:
        try:
            # Read the uploaded ZIP file
            with zipfile.ZipFile(uploaded_file, 'r') as zip_ref:
                # List all the file/folder names in the ZIP file
                zip_contents = zip_ref.namelist()
                
                # Filter out only directories (folders) from the list
                folders = set()
                for item in zip_contents:
                    if item.endswith('/'):  # Folders end with '/'
                        folders.add(item)
                
                # Display the folder names
                if folders:
                    st.write("Folders in the ZIP file:")
                    for folder in sorted(folders):
                        st.write(folder)
                else:
                    st.write("No folders found in the ZIP file.")
        
        except zipfile.BadZipFile:
            st.error("The file uploaded is not a valid ZIP file.")
        except Exception as e:
            st.error(f"An error occurred: {e}")

# Main app code
def calibrate_page():
    st.title("ZIP File Folder Extractor")

    st.write("Upload a ZIP file, and this tool will list all the folder names inside it.")

    handle_zip_upload()  # Call the function to handle ZIP upload and extract folder names

