# app.py
import streamlit as st
from gaussian_fitting import upload_and_plot
from calibrate import calibrate_page
from generate_input_files import generate_input_dat_files_app
# Sidebar for navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Fit Gaussians to Data", "Calibrate", "Generate Input Files", "Process Outputs"])

# Check if the page has changed, to ensure no redundant loading
if "page" not in st.session_state or st.session_state["page"] != page:
    st.session_state.clear()  # Clear session state for fresh start

# Store the current page in session state
st.session_state["page"] = page

# Home Page
if page == "Home":
    st.title("👩🏻‍🔬 Barran Group IM-MS Processing Tools")
    st.subheader("← Navigate to the tool you need from the sidebar.")
    st.markdown('''This is a work in progress, so there are a few things to note:  
    1. I haven't validated this - sanity check your results, don't trust the outputs for important things.  
    2. If the site is down it is because I'm working on it, either find a previous version on my github or just ask me and I'll sort it!  
    3. I encourage you to make improvements to this - please make a new branch and do your thing, then open a pull request and I will likely accept.
    Other than that, I hope you enjoy. Let me know if you need any help (ana.bathalen@manchester.ac.uk or slack)!''')

# Fit Gaussians to Data Page
elif page == "Fit Gaussians to Data":
    st.title("Fit Gaussians to Data")
    upload_and_plot()

# Calibrate Page
elif page == "Calibrate":
    calibrate_page()

# Generate Input Files Page
elif page == "Generate Input Files":
    st.title("Get Input Files")
    generate_input_dat_files_app()
