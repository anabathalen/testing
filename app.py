# app.py
import streamlit as st
from gaussian_fitting import upload_and_plot
from calibrate import calibrate_page
from generate_input_files import generate_input_dat_files_app
from process_outputs import process_outputs_page
from calibrate_drift_files import calibrate_drift_files_page


# Sidebar for navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Fit Gaussians to Data", "Calibrate", "Generate Input Files", "Process Outputs", "Calibrate Drift Files"])

# Check if the page has changed, to ensure no redundant loading
if "page" not in st.session_state or st.session_state["page"] != page:
    st.session_state.clear()  # Clear session state for fresh start

# Store the current page in session state
st.session_state["page"] = page

# Home Page
if page == "Home":
    st.title("👩🏻‍🔬 Barran Group IM-MS Processing Tools")
    st.subheader("← Navigate to the tool you need from the sidebar.")
    st.markdown('''**This** is a work in progress, so there are a few things to note:  
    1. I haven't validated this - sanity check your results, don't trust the outputs for important things.  
    2. If the site is down it is because I'm working on it, either find a previous version on my github or just ask me and I'll sort it!  
    3. I encourage you to make improvements to this - please make a new branch and do your thing, then open a pull request and I will likely accept.  
    Other than that, I hope you enjoy. Let me know if you need any help (ana.bathalen@manchester.ac.uk or slack)!''')

# Fit Gaussians to Data Page
elif page == "Fit Gaussians to Data":
    st.title("Fit Gaussians to Data")
    st.subheader("Takes x, y data and helps you fit gaussians and plot a pretty graph.")
    st.markdown('''This is a super simple tool - upload your ATD/CCSD as a csv file with headings 'x' and 'y' for drift time / CCS and intensity respectively. When you do this, you'll see your data visualised with starting guesses at the peak maxima. Input the number of gaussians you see and their approximate positions. The script uses the data +/- 5% of the inputted peak position to fit gaussians to each peak. This will work well if you have relatively well-resolved peaks - if you don't, then this is not the tool for you. The purpose of this method is to avoid overfitting small features.''')
    upload_and_plot()

# Calibrate Page
elif page == "Calibrate":
    st.title("Calibrate (with IMSCal)")
    st.subheader("Takes ATDs of Calibrants and Returns IMSCal Reference File AND a CSV for the Old Method.")
    st.markdown('''For now this tool doesn't allow you to bypass the tedious MassLynx stuff - in future it will use TWIMExtract. To use this tool, make a folder for each calibrant titled the name of the calibrant (see Bush Database for calibrant names). In each folder, paste the ATD for each charge state in turn into a text file called x.txt where x is the charge state. Zip these folders and upload it here.''')
    calibrate_page()

# Generate Input Files Page
elif page == "Generate Input Files":
    st.title("Get Input Files")
    generate_input_dat_files_app()

# Process Outputs Page
elif page == "Process Outputs":
    st.title("Process Output Files")
    process_outputs_page()

elif page == "Calibrate Drift Files":
    calibrate_drift_files_page()

