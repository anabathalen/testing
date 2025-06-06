
# app.py
import streamlit as st
from gaussian_fitting import upload_and_plot
from calibrate import calibrate_page
from generate_input_files import generate_input_dat_files_app
from process_outputs import process_outputs_page
from calibrate_drift_files import calibrate_drift_files_page
from calibrate_aims_from_twimextract import twim_extract_page
from process_plot_ims import plot_and_scale_page
from fit_ccs_traces import fit_ccs_traces_page


# Sidebar for navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Fit Gaussians to Data", "Fit CCS Traces", "Calibrate", "Generate Input Files", "Process Outputs", "Calibrate Drift Files", "Calibrate CIU", "Process/Plot IMS"])

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
    st.title("Get Input Files (for IMSCal)")
    st.subheader("Takes ATDs of Samples and Returns Input Files for IMSCal.")
    st.markdown('''Make a zipped folder (like you did for the calibration) containing a folder for each protein titled the name of the protein. Each folder should contain text files called x.txt where x is the charge state.''')
    generate_input_dat_files_app()

# Process Outputs Page
elif page == "Process Outputs":
    st.title("Process Output Files (from IMSCal)")
    st.subheader("Combine IMSCal Outputs to Single CSV Per Protein.")
    st.markdown('''IMSCal will have populated your protein folders containing input files with corresponding calibrated output files for each charge state - this tool will compile all that information into one csv file per protein. The csv files you generate here allow you to calibrate any other ATDs obtained on the same day for the same protein, which will be useful e.g. for aIMS.''')
    process_outputs_page()

elif page == "Calibrate Drift Files":
    calibrate_drift_files_page()

elif page == "Calibrate CIU":
    twim_extract_page()

elif page == "Process/Plot IMS":
    st.title("Process/Plot IMS")
    plot_and_scale_page()
    
elif page == "Fit CCS Traces":
    st.title("Fit CCS Traces")
    fit_ccs_traces_page()

