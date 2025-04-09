# app.py
import streamlit as st
from gaussian_fitting import upload_and_fit_gaussians
from calibrate import calibrate_page

# Sidebar for navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Fit Gaussians to Data", "Calibrate"])

# Check if the page has changed, to ensure no redundant loading
if "page" not in st.session_state or st.session_state["page"] != page:
    st.session_state.clear()  # Clear session state for fresh start

# Store the current page in session state
st.session_state["page"] = page

# Home Page
if page == "Home":
    st.title("👩🏻‍🔬 Barran Group IM-MS Processing Tools")
    st.subheader("←←← Navigate to the tool you need from the sidebar.")

# Fit Gaussians to Data Page
elif page == "Fit Gaussians to Data":
    st.title("Fit Gaussians to Data")
    upload_and_fit_gaussians()

# Calibrate Page
elif page == "Calibrate":
    calibrate_page()
