import streamlit as st
import pandas as pd
import numpy as np
from io import StringIO

def twim_extract_page():
    # Upload TWIM extract CSV file (with drift time and intensity data)
    twim_extract_file = st.file_uploader("Upload the TWIM Extract CSV file", type="csv")
    
    # Upload the calibration CSV file
    calibration_file = st.file_uploader("Upload the calibration CSV file", type="csv")
    
    if twim_extract_file and calibration_file:
        # Read the TWIM extract data
        twim_df = pd.read_csv(twim_extract_file, header=None)
        
        # Explicitly set the first column as "Drift Time" and use the rest of the columns' headers
        twim_df.columns = ["Drift Time"] + [str(i) for i in range(1, len(twim_df.columns))]
        
        # Display the first few rows for user confirmation
        st.write("Uploaded TWIM Extract Data:")
        st.dataframe(twim_df.head())

        # Read the calibration data
        cal_df = pd.read_csv(calibration_file)
        
        # Display the calibration data
        st.write("Uploaded Calibration Data:")
        st.dataframe(cal_df.head())
        
        # Check if the 'Z' column exists in the calibration file
        if 'Z' not in cal_df.columns:
            st.error("Calibration data must include a 'Z' column for charge state. Please check the file.")
            return

        # Ask for the data type (Synapt or Cyclic)
        data_type = st.radio("Is your data from a Synapt or Cyclic instrument?", ["Synapt", "Cyclic"])
        
        # Ask for the charge state (Z)
        charge_state = st.number_input("Enter the charge state of the protein (Z)", min_value=1, max_value=10, step=1)

        # If the data is Cyclic, ask for the injection time
        inject_time = None
        if data_type == "Cyclic":
            inject_time = st.number_input("Enter the injection time (ms)", min_value=0.0, value=0.0, step=0.1)

        # Process the TWIM extract data
        if st.button("Process Data"):
            # Adjust drift time for Cyclic data (if needed)
            if inject_time is not None and data_type == "Cyclic":
                twim_df["Drift Time"] = twim_df["Drift Time"] - inject_time  # Subtract the injection time from the drift time

            # Extract calibration data for the given charge state (Z)
            cal_data = cal_df[cal_df["Z"] == charge_state]

            if cal_data.empty:
                st.error(f"No calibration data found for charge state {charge_state}")
                return
            
            # Ensure the calibration data has the necessary columns
            if "Drift" not in cal_data.columns or "CCS" not in cal_data.columns:
                st.error("Calibration data must include 'Drift' and 'CCS' columns.")
                return

            # Filter out calibration data where the standard deviation of CCS is greater than 10% of the CCS value
            cal_data['CCS Std.Dev.'] = cal_data['CCS Std.Dev.'].fillna(0)  # Handle missing std deviation values
            cal_data = cal_data[cal_data['CCS Std.Dev.'] <= 0.1 * cal_data['CCS']]  # Keep only valid data
            
            # Now, calibrate the TWIM extract data
            calibrated_data = []

            # The first column of the TWIM extract file contains drift time
            drift_times = twim_df["Drift Time"]

            # The remaining columns represent the intensity at various collision energies
            collision_voltages = twim_df.columns[1:]  # Starting from the second column (column indices 1 to N)

            for idx, drift_time in enumerate(drift_times):
                intensities = twim_df.iloc[idx, 1:].values  # All columns after the first one are intensities
                
                # Find the closest drift time in the calibration data
                if pd.isna(drift_time):
                    # Skip if the drift time is NaN
                    continue
                
                # Round drift time for precision matching
                drift_time_rounded = round(drift_time, 4)  # Adjust to the desired decimal precision

                # Multiply the calibration drift times by 1000 (to match TWIM drift times in ms)
                cal_data["Drift (ms)"] = cal_data["Drift"] * 1000

                # Find the closest matching drift time
                closest_drift_idx = (cal_data["Drift (ms)"] - drift_time_rounded).abs().idxmin()
                
                ccs_value = cal_data.loc[closest_drift_idx, "CCS"]
                
                # Store the calibrated data
                for col_idx, intensity in enumerate(intensities):
                    collision_voltage = collision_voltages[col_idx]  # The collision voltage is in the column header
                    calibrated_data.append([ccs_value, drift_time, collision_voltage, intensity])

            # Create a DataFrame from the calibrated data
            calibrated_df = pd.DataFrame(calibrated_data, columns=["CCS", "Drift Time", "Collision Voltage", "Intensity"])

            # Display the calibrated data
            st.write("Calibrated Data:")
            st.dataframe(calibrated_df.head())

            # Provide a download button for the calibrated CSV
            csv = calibrated_df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="Download Calibrated Data",
                data=csv,
                file_name="calibrated_twim_extract.csv",
                mime="text/csv"
            )
