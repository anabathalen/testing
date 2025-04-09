import pandas as pd
import os
import streamlit as st

# Function to read bush.csv and get protein data
def read_bush_csv():
    calibrant_file_path = os.path.join(os.path.dirname(__file__), 'bush.csv')
    if os.path.exists(calibrant_file_path):
        bush_df = pd.read_csv(calibrant_file_path)
    else:
        st.error(f"'{calibrant_file_path}' not found. Make sure 'bush.csv' is in the same directory as the script.")
        bush_df = pd.DataFrame()  # Empty DataFrame if the file isn't found
    return bush_df

# Function to generate the .dat file
def generate_dat_file(bush_df, velocity, voltage, pressure, length):
    # Choose the correct column based on user input for calibrant type
    calibrant_type = st.selectbox("Select Calibrant Type", options=["He", "N2"])
    calibrant_column = 'CCS_he' if calibrant_type == 'He' else 'CCS_n2'

    # Open the file for writing
    dat_file_path = os.path.join(os.path.dirname(__file__), 'calibration_data.dat')
    with open(dat_file_path, 'w') as f:
        # Write metadata information
        f.write(f"# length {length}\n")
        f.write(f"# velocity {velocity}\n")
        f.write(f"# voltage {voltage}\n")
        f.write(f"# pressure {pressure}\n")

        # Iterate through the protein data and write to the .dat file
        for _, row in bush_df.iterrows():
            protein = row['protein']
            charge_states = range(13, 21)  # Example range for charge states (adjust if needed)
            
            for charge_state in charge_states:
                # Create the protein name (e.g., myoglobin_13)
                protein_name = f"{protein}_{charge_state}"
                
                # Get the mass and CCS value (calibrant data)
                mass = row['mass']
                reference_value = row[calibrant_column] * 100  # Multiply CCS value by 100
                
                # For now, we can use some estimated drift time or a placeholder value (adjust later)
                drift_time = reference_value / 100  # Placeholder drift time (you'll replace this with actual values)

                # Write the protein data line in the specified format
                f.write(f"{protein_name} {mass} {charge_state} {reference_value} {drift_time}\n")

    st.success(f"Calibration data .dat file generated successfully: {dat_file_path}")

# Streamlit UI for calibration page
def calibrate_page():
    st.title("Calibration Data Generator")

    # Step 1: User inputs for calibration parameters
    velocity = st.number_input("Enter velocity (e.g., 281)", min_value=1.0)
    voltage = st.number_input("Enter voltage (e.g., 20)", min_value=1.0)
    pressure = st.number_input("Enter pressure (e.g., 1.63)", min_value=0.1)
    length = st.number_input("Enter length (e.g., 0.980)", min_value=0.1)

    # Step 2: Load the bush.csv file
    bush_df = read_bush_csv()

    # Step 3: Generate the .dat file upon button press
    if st.button("Generate .dat File"):
        generate_dat_file(bush_df, velocity, voltage, pressure, length)

# Run the Streamlit app
if __name__ == "__main__":
    calibrate_page()

