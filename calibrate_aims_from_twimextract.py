import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from io import BytesIO

def twim_extract_page():
    st.title("TWIM CCS Calibration and CIU Heatmap Generator")

    # Upload files
    twim_extract_file = st.file_uploader("Upload the TWIM Extract CSV file", type="csv")
    calibration_file = st.file_uploader("Upload the calibration CSV file", type="csv")
    
    if twim_extract_file and calibration_file:
        # Read TWIM extract data
        twim_df = pd.read_csv(twim_extract_file, header=None)
        twim_df.columns = ["Drift Time"] + [str(i) for i in range(1, len(twim_df.columns))]
        st.write("Uploaded TWIM Extract Data:")
        st.dataframe(twim_df.head())

        # Read calibration data
        cal_df = pd.read_csv(calibration_file)
        st.write("Uploaded Calibration Data:")
        st.dataframe(cal_df.head())

        if 'Z' not in cal_df.columns:
            st.error("Calibration data must include a 'Z' column for charge state.")
            return

        data_type = st.radio("Is your data from a Synapt or Cyclic instrument?", ["Synapt", "Cyclic"])
        charge_state = st.number_input("Enter the charge state of the protein (Z)", min_value=1, max_value=10, step=1)
        inject_time = None

        if data_type == "Cyclic":
            inject_time = st.number_input("Enter the injection time (ms)", min_value=0.0, value=0.0, step=0.1)

        # Button to process data
        if st.button("Process Data"):
            if inject_time is not None and data_type == "Cyclic":
                twim_df["Drift Time"] = twim_df["Drift Time"] - inject_time

            cal_data = cal_df[cal_df["Z"] == charge_state]
            if cal_data.empty:
                st.error(f"No calibration data found for charge state {charge_state}")
                return

            if "Drift" not in cal_data.columns or "CCS" not in cal_data.columns:
                st.error("Calibration data must include 'Drift' and 'CCS' columns.")
                return

            cal_data["CCS Std.Dev."] = cal_data["CCS Std.Dev."].fillna(0)
            cal_data = cal_data[cal_data["CCS Std.Dev."] <= 0.1 * cal_data["CCS"]]
            cal_data["Drift (ms)"] = cal_data["Drift"] * 1000

            calibrated_data = []

            drift_times = twim_df["Drift Time"]
            collision_voltages = twim_df.columns[1:]

            for idx, drift_time in enumerate(drift_times):
                intensities = twim_df.iloc[idx, 1:].values
                if pd.isna(drift_time):
                    continue
                drift_time_rounded = round(drift_time, 4)
                closest_idx = (cal_data["Drift (ms)"] - drift_time_rounded).abs().idxmin()
                ccs_value = cal_data.loc[closest_idx, "CCS"]

                for col_idx, intensity in enumerate(intensities):
                    cv = collision_voltages[col_idx]
                    calibrated_data.append([ccs_value, drift_time, float(cv), intensity])

            # Convert calibrated data to NumPy array
            calibrated_array = np.array(calibrated_data)

            # Store the result in session state
            st.session_state["calibrated_array"] = calibrated_array

    # If processed data exists, allow customization and visualization
    if "calibrated_array" in st.session_state:
        calibrated_array = st.session_state["calibrated_array"]

        # Customization section
        st.header("📊 CIU Heatmap Customization")

        color_map = st.selectbox("Color Map", ["viridis", "plasma", "inferno", "cividis", "coolwarm", "magma"])
        font_size = st.slider("Font Size", 8, 24, 12, 1)
        figure_size = st.slider("Figure Size (inches)", 5, 15, 10, 1)
        dpi = st.slider("Figure Resolution (DPI)", 100, 1000, 300, 50)

        x_min, x_max = st.slider("Crop Collision Voltage Range",
            float(calibrated_array[:, 2].min()),
            float(calibrated_array[:, 2].max()),
            (float(calibrated_array[:, 2].min()), float(calibrated_array[:, 2].max()))
        )
        y_min, y_max = st.slider("Crop CCS Range",
            float(calibrated_array[:, 0].min()),
            float(calibrated_array[:, 0].max()),
            (float(calibrated_array[:, 0].min()), float(calibrated_array[:, 0].max()))
        )

        # Create heatmap grid with CCS and Collision Voltage as axes
        grid_x = np.linspace(x_min, x_max, num=100)
        grid_y = np.linspace(y_min, y_max, num=100)

        # Create a meshgrid for the x and y axes
        X, Y = np.meshgrid(grid_x, grid_y)

        # Interpolate intensities over the meshgrid
        from scipy.interpolate import griddata
        Z = griddata(
            (calibrated_array[:, 2], calibrated_array[:, 0]),  # Points
            calibrated_array[:, 3],  # Intensity
            (X, Y),  # Grid
            method='cubic'  # Interpolation method
        )

        # Plot the heatmap
        fig, ax = plt.subplots(figsize=(figure_size, figure_size), dpi=dpi)
        c = ax.pcolormesh(X, Y, Z, cmap=color_map, shading='auto')

        # Customize the plot with labels and ticks
        ax.set_xlabel("Collision Voltage", fontsize=font_size)
        ax.set_ylabel("CCS", fontsize=font_size)
        ax.tick_params(labelsize=font_size)

        # Add color bar
        fig.colorbar(c, ax=ax)

        # Add dashed lines based on user input for x-values and y-values
        num_x_labels = st.slider("How many x-values to label (0-5)?", 0, 5, 0)
        x_values = []
        x_labels = []
        for i in range(num_x_labels):
            value = st.number_input(f"Enter x-value {i+1}", min_value=float(x_min), max_value=float(x_max))
            label = st.text_input(f"Enter label for x-value {i+1}")
            x_values.append(value)
            x_labels.append(label)

        num_y_labels = st.slider("How many y-values to label (0-5)?", 0, 5, 0)
        y_values = []
        y_labels = []
        for i in range(num_y_labels):
            value = st.number_input(f"Enter y-value {i+1}", min_value=float(y_min), max_value=float(y_max))
            label = st.text_input(f"Enter label for y-value {i+1}")
            y_values.append(value)
            y_labels.append(label)

        # Plot the dashed lines and labels based on user input
        for i in range(num_x_labels):
            ax.axvline(x=x_values[i], color='white', linestyle='--', linewidth=1)
            ax.text(x=x_values[i], y=y_max, s=x_labels[i], color='white', va='bottom', ha='center', fontsize=font_size)

        for i in range(num_y_labels):
            ax.axhline(y=y_values[i], color='white', linestyle='--', linewidth=1)
            ax.text(x=x_min, y=y_values[i], s=y_labels[i], color='white', va='center', ha='left', fontsize=font_size)

        plt.tight_layout()
        st.pyplot(fig)

        # Download buttons
        csv = pd.DataFrame(calibrated_array, columns=["CCS", "Drift Time", "Collision Voltage", "Intensity"]).to_csv(index=False).encode('utf-8')
        st.download_button("Download Calibrated CSV", data=csv, file_name="calibrated_twim_extract.csv", mime="text/csv")
        img = BytesIO()
        fig.savefig(img, format='png', bbox_inches="tight")
        img.seek(0)
        st.download_button("Download CIU Heatmap Image", data=img, file_name="ciu_heatmap.png", mime="image/png")


# Run the app
if __name__ == "__main__":
    twim_extract_page()
