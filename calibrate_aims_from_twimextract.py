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
        st.dataframe(calibrated_df.head())

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

            calibrated_df = pd.DataFrame(calibrated_data, columns=["CCS", "Drift Time", "Collision Voltage", "Intensity"])
            calibrated_df["Collision Voltage"] = pd.to_numeric(calibrated_df["Collision Voltage"], errors='coerce')
            calibrated_df["CCS"] = pd.to_numeric(calibrated_df["CCS"], errors='coerce')
            calibrated_df = calibrated_df.sort_values(by=["Collision Voltage", "CCS"])

            # Store the result in session state
            st.session_state["calibrated_df"] = calibrated_df

    # If processed data exists, allow customization and visualization
    if "calibrated_df" in st.session_state:
        calibrated_df = st.session_state["calibrated_df"]

        st.write("Calibrated Data:")
        st.dataframe(calibrated_df.head())

        # Customization section
        st.header("📊 CIU Heatmap Customization")

        color_map = st.selectbox("Color Map", ["viridis", "plasma", "inferno", "cividis", "coolwarm", "magma"])
        font_size = st.slider("Font Size", 8, 24, 12, 1)
        figure_size = st.slider("Figure Size (inches)", 5, 15, 10, 1)
        dpi = st.slider("Figure Resolution (DPI)", 100, 1000, 300, 50)

        x_min, x_max = st.slider("Crop Collision Voltage Range",
            float(calibrated_df["Collision Voltage"].min()),
            float(calibrated_df["Collision Voltage"].max()),
            (float(calibrated_df["Collision Voltage"].min()), float(calibrated_df["Collision Voltage"].max()))
        )
        y_min, y_max = st.slider("Crop CCS Range",
            float(calibrated_df["CCS"].min()),
            float(calibrated_df["CCS"].max()),
            (float(calibrated_df["CCS"].min()), float(calibrated_df["CCS"].max()))
        )

        # Prepare pivot table
        heatmap_data = calibrated_df.pivot_table(
            index="CCS",
            columns="Collision Voltage",
            values="Intensity",
            aggfunc="mean"
        )

        # Crop heatmap data based on user-selected ranges
        heatmap_data = heatmap_data.loc[
            (heatmap_data.index >= y_min) & (heatmap_data.index <= y_max),
            (heatmap_data.columns >= x_min) & (heatmap_data.columns <= x_max)
        ]

        fig, ax = plt.subplots(figsize=(figure_size, figure_size), dpi=dpi)

        # Generate the heatmap
        sns.heatmap(
            heatmap_data.sort_index(ascending=False),  # Sort index to get CCS in descending order
            cmap=color_map,
            ax=ax,
            cbar=False,  # Remove the color bar
            xticklabels=True,
            yticklabels=True,
            linewidths=0.5,  # Optional: add some line widths for better separation
            linecolor='white'  # Optional: set the line color between blocks for separation
        )

        # Set axis labels and adjust tick marks
        ax.set_xlabel("Collision Voltage", fontsize=font_size)
        ax.set_ylabel("CCS", fontsize=font_size)
        ax.tick_params(labelsize=font_size)

        # Adjust axes to have sensible tick marks (rounded to nearest 100 or 200)
        ax.set_xticks(np.linspace(x_min, x_max, num=10))  # 10 ticks along x-axis
        ax.set_xticklabels(np.round(np.linspace(x_min, x_max, num=10), 0))  # Round values
        ax.set_yticks(np.linspace(y_min, y_max, num=10))  # 10 ticks along y-axis
        ax.set_yticklabels(np.round(np.linspace(y_min, y_max, num=10), 0))  # Round values

        # User input for x-value labeling (0-5)
        num_x_labels = st.slider("How many x-values to label (0-5)?", 0, 5, 0)
        x_values = []
        x_labels = []
        for i in range(num_x_labels):
            value = st.number_input(f"Enter x-value {i+1}", min_value=float(x_min), max_value=float(x_max))
            label = st.text_input(f"Enter label for x-value {i+1}")
            x_values.append(value)
            x_labels.append(label)

        # User input for y-value labeling (0-5)
        num_y_labels = st.slider("How many y-values to label (0-5)?", 0, 5, 0)
        y_values = []
        y_labels = []
        for i in range(num_y_labels):
            value = st.number_input(f"Enter y-value {i+1}", min_value=float(y_min), max_value=float(y_max))
            label = st.text_input(f"Enter label for y-value {i+1}")
            y_values.append(value)
            y_labels.append(label)

        # Plot the dashed lines and add labels based on user input
        for i in range(num_x_labels):
            x_pos = np.argmin(np.abs(np.array(heatmap_data.columns) - x_values[i]))  # Find nearest x value
            ax.axvline(x=x_pos, color='white', linestyle='--', linewidth=1)
            ax.text(x=x_pos, y=heatmap_data.index[0], s=x_labels[i], color='white', va='bottom', ha='center', fontsize=font_size)

        for i in range(num_y_labels):
            y_pos = np.argmin(np.abs(np.array(heatmap_data.index) - y_values[i]))  # Find nearest y value
            ax.axhline(y=y_pos, color='white', linestyle='--', linewidth=1)
            ax.text(x=heatmap_data.columns[0], y=y_pos, s=y_labels[i], color='white', va='center', ha='left', fontsize=font_size)

        plt.tight_layout()
        st.pyplot(fig)

        # Download buttons
        csv = calibrated_df.to_csv(index=False).encode('utf-8')
        st.download_button("Download Calibrated CSV", data=csv, file_name="calibrated_twim_extract.csv", mime="text/csv")

        buf = BytesIO()
        fig.savefig(buf, format="png", dpi=dpi)
        st.download_button("Download Heatmap PNG", data=buf.getvalue(), file_name="ciu_heatmap.png", mime="image/png")


# Run the app
if __name__ == "__main__":
    twim_extract_page()
