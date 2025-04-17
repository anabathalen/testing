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

        ccs_labels = st.multiselect("Label specific CCS values", options=np.round(calibrated_df["CCS"].unique(), 2).tolist())
        cv_labels = st.multiselect("Label specific Collision Voltages", options=np.round(calibrated_df["Collision Voltage"].unique(), 2).tolist())

        # Prepare pivot table
        heatmap_data = calibrated_df.pivot_table(
            index="CCS",
            columns="Collision Voltage",
            values="Intensity",
            aggfunc="mean"
        )

        heatmap_data = heatmap_data.loc[
            (heatmap_data.index >= y_min) & (heatmap_data.index <= y_max),
            (heatmap_data.columns >= x_min) & (heatmap_data.columns <= x_max)
        ]

        fig, ax = plt.subplots(figsize=(figure_size, figure_size), dpi=dpi)
        sns.heatmap(
            heatmap_data.sort_index(ascending=False),
            cmap=color_map,
            ax=ax,
            cbar=True,
            xticklabels=True,
            yticklabels=True
        )

        ax.set_xlabel("Collision Voltage", fontsize=font_size)
        ax.set_ylabel("CCS", fontsize=font_size)
        ax.set_title("CIU Heatmap", fontsize=font_size + 2)
        ax.tick_params(labelsize=font_size)

        # Add optional label lines
        for label in ccs_labels:
            if label in heatmap_data.index:
                ax.axhline(y=heatmap_data.index.get_loc(label), color='white', linestyle='--', linewidth=1)
        for label in cv_labels:
            if label in heatmap_data.columns:
                ax.axvline(x=heatmap_data.columns.get_loc(label), color='white', linestyle='--', linewidth=1)

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
