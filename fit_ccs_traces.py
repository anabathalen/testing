# fit_ccs_traces.py
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

def gaussian(x, mean, std, amplitude):
    return amplitude * norm.pdf(x, mean, std)

def fit_ccs_traces_page():
    st.title("Fit CCS Traces with Custom Gaussians")

    uploaded_file = st.file_uploader("Upload a CSV with columns: CCS, Intensity, Charge", type="csv")

    if uploaded_file:
        df = pd.read_csv(uploaded_file)
        if not {"CCS", "Intensity", "Charge"}.issubset(df.columns):
            st.error("Your file must contain 'CCS', 'Intensity', and 'Charge' columns.")
            return

        df = df.dropna(subset=["CCS", "Intensity", "Charge"])
        charges = sorted(df["Charge"].unique())

        fit_mode = st.radio("Choose fit mode:", ["Sum first, then fit", "Fit individual charge states, then sum"])

        # User-defined x range
        st.subheader("Define X-axis (CCS) Range")
        default_min = float(df["CCS"].min())
        default_max = float(df["CCS"].max())
        x_min = float(st.text_input("Minimum CCS (x-axis)", value=str(round(default_min, 2))))
        x_max = float(st.text_input("Maximum CCS (x-axis)", value=str(round(default_max, 2))))
        if x_min >= x_max:
            st.error("Minimum CCS must be less than maximum CCS.")
            return

        ccs_grid = np.linspace(x_min, x_max, 1000)
        interpolated = {}

        for charge, group in df.groupby("Charge"):
            group_sorted = group.sort_values("CCS")
            interp = np.interp(ccs_grid, group_sorted["CCS"], group_sorted["Intensity"], left=0, right=0)
            interpolated[charge] = interp

        if fit_mode == "Sum first, then fit":
            target_y = np.sum(list(interpolated.values()), axis=0)
        else:
            target_y = np.sum([interpolated[c] for c in charges], axis=0)

        target_x = ccs_grid

        st.subheader("Specify Number of Gaussian Peaks")
        num_peaks = st.number_input("Number of Gaussians", min_value=1, max_value=10, value=2, step=1)

        st.markdown("Enter parameters for each Gaussian (mean, std dev, amplitude):")

        peaks = []
        for i in range(num_peaks):
            st.markdown(f"**Gaussian {i+1}**")
            col1, col2, col3 = st.columns(3)
            with col1:
                mean = st.text_input(f"Mean {i+1}", value=str(x_min + (i+1)*(x_max - x_min)/(num_peaks + 1)))
            with col2:
                std = st.text_input(f"Std Dev {i+1}", value="3.0")
            with col3:
                amp = st.text_input(f"Amplitude {i+1}", value=str(np.max(target_y)/2))

            try:
                peaks.append((float(mean), float(std), float(amp)))
            except ValueError:
                st.error(f"Could not parse Gaussian {i+1} parameters. Ensure all values are numeric.")
                return

        # Plotting
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(target_x, target_y, label="Original", color="black")

        combined_fit = np.zeros_like(target_x)
        for i, (mean, std, amp) in enumerate(peaks):
            gauss = gaussian(target_x, mean, std, amp)
            combined_fit += gauss
            ax.plot(target_x, gauss, "--", label=f"Gaussian {i+1}")

        ax.plot(target_x, combined_fit, color="red", linewidth=2, label="Total Fit")
        ax.set_xlim(x_min, x_max)
        ax.set_xlabel("CCS")
        ax.set_ylabel("Intensity")
        ax.legend()
        st.pyplot(fig)
