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

    # Upload CSV
    uploaded_file = st.file_uploader("Upload a CSV with columns: CCS, Intensity, Charge", type="csv")

    if uploaded_file:
        df = pd.read_csv(uploaded_file)
        if not {"CCS", "Intensity", "Charge"}.issubset(df.columns):
            st.error("Your file must contain 'CCS', 'Intensity', and 'Charge' columns.")
            return

        df = df.dropna(subset=["CCS", "Intensity", "Charge"])
        charges = sorted(df["Charge"].unique())
        fit_mode = st.radio("Choose fit mode:", ["Sum first, then fit", "Fit individual charge states, then sum"])

        ccs_grid = np.linspace(df["CCS"].min(), df["CCS"].max(), 1000)
        interpolated = {}

        # Create traces
        for charge, group in df.groupby("Charge"):
            group_sorted = group.sort_values("CCS")
            interp = np.interp(ccs_grid, group_sorted["CCS"], group_sorted["Intensity"], left=0, right=0)
            interpolated[charge] = interp

        if fit_mode == "Sum first, then fit":
            summed = np.sum(list(interpolated.values()), axis=0)
            target_x = ccs_grid
            target_y = summed
        else:
            fits = []
            for charge in charges:
                fits.append(interpolated[charge])
            target_x = ccs_grid
            target_y = np.sum(fits, axis=0)

        # Ask user how many Gaussians
        num_peaks = st.slider("Number of Gaussian peaks to fit", 1, 6, 2)

        st.markdown("Adjust each Gaussian peak's parameters:")

        peaks = []
        for i in range(num_peaks):
            st.subheader(f"Gaussian {i+1}")
            col1, col2, col3 = st.columns(3)
            with col1:
                mean = st.slider(f"Mean {i+1}", float(target_x.min()), float(target_x.max()), float(target_x.min() + (i+1) * 10), step=0.1)
            with col2:
                std = st.slider(f"Std Dev {i+1}", 0.1, (target_x.max() - target_x.min()) / 2, 3.0, step=0.1)
            with col3:
                amplitude = st.slider(f"Amplitude {i+1}", 0.0, max(target_y)*1.5, max(target_y)/2, step=0.1)
            peaks.append((mean, std, amplitude))

        # Plotting
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(target_x, target_y, label="Original", color="black")

        combined_fit = np.zeros_like(target_x)
        for i, (mean, std, amp) in enumerate(peaks):
            gauss = gaussian(target_x, mean, std, amp)
            combined_fit += gauss
            ax.plot(target_x, gauss, linestyle="--", label=f"Gaussian {i+1}")

        ax.plot(target_x, combined_fit, color="red", linewidth=2, label="Total Fit")
        ax.legend()
        ax.set_xlabel("CCS")
        ax.set_ylabel("Intensity")
        st.pyplot(fig)
