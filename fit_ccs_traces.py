import streamlit as st
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.pyplot as plt

def gaussian(x, mean, std, amplitude):
    return amplitude * norm.pdf(x, mean, std)

def sum_of_gaussians(x, *params):
    n = len(params) // 3
    y = np.zeros_like(x)
    for i in range(n):
        m, s, a = params[i*3:i*3+3]
        y += gaussian(x, m, s, a)
    return y

def fit_ccs_traces_page():
    st.title("Fit CCS Traces with Gaussian Models")
    st.markdown("Upload a CSV with columns: `CCS`, `Intensity`, and optionally `Charge`. You can fit all data summed or each charge state separately.")

    uploaded_file = st.file_uploader("Upload CSV", type=["csv"])
    if uploaded_file is None:
        return

    df = pd.read_csv(uploaded_file)
    if not {"CCS", "Intensity"}.issubset(df.columns):
        st.error("CSV must contain 'CCS' and 'Intensity' columns.")
        return

    fit_mode = st.radio("Choose fitting mode", ["Fit Summed Trace", "Fit Each Charge Separately"])

    x_min, x_max = st.slider("Select x-axis range (CCS)", float(df["CCS"].min()), float(df["CCS"].max()), (float(df["CCS"].min()), float(df["CCS"].max())))

    x_fit = np.linspace(x_min, x_max, 1000)

    if fit_mode == "Fit Each Charge Separately" and "Charge" not in df.columns:
        st.warning("No 'Charge' column found. Defaulting to summed trace.")
        fit_mode = "Fit Summed Trace"

    data_groups = [("Summed", df)] if fit_mode == "Fit Summed Trace" else list(df.groupby("Charge"))

    for group_label, group_df in data_groups:
        subset = group_df[(group_df["CCS"] >= x_min) & (group_df["CCS"] <= x_max)]
        x = subset["CCS"].values
        y = subset["Intensity"].values
        y = y / y.max()  # Normalize

        st.markdown(f"### {group_label}")
        num_peaks = st.number_input(f"Number of Gaussians for {group_label}", min_value=1, max_value=10, value=2, key=f"npeaks_{group_label}")
        peak_params = []

        for i in range(num_peaks):
            m = st.text_input(f"Mean of peak {i+1}", value=str(round(x.min() + (i+1)*(x.max()-x.min())/(num_peaks+1), 2)), key=f"mean_{i}_{group_label}")
            s = st.text_input(f"Std dev of peak {i+1}", value="5.0", key=f"std_{i}_{group_label}")
            a = st.text_input(f"Amplitude of peak {i+1}", value="1.0", key=f"amp_{i}_{group_label}")
            peak_params += [float(m), float(s), float(a)]

        autofit = st.button(f"Auto-Fit Gaussians for {group_label}")

        y_fit = sum_of_gaussians(x_fit, *peak_params)

        if autofit:
            bounds_lower = [p * 0.95 for p in peak_params]
            bounds_upper = [p * 1.05 for p in peak_params]
            try:
                popt, _ = curve_fit(sum_of_gaussians, x, y, p0=peak_params, bounds=(bounds_lower, bounds_upper))
                st.success("Optimized fit successful.")
                y_fit = sum_of_gaussians(x_fit, *popt)
                st.text("Fitted parameters:\n" + str(np.round(popt, 3)))
            except Exception as e:
                st.error(f"Fit failed: {e}")

        fig, ax = plt.subplots()
        ax.plot(x, y, label="Data")
        ax.plot(x_fit, y_fit, label="Fit", linestyle="--")
        ax.set_title(f"Gaussian Fit: {group_label}")
        ax.set_xlabel("CCS")
        ax.set_ylabel("Normalized Intensity")
        ax.legend()
        st.pyplot(fig)
