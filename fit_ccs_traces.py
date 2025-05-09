# fit_ccs_traces.py

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from io import StringIO
import base64

PROTON_MASS = 1.007276

def gaussian(x, mu, sigma, amp):
    return amp * np.exp(-(x - mu)**2 / (2 * sigma**2))

def sum_of_gaussians(x, *params):
    n = len(params) // 3
    return sum(gaussian(x, params[3*i], params[3*i+1], params[3*i+2]) for i in range(n))

def optimize_fit(x_data, y_data, initial_params):
    try:
        bounds_lower = [p * 0.95 for p in initial_params]
        bounds_upper = [p * 1.05 for p in initial_params]
        popt, _ = curve_fit(sum_of_gaussians, x_data, y_data, p0=initial_params, bounds=(bounds_lower, bounds_upper))
        return popt
    except Exception as e:
        st.warning(f"Fit failed: {e}")
        return initial_params

def download_link(df, filename, label):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">{label}</a>'
    return href

def fit_ccs_traces_page():
    st.header("Fit CCS Traces")

    cal_file = st.file_uploader("Upload calibrated CSV", type="csv")
    ms_file = st.file_uploader("Upload mass spectrum TXT", type="txt")
    protein_mass = st.number_input("Enter protein mass (Da)", min_value=0.0, step=1.0)

    x_min = st.number_input("X-axis min (CCS)", min_value=0.0, value=1000.0, step=10.0)
    x_max = st.number_input("X-axis max (CCS)", min_value=0.0, value=3000.0, step=10.0)

    if cal_file and ms_file and protein_mass > 0:
        cal_df = pd.read_csv(cal_file)
        if not {'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'}.issubset(cal_df.columns):
            st.error("Calibrated file must include 'Charge', 'CCS', 'CCS Std.Dev.', and 'Intensity'")
            return

        cal_df = cal_df[cal_df['CCS Std.Dev.'] < 0.5 * cal_df['CCS']].copy()
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df.dropna(inplace=True)

        all_charges = sorted(cal_df["Charge"].unique())
        selected_charges = st.multiselect("Select charge states", all_charges, default=all_charges)
        cal_df = cal_df[cal_df["Charge"].isin(selected_charges)]

        scale_factors = {}
        for z in selected_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_range = ms_df[(ms_df["m/z"] >= mz * 0.98) & (ms_df["m/z"] <= mz * 1.02)]
            scale_factors[z] = mz_range["Intensity"].sum()

        cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        x_common = np.linspace(x_min, x_max, 1000)
        all_fit_data = []
        summed_y = np.zeros_like(x_common)

        for z in selected_charges:
            subset = cal_df[cal_df["Charge"] == z]
            f = interp1d(subset["CCS"], subset["Scaled Intensity"], kind="linear", bounds_error=False, fill_value=0)
            y_interp = f(x_common)
            summed_y += y_interp

            st.subheader(f"Charge State {z}")
            st.line_chart(pd.DataFrame({"CCS": x_common, "Scaled Intensity": y_interp}))

            n_peaks = st.number_input(f"# of Gaussians for {z}", min_value=1, max_value=5, value=1, key=f"num_{z}")
            initial_params = []

            for i in range(n_peaks):
                mu = st.number_input(f"Mean (mu) Peak {i+1} [z={z}]", value=float(x_common[np.argmax(y_interp)]), key=f"mu_{z}_{i}")
                sigma = st.number_input(f"Std Dev (sigma) Peak {i+1} [z={z}]", value=30.0, key=f"sigma_{z}_{i}")
                amp = st.number_input(f"Amplitude Peak {i+1} [z={z}]", value=max(y_interp), key=f"amp_{z}_{i}")
                initial_params.extend([mu, sigma, amp])

            if st.button(f"Optimize Fit for Charge {z}"):
                fitted_params = optimize_fit(x_common, y_interp, initial_params)
            else:
                fitted_params = initial_params

            fit = sum_of_gaussians(x_common, *fitted_params)

            fig, ax = plt.subplots()
            ax.plot(x_common, y_interp, label="Raw", linestyle="--")
            ax.plot(x_common, fit, label="Fitted")
            ax.set_title(f"Charge {z} Fit")
            ax.legend()
            st.pyplot(fig)

            fit_df = pd.DataFrame({"CCS": x_common, f"Fit z={z}": fit})
            all_fit_data.append(fit_df)
            st.markdown(download_link(fit_df, f"fit_charge_{z}.csv", f"Download Fit CSV for z={z}"), unsafe_allow_html=True)

        st.subheader("Summed Fit")
        n_sum = st.number_input("# Gaussians for Summed Fit", min_value=1, max_value=5, value=1, key="sum_gauss")
        sum_init = []

        for i in range(n_sum):
            mu = st.number_input(f"Summed Mean {i+1}", value=float(x_common[np.argmax(summed_y)]), key=f"sum_mu_{i}")
            sigma = st.number_input(f"Summed StdDev {i+1}", value=30.0, key=f"sum_sigma_{i}")
            amp = st.number_input(f"Summed Amplitude {i+1}", value=max(summed_y), key=f"sum_amp_{i}")
            sum_init.extend([mu, sigma, amp])

        if st.button("Optimize Summed Fit"):
            summed_params = optimize_fit(x_common, summed_y, sum_init)
        else:
            summed_params = sum_init

        summed_fit = sum_of_gaussians(x_common, *summed_params)
        fig2, ax2 = plt.subplots()
        ax2.plot(x_common, summed_y, label="Raw Summed", linestyle="--")
        ax2.plot(x_common, summed_fit, label="Fitted Summed")
        ax2.set_title("Summed Fit")
        ax2.legend()
        st.pyplot(fig2)

        summed_df = pd.DataFrame({"CCS": x_common, "Summed Fit": summed_fit})
        st.markdown(download_link(summed_df, "summed_fit.csv", "Download Summed Fit CSV"), unsafe_allow_html=True)
