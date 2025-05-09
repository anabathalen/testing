import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import io

PROTON_MASS = 1.007276

# Define Gaussian function
def gaussian(x, amp, mu, sigma):
    return amp * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

def multi_gaussian(x, *params):
    return sum(gaussian(x, *params[i:i+3]) for i in range(0, len(params), 3))

def fit_gaussians(x, y, initial_params):
    bounds = ([], [])
    for i in range(0, len(initial_params), 3):
        amp, mu, sigma = initial_params[i:i+3]
        bounds[0].extend([amp * 0.95, mu * 0.95, sigma * 0.95])
        bounds[1].extend([amp * 1.05 + 1e-9, mu * 1.05 + 1e-9, sigma * 1.05 + 1e-9])
    try:
        popt, _ = curve_fit(multi_gaussian, x, y, p0=initial_params, bounds=bounds)
        return popt
    except Exception as e:
        st.error(f"Fit failed: {e}")
        return initial_params

def fit_ccs_traces_page():
    st.header("Fit CCS Traces")

    cal_file = st.file_uploader("Upload a calibrated CSV file for a protein", type="csv")
    ms_file = st.file_uploader("Upload the mass spectrum TXT file (no headers)", type="txt")
    protein_mass = st.number_input("Enter the protein mass (Da)", min_value=0.0, step=1.0)
    x_start = st.number_input("X-axis start (CCS)", value=1000.0)
    x_end = st.number_input("X-axis end (CCS)", value=3000.0)

    if cal_file and ms_file and protein_mass > 0:
        cal_df = pd.read_csv(cal_file)
        if not {'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'}.issubset(cal_df.columns):
            st.error("The CSV file must contain 'Charge', 'CCS', 'CCS Std.Dev.', and 'Intensity' columns.")
            return

        cal_df = cal_df[cal_df['CCS Std.Dev.'] < 0.5 * cal_df['CCS']].copy()
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df.dropna(inplace=True)

        all_charges = sorted(cal_df["Charge"].unique())
        selected_charges = st.multiselect("Select charge states to include", all_charges, default=all_charges)
        cal_df = cal_df[cal_df["Charge"].isin(selected_charges)]

        scale_factors = {}
        for z in selected_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min, mz_max = mz * 0.98, mz * 1.02
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[z] = intensity_sum

        cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        x_fit = np.linspace(x_start, x_end, 1000)
        fit_results = []
        all_traces = []
        all_params = []

        st.subheader("Fit Individual Charge States")
        for charge in selected_charges:
            sub_df = cal_df[cal_df["Charge"] == charge]
            x = sub_df["CCS"].values
            y = sub_df["Scaled Intensity"].values
            y_interp = np.interp(x_fit, x, y)
            all_traces.append(y_interp)

            st.markdown(f"#### Charge State {charge}")
            num_peaks = st.number_input(f"Number of peaks for charge {charge}", min_value=1, max_value=5, value=1, key=f"num_peaks_{charge}")

            init_params = []
            for i in range(num_peaks):
                amp = st.text_input(f"Amp {i+1} for charge {charge}", value=str(max(y)), key=f"amp_{charge}_{i}")
                mu = st.text_input(f"Mu {i+1} for charge {charge}", value=str(x[np.argmax(y)]), key=f"mu_{charge}_{i}")
                sigma = st.text_input(f"Sigma {i+1} for charge {charge}", value="10.0", key=f"sigma_{charge}_{i}")
                init_params.extend([float(amp), float(mu), float(sigma)])

            if st.button(f"Optimize Fit for Charge {charge}"):
                optimized_params = fit_gaussians(x_fit, y_interp, init_params)
            else:
                optimized_params = init_params

            fitted_y = multi_gaussian(x_fit, *optimized_params)
            all_params.extend([[charge] + optimized_params[i:i+3] for i in range(0, len(optimized_params), 3)])

            fig, ax = plt.subplots()
            ax.plot(x_fit, y_interp, label="Raw Data")
            ax.plot(x_fit, fitted_y, label="Fitted Curve")
            ax.legend()
            st.pyplot(fig)

        st.subheader("Fit Summed Trace")
        summed_trace = np.sum(all_traces, axis=0)

        show_summed = st.checkbox("Enable summed trace fitting")
        if show_summed:
            num_peaks_sum = st.number_input("Number of peaks for summed trace", min_value=1, max_value=5, value=1, key="num_peaks_sum")
            init_params_sum = []
            for i in range(num_peaks_sum):
                amp = st.text_input(f"Summed Amp {i+1}", value=str(max(summed_trace)), key=f"amp_sum_{i}")
                mu = st.text_input(f"Summed Mu {i+1}", value=str(x_fit[np.argmax(summed_trace)]), key=f"mu_sum_{i}")
                sigma = st.text_input(f"Summed Sigma {i+1}", value="10.0", key=f"sigma_sum_{i}")
                init_params_sum.extend([float(amp), float(mu), float(sigma)])

            if st.button("Optimize Summed Fit"):
                optimized_sum_params = fit_gaussians(x_fit, summed_trace, init_params_sum)
            else:
                optimized_sum_params = init_params_sum

            fitted_summed = multi_gaussian(x_fit, *optimized_sum_params)
            fig2, ax2 = plt.subplots()
            ax2.plot(x_fit, summed_trace, label="Summed Raw")
            ax2.plot(x_fit, fitted_summed, label="Summed Fit")
            ax2.legend()
            st.pyplot(fig2)

        # Save all charge state fit params to CSV
        if all_params:
            fit_df = pd.DataFrame(all_params, columns=["Charge", "Amplitude", "Mu", "Sigma"])
            csv_data = fit_df.to_csv(index=False).encode("utf-8")
            st.download_button("Download All Gaussian Parameters (CSV)", csv_data, file_name="gaussian_parameters.csv", mime="text/csv")
