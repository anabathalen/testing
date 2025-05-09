import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from io import StringIO
import base64

# Constants
PROTON_MASS = 1.007276

# Gaussian function
def gaussian(x, amp, mean, std):
    return amp * np.exp(-(x - mean) ** 2 / (2 * std ** 2))

def multi_gaussian(x, *params):
    return sum(gaussian(x, params[i], params[i+1], params[i+2]) for i in range(0, len(params), 3))

# Fit a multi-Gaussian model
def fit_gaussians(x_data, y_data, initial_params, bounds):
    try:
        popt, _ = curve_fit(multi_gaussian, x_data, y_data, p0=initial_params, bounds=bounds, maxfev=10000)
        return popt
    except Exception as e:
        st.warning(f"Fit failed: {e}")
        return initial_params

def download_csv_button(csv_str, filename, label):
    b64 = base64.b64encode(csv_str.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">{label}</a>'
    st.markdown(href, unsafe_allow_html=True)

def fit_ccs_traces_page():
    st.header("Fit CCS Traces")

    cal_file = st.file_uploader("Upload a calibrated CSV file for a protein", type="csv")
    ms_file = st.file_uploader("Upload the mass spectrum TXT file (no headers)", type="txt")
    protein_mass = st.number_input("Enter the protein mass (Da)", min_value=0.0, step=1.0)
    x_min = st.number_input("X-axis minimum (CCS)", value=500.0)
    x_max = st.number_input("X-axis maximum (CCS)", value=3000.0)

    if cal_file and ms_file and protein_mass > 0:
        cal_df = pd.read_csv(cal_file)
        if not {'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'}.issubset(cal_df.columns):
            st.error("CSV must contain 'Charge', 'CCS', 'CCS Std.Dev.', and 'Intensity'.")
            return

        cal_df = cal_df[cal_df['CCS Std.Dev.'] < 0.5 * cal_df['CCS']].copy()
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df.dropna(inplace=True)

        all_charges = sorted(cal_df["Charge"].unique())
        selected_charges = st.multiselect("Select charge states to include", all_charges, default=all_charges)
        cal_df = cal_df[cal_df["Charge"].isin(selected_charges)]

        # Scaling
        scale_factors = {}
        for z in selected_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.98
            mz_max = mz * 1.02
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[z] = intensity_sum

        cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        x_common = np.linspace(x_min, x_max, 1000)
        all_fits = []
        fit_results = []

        for z in selected_charges:
            st.subheader(f"Charge State {z}")
            charge_df = cal_df[cal_df['Charge'] == z]

            if charge_df.empty:
                continue

            interp_func = interp1d(charge_df['CCS'], charge_df['Scaled Intensity'], bounds_error=False, fill_value=0)
            y_interp = interp_func(x_common)

            num_peaks = st.number_input(f"Number of Gaussians for {z}", min_value=1, max_value=5, value=1, key=f"num_peaks_{z}")
            initial_params = []

            for i in range(num_peaks):
                amp = st.text_input(f"Amp {i+1} (z={z})", value="1000.0", key=f"amp_{z}_{i}")
                mean = st.text_input(f"Mean {i+1} (z={z})", value="1500.0", key=f"mean_{z}_{i}")
                std = st.text_input(f"Std {i+1} (z={z})", value="50.0", key=f"std_{z}_{i}")
                try:
                    initial_params.extend([float(amp), float(mean), float(std)])
                except ValueError:
                    st.warning("Invalid Gaussian parameter. Please enter numbers.")

            bounds_lower = [0.95 * p for p in initial_params]
            bounds_upper = [1.05 * p for p in initial_params]

            col1, col2 = st.columns([2, 1])
            with col1:
                fig, ax = plt.subplots()
                ax.plot(x_common, y_interp, label="Raw Data", linestyle="--")
                y_fit = multi_gaussian(x_common, *initial_params)
                ax.plot(x_common, y_fit, label="Initial Fit")
                ax.set_title(f"Gaussian Fit for Charge {z}")
                ax.set_xlabel("CCS")
                ax.set_ylabel("Scaled Intensity")
                ax.legend()
                st.pyplot(fig)

            with col2:
                if st.button(f"Optimize Fit (z={z})"):
                    fitted_params = fit_gaussians(x_common, y_interp, initial_params, bounds=(bounds_lower, bounds_upper))
                    y_fit = multi_gaussian(x_common, *fitted_params)

                    fig_opt, ax_opt = plt.subplots()
                    ax_opt.plot(x_common, y_interp, label="Raw Data", linestyle="--")
                    ax_opt.plot(x_common, y_fit, label="Optimized Fit")
                    ax_opt.set_title(f"Optimized Gaussian Fit for Charge {z}")
                    ax_opt.set_xlabel("CCS")
                    ax_opt.set_ylabel("Scaled Intensity")
                    ax_opt.legend()
                    st.pyplot(fig_opt)

                    for i in range(0, len(fitted_params), 3):
                        fit_results.append({
                            "Charge": z,
                            "Peak": i // 3 + 1,
                            "Amplitude": fitted_params[i],
                            "Mean": fitted_params[i + 1],
                            "StdDev": fitted_params[i + 2]
                        })

                    all_fits.append(y_fit)
                else:
                    all_fits.append(y_fit)

        if all_fits:
            st.markdown("### Summed Trace and Fit")
            summed_trace = np.sum(all_fits, axis=0)
            fig_sum, ax_sum = plt.subplots()
            ax_sum.plot(x_common, summed_trace, label="Summed Fit")
            ax_sum.set_title("Summed Gaussian Fit")
            ax_sum.set_xlabel("CCS")
            ax_sum.set_ylabel("Intensity")
            ax_sum.legend()
            st.pyplot(fig_sum)

        if fit_results:
            results_df = pd.DataFrame(fit_results)
            st.dataframe(results_df)
            csv = results_df.to_csv(index=False)
            download_csv_button(csv, "gaussian_fit_parameters.csv", "Download Fit Parameters CSV")

