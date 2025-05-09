import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from io import StringIO

PROTON_MASS = 1.007276

def gaussian(x, mean, stddev, amplitude):
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

def fit_ccs_traces_page():
    st.header("Fit CCS Traces")

    cal_file = st.file_uploader("Upload calibrated CCS CSV", type="csv")
    ms_file = st.file_uploader("Upload mass spectrum TXT file", type="txt")
    protein_mass = st.number_input("Enter protein mass (Da)", min_value=0.0, step=1.0)

    x_min = st.number_input("X-axis min (CCS)", value=1000.0)
    x_max = st.number_input("X-axis max (CCS)", value=3000.0)

    if cal_file and ms_file and protein_mass > 0:
        cal_df = pd.read_csv(cal_file)
        required_cols = {'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'}
        if not required_cols.issubset(cal_df.columns):
            st.error("CSV must contain: 'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'")
            return

        cal_df = cal_df[cal_df['CCS Std.Dev.'] < 0.5 * cal_df['CCS']].copy()
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df.dropna(inplace=True)

        all_charges = sorted(cal_df["Charge"].unique())
        selected_charges = st.multiselect("Select charge states to include", all_charges, default=all_charges)
        cal_df = cal_df[cal_df["Charge"].isin(selected_charges)]

        # Scale intensities from MS
        scale_factors = {}
        for z in selected_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.98
            mz_max = mz * 1.02
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[z] = intensity_sum

        cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        # Create shared x-axis
        x_fit = np.linspace(x_min, x_max, 1000)
        all_traces = []
        all_fit_dfs = []
        raw_trace_dfs = []

        st.markdown("### Gaussian Fit Parameters")
        for charge in selected_charges:
            st.subheader(f"Charge State {charge}")
            sub_df = cal_df[cal_df["Charge"] == charge]
            interp_y = np.interp(x_fit, sub_df["CCS"], sub_df["Scaled Intensity"])
            all_traces.append(interp_y)
            raw_trace_dfs.append(pd.DataFrame({"CCS": x_fit, "Raw Intensity": interp_y}))

            n_peaks = st.number_input(f"Number of Gaussians for charge {charge}", min_value=1, max_value=5, value=1, key=f"npeaks_{charge}")
            fit = np.zeros_like(x_fit)
            fit_params = []
            for i in range(n_peaks):
                mean = st.text_input(f"Peak {i+1} mean (charge {charge})", value=str(x_fit[np.argmax(interp_y)]), key=f"mean_{charge}_{i}")
                std = st.text_input(f"Peak {i+1} stddev (charge {charge})", value="20", key=f"std_{charge}_{i}")
                amp = st.text_input(f"Peak {i+1} amplitude (charge {charge})", value="1.0", key=f"amp_{charge}_{i}")
                try:
                    g = gaussian(x_fit, float(mean), float(std), float(amp))
                    fit += g
                    fit_params.append({"Mean": float(mean), "StdDev": float(std), "Amplitude": float(amp)})
                except ValueError:
                    st.error("Please enter valid numerical values for Gaussian parameters.")

            fig, ax = plt.subplots()
            ax.plot(x_fit, interp_y, label="Raw Data", linestyle=":")
            ax.plot(x_fit, fit, label="Fit")
            ax.set_title(f"Charge {charge} Fit")
            ax.set_xlabel("CCS")
            ax.set_ylabel("Intensity")
            ax.legend()
            st.pyplot(fig)

            fit_df = pd.DataFrame({"CCS": x_fit, "Fitted Intensity": fit})
            all_fit_dfs.append((charge, fit_df, pd.DataFrame(fit_params)))

        summed_trace = np.sum(all_traces, axis=0)
        st.markdown("### Summed Trace Fit")

        n_sum_peaks = st.number_input("Number of Gaussians for summed trace", min_value=1, max_value=5, value=1, key="sum_npeaks")
        summed_fit = np.zeros_like(x_fit)
        sum_params = []
        for i in range(n_sum_peaks):
            mean = st.text_input(f"Summed Peak {i+1} mean", value=str(x_fit[np.argmax(summed_trace)]), key=f"sum_mean_{i}")
            std = st.text_input(f"Summed Peak {i+1} stddev", value="20", key=f"sum_std_{i}")
            amp = st.text_input(f"Summed Peak {i+1} amplitude", value="1.0", key=f"sum_amp_{i}")
            try:
                g = gaussian(x_fit, float(mean), float(std), float(amp))
                summed_fit += g
                sum_params.append({"Mean": float(mean), "StdDev": float(std), "Amplitude": float(amp)})
            except ValueError:
                st.error("Invalid Gaussian parameter in summed fit.")

        fig2, ax2 = plt.subplots()
        ax2.plot(x_fit, summed_trace, label="Summed Raw Data", linestyle=":")
        ax2.plot(x_fit, summed_fit, label="Summed Fit")
        ax2.set_title("Summed Fit")
        ax2.set_xlabel("CCS")
        ax2.set_ylabel("Intensity")
        ax2.legend()
        st.pyplot(fig2)

        st.markdown("### Download Fitted Data")
        for charge, fit_df, params_df in all_fit_dfs:
            st.download_button(
                label=f"Download Charge {charge} Fitted Curve CSV",
                data=fit_df.to_csv(index=False).encode(),
                file_name=f"fit_charge_{charge}.csv",
                mime="text/csv"
            )
            st.download_button(
                label=f"Download Charge {charge} Gaussian Parameters CSV",
                data=params_df.to_csv(index=False).encode(),
                file_name=f"params_charge_{charge}.csv",
                mime="text/csv"
            )

        summed_df = pd.DataFrame({"CCS": x_fit, "Fitted Intensity": summed_fit})
        summed_param_df = pd.DataFrame(sum_params)

        st.download_button(
            label="Download Summed Fit CSV",
            data=summed_df.to_csv(index=False).encode(),
            file_name="summed_fit.csv",
            mime="text/csv"
        )
        st.download_button(
            label="Download Summed Gaussian Parameters CSV",
            data=summed_param_df.to_csv(index=False).encode(),
            file_name="summed_gaussians.csv",
            mime="text/csv"
        )

