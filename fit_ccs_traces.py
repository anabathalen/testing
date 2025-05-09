import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import norm
from io import StringIO

def fit_ccs_traces_page():
    st.header("Fit CCS Traces")

    cal_file = st.file_uploader("Upload a calibrated CSV file for a protein", type="csv")
    ms_file = st.file_uploader("Upload the mass spectrum TXT file (no headers)", type="txt")
    protein_mass = st.number_input("Enter the protein mass (Da)", min_value=0.0, step=1.0)

    if cal_file and ms_file and protein_mass > 0:
        cal_df = pd.read_csv(cal_file)
        if not {'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'}.issubset(cal_df.columns):
            st.error("CSV must contain 'Charge', 'CCS', 'CCS Std.Dev.', and 'Intensity'.")
            return

        cal_df = cal_df[cal_df['CCS Std.Dev.'] < 0.5 * cal_df['CCS']].copy()
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df.dropna(inplace=True)

        PROTON_MASS = 1.007276
        all_charges = sorted(cal_df["Charge"].unique())
        selected_charges = st.multiselect("Select charge states to include", all_charges, default=all_charges)
        cal_df = cal_df[cal_df["Charge"].isin(selected_charges)]

        # Scaling logic
        scale_factors = {}
        for z in selected_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.98
            mz_max = mz * 1.02
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[z] = intensity_sum

        cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        # Fit mode
        fit_mode = st.radio("Fit mode", ["Sum trace then fit", "Fit each charge then sum"])

        all_traces = []
        gaussian_params = []
        common_ccs = np.linspace(cal_df["CCS"].min(), cal_df["CCS"].max(), 1000)

        fig, ax = plt.subplots()
        ax.set_title("CCS Traces and Gaussian Fits")
        ax.set_xlabel("CCS")
        ax.set_ylabel("Scaled Intensity")

        for z in selected_charges:
            sub = cal_df[cal_df["Charge"] == z]
            if sub.empty:
                continue
            interp_func = interp1d(sub["CCS"], sub["Scaled Intensity"], bounds_error=False, fill_value=0)
            y_interp = interp_func(common_ccs)
            all_traces.append(y_interp)

            if fit_mode == "Fit each charge then sum":
                st.subheader(f"Charge {z} Gaussian Fit")
                num_peaks = st.slider(f"Number of Gaussians for charge {z}", 1, 5, 1)
                y_fit = np.zeros_like(common_ccs)
                for i in range(num_peaks):
                    mean = st.slider(f"Mean {i+1} (z={z})", float(common_ccs.min()), float(common_ccs.max()), float(sub["CCS"].mean()))
                    std = st.slider(f"Std. Dev. {i+1} (z={z})", 1.0, 30.0, 10.0)
                    amp = st.slider(f"Amplitude {i+1} (z={z})", 0.0, float(y_interp.max()) * 2, float(y_interp.max()))
                    y_fit += amp * norm.pdf(common_ccs, mean, std)
                    gaussian_params.append({"Charge": z, "Peak": i+1, "Mean": mean, "StdDev": std, "Amplitude": amp})
                ax.plot(common_ccs, y_interp, label=f"Raw z={z}", alpha=0.3)
                ax.plot(common_ccs, y_fit, label=f"Fit z={z}", linestyle="--")

        if fit_mode == "Sum trace then fit":
            summed_trace = np.sum(all_traces, axis=0)
            st.subheader("Summed Gaussian Fit")
            num_peaks = st.slider("Number of Gaussians for summed trace", 1, 5, 1)
            y_fit = np.zeros_like(common_ccs)
            for i in range(num_peaks):
                mean = st.slider(f"Mean {i+1}", float(common_ccs.min()), float(common_ccs.max()), float(common_ccs.mean()))
                std = st.slider(f"Std. Dev. {i+1}", 1.0, 30.0, 10.0)
                amp = st.slider(f"Amplitude {i+1}", 0.0, float(summed_trace.max()) * 2, float(summed_trace.max()))
                y_fit += amp * norm.pdf(common_ccs, mean, std)
                gaussian_params.append({"Charge": "summed", "Peak": i+1, "Mean": mean, "StdDev": std, "Amplitude": amp})
            ax.plot(common_ccs, summed_trace, label="Summed Trace", alpha=0.3)
            ax.plot(common_ccs, y_fit, label="Summed Fit", linestyle="-")

        ax.legend()
        st.pyplot(fig)

        # Download CSVs
        if gaussian_params:
            params_df = pd.DataFrame(gaussian_params)
            csv_gauss = params_df.to_csv(index=False).encode("utf-8")
            st.download_button("Download Gaussian Parameters CSV", csv_gauss, "gaussian_parameters.csv", "text/csv")

        if fit_mode == "Sum trace then fit":
            trace_df = pd.DataFrame({"CCS": common_ccs, "Summed Intensity": y_fit})
            csv_trace = trace_df.to_csv(index=False).encode("utf-8")
            st.download_button("Download Summed Fitted Trace CSV", csv_trace, "summed_fitted_trace.csv", "text/csv")
        else:
            for z in selected_charges:
                trace_df = pd.DataFrame({"CCS": common_ccs})
                trace_df[f"Charge {z}"] = interp1d(cal_df[cal_df["Charge"] == z]["CCS"],
                                                   cal_df[cal_df["Charge"] == z]["Scaled Intensity"],
                                                   bounds_error=False, fill_value=0)(common_ccs)
                csv_trace = trace_df.to_csv(index=False).encode("utf-8")
                st.download_button(f"Download Fitted Trace CSV for Charge {z}", csv_trace, f"fitted_trace_charge_{z}.csv", "text/csv")
