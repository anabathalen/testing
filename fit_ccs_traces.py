import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from io import BytesIO

PROTON_MASS = 1.007276

def gaussian(x, mu, sigma, amplitude):
    return amplitude * np.exp(-(x - mu)**2 / (2 * sigma**2))

def generate_gaussian_sum(x, params):
    return sum(gaussian(x, mu, sigma, amp) for mu, sigma, amp in params)

def fit_ccs_traces_page():
    st.header("Gaussian Fitting of CCS Traces")

    cal_file = st.file_uploader("Upload a calibrated CSV file for a protein", type="csv")
    ms_file = st.file_uploader("Upload the mass spectrum TXT file (no headers)", type="txt")
    protein_mass = st.number_input("Enter the protein mass (Da)", min_value=0.0, step=1.0)

    x_min = st.number_input("X-axis minimum (CCS)", value=0.0)
    x_max = st.number_input("X-axis maximum (CCS)", value=3000.0)

    fit_mode = st.radio("Fit mode", ["Fit individual charge states before summing", "Fit the summed trace"])

    if cal_file and ms_file and protein_mass > 0 and x_max > x_min:
        cal_df = pd.read_csv(cal_file)
        if not {'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'}.issubset(cal_df.columns):
            st.error("CSV file must contain 'Charge', 'CCS', 'CCS Std.Dev.', and 'Intensity'.")
            return

        cal_df = cal_df[cal_df['CCS Std.Dev.'] < 0.5 * cal_df['CCS']].copy()
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df.dropna(inplace=True)

        all_charges = sorted(cal_df["Charge"].unique())
        selected_charges = st.multiselect("Select charge states to include", all_charges, default=all_charges)

        cal_df = cal_df[cal_df["Charge"].isin(selected_charges)]

        # Calculate scale factors
        scale_factors = {}
        for z in selected_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.98
            mz_max = mz * 1.02
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[z] = intensity_sum

        cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        x_fit = np.linspace(x_min, x_max, 1000)

        output_dfs = []
        output_labels = []

        if fit_mode == "Fit individual charge states before summing":
            summed_trace = np.zeros_like(x_fit)

            for z in selected_charges:
                st.subheader(f"Charge State {z}")
                n_peaks = st.number_input(f"Number of Gaussians for {z}", min_value=1, max_value=5, value=1, key=f"peaks_{z}")
                params = []
                for i in range(n_peaks):
                    mu = st.text_input(f"Mean of Gaussian {i+1} (z={z})", value="1000", key=f"mu_{z}_{i}")
                    sigma = st.text_input(f"Sigma of Gaussian {i+1} (z={z})", value="20", key=f"sigma_{z}_{i}")
                    amp = st.text_input(f"Amplitude of Gaussian {i+1} (z={z})", value="1", key=f"amp_{z}_{i}")
                    try:
                        params.append((float(mu), float(sigma), float(amp)))
                    except ValueError:
                        st.warning("Invalid parameter input. Using default.")
                        params.append((1000, 20, 1))

                y_fit = generate_gaussian_sum(x_fit, params)
                summed_trace += y_fit

                df = pd.DataFrame({"CCS": x_fit, f"Fitted Intensity z={z}": y_fit})
                output_dfs.append(df)
                output_labels.append(f"fit_z{z}.csv")

                fig, ax = plt.subplots()
                ax.plot(x_fit, y_fit, label=f"Charge {z} Fit")
                ax.set_xlim(x_min, x_max)
                ax.set_title(f"Gaussian Fit for Charge State {z}")
                ax.set_xlabel("CCS")
                ax.set_ylabel("Fitted Intensity")
                st.pyplot(fig)

            # Show summed trace
            fig2, ax2 = plt.subplots()
            ax2.plot(x_fit, summed_trace, label="Summed Fit")
            ax2.set_xlim(x_min, x_max)
            ax2.set_title("Summed Gaussian Fit")
            ax2.set_xlabel("CCS")
            ax2.set_ylabel("Fitted Intensity")
            st.pyplot(fig2)

            df_sum = pd.DataFrame({"CCS": x_fit, "Summed Fitted Intensity": summed_trace})
            output_dfs.append(df_sum)
            output_labels.append("summed_fit.csv")

        else:
            st.subheader("Summed Trace Fit")
            trace = np.zeros_like(x_fit)
            for z in selected_charges:
                sub_df = cal_df[cal_df['Charge'] == z]
                trace_z = np.interp(x_fit, sub_df['CCS'], sub_df['Scaled Intensity'], left=0, right=0)
                trace += trace_z

            n_peaks = st.number_input("Number of Gaussians for summed fit", min_value=1, max_value=5, value=1, key="summed")
            params = []
            for i in range(n_peaks):
                mu = st.text_input(f"Mean of Gaussian {i+1} (summed)", value="1000", key=f"mu_sum_{i}")
                sigma = st.text_input(f"Sigma of Gaussian {i+1} (summed)", value="20", key=f"sigma_sum_{i}")
                amp = st.text_input(f"Amplitude of Gaussian {i+1} (summed)", value="1", key=f"amp_sum_{i}")
                try:
                    params.append((float(mu), float(sigma), float(amp)))
                except ValueError:
                    st.warning("Invalid parameter input. Using default.")
                    params.append((1000, 20, 1))

            y_fit = generate_gaussian_sum(x_fit, params)

            fig, ax = plt.subplots()
            ax.plot(x_fit, trace, label="Observed", linestyle=":")
            ax.plot(x_fit, y_fit, label="Fitted", linestyle="--")
            ax.set_xlim(x_min, x_max)
            ax.set_title("Summed Gaussian Fit")
            ax.set_xlabel("CCS")
            ax.set_ylabel("Intensity")
            ax.legend()
            st.pyplot(fig)

            df = pd.DataFrame({"CCS": x_fit, "Observed": trace, "Fitted": y_fit})
            output_dfs.append(df)
            output_labels.append("summed_fit.csv")

        for df, label in zip(output_dfs, output_labels):
            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button(f"Download {label}", data=csv, file_name=label, mime="text/csv")

