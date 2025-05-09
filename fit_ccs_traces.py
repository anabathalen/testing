import streamlit as st
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO

# Gaussian and sum of Gaussians functions for fitting
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

    # Load the CSV containing CCS and Intensity data
    df = pd.read_csv(uploaded_file)
    if not {"CCS", "Intensity"}.issubset(df.columns):
        st.error("CSV must contain 'CCS' and 'Intensity' columns.")
        return

    # Scaling and interpolation setup
    st.markdown("### Scaling and Interpolation Settings")
    protein_mass = st.number_input("Enter the protein mass (Da)", min_value=0.0, step=1.0)
    cal_file = st.file_uploader("Upload a calibrated CSV file", type="csv")
    ms_file = st.file_uploader("Upload the mass spectrum TXT file", type="txt")

    if cal_file and ms_file and protein_mass > 0:
        cal_df = pd.read_csv(cal_file)
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])

        PROTON_MASS = 1.007276
        all_charges = sorted(cal_df["Charge"].unique())
        selected_charges = st.multiselect("Select charge states to include", all_charges, default=all_charges)

        cal_df = cal_df[cal_df["Charge"].isin(selected_charges)]
        scale_factors = {}
        for z in selected_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.98
            mz_max = mz * 1.02
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[z] = intensity_sum

        cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        st.subheader("Scaled Calibrated Data")
        st.dataframe(cal_df)

        # Gaussian fitting interface
        st.subheader("Gaussian Fitting Interface")
        fit_mode = st.radio("Choose fitting mode", ["Fit Summed Trace", "Fit Each Charge Separately"])
        x_min, x_max = st.slider("Select x-axis range (CCS)", float(df["CCS"].min()), float(df["CCS"].max()), (float(df["CCS"].min()), float(df["CCS"].max())))
        x_fit = np.linspace(x_min, x_max, 1000)

        if fit_mode == "Fit Each Charge Separately" and "Charge" not in df.columns:
            st.warning("No 'Charge' column found. Defaulting to summed trace.")
            fit_mode = "Fit Summed Trace"

        # Grouping by charge state
        data_groups = [("Summed", cal_df)] if fit_mode == "Fit Summed Trace" else list(cal_df.groupby("Charge"))

        all_fit_params = []
        for group_label, group_df in data_groups:
            subset = group_df[(group_df["CCS"] >= x_min) & (group_df["CCS"] <= x_max)]
            x = subset["CCS"].values
            y = subset["Scaled Intensity"].values
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
                    all_fit_params.append(popt)
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

            # Saving the Gaussian fit CSV
            fit_df = pd.DataFrame(np.round(popt, 3), columns=["Mean", "Std Dev", "Amplitude"])
            fit_csv = fit_df.to_csv(index=False).encode("utf-8")
            st.download_button("Download Gaussian Fit CSV", data=fit_csv, file_name=f"{group_label}_gaussian_fit.csv", mime="text/csv")

        # Saving fitted CCS vs Intensity data CSV
        fit_data = pd.DataFrame({"CCS": x_fit, "Fitted Intensity": y_fit})
        fit_data_csv = fit_data.to_csv(index=False).encode("utf-8")
        st.download_button("Download CCS vs Fitted Intensity CSV", data=fit_data_csv, file_name="ccs_vs_fitted_intensity.csv", mime="text/csv")
