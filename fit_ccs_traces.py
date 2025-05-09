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
    
    # Upload the CCS data
    uploaded_file = st.file_uploader("Upload CCS Data CSV", type=["csv"])
    if uploaded_file is None:
        return

    df = pd.read_csv(uploaded_file)
    if not {"CCS", "Intensity"}.issubset(df.columns):
        st.error("CSV must contain 'CCS' and 'Intensity' columns.")
        return

    # Upload the mass spectrum data
    ms_file = st.file_uploader("Upload Mass Spectrum Data CSV", type=["csv"])
    if ms_file is None:
        st.error("You need to upload a mass spectrum file to continue.")
        return
    
    ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
    ms_df.dropna(inplace=True)
    
    fit_mode = st.radio("Choose fitting mode", ["Fit Summed Trace", "Fit Each Charge Separately"])

    x_min, x_max = st.slider("Select x-axis range (CCS)", float(df["CCS"].min()), float(df["CCS"].max()), (float(df["CCS"].min()), float(df["CCS"].max())))

    x_fit = np.linspace(x_min, x_max, 1000)

    if fit_mode == "Fit Each Charge Separately" and "Charge" not in df.columns:
        st.warning("No 'Charge' column found. Defaulting to summed trace.")
        fit_mode = "Fit Summed Trace"

    data_groups = [("Summed", df)] if fit_mode == "Fit Summed Trace" else list(df.groupby("Charge"))

    scale_factors = {}

    all_traces = []  # For summed trace fitting

    for group_label, group_df in data_groups:
        subset = group_df[(group_df["CCS"] >= x_min) & (group_df["CCS"] <= x_max)]
        x = subset["CCS"].values
        y = subset["Intensity"].values
        y = y / y.max()  # Normalize the intensities

        if fit_mode == "Fit Each Charge Separately":
            charge = group_label
            PROTON_MASS = 1.007276
            protein_mass = st.number_input("Enter Protein Mass (Da)", min_value=0.0, step=0.1, value=10000.0)

            mz = (protein_mass + charge * PROTON_MASS) / charge
            mz_min = mz * 0.98
            mz_max = mz * 1.02

            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[charge] = intensity_sum

            # Apply the scaling factor
            scale_factor = scale_factors[charge]
            y_scaled = y * scale_factor

            st.markdown(f"### Charge {charge} (Scaled)")
            num_peaks = st.number_input(f"Number of Gaussians for Charge {charge}", min_value=1, max_value=10, value=2, key=f"npeaks_{charge}")
            peak_params = []

            for i in range(num_peaks):
                m = st.text_input(f"Mean of peak {i+1} for Charge {charge}", value=str(round(x.min() + (i+1)*(x.max()-x.min())/(num_peaks+1), 2)), key=f"mean_{i}_{charge}")
                s = st.text_input(f"Std dev of peak {i+1} for Charge {charge}", value="5.0", key=f"std_{i}_{charge}")
                a = st.text_input(f"Amplitude of peak {i+1} for Charge {charge}", value="1.0", key=f"amp_{i}_{charge}")
                peak_params += [float(m), float(s), float(a)]

            autofit = st.button(f"Auto-Fit Gaussians for Charge {charge}")

            y_fit = sum_of_gaussians(x_fit, *peak_params)

            if autofit:
                bounds_lower = [p * 0.95 for p in peak_params]
                bounds_upper = [p * 1.05 for p in peak_params]
                try:
                    popt, _ = curve_fit(sum_of_gaussians, x, y_scaled, p0=peak_params, bounds=(bounds_lower, bounds_upper))
                    st.success("Optimized fit successful.")
                    y_fit = sum_of_gaussians(x_fit, *popt)
                    st.text("Fitted parameters:\n" + str(np.round(popt, 3)))
                except Exception as e:
                    st.error(f"Fit failed: {e}")

            fig, ax = plt.subplots()
            ax.plot(x, y_scaled, label="Scaled Data")
            ax.plot(x_fit, y_fit, label="Fit", linestyle="--")
            ax.set_title(f"Gaussian Fit: Charge {charge}")
            ax.set_xlabel("CCS")
            ax.set_ylabel("Scaled Intensity")
            ax.legend()
            st.pyplot(fig)

            if fit_mode == "Fit Summed Trace":
                all_traces.append(y_fit)

    if fit_mode == "Fit Summed Trace":
        # Make sure summed_trace is the same length as x_fit
        summed_trace = np.sum(all_traces, axis=0)
        summed_trace = np.interp(x_fit, x, summed_trace)  # Interpolate to match x_fit
        st.markdown("### Summed Gaussian Fit")
        fig2, ax2 = plt.subplots()
        ax2.plot(x_fit, summed_trace, label="Summed Fit", linestyle="--")
        ax2.set_title("Summed Gaussian Fit")
        ax2.set_xlabel("CCS")
        ax2.set_ylabel("Intensity")
        ax2.legend()
        st.pyplot(fig2)

        # Save CSV with summed data
        summed_data = pd.DataFrame({"CCS": x_fit, "Fitted Intensity": summed_trace})
        csv_output = summed_data.to_csv(index=False).encode("utf-8")
        st.download_button("Download Summed Fitted Data", data=csv_output, file_name="summed_fitted_data.csv", mime="text/csv", key="summed_fitted_data")

