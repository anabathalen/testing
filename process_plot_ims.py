import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO

def plot_and_scale_page():
    st.title("Plot and Scale Calibrated Data")

    cal_file = st.file_uploader("Upload a calibrated CSV file for a protein", type="csv")
    ms_file = st.file_uploader("Upload the mass spectrum TXT file (no headers)", type="txt")
    protein_mass = st.number_input("Enter the protein mass (Da)", min_value=0.0, step=1.0)

    if cal_file and ms_file and protein_mass > 0:
        cal_df = pd.read_csv(cal_file)
        if not {'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'}.issubset(cal_df.columns):
            st.error("The CSV file must contain 'Charge', 'CCS', 'CCS Std.Dev.', and 'Intensity' columns.")
            return

        cal_df = cal_df[cal_df['CCS Std.Dev.'] < 0.5 * cal_df['CCS']].copy()
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df.dropna(inplace=True)

        PROTON_MASS = 1.007276
        all_charges = sorted(cal_df["Charge"].unique())
        selected_charges = st.multiselect("Select charge states to include", all_charges, default=all_charges)

        cal_df = cal_df[cal_df["Charge"].isin(selected_charges)]

        scale_factors = {}
        for z in selected_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[z] = intensity_sum

        cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        st.subheader("Scaled Calibrated Data")
        st.dataframe(cal_df)

        csv_output = cal_df.to_csv(index=False).encode("utf-8")
        st.download_button("Download Scaled CSV", data=csv_output, file_name="scaled_calibrated_data.csv", mime="text/csv", key="csv_download")

        st.subheader("Plot Options")
        palette_choice = st.selectbox("Choose a color palette", list(sns.palettes.SEABORN_PALETTES.keys()))
        fig_width = st.slider("Figure width", min_value=2, max_value=20, value=6)
        fig_height = st.slider("Figure height", min_value=2, max_value=20, value=4)
        fig_dpi = st.slider("Figure DPI", min_value=100, max_value=1000, value=300)
        font_size = st.slider("Font size", min_value=5, max_value=24, value=12)
        line_thickness = st.slider("Line thickness", min_value=0.1, max_value=5.0, value=1, step=0.1)
        plot_mode = st.radio("Display Mode", ["Summed", "Stacked"])
        use_scaled = st.radio("Use Scaled or Unscaled Intensities?", ["Scaled", "Unscaled"]) == "Scaled"

        ccs_min_input = st.number_input("CCS x-axis min", value=float(np.floor(cal_df["CCS"].min())))
        ccs_max_input = st.number_input("CCS x-axis max", value=float(np.ceil(cal_df["CCS"].max())))
        ccs_grid = np.arange(ccs_min_input, ccs_max_input + 1, 1.0)

        ccs_label_value = st.number_input("Optional CCS label position (leave blank if unused)", value=0.0, step=1.0, format="%.1f")

        palette = sns.color_palette(palette_choice, n_colors=len(selected_charges))

        # === Mass Spectrum Plot ===
        st.subheader("Mass Spectrum with Charge State Integration Regions")
        fig1, ax1 = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)
        ax1.plot(ms_df["m/z"], ms_df["Intensity"], color="gray", label="Mass Spectrum")

        for i, z in enumerate(selected_charges):
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.98
            mz_max = mz * 1.02
            region = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]
            ax1.fill_between(region["m/z"], region["Intensity"], color=palette[i], alpha=0.5, label=f"{z}+")

        ax1.set_xlabel("m/z", fontsize=font_size)
        ax1.set_ylabel("")
        ax1.set_yticks([])
        ax1.set_title("Mass Spectrum with Integration Windows", fontsize=font_size)
        ax1.legend(fontsize=font_size, frameon=False)

        for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            label.set_fontsize(font_size)
        for spine in ax1.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)

        st.pyplot(fig1)

        # === CCS Plot ===
        st.subheader("Scaled Intensity vs CCS")
        fig2, ax2 = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)
        interpolated_traces = []
        max_y_value = 0

        if plot_mode == "Summed":
            for i, (charge, group) in enumerate(cal_df.groupby("Charge")):
                group_sorted = group.sort_values("CCS")
                y_values = group_sorted["Scaled Intensity"] if use_scaled else group_sorted["Intensity"]
                interp = np.interp(ccs_grid, group_sorted["CCS"], y_values, left=0, right=0)
                interpolated_traces.append(interp)
                ax2.plot(ccs_grid, interp, color=palette[i], label=f"{int(charge)}+", linewidth=line_thickness)
                ax2.fill_between(ccs_grid, 0, interp, color=palette[i], alpha=0.3)
                max_y_value = max(max_y_value, interp.max())

            total_trace = np.sum(interpolated_traces, axis=0)
            ax2.plot(ccs_grid, total_trace, color="black", linewidth=line_thickness + 0.5, label="Total (Interpolated)")
            ax2.legend(fontsize=font_size, frameon=False)

        elif plot_mode == "Stacked":
            offset_unit = 1.0 / len(selected_charges)
            base_max = 0

            # Calculate base max intensity
            for i, (charge, group) in enumerate(cal_df.groupby("Charge")):
                y_values = group["Scaled Intensity"] if use_scaled else group["Intensity"]
                base_max = max(base_max, y_values.max())

            for i, (charge, group) in enumerate(cal_df.groupby("Charge")):
                group_sorted = group.sort_values("CCS")
                y_values = group_sorted["Scaled Intensity"] if use_scaled else group_sorted["Intensity"]
                interp = np.interp(ccs_grid, group_sorted["CCS"], y_values, left=0, right=0)

                if not use_scaled:
                    interp = interp / interp.max() if interp.max() > 0 else interp

                offset = i * offset_unit * base_max
                offset_interp = interp + offset

                ax2.plot(ccs_grid, offset_interp, color=palette[i], linewidth=line_thickness)
                ax2.fill_between(ccs_grid, offset, offset_interp, color=palette[i], alpha=0.3)

                label_x = ccs_grid[0] - (ccs_max_input - ccs_min_input) * 0.1
                label_y = offset + offset_interp.max() * 0.1
                ax2.text(label_x, label_y, f"{int(charge)}+", fontsize=font_size,
                         verticalalignment="bottom", horizontalalignment="left", color=palette[i])

            max_y_value = (len(selected_charges)) * offset_unit * base_max

        if ccs_label_value > 0:
            ax2.axvline(ccs_label_value, color="black", linewidth=1.0, linestyle="--")
            ax2.text(ccs_label_value - 0.5, max_y_value * 0.95, "txt", color="black", fontsize=font_size,
                     horizontalalignment="right", verticalalignment="top")

        ax2.set_xlim([ccs_min_input, ccs_max_input])
        ax2.set_xlabel("CCS (Å²)", fontsize=font_size)
        ax2.set_ylabel("")
        ax2.set_yticks([])
        ax2.grid(False)

        for label in ax2.get_xticklabels():
            label.set_fontsize(font_size)
        for spine in ax2.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)

        st.pyplot(fig2)

        fig_buffer = BytesIO()
        fig2.savefig(fig_buffer, format='png', dpi=fig_dpi, bbox_inches='tight')
        fig_buffer.seek(0)
        st.download_button("Download CCS Plot as PNG", data=fig_buffer, file_name="ccs_plot.png", mime="image/png", key="ccs_download")

