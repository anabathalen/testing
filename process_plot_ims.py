import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO


def plot_and_scale_page():
    st.title("Plot and Scale Calibrated Data")

    # === File Inputs ===
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
        unique_charges = sorted(cal_df["Charge"].unique())

        # === Integration for scaling ===
        scale_factors = {}
        for z in unique_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min, mz_max = mz * 0.99, mz * 1.01
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[z] = intensity_sum

        cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        # === Scaled vs. Unscaled Toggle ===
        use_scaled = st.radio("Use scaled or unscaled intensities?", ["Scaled", "Unscaled"]) == "Scaled"
        cal_df["Plot Intensity"] = cal_df["Scaled Intensity"] if use_scaled else cal_df["Intensity"]

        st.subheader("Processed Calibrated Data")
        st.dataframe(cal_df)

        csv_output = cal_df.to_csv(index=False).encode("utf-8")
        st.download_button("Download CSV", data=csv_output, file_name="processed_calibrated_data.csv", mime="text/csv")

        # === Plot Controls ===
        st.subheader("Plot Options")
        palette_choice = st.selectbox("Color palette", list(sns.palettes.SEABORN_PALETTES.keys()))
        fig_width = st.slider("Figure width", min_value=3, max_value=20, value=10)
        fig_height = st.slider("Figure height", min_value=3, max_value=15, value=6)
        fig_dpi = st.slider("Figure DPI (resolution)", min_value=100, max_value=1000, value=300)
        font_size = st.slider("Font size", min_value=5, max_value=24, value=12)
        line_thickness = st.slider("Line thickness", min_value=0.5, max_value=5.0, value=2.0, step=0.1)
        plot_mode = st.radio("Display Mode", ["Summed", "Stacked"])

        # === CCS range (manual override allowed) ===
        st.subheader("Set CCS X-Axis Range")
        user_ccs_min = st.number_input("Minimum CCS (Å²)", value=float(np.floor(cal_df["CCS"].min())))
        user_ccs_max = st.number_input("Maximum CCS (Å²)", value=float(np.ceil(cal_df["CCS"].max())))
        if user_ccs_min >= user_ccs_max:
            st.warning("Minimum CCS must be less than maximum CCS.")
            return
        ccs_range = (user_ccs_min, user_ccs_max)
        ccs_grid = np.arange(user_ccs_min, user_ccs_max + 1, 1.0)

        palette = sns.color_palette(palette_choice, n_colors=len(unique_charges))

        # === Mass Spectrum Plot ===
        st.subheader("Mass Spectrum with Integration Regions")
        fig1, ax1 = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)
        ax1.plot(ms_df["m/z"], ms_df["Intensity"], color="gray", label="Mass Spectrum")

        for i, z in enumerate(unique_charges):
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min, mz_max = mz * 0.99, mz * 1.01
            region = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]
            ax1.fill_between(region["m/z"], region["Intensity"], color=palette[i], alpha=0.5, label=f"{z}+")

        ax1.set_xlabel("m/z", fontsize=font_size)
        ax1.set_ylabel("Intensity", fontsize=font_size)
        ax1.set_title("Mass Spectrum with Integration Windows", fontsize=font_size)
        ax1.legend(fontsize=font_size)
        for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            label.set_fontsize(font_size)
        for spine in ax1.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)
        st.pyplot(fig1)

        # === CCS Plot ===
        st.subheader("CCS Plot")

        fig2, ax2 = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)

        if plot_mode == "Summed":
            interpolated_traces = []
            for i, (charge, group) in enumerate(cal_df.groupby("Charge")):
                group_sorted = group.sort_values("CCS")
                interp = np.interp(ccs_grid, group_sorted["CCS"], group_sorted["Plot Intensity"], left=0, right=0)
                interpolated_traces.append(interp)
                ax2.plot(ccs_grid, interp, color=palette[i], linewidth=line_thickness, label=f"{int(charge)}+")
                ax2.fill_between(ccs_grid, 0, interp, color=palette[i], alpha=0.3)

            total_trace = np.sum(interpolated_traces, axis=0)
            ax2.plot(ccs_grid, total_trace, color="black", linewidth=line_thickness, label="Summed")
            ax2.legend(fontsize=font_size, frameon=False)

        elif plot_mode == "Stacked":
            vertical_offset = 1.2
            for i, (charge, group) in enumerate(cal_df.groupby("Charge")):
                group_sorted = group.sort_values("CCS")
                interp = np.interp(ccs_grid, group_sorted["CCS"], group_sorted["Plot Intensity"], left=0, right=0)
                if not use_scaled and interp.max() > 0:
                    interp = interp / interp.max()
                offset_interp = interp + i * vertical_offset
                ax2.plot(ccs_grid, offset_interp, color="black", linewidth=line_thickness)
                ax2.fill_between(ccs_grid, i * vertical_offset, offset_interp, color=palette[i], alpha=0.4)

                label_x = ccs_grid[0] + 1
                label_y = i * vertical_offset + 0.1
                ax2.text(label_x, label_y, f"{int(charge)}+", fontsize=font_size,
                         verticalalignment="bottom", horizontalalignment="left", color=palette[i])

        ax2.set_xlim(ccs_range)
        ax2.set_xlabel("CCS (Å²)", fontsize=font_size)
        ax2.set_yticks([])
        ax2.set_ylabel("")
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
        st.download_button("Download CCS Plot as PNG", data=fig_buffer, file_name="ccs_plot.png", mime="image/png")
