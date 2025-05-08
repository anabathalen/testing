import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO


def plot_and_scale_page():
    st.title("Plot and Scale Calibrated Data")

    # === File Upload ===
    cal_file = st.file_uploader("Upload a calibrated CSV file for a protein", type="csv")
    ms_file = st.file_uploader("Upload the mass spectrum TXT file (no headers)", type="txt")
    protein_mass = st.number_input("Enter the protein mass (Da)", min_value=0.0, step=1.0)

    if cal_file and ms_file and protein_mass > 0:
        # === Load and Validate Data ===
        cal_df = pd.read_csv(cal_file)
        if not {'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'}.issubset(cal_df.columns):
            st.error("The CSV file must contain 'Charge', 'CCS', 'CCS Std.Dev.', and 'Intensity' columns.")
            return

        cal_df = cal_df[cal_df['CCS Std.Dev.'] < 0.5 * cal_df['CCS']].copy()
        ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
        ms_df.dropna(inplace=True)

        PROTON_MASS = 1.007276
        unique_charges = sorted(cal_df["Charge"].unique())

        # === Calculate Scale Factors ===
        scale_factors = {}
        for z in unique_charges:
            mz = (protein_mass + z * PROTON_MASS) / z
            mz_min = mz * 0.99
            mz_max = mz * 1.01
            intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
            scale_factors[z] = intensity_sum

        cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
        cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]

        # === Let user choose scaled vs unscaled ===
        use_scaled = st.radio("Use Scaled or Unscaled Intensities?", ["Scaled", "Unscaled"], index=0)
        intensity_column = "Scaled Intensity" if use_scaled == "Scaled" else "Intensity"

        st.subheader("Data Preview")
        st.dataframe(cal_df)

        csv_output = cal_df.to_csv(index=False).encode("utf-8")
        st.download_button("Download CSV", data=csv_output, file_name="calibrated_data.csv", mime="text/csv")

        # === Plot Controls ===
        st.subheader("Plot Options")
        palette_choice = st.selectbox("Choose a color palette", list(sns.palettes.SEABORN_PALETTES.keys()))
        fig_width = st.slider("Figure width", 4, 20, 10)
        fig_height = st.slider("Figure height", 3, 15, 6)
        fig_dpi = st.slider("Figure DPI", 72, 300, 100)
        font_size = st.slider("Font size", 8, 24, 12)
        plot_mode = st.radio("Display Mode", ["Summed", "Stacked"])

        # === Additional Plot Settings ===
        selected_charges = st.multiselect(
            "Select charge states to display",
            options=unique_charges,
            default=unique_charges,
            format_func=lambda x: f"{x}+"
        )

        col1, col2 = st.columns(2)
        with col1:
            x_min = st.number_input("CCS x-axis min (Å²)", value=float(np.floor(cal_df["CCS"].min())))
        with col2:
            x_max = st.number_input("CCS x-axis max (Å²)", value=float(np.ceil(cal_df["CCS"].max())))

        ref_line_value = st.number_input("Reference CCS value (Å²)", value=0.0)
        ref_line_label = st.text_input("Reference label (e.g., AlphaFold2)", value="AlphaFold2")

        palette = sns.color_palette(palette_choice, n_colors=len(unique_charges))

        # === CCS Plot ===
        st.subheader("Scaled Intensity vs CCS")
        ccs_grid = np.arange(x_min, x_max + 1, 1.0)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)

        interpolated_traces = []
        filtered_df = cal_df[cal_df["Charge"].isin(selected_charges)]
        charges_to_plot = sorted(filtered_df["Charge"].unique())

        if plot_mode == "Summed":
            for i, (charge, group) in enumerate(filtered_df.groupby("Charge")):
                group_sorted = group.sort_values("CCS")
                interp = np.interp(ccs_grid, group_sorted["CCS"], group_sorted[intensity_column], left=0, right=0)
                interpolated_traces.append(interp)
                ax.plot(ccs_grid, interp, color=palette[i], label=f"{int(charge)}+")
                ax.fill_between(ccs_grid, 0, interp, color=palette[i], alpha=0.3)

            total_trace = np.sum(interpolated_traces, axis=0)
            ax.plot(ccs_grid, total_trace, color="black", linewidth=2.0, label="Total (Interpolated)")
            ax.legend(fontsize=font_size)

        elif plot_mode == "Stacked":
            vertical_offset = 1.2
            for i, (charge, group) in enumerate(filtered_df.groupby("Charge")):
                group_sorted = group.sort_values("CCS")
                interp = np.interp(ccs_grid, group_sorted["CCS"], group_sorted[intensity_column], left=0, right=0)
                norm_interp = interp / interp.max() if interp.max() > 0 else interp
                offset_interp = norm_interp + i * vertical_offset
                ax.plot(ccs_grid, offset_interp, color=palette[i])
                ax.fill_between(ccs_grid, i * vertical_offset, offset_interp, color=palette[i], alpha=0.3)
                # Shifted label
                label_x = ccs_grid[0] + 20  # right shift
                label_y = i * vertical_offset + 0.1  # slight raise
                ax.text(label_x, label_y, f"{int(charge)}+", fontsize=font_size, color=palette[i], va="bottom")

        # Reference Line
        if ref_line_value > 0:
            ax.axvline(ref_line_value, color='black', linestyle='--', linewidth=1.2)
            ax.text(ref_line_value - 5, ax.get_ylim()[1], ref_line_label,
                    color='black', fontsize=font_size, rotation=90,
                    verticalalignment='top', horizontalalignment='right')

        ax.set_xlim(x_min, x_max)
        ax.set_xlabel("CCS (Å²)", fontsize=font_size)
        ax.set_yticks([])
        ax.set_ylabel("")
        ax.set_title("Scaled Intensity vs CCS", fontsize=font_size)
        ax.grid(False)

        for label in ax.get_xticklabels():
            label.set_fontsize(font_size)
        for spine in ax.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)

        st.pyplot(fig)

        # Download PNG
        fig_buffer = BytesIO()
        fig.savefig(fig_buffer, format='png', dpi=fig_dpi, bbox_inches='tight')
        fig_buffer.seek(0)
        st.download_button("Download CCS Plot as PNG", data=fig_buffer, file_name="ccs_plot.png", mime="image/png")

