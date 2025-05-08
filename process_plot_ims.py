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
        unique_charges = sorted(cal_df["Charge"].unique())

        # === Scale Factors ===
        scale_factors = {}
        for z in unique_charges:
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

        # === Plotting Controls ===
        st.subheader("Plot Options")
        use_scaled = st.radio("Choose intensity type", ["Scaled", "Unscaled"])
        palette_choice = st.selectbox("Choose a color palette", list(sns.palettes.SEABORN_PALETTES.keys()))
        fig_width = st.slider("Figure width", min_value=2, max_value=20, value=10)
        fig_height = st.slider("Figure height", min_value=2, max_value=20, value=6)
        fig_dpi = st.slider("Figure DPI", min_value=100, max_value=1000, value=300)
        font_size = st.slider("Font size", min_value=8, max_value=24, value=12)
        line_thickness = st.slider("Line thickness", min_value=0.5, max_value=5.0, value=2.0, step=0.1)
        plot_mode = st.radio("Display Mode", ["Summed", "Stacked"])

        # === Charge Filter ===
        min_charge = int(cal_df["Charge"].min())
        max_charge = int(cal_df["Charge"].max())
        selected_charges = st.slider("Select charge state range", min_value=min_charge, max_value=max_charge, value=(min_charge, max_charge))
        filtered_df = cal_df[(cal_df["Charge"] >= selected_charges[0]) & (cal_df["Charge"] <= selected_charges[1])]
        filtered_charges = sorted(filtered_df["Charge"].unique())

        # === x-axis filter ===
        ccs_min_all = float(np.floor(cal_df["CCS"].min()))
        ccs_max_all = float(np.ceil(cal_df["CCS"].max()))
        ccs_range = st.slider("Select CCS x-axis range", min_value=ccs_min_all, max_value=ccs_max_all, value=(ccs_min_all, ccs_max_all))
        ccs_grid = np.arange(ccs_range[0], ccs_range[1] + 1, 1.0)

        intensity_column = "Scaled Intensity" if use_scaled == "Scaled" else "Intensity"
        palette = sns.color_palette(palette_choice, n_colors=len(filtered_charges))

        # === CCS Plot ===
        st.subheader("Scaled Intensity vs CCS")
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)
        interpolated_traces = []

        if plot_mode == "Summed":
            for i, (charge, group) in enumerate(filtered_df.groupby("Charge")):
                group_sorted = group.sort_values("CCS")
                interp = np.interp(ccs_grid, group_sorted["CCS"], group_sorted[intensity_column], left=0, right=0)
                interpolated_traces.append(interp)
                ax.plot(ccs_grid, interp, color=palette[i], label=f"{int(charge)}+", linewidth=line_thickness)
                ax.fill_between(ccs_grid, 0, interp, color=palette[i], alpha=0.3)

            total_trace = np.sum(interpolated_traces, axis=0)
            ax.plot(ccs_grid, total_trace, color="black", linewidth=line_thickness, label="Total (Interpolated)")
            ax.legend(fontsize=font_size, frameon=False)

        elif plot_mode == "Stacked":
            vertical_offset = 1.2
            for i, (charge, group) in enumerate(filtered_df.groupby("Charge")):
                group_sorted = group.sort_values("CCS")
                interp = np.interp(ccs_grid, group_sorted["CCS"], group_sorted[intensity_column], left=0, right=0)

                if use_scaled == "Unscaled" and interp.max() > 0:
                    interp = interp / interp.max()

                offset_interp = interp + i * vertical_offset
                ax.plot(ccs_grid, offset_interp, color=palette[i], linewidth=line_thickness)
                ax.fill_between(ccs_grid, i * vertical_offset, offset_interp, color=palette[i], alpha=0.3)

                label_x = ccs_grid[0] + 100
                label_y = i * vertical_offset + 0.1
                ax.text(label_x, label_y, f"{int(charge)}+", fontsize=font_size, color=palette[i], va="bottom")

        ax.set_xlim(ccs_range)
        ax.set_xlabel("CCS (Å²)", fontsize=font_size)
        ax.set_yticks([] if plot_mode == "Stacked" else None)
        ax.set_ylabel("" if plot_mode == "Stacked" else "Intensity", fontsize=font_size)
        ax.set_title("Scaled Intensity vs CCS", fontsize=font_size)
        ax.grid(False)

        for label in ax.get_xticklabels():
            label.set_fontsize(font_size)
        for spine in ax.spines.values():
            spine.set_edgecolor("black")
            spine.set_linewidth(1.5)

        st.pyplot(fig)

        fig_buffer = BytesIO()
        fig.savefig(fig_buffer, format='png', dpi=fig_dpi, bbox_inches='tight')
        fig_buffer.seek(0)
        st.download_button("Download CCS Plot as PNG", data=fig_buffer, file_name="ccs_plot.png", mime="image/png", key="ccs_download")

        fig_buffer.seek(0)
        st.download_button("Download CCS Plot as PNG", data=fig_buffer, file_name="ccs_plot.png", mime="image/png")

