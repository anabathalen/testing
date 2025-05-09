import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import io
import base64

# Constants
PROTON_MASS = 1.007276

# Set page configuration
st.set_page_config(page_title="CCS Trace Fitting", layout="wide")

# Define Gaussian function
def gaussian(x, amp, mu, sigma):
    """Single Gaussian function"""
    return amp * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

def multi_gaussian(x, *params):
    """Sum of multiple Gaussian functions"""
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        amp, mu, sigma = params[i:i+3]
        y += gaussian(x, amp, mu, sigma)
    return y

def fit_gaussians(x, y, initial_params):
    """Fit multiple gaussians to data using the same core algorithm"""
    # Create reasonable bounds that allow parameters to vary
    lower_bounds = []
    upper_bounds = []
    
    for i in range(0, len(initial_params), 3):
        amp, mu, sigma = initial_params[i:i+3]
        # Amplitude: can vary but stay positive
        lower_bounds.append(max(0.0, amp * 0.1))
        upper_bounds.append(amp * 10.0)
        
        # Mean: can move around but not too far
        lower_bounds.append(max(mu * 0.8, 0))
        upper_bounds.append(mu * 1.2)
        
        # Sigma: can vary but stay positive
        lower_bounds.append(max(1.0, sigma * 0.1))
        upper_bounds.append(sigma * 10.0)
    
    bounds = (lower_bounds, upper_bounds)
    
    try:
        popt, pcov = curve_fit(
            multi_gaussian, x, y, 
            p0=initial_params, 
            bounds=bounds,
            maxfev=5000  # Increase maximum function evaluations
        )
        
        # Calculate goodness of fit (R²)
        fitted_y = multi_gaussian(x, *popt)
        residuals = y - fitted_y
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        
        return popt, r_squared, fitted_y
    
    except Exception as e:
        st.error(f"Fit failed: {str(e)}")
        return initial_params, 0.0, multi_gaussian(x, *initial_params)

def generate_csv_download_link(dataframe, filename):
    """Generate a download link for a CSV file"""
    csv = dataframe.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">{filename}</a>'
    return href

def generate_excel_download_link(fit_results, trace_data, filename="ccs_fit_results.xlsx"):
    """Generate a download link for Excel file with all fit results and traces"""
    output = io.BytesIO()
    
    # Create Excel writer
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        # Write fit parameters
        fit_results.to_excel(writer, sheet_name='Fit Parameters', index=False)
        
        # Write trace data
        trace_data.to_excel(writer, sheet_name='Trace Data', index=False)
        
        # Get the xlsxwriter workbook and worksheet objects
        workbook = writer.book
        
        # Add chart for fit parameters on its own sheet
        chart_sheet = workbook.add_worksheet('Charts')
        chart = workbook.add_chart({'type': 'scatter'})
        
        # Configure series for each peak
        for i in range(len(fit_results['Charge'].unique())):
            chart.add_series({
                'name': f'Charge {fit_results["Charge"].unique()[i]}',
                'categories': ['Fit Parameters', 1, 2, len(fit_results), 2],  # mu values
                'values': ['Fit Parameters', 1, 1, len(fit_results), 1],      # amp values
                'marker': {'type': 'circle', 'size': 8}
            })
            
        chart.set_title({'name': 'Peak Parameters by Charge State'})
        chart.set_x_axis({'name': 'CCS (Å²)'})
        chart.set_y_axis({'name': 'Amplitude'})
        chart_sheet.insert_chart('A1', chart, {'x_scale': 1.5, 'y_scale': 1.5})
    
    # Create download link
    b64 = base64.b64encode(output.getvalue()).decode()
    href = f'<a href="data:application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;base64,{b64}" download="{filename}">{filename}</a>'
    return href

def generate_figure_download_link(fig, filename="ccs_fit_figure.png"):
    """Generate a download link for a figure"""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    buf.seek(0)
    b64 = base64.b64encode(buf.getvalue()).decode()
    href = f'<a href="data:image/png;base64,{b64}" download="{filename}">{filename}</a>'
    return href

def fit_ccs_traces_page():
    st.title("CCS Trace Fitting")
    
    # Instructions in a collapsible section
    with st.expander("Instructions"):
        st.write("""
        1. Upload your calibrated CSV file with columns: 'Charge', 'CCS', 'CCS Std.Dev.', and 'Intensity'
        2. Upload your mass spectrum TXT file (two columns, no headers)
        3. Enter the protein mass in Daltons
        4. Select the CCS range and charge states to analyze
        5. Set peak parameters for each charge state
        6. Click "Fit All Charge States" to process everything at once
        7. Download your results in Origin-compatible format
        """)
    
    # File uploads in two columns
    col1, col2 = st.columns(2)
    with col1:
        cal_file = st.file_uploader("Upload calibrated CCS CSV file", type="csv")
    
    with col2:
        ms_file = st.file_uploader("Upload mass spectrum TXT file", type="txt")
    
    # Parameters in three columns
    col1, col2, col3 = st.columns(3)
    with col1:
        protein_mass = st.number_input("Protein mass (Da)", min_value=1000.0, value=20000.0, step=1000.0)
    
    with col2:
        x_start = st.number_input("CCS range start (Å²)", min_value=500.0, value=1000.0, step=100.0)
    
    with col3:
        x_end = st.number_input("CCS range end (Å²)", min_value=1000.0, value=3000.0, step=100.0)
    
    # Display options in two columns
    col1, col2 = st.columns(2)
    with col1:
        show_components = st.checkbox("Show individual Gaussian components", value=True)
    
    with col2:
        show_residuals = st.checkbox("Show residuals plot", value=True)
    
    # Continue only if files are uploaded
    if not cal_file or not ms_file:
        st.warning("Please upload both required files to continue.")
        return
    
    # Load and process data
    try:
        cal_df = pd.read_csv(cal_file)
        required_columns = {'Charge', 'CCS', 'CCS Std.Dev.', 'Intensity'}
        
        if not required_columns.issubset(cal_df.columns):
            st.error(f"Missing required columns in CSV file. Found: {list(cal_df.columns)}, Need: {list(required_columns)}")
            return
        
        # Filter out unreasonable standard deviations
        cal_df = cal_df[cal_df['CCS Std.Dev.'] < 0.5 * cal_df['CCS']].copy()
        
        # Load mass spectrum data
        try:
            ms_df = pd.read_csv(ms_file, sep="\t", header=None, names=["m/z", "Intensity"])
            ms_df.dropna(inplace=True)
        except Exception as e:
            st.error(f"Error reading mass spectrum file: {str(e)}")
            st.info("Make sure the file is tab-separated with two columns (m/z and intensity).")
            return
            
    except Exception as e:
        st.error(f"Error processing input files: {str(e)}")
        return
    
    # Get and select available charge states
    all_charges = sorted(cal_df["Charge"].unique())
    
    if not all_charges:
        st.error("No charge states found in the data.")
        return
    
    selected_charges = st.multiselect(
        "Select charge states to include", 
        all_charges, 
        default=all_charges
    )
    
    if not selected_charges:
        st.warning("Please select at least one charge state.")
        return
    
    # Filter data by selected charges
    cal_df = cal_df[cal_df["Charge"].isin(selected_charges)]
    
    # Calculate scaling factors based on MS data
    scale_factors = {}
    for z in selected_charges:
        mz = (protein_mass + z * PROTON_MASS) / z
        window_pct = 0.01  # 1% window
        mz_min, mz_max = mz * (1 - window_pct), mz * (1 + window_pct)
        intensity_sum = ms_df[(ms_df["m/z"] >= mz_min) & (ms_df["m/z"] <= mz_max)]["Intensity"].sum()
        
        # If no intensity found, use a small value to avoid zeros
        if intensity_sum <= 0:
            intensity_sum = 1.0
            
        scale_factors[z] = intensity_sum
    
    # Apply scaling factors
    cal_df["Scale Factor"] = cal_df["Charge"].map(scale_factors)
    cal_df["Scaled Intensity"] = cal_df["Intensity"] * cal_df["Scale Factor"]
    
    # Create fit grid and containers for all results
    x_fit = np.linspace(x_start, x_end, 1000)
    all_traces = []
    all_params = []
    all_fitted_curves = []
    all_r_squared = []
    
    # Peak parameters setup
    st.write("---")
    st.subheader("Peak Configuration")
    
    # Create a dictionary to store peak configurations for each charge state
    peak_configs = {}
    
    # Use tabs for individual vs. batch fitting
    tab1, tab2 = st.tabs(["Individual Charge States", "Batch Fitting"])
    
    with tab1:
        # Create columns for charge state parameters - 3 charge states per row
        charge_groups = [selected_charges[i:i+3] for i in range(0, len(selected_charges), 3)]
        
        for group in charge_groups:
            cols = st.columns(len(group))
            for i, charge in enumerate(group):
                with cols[i]:
                    st.write(f"**Charge {charge}**")
                    num_peaks = st.number_input(
                        f"Number of peaks", 
                        min_value=1, 
                        max_value=5, 
                        value=1, 
                        key=f"num_peaks_{charge}"
                    )
                    peak_configs[charge] = {"num_peaks": num_peaks, "params": []}
                    
                    # Create fields for initial parameters of each peak
                    for j in range(num_peaks):
                        st.write(f"Peak {j+1}")
                        
                        # Extract data for this charge state to get reasonable defaults
                        sub_df = cal_df[cal_df["Charge"] == charge]
                        x_raw = sub_df["CCS"].values
                        y_raw = sub_df["Scaled Intensity"].values
                        
                        if len(x_raw) > 0:
                            max_idx = np.argmax(y_raw)
                            default_amp = float(y_raw[max_idx])
                            default_mu = float(x_raw[max_idx])
                        else:
                            default_amp = 1000.0
                            default_mu = (x_start + x_end) / 2
                        
                        amp = st.number_input(
                            "Amplitude", 
                            value=default_amp,
                            format="%.1f", 
                            key=f"amp_{charge}_{j}"
                        )
                        
                        mu = st.number_input(
                            "Center (CCS)", 
                            value=default_mu,
                            step=10.0, 
                            key=f"mu_{charge}_{j}"
                        )
                        
                        sigma = st.number_input(
                            "Width (σ)", 
                            value=30.0,
                            min_value=1.0, 
                            step=5.0, 
                            key=f"sigma_{charge}_{j}"
                        )
                        
                        peak_configs[charge]["params"].extend([amp, mu, sigma])
        
        # Individual fitting buttons
        for charge in selected_charges:
            st.write("---")
            st.subheader(f"Fit Charge State {charge}")
            
            # Extract data for this charge state
            sub_df = cal_df[cal_df["Charge"] == charge]
            x_raw = sub_df["CCS"].values
            y_raw = sub_df["Scaled Intensity"].values
            
            if len(x_raw) < 3:
                st.warning(f"Not enough data points for charge state {charge}. Skipping.")
                continue
                
            # Interpolate to fit grid
            y_interp = np.interp(x_fit, x_raw, y_raw)
            
            if st.button(f"Fit Charge {charge}"):
                with st.spinner(f"Fitting charge state {charge}..."):
                    # Get initial parameters for this charge
                    init_params = peak_configs[charge]["params"]
                    
                    # Perform the fit
                    optimized_params, r_squared, fitted_y = fit_gaussians(x_fit, y_interp, init_params)
                    
                    # Store results
                    for j in range(0, len(optimized_params), 3):
                        peak_num = j // 3 + 1
                        amp, mu, sigma = optimized_params[j:j+3]
                        all_params.append([charge, peak_num, amp, mu, sigma])
                    
                    all_r_squared.append([charge, r_squared])
                    all_traces.append(y_interp)
                    all_fitted_curves.append(fitted_y)
                    
                    # Plot results
                    fig, axs = plt.subplots(2 if show_residuals else 1, 1, 
                                           figsize=(10, 8 if show_residuals else 5),
                                           gridspec_kw={'height_ratios': [3, 1]} if show_residuals else None)
                    
                    ax = axs[0] if show_residuals else axs
                    
                    # Plot raw data
                    ax.scatter(x_raw, y_raw, color='blue', alpha=0.7, s=30, label="Raw Data")
                    ax.plot(x_fit, y_interp, color='blue', alpha=0.5, label="Interpolated")
                    
                    # Plot fitted curve
                    ax.plot(x_fit, fitted_y, color='red', linewidth=2, label=f"Fitted (R² = {r_squared:.4f})")
                    
                    # Plot individual components if requested
                    if show_components:
                        for j in range(0, len(optimized_params), 3):
                            amp, mu, sigma = optimized_params[j:j+3]
                            component = gaussian(x_fit, amp, mu, sigma)
                            ax.plot(x_fit, component, '--', alpha=0.7, label=f"Peak {j//3+1}: CCS={mu:.1f}, σ={sigma:.1f}")
                    
                    ax.set_title(f"Charge State {charge} - Fitted Curve")
                    ax.set_xlabel("CCS (Å²)")
                    ax.set_ylabel("Scaled Intensity")
                    ax.legend(loc='upper right')
                    ax.grid(alpha=0.3)
                    
                    # Add residuals plot if requested
                    if show_residuals:
                        residuals = y_interp - fitted_y
                        axs[1].plot(x_fit, residuals, 'k-')
                        axs[1].axhline(y=0, color='r', linestyle='-', alpha=0.3)
                        axs[1].set_xlabel("CCS (Å²)")
                        axs[1].set_ylabel("Residuals")
                        axs[1].grid(alpha=0.3)
                    
                    plt.tight_layout()
                    st.pyplot(fig)
                    
                    # Display fit parameters
                    results_df = pd.DataFrame(
                        [[j//3+1, optimized_params[j], optimized_params[j+1], optimized_params[j+2]] 
                         for j in range(0, len(optimized_params), 3)],
                        columns=["Peak", "Amplitude", "Center (CCS)", "Width (σ)"]
                    )
                    st.write("Fitted parameters:")
                    st.dataframe(results_df)
                    
                    # Add download buttons for this specific fit
                    fit_csv = results_df.copy()
                    fit_csv.insert(0, "Charge", charge)
                    
                    st.markdown(
                        generate_csv_download_link(fit_csv, f"charge_{charge}_fit_params.csv"),
                        unsafe_allow_html=True
                    )
                    
                    # Raw data with fit as CSV
                    trace_df = pd.DataFrame({
                        "CCS": x_fit,
                        "Raw": y_interp,
                        "Fitted": fitted_y,
                        "Residuals": y_interp - fitted_y
                    })
                    
                    # Add individual components
                    for j in range(0, len(optimized_params), 3):
                        peak_num = j // 3 + 1
                        amp, mu, sigma = optimized_params[j:j+3]
                        trace_df[f"Peak_{peak_num}"] = gaussian(x_fit, amp, mu, sigma)
                    
                    st.markdown(
                        generate_csv_download_link(trace_df, f"charge_{charge}_trace_data.csv"),
                        unsafe_allow_html=True
                    )
                    
                    # Export figure
                    st.markdown(
                        generate_figure_download_link(fig, f"charge_{charge}_fit.png"),
                        unsafe_allow_html=True
                    )
    
    with tab2:
        st.subheader("Batch Process All Charge States")
        
        # Button to fit all charge states at once
        if st.button("Fit All Charge States"):
            batch_results = []
            batch_traces = {}
            all_components = {}
            
            # Store CCS values for traces
            batch_traces["CCS"] = x_fit
            
            # Create progress bar
            progress_bar = st.progress(0)
            
            for i, charge in enumerate(selected_charges):
                progress = (i / len(selected_charges))
                progress_bar.progress(progress)
                
                # Extract data for this charge state
                sub_df = cal_df[cal_df["Charge"] == charge]
                x_raw = sub_df["CCS"].values
                y_raw = sub_df["Scaled Intensity"].values
                
                if len(x_raw) < 3:
                    st.warning(f"Not enough data points for charge state {charge}. Skipping.")
                    continue
                    
                # Interpolate to fit grid
                y_interp = np.interp(x_fit, x_raw, y_raw)
                batch_traces[f"Raw_{charge}"] = y_interp
                
                # Get initial parameters for this charge
                init_params = peak_configs[charge]["params"]
                
                # Perform the fit
                optimized_params, r_squared, fitted_y = fit_gaussians(x_fit, y_interp, init_params)
                batch_traces[f"Fit_{charge}"] = fitted_y
                batch_traces[f"Residual_{charge}"] = y_interp - fitted_y
                
                # Store components
                all_components[charge] = []
                for j in range(0, len(optimized_params), 3):
                    peak_num = j // 3 + 1
                    amp, mu, sigma = optimized_params[j:j+3]
                    
                    # Store parameter
                    batch_results.append([charge, peak_num, amp, mu, sigma, r_squared])
                    
                    # Store component trace
                    component = gaussian(x_fit, amp, mu, sigma)
                    batch_traces[f"Peak_{charge}_{peak_num}"] = component
                    all_components[charge].append(component)
            
            progress_bar.progress(1.0)
            
            # Create combined results DataFrame
            results_df = pd.DataFrame(
                batch_results,
                columns=["Charge", "Peak", "Amplitude", "Center (CCS)", "Width (σ)", "R²"]
            )
            
            # Sort by charge and peak number
            results_df = results_df.sort_values(["Charge", "Peak"])
            
            # Create trace DataFrame
            trace_df = pd.DataFrame(batch_traces)
            
            # Show combined results
            st.subheader("Fit Results for All Charge States")
            st.dataframe(results_df)
            
            # Plot all charge states
            st.subheader("Visualization of All Fits")
            
            # Create tabs for different visualization types
            vis_tab1, vis_tab2, vis_tab3 = st.tabs(["Individual Fits", "Summed Trace", "Peak Summary"])
            
            with vis_tab1:
                for charge in selected_charges:
                    # Skip charges without enough data
                    if f"Raw_{charge}" not in batch_traces:
                        continue
                        
                    fig, ax = plt.subplots(figsize=(10, 5))
                    
                    # Plot raw data and fit
                    ax.plot(x_fit, batch_traces[f"Raw_{charge}"], 'b-', label=f"Raw Charge {charge}")
                    ax.plot(x_fit, batch_traces[f"Fit_{charge}"], 'r-', label=f"Fitted Curve")
                    
                    # Plot components
                    if show_components and charge in all_components:
                        for i, component in enumerate(all_components[charge]):
                            ax.plot(x_fit, component, '--', alpha=0.7, label=f"Peak {i+1}")
                    
                    ax.set_title(f"Charge State {charge}")
                    ax.set_xlabel("CCS (Å²)")
                    ax.set_ylabel("Scaled Intensity")
                    ax.legend()
                    ax.grid(alpha=0.3)
                    
                    st.pyplot(fig)
            
            with vis_tab2:
                # Create summed trace
                summed_raw = np.zeros_like(x_fit)
                summed_fit = np.zeros_like(x_fit)
                
                for charge in selected_charges:
                    if f"Raw_{charge}" in batch_traces:
                        summed_raw += batch_traces[f"Raw_{charge}"]
                        summed_fit += batch_traces[f"Fit_{charge}"]
                
                fig, ax = plt.subplots(figsize=(10, 5))
                ax.plot(x_fit, summed_raw, 'b-', label="Summed Raw Data")
                ax.plot(x_fit, summed_fit, 'r-', label="Summed Fits")
                ax.set_title("Summed Trace Across All Charge States")
                ax.set_xlabel("CCS (Å²)")
                ax.set_ylabel("Scaled Intensity")
                ax.legend()
                ax.grid(alpha=0.3)
                
                st.pyplot(fig)
                
                # Add summed traces to trace dataframe
                trace_df["Summed_Raw"] = summed_raw
                trace_df["Summed_Fit"] = summed_fit
            
            with vis_tab3:
                # Create scatter plot of peak centers and amplitudes
                fig, ax = plt.subplots(figsize=(10, 5))
                
                for charge in results_df["Charge"].unique():
                    charge_data = results_df[results_df["Charge"] == charge]
                    ax.scatter(
                        charge_data["Center (CCS)"], 
                        charge_data["Amplitude"],
                        label=f"Charge {charge}",
                        s=50
                    )
                
                ax.set_title("Peak Parameters Summary")
                ax.set_xlabel("Peak Center (CCS)")
                ax.set_ylabel("Peak Amplitude")
                ax.legend()
                ax.grid(alpha=0.3)
                
                st.pyplot(fig)
            
            # Download options
            st.subheader("Download Results")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.markdown(
                    generate_csv_download_link(results_df, "all_fit_parameters.csv"),
                    unsafe_allow_html=True
                )
            
            with col2:
                st.markdown(
                    generate_csv_download_link(trace_df, "all_trace_data.csv"),
                    unsafe_allow_html=True
                )
            
            with col3:
                st.markdown(
                    generate_excel_download_link(results_df, trace_df, "ccs_fit_results.xlsx"),
                    unsafe_allow_html=True
                )

# Run the application
fit_ccs_traces_page()
