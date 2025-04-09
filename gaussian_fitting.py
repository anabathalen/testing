# gaussian_fitting.py
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import seaborn as sns
import io

# Gaussian function definition
def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-(x - mean)**2 / (2 * stddev**2))

# Find major local maxima (points where the y values go down on both sides 5 times) and intensity >= 0.1 * max intensity
def find_major_local_maxima(x, y, window_size=5, intensity_threshold=0.1):
    maxima_indices = []
    max_intensity = max(y)
    for i in range(window_size, len(y) - window_size):
        if y[i] == max(y[i - window_size:i + window_size + 1]) and y[i] >= intensity_threshold * max_intensity:
            maxima_indices.append(i)
    return maxima_indices

# Main function to upload data, detect peaks, and fit Gaussians
def upload_and_plot():
    # Allow the user to upload a CSV file
    uploaded_file = st.file_uploader("Upload your csv file here:", type="csv")
    
    if uploaded_file is not None:
        
        # Read the CSV into a pandas DataFrame
        df = pd.read_csv(uploaded_file)
        
        # Show the uploaded data as a DataFrame for reference
        st.write("Your data:")
        st.dataframe(df)
        
        # Check if the DataFrame contains 'x' and 'y' columns
        if 'x' in df.columns and 'y' in df.columns:
            # Find and label major local maxima
            maxima_indices = find_major_local_maxima(df['x'], df['y'])
            maxima_x = df['x'].iloc[maxima_indices]
            maxima_y = df['y'].iloc[maxima_indices]

            # Plot the raw data with labeled major maxima
            fig, ax = plt.subplots()
            ax.plot(df['x'], df['y'], label='Data', color='black', alpha=1.0, linewidth=1)  # Data as line
            ax.scatter(maxima_x, maxima_y, color='red', label='Major Local Maxima', zorder=5)

            # Annotate the major local maxima with x values
            for i, x_val in enumerate(maxima_x):
                # Adjust the annotation position to ensure it stays within the plot
                y_val = maxima_y.iloc[i]
                offset = 10
                # Ensure the label doesn't go outside the plot boundaries
                if x_val < min(df['x']) + 0.1 * (max(df['x']) - min(df['x'])):
                    offset = -10  # Move the label slightly to the right if too close to left boundary
                elif x_val > max(df['x']) - 0.1 * (max(df['x']) - min(df['x'])):
                    offset = -10  # Move the label slightly to the left if too close to right boundary
                
                ax.annotate(f'{x_val:.2f}', (x_val, y_val), textcoords="offset points", xytext=(0, offset), ha='center')

            ax.set_xlabel("Drift Time (Bins)")
            ax.set_ylabel("Intensity")
            ax.legend()
            st.pyplot(fig)

            # Ask for initial guesses for the Gaussian means (peak x-values)
            st.write("Now that the major maxima have been identified, please input the number of Gaussian peaks and their positions.")
            num_gaussians = st.number_input("How many Gaussians would you like to fit to the data?", min_value=1, max_value=10, value=1)

            # Ask for the Gaussian center guesses (the means of the Gaussians)
            peaks = []
            for i in range(num_gaussians):
                peak_guess = st.number_input(f"Enter the initial guess for the {i+1}th peak (mean of Gaussian {i+1}):", value=float(df['x'].median()))
                peaks.append(peak_guess)

            # Once the user has entered all the Gaussian centers, allow the user to customize the plot
            with st.expander("Plot Customization Options", expanded=True):
                # First, figure-related customization (size, DPI)
                st.header("Figure Customization")
                fig_width = st.slider("Figure Width (inches)", min_value=2, max_value=15, value=8)
                fig_height = st.slider("Figure Height (inches)", min_value=2, max_value=15, value=6)
                dpi = st.slider("Select DPI", min_value=50, max_value=1000, value=300)

                # Then, aesthetic-related customization (font, color, etc.)
                st.header("Aesthetic Customization")
                font_size = st.slider("Font Size", min_value=6, max_value=20, value=12)
                x_label = st.text_input("Enter X-axis Label", "Drift Time (Bins)")
                color_palette = st.selectbox("Choose a Color Palette", options=["Set1", "Set2", "Paired", "Pastel1", "Dark2"])
                line_width = st.slider("Line Width for Data Plot", min_value=0, max_value=5, value=1)

                # Allow the user to proceed once all inputs are filled
                if st.button("Generate Plot and Fit Gaussians"):
                    # Prepare for Gaussian fitting
                    x_data = df['x']
                    y_data = df['y']

                    # Prepare the initial parameters for curve fitting
                    initial_guess = []
                    for peak in peaks:
                        # Initial guess: amplitude (max y), mean (user input), and standard deviation (arbitrary, set to 1)
                        initial_guess += [max(y_data), peak, 100]

                    # Fitting function to include only a local region (x-10 to x+10)
                    def multi_gaussian(x, *params):
                        result = np.zeros_like(x)
                        for i in range(num_gaussians):
                            amp, mean, stddev = params[3*i:3*(i+1)]
                            result += gaussian(x, amp, mean, stddev)
                        return result

                    # Once customization is done, fit the Gaussians and plot the result
                    fig, ax = plt.subplots()
                    ax.plot(df['x'], df['y'], label='Data', color='black', alpha=1.0, linewidth=line_width)  # Data as line

                    # Get the selected color palette
                    colors = sns.color_palette(color_palette, n_colors=num_gaussians)

                    # Create a high-resolution x-axis (full range) for the fit
                    x_full = np.linspace(min(x_data), max(x_data), 1000)

                    # Loop through each peak guess to perform the fitting and plot the results
                    for i, peak in enumerate(peaks):
                        # Define the local region from peak - 10 to peak + 10
                        x_range_min = peak - peak*0.05
                        x_range_max = peak + peak*0.05

                        # Get the subset of data within the x-range [peak-10, peak+10]
                        mask = (x_data >= x_range_min) & (x_data <= x_range_max)
                        x_local = x_data[mask]
                        y_local = y_data[mask]

                        # Ensure we have enough points for fitting
                        if len(x_local) < 3:
                            continue  # Skip this peak if there's not enough data

                        # Initial guess for this local region
                        local_guess = [max(y_local), peak, 1]  # Amplitude, Mean (the peak), Stddev

                        # Fit the Gaussians for this region only
                        try:
                            # Fit the Gaussian to the local data
                            popt, _ = curve_fit(gaussian, x_local, y_local, p0=local_guess)
                            
                            # Extract parameters from the fitting result
                            amp, mean, stddev = popt

                            # Generate y values across the full x_full range to plot the Gaussian
                            y_fit = gaussian(x_full, amp, mean, stddev)

                            # Plot the fitted Gaussian across the full range with a transparent fill
                            ax.fill_between(x_full, y_fit, color=colors[i], alpha=0.3, label=f'mean = {mean:.2f}')
                            
                        except Exception as e:
                            continue  # Skip this peak if fitting fails

                    # Update plot aesthetics based on user settings
                    ax.set_xlabel(x_label, fontsize=font_size)
                    ax.tick_params(axis='y', labelleft=False, left=False, right=False)
                    ax.set_ylabel("", fontsize=font_size)

                    # Update X-axis tick label font size
                    ax.tick_params(axis='x', labelsize=font_size)

                    # Remove grey line around the legend
                    ax.legend(fontsize=font_size, frameon=False)

                    # Adjust figure size
                    fig.set_size_inches(fig_width, fig_height)
                    plt.rcParams.update({'font.size': font_size})  # Update font size globally

                    # Show the plot to the user
                    st.pyplot(fig)

                    # Allow user to download the customized plot
                    buf = io.BytesIO()
                    fig.savefig(buf, format="png", dpi=dpi)
                    buf.seek(0)

                    st.download_button(
                        label="Download Customized Plot as PNG",
                        data=buf,
                        file_name="customized_gaussian_plot.png",
                        mime="image/png"
                    )

