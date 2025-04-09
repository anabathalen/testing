import pandas as pd
import os

def read_input_file(input_filename):
    """Read and return data from the input_X.txt file."""
    data = pd.read_csv(input_filename, sep=' ', header=None, names=['Index', 'Mass', 'Z', 'DriftTime'])
    return data

def read_output_file(output_filename):
    """Read and extract the calibrated data section from output_X.txt."""
    with open(output_filename, 'r') as file:
        lines = file.readlines()

    # Find the starting line of the 'DIAGNOSTICS' and 'CALIBRATED DATA' sections
    start_index_diagnostics = lines.index('[DIAGNOSTICS]\n') + 2  # Skip header lines
    start_index_calibrated = lines.index('[CALIBRATED DATA]\n') + 2  # Skip header lines
    diagnostics_data = []
    calibrated_data = []

    # Extract Diagnostics Data (Protein Information)
    for line in lines[start_index_diagnostics:start_index_calibrated-2]:
        fields = line.strip().split(',')
        diagnostics_data.append(fields)

    # Extract Calibrated Data
    for line in lines[start_index_calibrated:]:
        if line.startswith('['):  # Stop when another section begins
            break
        fields = line.strip().split(',')
        calibrated_data.append(fields)

    # Create DataFrames from the extracted data
    diagnostics_df = pd.DataFrame(diagnostics_data, columns=['ID', 'Mass', 'Charge', 'Mobility', 'Alpha', 'Gamma', 'Model_Vel', 'Exp_Vel', 'Error'])
    calibrated_df = pd.DataFrame(calibrated_data, columns=['ID', 'Mass', 'Charge', 'DriftTime', 'CCS', 'CCS_StdDev'])
    
    # Merge the diagnostics and calibrated data on 'ID'
    merged_output_df = pd.merge(calibrated_df, diagnostics_df[['ID', 'Mass', 'Charge']], on='ID', how='left')
    
    return merged_output_df

def merge_data(input_df, output_df):
    """Merge the input data with the output calibrated data on Drift Time and Mass."""
    merged_df = pd.merge(input_df, output_df, on='DriftTime', how='left')

    # Handle cases with multiple drift times or masses by dropping duplicates
    merged_df = merged_df.drop_duplicates(subset=['Mass', 'DriftTime'])

    return merged_df

def process_data(input_filename, output_filename):
    """Read the input and output files, merge them, and return the result."""
    input_df = read_input_file(input_filename)
    output_df = read_output_file(output_filename)

    # Merge input data with output data based on DriftTime
    merged_df = merge_data(input_df, output_df)

    return merged_df

def combine_all_files(input_filenames, output_filenames, output_csv_path):
    """Process and combine data from multiple input/output files into a single CSV."""
    combined_df = pd.DataFrame()

    # Loop over each pair of input and output files
    for input_filename, output_filename in zip(input_filenames, output_filenames):
        # Process each file
        merged_data = process_data(input_filename, output_filename)
        
        # Append the result to the combined dataframe
        combined_df = pd.concat([combined_df, merged_data], ignore_index=True)

    # Save the combined data to a CSV file
    combined_df.to_csv(output_csv_path, index=False)
    print(f"Saved combined data to {output_csv_path}")

# Example usage
input_filenames = ['input_5.txt', 'input_6.txt']  # List of input files (adjust as needed)
output_filenames = ['output_5.txt', 'output_6.txt']  # List of corresponding output files
output_csv_path = 'combined_output.csv'  # Path to save the combined CSV

# Combine all data and save to a CSV
combine_all_files(input_filenames, output_filenames, output_csv_path)

