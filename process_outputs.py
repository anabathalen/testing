import os
import pandas as pd

def process_data(input_filename, output_filename):
    # Read the input file
    input_df = pd.read_csv(input_filename, sep=' ', header=None, names=['Index', 'Mass', 'Z', 'DriftTime'])
    
    # Read the output file
    with open(output_filename, 'r') as file:
        lines = file.readlines()
    
    # Extract the protein name from the output file
    protein_name = lines[1].split(',')[0]
    
    # Extract calibrated data from the output file (from the [CALIBRATED DATA] section)
    start_idx = lines.index('[CALIBRATED DATA]\n') + 1
    calibrated_data = []
    
    for line in lines[start_idx:]:
        if line.strip():  # Ignore empty lines
            data = line.strip().split(',')
            calibrated_data.append([data[0], float(data[1]), int(data[2]), float(data[3]), float(data[4]), float(data[5])])

    # Convert the calibrated data to a DataFrame
    calibrated_df = pd.DataFrame(calibrated_data, columns=['ID', 'Mass', 'Z', 'Drift', 'CCS', 'CCS_StdDev'])

    # Merge the input and calibrated data on 'Mass' and 'DriftTime'
    merged_data = pd.merge(input_df, calibrated_df, on='Mass', how='left')

    # Add the protein name
    merged_data['Protein'] = protein_name

    return merged_data

def combine_all_files(input_dir, output_dir, output_csv_path):
    # Initialize an empty DataFrame to collect all data
    all_data = []

    # Walk through the input and output directories
    for protein_folder in os.listdir(input_dir):
        protein_folder_path = os.path.join(input_dir, protein_folder)
        
        if os.path.isdir(protein_folder_path):  # Make sure it's a directory
            for filename in os.listdir(protein_folder_path):
                if filename.startswith('input_') and filename.endswith('.txt'):
                    charge_state = filename.split('_')[1].split('.')[0]
                    
                    # Construct corresponding output file path
                    output_filename = f'output_{charge_state}.txt'
                    output_file_path = os.path.join(output_dir, protein_folder, output_filename)
                    input_filename = os.path.join(protein_folder_path, filename)
                    
                    # Ensure both input and output files exist
                    if os.path.exists(input_filename) and os.path.exists(output_file_path):
                        print(f"Processing {input_filename} and {output_file_path}")
                        # Process the data and merge
                        merged_data = process_data(input_filename, output_file_path)
                        all_data.append(merged_data)

    # Combine all the data into a single DataFrame
    final_data = pd.concat(all_data, ignore_index=True)

    # Save the combined data to a CSV file
    final_data.to_csv(output_csv_path, index=False)
    print(f"Combined data saved to {output_csv_path}")

# Example usage
input_dir = '/path/to/input/directory'  # Path to the directory containing protein folders
output_dir = '/path/to/output/directory'  # Path to the directory containing output files
output_csv_path = 'combined_data.csv'  # Path to save the final combined CSV

combine_all_files(input_dir, output_dir, output_csv_path)
