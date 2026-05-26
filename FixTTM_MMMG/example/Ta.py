import os
import sys

def analyze_chunk_file(file_path, output_file):
    """
    Parses a LAMMPS chunk/time output file, calculates the max and mean values 
    for each timestep block, and saves the summary into the specified output file.
    """
    if not os.path.exists(file_path):
        print(f"Warning: File '{file_path}' not found. Skipping...")
        return False

    print(f"\nAnalyzing file: {file_path}")
    print(f"{'File_Suffix_n':<15} | {'LAMMPS_Timestep':<15} | {'Max_Temp':<12} | {'Mean_Temp':<12}")
    print("-" * 62)

    with open(file_path, 'r') as f:
        lines = f.readlines()

    line_idx = 0
    total_lines = len(lines)
    
    # List to store our summary rows for file output
    summary_data = []
    block_counter = 1  # Tracks the sequential suffix index (1, 2, 3...)

    # Loop through the entire file line by line
    while line_idx < total_lines:
        line = lines[line_idx].strip()

        # Skip comment lines and empty lines
        if not line or line.startswith('#'):
            line_idx += 1
            continue

        parts = line.split()
        
        # Detect block header line (e.g., "0 1000" or "1000 1000")
        if len(parts) == 2:
            try:
                timestep = int(parts[0])
                num_chunks = int(parts[1])
                
                temperatures = []
                line_idx += 1  # Advance to the first row of data inside the block
                chunks_read = 0

                # Read exactly the number of chunks specified by the header
                while chunks_read < num_chunks and line_idx < total_lines:
                    data_line = lines[line_idx].strip()
                    
                    # Ignore unexpected empty lines or comments inside a block
                    if not data_line or data_line.startswith('#'):
                        line_idx += 1
                        continue

                    data_parts = data_line.split()
                    if len(data_parts) >= 2:
                        try:
                            # Column 0 is the row index/chunk ID, Column 1 is the temperature value
                            temp_val = float(data_parts[1])
                            temperatures.append(temp_val)
                            chunks_read += 1
                        except ValueError:
                            pass
                    
                    line_idx += 1

                # Calculate metrics for the completed block
                if temperatures:
                    max_temp = max(temperatures)
                    mean_temp = sum(temperatures) / len(temperatures)
                    
                    # Print progress directly to the screen
                    print(f"{block_counter:<15} | {timestep:<15} | {max_temp:<12.6f} | {mean_temp:<12.6f}")
                    
                    # Save structured text line for our output file array
                    summary_data.append(f"{block_counter} {timestep} {max_temp:.6f} {mean_temp:.6f}\n")
                    block_counter += 1
                else:
                    print(f"Warning: Timestep {timestep} contained no valid numeric values.")
                
                # Continue loop without incrementing line_idx since inner loop already advanced it
                continue

            except ValueError:
                # Line has 2 columns but they are not integers (not a header line)
                line_idx += 1
                continue
        else:
            line_idx += 1

    # Write data out to the specified output file if data was found
    if summary_data:
        with open(output_file, 'w') as out_f:
            # Write the header line specified by your prompt
            out_f.write("#File_Suffix_n,LAMMPS_Timestep,Max_Temp,Mean_Temp\n")
            out_f.writelines(summary_data)
        print(f"\nSuccess: Summary metrics compiled and saved to '{output_file}'.")
        return True
    else:
        print(f"\nWarning: No valid snapshot datasets were discovered in '{file_path}' to generate a summary.")
        return False

if __name__ == "__main__":
    # Define the files to process
    files_to_process = [
        ("CIT_filtered.out", "Ta_filtered_summary.txt"),
        ("CIT.out", "Ta_summary.txt")
    ]
    
    # Check if command line arguments are provided
    if len(sys.argv) > 1:
        # If arguments provided, use the first as input file and second as output file (if given)
        input_file = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else "Ta_summary.txt"
        analyze_chunk_file(input_file, output_file)
    else:
        # Process both default files
        print("Processing default files...")
        for input_file, output_file in files_to_process:
            analyze_chunk_file(input_file, output_file)
        print("\nAll files processed.")
