import glob
import re
import os

def extract_step_metrics(file_path):
    """
    Parses a single CET.out.n file to extract the LAMMPS timestep
    and calculate max/mean electronic temperatures.
    """
    temperatures = []
    lammps_step = None
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            # Read LAMMPS execution step from the header comment
            if line.startswith("#"):
                step_match = re.search(r'at step (\d+)', line)
                if step_match:
                    lammps_step = int(step_match.group(1))
                continue
            
            # Extract temperature from the 4th column
            parts = line.split()
            if len(parts) == 4:
                try:
                    temperatures.append(float(parts[3]))
                except ValueError:
                    continue
                    
    if not temperatures:
        return None
        
    return lammps_step, max(temperatures), sum(temperatures) / len(temperatures)

def main():
    # 1. Dynamically find all files matching 'CET.out.*' in the current folder
    file_pattern = "CET.out.*"
    found_files = glob.glob(file_pattern)
    
    valid_files_data = []
    
    # 2. Extract the suffix 'n' using regular expressions to ensure proper numerical sorting
    for f_path in found_files:
        filename = os.path.basename(f_path)
        match = re.search(r'CET\.out\.(\d+)$', filename)
        if match:
            n_val = int(match.group(1))
            valid_files_data.append((n_val, f_path))
            
    if not valid_files_data:
        print("No files matching 'CET.out.n' were found in the current directory.")
        return
        
    # Sort files numerically based on 'n' (handles 10 before 100, etc.)
    valid_files_data.sort(key=lambda x: x[0])
    
    # 3. Calculate sequence traits
    n_list = [item[0] for item in valid_files_data]
    total_files = len(n_list)
    
    # Determine the stride/step between 'n' suffixes if there are at least two files
    if total_files > 1:
        n_stride = n_list[1] - n_list[0]
    else:
        n_stride = 0 # Only one file present

    print(f"Found {total_files} matching files.")
    print(f"File suffix sequence 'n' starts at {n_list[0]}, ends at {n_list[-1]} (Stride = {n_stride})")
    print("\n" + "="*70)
    print(f"{'File Name':<16} | {'Suffix (n)':<12} | {'LAMMPS Step':<12} | {'Max T (K)':<11} | {'Mean T (K)':<11}")
    print("="*70)
    
    summary_data = []
    
    # 4. Process the dynamically discovered list
    for n_val, file_path in valid_files_data:
        metrics = extract_step_metrics(file_path)
        if metrics:
            lammps_step, t_max, t_mean = metrics
            # If header didn't specify a step, fallback gracefully to the suffix
            display_step = lammps_step if lammps_step is not None else "N/A"
            
            print(f"{os.path.basename(file_path):<16} | {n_val:<12} | {display_step:<12} | {t_max:<11.2f} | {t_mean:<11.2f}")
            summary_data.append((n_val, lammps_step, t_max, t_mean))

    # 5. Save results to log
    output_log = "Te_summary.txt"
    with open(output_log, "w") as log:
        log.write("#File_Suffix_n,LAMMPS_Timestep,Max_Temp,Mean_Temp\n")
        for n_val, lammps_step, t_max, t_mean in summary_data:
            step_str = str(lammps_step) if lammps_step is not None else ""
            log.write(f"{n_val} {step_str} {t_max:.6f} {t_mean:.6f}\n")
            
    print("="*70)
    print(f"Analysis complete. Metrics written to '{output_log}'")

if __name__ == "__main__":
    main()
