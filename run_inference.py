import os
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Run inference for protein structure prediction.")
parser.add_argument("msas_dir", help="Directory containing MSA files.")
parser.add_argument("output_dir", help="Directory to store the output PDB files.")
parser.add_argument("prefix", help="Prefix for input and output files (e.g., '1ake').")
parser.add_argument("--af2_run_path", default="scripts/RunAF2.py", help="Path to the AF2 inference script (default: scripts/RunAF2.py).")
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Process each MSA file, default to 100 files
for index in range(100):
    msa_filename = os.path.join(args.msas_dir, f"{args.prefix}_{index:03}.a3m")
    pdb_filename = os.path.join(args.output_dir, f"{args.prefix}_{index:03}.pdb")
    
    if os.path.exists(pdb_filename):
        print(f"PDB for {args.prefix} index {index} already exists, skipping inference.")
    else:
        print(f"Running inference for {args.prefix} index {index}.")
        cmd = (
            f"python {args.af2_run_path} {msa_filename} "
            f"--model_num 3 --recycles 1 --output_dir {args.output_dir}"
        )
        os.system(cmd)
