import os
import re
from pathlib import Path
import argparse

# list of files and names of variables pointing to external data so that the script can check
file_with_external_data_list = [("simulations/0D/ionization/0D_ionization_1neutralspecies_es.jl", "cs_n_e_filepath"),
                                ("simulations/0D/ionization/0D_ionization_1neutralspecies_nnls_es.jl", "cs_n_e_filepath"),]

def extract_key_value(filepath, key):
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip().startswith(f"{key} = "):
                value = line.split('=')[1].strip().strip('"').strip("'")
                return value
    return None

def find_and_replace_nt(file_path):
    # Read the file
    with open(file_path, 'r') as f:
        content = f.readlines()

    original_nt = None
    modified_lines = []
    for line in content:
        if not line.lstrip().startswith('#'):
            match = re.search(r'^(.*?)(const\s+)?n_t\s*=\s*(\d+)', line)
            if match:
                original_nt = match.group(3)
                prefix = match.group(1)
                const_prefix = match.group(2) or ""
                # Replace the value with 10, preserving const and whitespace
                modified_line = f"{prefix}{const_prefix}n_t = 10\n"
                modified_lines.append(modified_line)
                continue
        modified_lines.append(line)

    if original_nt is None:
        raise ValueError(f"Could not find uncommented 'n_t = NUMBER' or 'const n_t = NUMBER' in {file_path}")

    # Write the modified content back to the file
    with open(file_path, 'w') as f:
        f.writelines(modified_lines)

    return original_nt

def restore_nt(file_path, original_nt):
    # Read the file
    with open(file_path, 'r') as f:
        content = f.readlines()

    restored_lines = []
    for line in content:
        if not line.lstrip().startswith('#'):
            match = re.search(r'^(.*?)(const\s+)?n_t\s*=\s*\d+', line)
            if match:
                prefix = match.group(1)
                const_prefix = match.group(2) or ""
                # Restore the original value, preserving const and whitespace
                restored_line = f"{prefix}{const_prefix}n_t = {original_nt}\n"
                restored_lines.append(restored_line)
                continue
        restored_lines.append(line)

    # Write the restored content back to the file
    with open(file_path, 'w') as f:
        f.writelines(restored_lines)

def ensure_path_exists(path):
    path = Path(path)
    if not path.exists():
        answer = input(f"Path '{path}' does not exist. Create it (scratch directory is already in gitignore)? (y/n) ")
        if answer.lower() == 'y':
            path.mkdir(parents=True, exist_ok=True)
            print(f"Created path: {path}")
        else:
            print(f"Path '{path}' was not created, quitting script execution.")
            exit()
    # else:
    #     print(f"Path '{path}' already exists.")
    return path.exists()


def main(logdir=None):
    for path in Path('simulations').rglob('*.jl'):
        ensure_path_exists("scratch/data")
        ensure_path_exists(logdir)
        print(f"Processing {path}")
        path2 = str(path).replace("/", "_")

        run_sim = True
        for file_kv in file_with_external_data_list:
            if str(path).endswith(file_kv[0]):
                value = extract_key_value(path, file_kv[1])
                if value is not None:
                    exists = os.path.exists(value)
                    if not exists:
                        answer = input(f"External file {value} required by {file_kv[0]} not found. Skip simulation (if no, error will be produced)? (y/n) ")
                        if answer.lower() == 'y':
                            run_sim = False
        
        if run_sim:
            # Find and replace n_t
            original_nt = find_and_replace_nt(path)

            # Run the Julia script
            print(f"Running {path} with n_t = 10")
            os.system(f"julia --project=. {path} > scratch/logs/out_{path2}.log 2> scratch/logs/error_{path2}.log")

            # Restore the original n_t value
            restore_nt(path, original_nt)
            print(f"Restored n_t = {original_nt} in {path}")
        else:
            print(f"Skipping {path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to run example simulations in simulations/")
    parser.add_argument("--logdir", type=str, help="Directory for log file output", default="scratch/logs")
    args = parser.parse_args()
    print(f"Will write logs to {args.logdir}")

    main(logdir=args.logdir)