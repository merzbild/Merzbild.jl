import os
import re
from pathlib import Path

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

for path in Path('simulations').rglob('*.jl'):
    print(f"Processing {path}")
    path2 = str(path).replace("/", "_")

    # Find and replace n_t
    original_nt = find_and_replace_nt(path)

    # Run the Julia script
    print(f"Running {path} with n_t = 10")
    os.system(f"julia --project=. {path} > scratch/logs/out_{path2}.log 2> scratch/logs/error_{path2}.log")

    # Restore the original n_t value
    restore_nt(path, original_nt)
    print(f"Restored n_t = {original_nt} in {path}")