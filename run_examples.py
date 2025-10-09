# script to test that all simulations in simulations directory are working - when breaking changes are made in the code
# requires external LXCat data for the ionization simulations
import os
from pathlib import Path

for path in Path('simulations').rglob('*.jl'):
    print(f"Running {path}")
    path2 = str(path).replace("/", "_")
    os.system(f"julia --project=. {path} > scratch/logs/out_{path2}.log 2> scratch/logs/error_{path2}.log")