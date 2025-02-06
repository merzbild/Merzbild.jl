# script to test that all simulations in simulations directory are working - when breaking changes are made in the code
# requires external LXCat data for the ionization simulations
import os
from pathlib import Path

for path in Path('simulations').rglob('*.jl'):
    print(f"Running {path}")
    os.system(f"julia --project=. {path}")