# Overview of scripts
This directory contains various utility scripts for fast plotting/post-processing of results, as well as scripts
to produce plots from papers (in the `reproducibility` sub-directory). For a description of scripts
for reproducing plots from a specific paper, navigate to the corresponding sub-directory of `reproducibility`
and take a look at the `README.md` there.

The Python scripts require Python3 and the `numpy`, `matplotlib`, and `netCDF4` packages.

## Time-averaged ionization rate and electron temperature
The script `compute_rate.py` computes time-averaged ionization rate and electron temperature from a 0D unsteady ionization simulation
with a single neutral and single electron species.
Example usage (if order of species in simulation was neutral/ion/electron, and timestep was 5e-14):
```bash
python3 scripts/compute_rate.py --filename scratch/data/ionization_Ar_400Tn_NNLSrate_approx_5full_69_es.nc --tstart 150000 --tend 500000 --nid 0 --ionid 1 --eid 2 --dt 5e-14
```

## Unsteady 0D simulations
The script `plot_0D_1species.py` plots a single macroscopic variable from multiple files over timesteps.
Example usage:
```bash
python3 scripts/plot_0D_1species.py --files test/data/bkw_20k_seed1234.nc test/data/bkw_vw_grid_seed1234.nc --propname=T --plotname plots/T_BKW_and_vw.png --labels "Fixed weight" "VW+merging"
```

## 1D simulations
The script `plot_1D.py` plots a single macroscopic variable from multiple files over cell index for single-species simulations.
All timesteps in the output file with index `>=startt` will be plotted (so if the file has only a single output, for example, the file
is the time-averaged output, set `startt` to 0). Optionally, output from a SPARTA simulation can also be plotted
(it is assumed that in the SPARTA output the cells are sorted, and they are the same ones as used in Merzbild.jl; the SPARTA output should begin at line 10, and line 9 is the header line).
Example usage:
```bash
python scripts/plot_1D.py --file scratch/data/couette_0.0005_50_500.0_300.0_1000.nc --propname T --startt 1 --plotname plots/couette_T.png --labels "Couette"
```

## Paper reproducibility scripts
Scripts to post-process simulation results for various papers can be found in the `reproducibility` subdirectory.
