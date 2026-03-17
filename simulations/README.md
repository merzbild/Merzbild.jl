# Overview of example simulations
This file provides a short overview of the different simulations in `simulations`. Specific simulation setups that were used
for results in publications are found at the end of this file in the "Paper reproducibility" section.
Unless noted, simulations use the NTC collision algorithm (or the variable-weight version thereof).
By default simulations output data into `scratch/data`, and `scratch` is included in `.gitignore`.

## 0D simulations

### Basic set-ups
* `0D/basic/basic_collisions_multispecies.jl` - initializes 2 species at different densities and temperature, performs collisions for multiple
time-steps
* `0D/basic/basic_sampling_and_io.jl` - samples particles, computes physical properties and writes them to disk
* `0D/basic/sample_and_merge.jl` - samples particles and merges them immediately, outputting various statistics

### BKW relaxation
* `0D/BKW/bkw_varweight_nnls.jl` - BKW simulation using variable-weight particles and grid-based merging
* `0D/BKW/bkw_varweight_nnls.jl` - BKW simulation using variable-weight particles and NNLS merging
* `0D/BKW/bkw_varweight_octree_swpm.jl` - BKW simulation using variable-weight particles and octree N:2 merging; uses the SWPM collision algorithm.
* `0D/BKW/bkw_varweight_octree_timing.jl` - BKW simulation using variable-weight particles and octree N:2 merging; uses `TimerOutputs.jl` for
timing of various parts of the simulation
* `0D/BKW/bkw_varweight_octree.jl` - BKW simulation using variable-weight particles and octree N:2 merging
* `0D/BKW/bkw.jl` - fixed-weight BKW simulation

### Ionization
* `0D/ionization/0D_ionization_1neutralspecies_es.jl` - single neutral species with elastic and ionizing neutral-electron collisions,
accelerated by a constant electric field. Can optionally use Event Splitting. Uses octree N:2 merging.
* `0D/ionization/0D_ionization_1neutralspecies_nnls_es.jl` - same as above, but uses NNLS merging (either moment-preserving,
moment-and-approximate-rate-preserving, or moment-and-exact-rate-preserving).

## 1D simulations
* `1D/couette_benchmarking.jl` - Couette flow simulation with fixed-weight particles; uses `TimerOutputs.jl` for
timing of various parts of the simulation
* `1D/couette_fp.jl` - Couette flow simulation with fixed-weight particles; uses the linear Fokker-Planck model for particle collisions
* `1D/couette_multithreaded_varweight_octree.jl` - multithreaded Couette flow simulation with variable weight particles and octree N:2 merging
* `1D/couette_multithreaded.jl` - multithreaded Couette flow simulation with fixed-weight particles
* `1D/couette_varweight_nnls.jl` - Couette flow simulation with variable weight particles and NNLS merging
* `1D/couette_varweight_octree_swpm.jl` - Couette flow simulation with variable weight particles and octree N:2 merging; uses the SWPM collision algorithm
* `1D/couette_varweight_octree.jl` -  Couette flow simulation with variable weight particles and octree N:2 merging
* `1D/couette_with_surface_quantities.jl` - Couette flow simulation with fixed-weight particles and surface properties output
* `1D/fourier_varweight_fw.jl` - Fourier flow simulation with fixed-weight particles
* `1D/fourier_varweight_nnls.jl` - Fourier flow simulation with variable weight particles and NNLS merging
* `1D/fourier_varweight_octree.jl` - Fourier flow simulation with variable weight particles and octree N:2 merging

# Paper reproducibility
The simulation files required to reproduce results from specific papers that made use of `Merzbild.jl` are listed below, along with
the specific version of `Merzbild.jl` used to perform the simulations.

## "Moment-preserving particle merging via non-negative least squares" (2026)
For the paper ["Moment-preserving particle merging via non-negative least squares"](TODO: link) by G. Oblapenko and M. Torrilhon,
the following simulation files were used. `Merzbild.jl` version `0.7.8`, Julia version `1.12`. The simulation files are kept up-to-date
and whilst the results are subject to change between `Merzbild.jl` and Julia versions (changes in random number generators,
order of operations affecting round-off errors, etc.), they should be runnable for each release.
Python scripts to postprocess results to produce plots
are available in `scripts/reproducibility/2026_nnls`: `process_bkw.py`, `process_ionization_compute_rates.py`,
`process_ionization.py`, `process_fourier.py`.
Details on executing the scripts can be found at the start of the script files and in the `scripts/reproducibility/2026_nnls/README.md`
file.
The parameter sets used for the paper are commented out (as running over all parameter sets takes a while); the simulation files
have notes on which commented-out-lines need to be uncommented to run the full simulations used to produce the data for the paper.

For the sampling test case:
* `0D/basic/sample_and_merge.jl` - for the simulations that sample particles and merge them (using two different sampling
strategies)

For the BKW test case (see the commented-out part on "multiple runs with ensembling if needed" for running with multiple ensembles):
* `0D/BKW/bkw_varweight_octree.jl` - for the variable-weight simulations using octree N:2 merging
* `0D/BKW/bkw_varweight_nnls.jl` - for the variable-weight simulations using NNLS merging

For the ionization test case (**when run over all ensembles, these produce VERY LARGE amounts of data, 100s of GBs**);
external data from [LXCat](https://us.lxcat.net/home/) is required for the cross-sections (IST-Lisbon database):
* `0D/ionization/0D_ionization_1neutralspecies_es.jl` - for the variable-weight simulations using octree N:2 merging
(uncomment lines below comment "Uncomment set-up below ..." to get the full set-up running over all parameters and ensembles; setting `paramset` to [12000, 6000] will produce results used as reference values)
* `0D/ionization/0D_ionization_1neutralspecies_nnls_es.jl` - for the variable-weight simulations using different versions of NNLS merging
(uncomment lines below comment "Uncomment set-up below ..." to get the full set-up running over all parameters and ensembles)

For the Fourier flow test case:
* `1D/fourier_fw.jl` - for the reference fixed-weight simulations
* `1D/fourier_varweight_octree.jl` - for the variable-weight simulations using octree N:2 merging (the commented `const params = ...` line
contains all the different target and threshold particle numbers used in the simulation)
* `1D/fourier_varweight_nnls.jl` - for the variable-weight simulations using NNLS merging

### Reference values
Reference values that can be expected to be produced by some of the simulations are provided below. `Merzbild.jl` version `0.7.8`, Julia version `1.12`.

Running `0D/basic/sample_and_merge.jl` with `n_t = 10000` and `params = [[4, 36]]` produces the following output (rounded to third decimal place):
```bash
NNLS merging: [4, 36]; sampling method: equal_weight
10000/10000
Npost = 35.0
weight ratio: 2591.674
weight std: 0.0309
weight log std: 1.427
f_tail(500.0): 0.261 -> 0.251
f_tail(750.0): 0.029 -> 0.030

Octree merging: [4, 36]; sampling method: equal_weight
10000/10000
Npost = 26.856
weight ratio: 40.835
weight std: 0.028
weight log std: 1.225
f_tail(500.0): 0.261 -> 0.292
f_tail(750.0): 0.029 -> 0.003

NNLS merging: [4, 36]; sampling method: weighted_samples
10000/10000
Npost = 35.0
weight ratio: 79954.614
weight std: 0.041
weight log std: 2.424
f_tail(500.0): 0.276 -> 0.249
f_tail(750.0): 0.031 -> 0.021

Octree merging: [4, 36]; sampling method: weighted_samples
10000/10000
Npost = 29.9759
weight ratio: 1.015e15
weight std: 0.043
weight log std: 6.320
f_tail(500.0): 0.276 -> 0.272
f_tail(750.0): 0.031 -> 0.018
```

Running `0D/ionization/0D_ionization_1neutralspecies_es.jl` with `paramset = [12000, 6000]`, `external_E_field_Tn = 400.0`, and `n_t = 500000`,
and then running `compute_rate.py` on the two results (with and without event splitting) will give the following output (rounded to second decimal place):
* `python3 scripts/compute_rate.py --filename scratch/data/ionization_Ar_400Tn_octree_mid_12000_to_6000_es.nc --tstart 150000 --tend 500000 --nid 0 --ionid 1 --eid 2 --dt 5e-14`: `4.46e-15 6.80`
* `python3 scripts/compute_rate.py --filename scratch/data/ionization_Ar_400Tn_octree_mid_12000_to_6000.nc --tstart 150000 --tend 500000 --nid 0 --ionid 1 --eid 2 --dt 5e-14`: `4.47e-15 6.82`