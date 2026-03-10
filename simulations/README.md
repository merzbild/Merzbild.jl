# Overview of example simulations
This file provides a short overview of the different simulations in `simulations`. Specific simulation setups that were used
for results in publications are found at the end of this file in the "Paper reproducibility" section.
Unless noted, simulations use the NTC collision algorithm (or the variable-weight version thereof).

## 0D simulations

### Basic set-ups
* `0D/basic/basic_collisions_multispecies.jl` - initializes 2 species at different densities and temperature, performs collisions for multiple
time-steps
* `0D/basic/basic_sampling_and_io.jl` - samples particles, computes physical properties and writes them to disk

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
the following simulation files were used. `Merzbild.jl` version `0.7.8`, Julia version `TODO`. The simulation files are kept up-to-date
and whilst the results are subject to change between `Merzbild.jl` and Julia versions (changes in random number generators,
order of operations affecting round-off errors, etc.), they should be runnable for each release.
Python scripts to postprocess results to produce plots
are available in `scripts/reproducibility/nnls_2026`: `process_bkw.py`, `process_ionization.py`, `process_fourier.py`.
Details on executing the scripts can be found at the start of the script files.

For the BKW test case (see the commented-out part on "multiple runs with ensembling if needed" for running with multiple ensembles):
* `0D/BKW/bkw.jl` - for the reference fixed-weight simulations
* `0D/BKW/bkw_varweight_octree.jl` - for the variable-weight simulations using octree N:2 merging
* `0D/BKW/bkw_varweight_nnls.jl` - for the variable-weight simulations using NNLS merging

For the ionization test case (**when run over all ensembles, these produce VERY LARGE amounts of data, 100s of GBs**);
external data from [LXCat](https://us.lxcat.net/home/) is required for the cross-sections (IST-Lisbon database):
* `0D/ionization/0D_ionization_1neutralspecies_es.jl` - for the variable-weight simulations using octree N:2 merging
* `0D/ionization/0D_ionization_1neutralspecies_nnls_es.jl` - for the variable-weight simulations using different versions of NNLS merging

For the Fourier flow test case:
* `1D/fourier_fw.jl` - for the reference fixed-weight simulations
* `1D/fourier_varweight_octree.jl` - for the variable-weight simulations using octree N:2 merging (the commented `const params = ...` line
contains all the different target and threshold particle numbers used in the simulation)
* `1D/fourier_varweight_nnls.jl` - for the variable-weight simulations using octree N:2 merging