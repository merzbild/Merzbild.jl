# Overview of example simulations
This file provides a short overview of the different simulations in `simulations`. Specific simulation setups that were used
for results in publications are found at the end of this file in the "Paper reproducibility" section.
Unless noted, simulations use the NTC collision algorithm (or the variable-weight version thereof).
By default simulations output data into `scratch/data`, and `scratch` is included in `.gitignore`.
Some simulations require external data that is not provided with `Merzbild.jl`, see the section on "External data" below.

The root directory of `Merzbild.jl` has a `run_examples.py` script; calling that will iterate over all files in `simulations`
and run them, replacing `n_t = VALUE` with `n_t = 10` and erroring otherwise; so this will run any simulation for 10 timesteps only
just to check that they runs without errors. After executing a simulation, the script resets `n_t` to the original value.
If a simulation does not have a line `n_t = VALUE` or `const n_t = VALUE`, the script will throw an error.
The script has an optional parameter, `logdir` (defaults to `scratch/logs`), console output and errors are written to `logdir`.
Simulation data is written to `scratch/data`. If these two directories do not exist,
the script will ask whether to create them. If a simulation file requires external data, the script will check whether the data exists -
if not, script will ask user whether to skip the simulation or proceed with an error.

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
These simulations require external data from LXCat (see below).
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

# External data
Some external data is not provided directly with the code, either due to data licensing requirements or data size.
Below is an overview of what data is required.

## 0D ionization
The 0D ionization simulations require LXCat data for argon-electron collision cross-sections from the IST-Lisbon database.
The path to the file is set in the `cs_n_e_filepath` variable; the data is expected in XML format. Since LXCat export to XML
is not entirely reliable, certain parts of an XML file (like missing names and ids) might need to be fixed manually. The file should have the following structure (other databases,
groups, and fields may be present, but the ones listed below are required; `DataX` and `DataY` hold the actual data, here they are left empty):
```xml
<?xml version="1.0" encoding="UTF-8"?>
<lxcat version="1.1">
<References>
    <Source>LXCat, www.lxcat.net. Generated on ---. All rights reserved.</Source>
    <Reference>- IST-Lisbon database, www.lxcat.net, retrieved on ---</Reference>
</References>
<Database name="IST-Lisbon database" id="IST-Lisbon">
    <Groups>
    <Group id="Ar">
        <Processes>
            <Process class="Scattering Cross Sections" type="Elastic">
                <Species>
                    <Reactant>e</Reactant>
                    <Reactant>Ar</Reactant>
                    <Product> E</Product>
                    <Product>Ar</Product>
                </Species>
                
                <Reaction>E + Ar -&gt; E + Ar</Reaction>
                
                <Parameters>
                    <mM>1.371000e-5</mM>
                    <Parameter>complete set</Parameter>
                </Parameters>
                
                <DataX type="Energy" units="eV" size="60"></DataX>
                <DataY type="Cross section" units="m2" size="60"></DataY>
            </Process>
            
            <Process class="Scattering Cross Sections" type="Ionization">
                <Species>
                    <Reactant>e</Reactant>
                    <Reactant>Ar</Reactant>
                    <Product> E</Product>
                    <Product>E</Product>
                    <Product>Ar+</Product>
                </Species>
                
                <Reaction>E + Ar -&gt; E + E + Ar+</Reaction>
                
                <Parameters>
                    <E units="eV">1.576000e+1</E>
                    <Parameter>complete set</Parameter>
                </Parameters>

                <DataX type="Energy" units="eV" size="28"></DataX>
                <DataY type="Cross section" units="m2" size="28"></DataY>
            </Process>
        </Processes>
        </Group>
    </Groups>
</Database>
</lxcat>
```

# Paper reproducibility
The simulation files required to reproduce results from specific papers that made use of `Merzbild.jl` are listed below, along with
the specific version of `Merzbild.jl` used to perform the simulations.

## "Moment-preserving particle merging via non-negative least squares" (2026)
For the paper ["Moment-preserving particle merging via non-negative least squares"](TODO: link) by G. Oblapenko and M. Torrilhon,
the following simulation files were used. `Merzbild.jl` version `0.7.8`, Julia version `1.12`. The simulation files are kept up-to-date
and whilst the results are subject to change between `Merzbild.jl` and Julia versions (changes in random number generators,
order of operations affecting round-off errors, etc.), they should be runnable for each release. By default all of the simulations
write data to `scratch/data`, therefore the directory should be created before running the scripts.
Python scripts to postprocess results to produce plots
are available in `scripts/reproducibility/2026_nnls`: `process_sample_and_merge.py`, `process_bkw.py`, `convert_ionization_data.py`,
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
external data from [LXCat](https://us.lxcat.net/home/) is required for the cross-sections (IST-Lisbon database) (see notes above on exact XML format required):
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

Running `0D/basic/sample_and_merge.jl` with `n_t = 10000` and `params = [[4, 36]]` produces the following output (rounded to third decimal place); the output is also written to files (by default `scratch/data/nnls_{sm}.log` and `scratch/data/octree_{sm}.log`,
where `sm` is either `equal_weight` or `weighted_samples`):
```bash
NNLS, equal_weight
10000/10000
Npost = 35.0
pre-merge weight ratio: 1.0
pre-merge weight std: 1.301e-18
pre-merge weight log std: 5.151e-14
weight ratio: 2591.674
weight std: 0.031
weight log std: 1.427
f_tail(500.0): 0.261 -> 0.251
f_tail(750.0): 0.029 -> 0.030

Octree N:2, equal_weight
10000/10000
Npost = 26.8555
pre-merge weight ratio: 1.0
pre-merge weight std: 1.301e-18
pre-merge weight log std: 5.151e-14
weight ratio: 40.835
weight std: 0.028
weight log std: 1.225
f_tail(500.0): 0.261 -> 0.292
f_tail(750.0): 0.029 -> 0.003

NNLS, weighted_samples
10000/10000
Npost = 35.0
pre-merge weight ratio: 7.533e18
pre-merge weight std: 0.011
pre-merge weight log std: 8.25
weight ratio: 79954.614
weight std: 0.041
weight log std: 2.424
f_tail(500.0): 0.276 -> 0.249
f_tail(750.0): 0.031 -> 0.021

Octree N:2, weighted_samples
10000/10000
Npost = 29.976
pre-merge weight ratio: 7.468e18
pre-merge weight std: 0.011
pre-merge weight log std: 8.248
weight ratio: 1.015e15
weight std: 0.043
weight log std: 6.32
f_tail(500.0): 0.276 -> 0.272
f_tail(750.0): 0.031 -> 0.018
```

Running `0D/ionization/0D_ionization_1neutralspecies_es.jl` with `paramset = [12000, 6000]`, `external_E_field_Tn = 400.0`, and `n_t = 500000`,
and then running `compute_rate.py` on the two results (with and without event splitting) will give the following output (rounded to second decimal place):
* `python3 scripts/compute_rate.py --filename scratch/data/ionization_Ar_400Tn_octree_mid_12000_to_6000_es.nc --tstart 150000 --tend 500000 --nid 0 --ionid 1 --eid 2 --dt 5e-14`: `4.46e-15 6.80`
* `python3 scripts/compute_rate.py --filename scratch/data/ionization_Ar_400Tn_octree_mid_12000_to_6000.nc --tstart 150000 --tend 500000 --nid 0 --ionid 1 --eid 2 --dt 5e-14`: `4.47e-15 6.82`