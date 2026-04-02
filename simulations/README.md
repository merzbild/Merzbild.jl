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
For reproducing results from a specific paper via running the simulation files listed above, please refer
to [Paper reproducibility](../PAPER_REPRODUCIBILITY.md).
