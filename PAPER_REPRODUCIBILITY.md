
# Paper reproducibility
Instruction on reproducing results from specific papers that made use of `Merzbild.jl` are listed below, along with
the specific version of `Merzbild.jl` used to perform the simulations.

## "Moment-preserving particle merging via non-negative least squares" (2026)
For the paper ["Moment-preserving particle merging via non-negative least squares"](TODO: link) by G. Oblapenko and M. Torrilhon,
the following setup was used: `Merzbild.jl` version `0.7.8`, Julia version `1.12`.
First, one should produce the data by running the various simulation files located in the `simulations` directory,
then the output data can be processed via use of scripts in the `scripts/reproducibility/2026_nnls` directory.
The data produced by `Merzbild.jl` can also be downloaded from Zenodo: [10.5281/zenodo.19352779](https://doi.org/10.5281/zenodo.19352779) (note that the ionization data
is already post-processed to save space).

### Producing results from Merzbild.jl
In order to produce the data for subsequent post-processing, one should
1. Edit the simulation files listed below to use the parameters listed in the paper (the files have these parameters
in commented-out lines)
2. Download LXCat data for argon-neutral electron cross-sections (elastic scattering & ionization) in XML format,
see subsection below on "LXCat data" for notes on expected format of the file
3. By default, the LXCat data is assumed to be located in `../../Data/cross_sections/Ar_IST_Lisbon.xml` (relative to
the root directory of `Merzbild.jl`); the output data is written by default to `scratch/data` (relative to the root
directory of `Merzbild.jl`), therefore the directory should be created before running the simulations, or the output
path in the simulations changed.

Specific simulation files used for the paper:
For the sampling test case:
* `simulations/0D/basic/sample_and_merge.jl` - for the simulations that sample particles and merge them (using two different sampling
strategies)

For the BKW test case (see the commented-out part on "multiple runs with ensembling if needed" for running with multiple ensembles):
* `simulations/0D/BKW/bkw_varweight_octree.jl` - for the variable-weight simulations using octree N:2 merging
* `simulations/0D/BKW/bkw_varweight_nnls.jl` - for the variable-weight simulations using NNLS merging

For the ionization test case (**when run over all ensembles, these produce VERY LARGE amounts of data, 100s of GBs**);
external data from [LXCat](https://us.lxcat.net/home/) is required for the cross-sections (IST-Lisbon database) (see notes above on exact XML format required):
* `simulations/0D/ionization/0D_ionization_1neutralspecies_es.jl` - for the variable-weight simulations using octree N:2 merging
(uncomment lines below comment "Uncomment set-up below ..." to get the full set-up running over all parameters and ensembles; setting `paramset` to [12000, 6000] will produce results used as reference values)
* `simulations/0D/ionization/0D_ionization_1neutralspecies_nnls_es.jl` - for the variable-weight simulations using different versions of NNLS merging
(uncomment lines below comment "Uncomment set-up below ..." to get the full set-up running over all parameters and ensembles)

For the Fourier flow test case:
* `simulations/1D/fourier_fw.jl` - for the reference fixed-weight simulations
* `simulations/1D/fourier_varweight_octree.jl` - for the variable-weight simulations using octree N:2 merging (the commented `const params = ...` line
contains all the different target and threshold particle numbers used in the simulation)
* `simulations/1D/fourier_varweight_nnls.jl` - for the variable-weight simulations using NNLS merging

### Post-processing of results
Once the simulations have finished, Python scripts that post-process the results and produce plots should be run;
these are available in `scripts/reproducibility/2026_nnls`: `process_sample_and_merge.py`, `process_bkw.py`, `convert_ionization_data.py`, `process_ionization.py`, `process_fourier.py`.
By default it is assumed that the outputs of the simulations are located in `scratch/data`, this can be adjusted by changing the value of the `pref` variable
at the start of the scripts. Change `savefigs` to `False` to turn off saving figures as PDFs.

The results were computed using `Merzbild.jl` version `0.7.8`, Julia version `1.12`. The Python
scripts require `numpy`, `scipy`, `matplotlib`, `netCDF4`.

#### Post-processing of sample-and-merge results
The numerical results produced by `sample_and_merge.jl` can be post-processed
with the `process_sample_and_merge.py` script.
The script loads the files automatically: output files are assumed to be located in
`scratch/data`, octree results are assumed to be named `octree_equalweight.log` and `octree_weighted.log`;
NNLS results are assumed to be named `nnls_equalweight.log` and `octree_weighted.log`.

The plotting parameters (font sizes, font families) are set at the top of the script after the imports.

The file produces 6 plots.

#### Post-processing of BKW results
The numerical results produced by `bkw_varweight_octree.jl` and `bkw_varweight_nnls.jl` can be post-processed
with the `process_bkw.py` script. The script loads the files automatically: output files are assumed to be located in
`scratch/data`, octree results are assumed to be named `bkw_octree_mid_{threshold}_{Ntarget}_{seed}.nc`, 
NNLS results are assumed to be named `bkw_nnls_$(n_full_up_to_total)full_$(threshold)_$(seed).nc`.

The parameters of the simulations (number of random seeds, timestep, etc.) are set at the top of the script after the imports;
then plotting parameters (font sizes, font families) are set. If `savefigs` is set to true, the produced figures will be saved as PDFs.

The file produces 3 plots.

#### Post-processing of 0D ionization results
The numerical results produced by `0D_ionization_1neutralspecies_es.jl` and `0D_ionization_1neutralspecies_nnls_es.jl`
can be processed by two scripts:

* `convert_ionization_data.py`: converts output of `Merzbild.jl` to leaner files containing only relevant data
(ionization rate coefficients, electron temperatures, number of electron particles). The parameters of the simulations
(number of random seeds, applied electric field strength, timestep, etc.) are set at the top of the script after the imports;
they are used to produce a `run_names_and_seeds` list: each element is a tuple, the first element of which
is the prefix of the output file to process, and the second is the number of random seeds the simulation was run with.
Uncomment different lines to process results produced by the octree or NNLS merging methods. Note that the rates in the output
are multiplied by 1e15. Depending on the number of simulations, this might take a while.
* `process_ionization.py`: takes converted files and produces 5 plots. The parameters of the simulations
are set at the top of the script after the imports; the field strengths are hard-coded (100 and 400 Tn).
Reference values are two `{field_strength: value}1 dictionaries (for ionization rate coefficient and electron temperature). 

#### Post-processing of 1D Fourier results
The numerical results produced by `fourier_fw.jl`, `fourier_varweight_nnls.jl` and `fourier_varweight_octree.jl`
can be processed with the `process_fourier.py` script. The script loads the files automatically: output files are assumed to be located in
`scratch/data`, octree results are assumed to be named `avg_Fourier_octree_ntc_0.005_1000_0.0_300.0_600.0_{threshold}_{Ntarget}_after500000.nc`,
NNLS results are assumed to be named `avg_Fourier_nnls_ntc_0.005_1000_0.0_300.0_600.0_{threshold}_{Ntarget}_after500000.nc`.

The parameters of the simulations (size of cells, etc.) are set at the top of the script after the imports;
then plotting parameters (font sizes, font families) are set. If `savefigs` is set to true, the produced figures will be saved as PDFs.

The file produces 5 plots.

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

### LXCat data
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