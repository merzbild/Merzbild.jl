# Scripts for post-processing of data for the paper "Moment-preserving particle merging via non-negative least squares" (2026)
The scripts in this directory are
for post-processing the numerical results for the paper ["Moment-preserving particle merging via non-negative least squares"](https://doi.org/10.48550/arXiv.2604.00668)
by G. Oblapenko and M. Torrilhon. The scripts require data produced by the simulations (see the [Paper reproducibility](../../../PAPER_REPRODUCIBILITY.md)).
By default it is assumed that the outputs of the simulations are located in `scratch/data`, this can be adjusted by changing the value of the `pref` variable
at the start of the scripts. Change `savefigs` to `False` to turn off saving figures as PDFs.

The results were computed using `Merzbild.jl` version `0.7.8`, Julia version `1.12`. The Python
scripts require `numpy`, `scipy`, `matplotlib`, `netCDF4`.

## Post-processing of sample-and-merge results

The numerical results produced by `sample_and_merge.jl` can be post-processed
with the `process_sample_and_merge.py` script.
The script loads the files automatically: output files are assumed to be located in
`scratch/data`, octree results are assumed to be named `octree_equalweight.log` and `octree_weighted.log`;
NNLS results are assumed to be named `nnls_equalweight.log` and `octree_weighted.log`.

The plotting parameters (font sizes, font families) are set at the top of the script after the imports.

The file produces 6 plots.

## Post-processing of BKW results
The numerical results produced by `bkw_varweight_octree.jl` and `bkw_varweight_nnls.jl` can be post-processed
with the `process_bkw.py` script. The script loads the files automatically: output files are assumed to be located in
`scratch/data`, octree results are assumed to be named `bkw_octree_mid_{threshold}_{Ntarget}_{seed}.nc`, 
NNLS results are assumed to be named `bkw_nnls_$(n_full_up_to_total)full_$(threshold)_$(seed).nc`.

The parameters of the simulations (number of random seeds, timestep, etc.) are set at the top of the script after the imports;
then plotting parameters (font sizes, font families) are set. If `savefigs` is set to true, the produced figures will be saved as PDFs.

The file produces 3 plots.

## Post-processing of 0D ionization results
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

## Post-processing of 1D Fourier results
The numerical results produced by `fourier_fw.jl`, `fourier_varweight_nnls.jl` and `fourier_varweight_octree.jl`
can be processed with the `process_fourier.py` script. The script loads the files automatically: output files are assumed to be located in
`scratch/data`, octree results are assumed to be named `avg_Fourier_octree_ntc_0.005_1000_0.0_300.0_600.0_{threshold}_{Ntarget}_after500000.nc`,
NNLS results are assumed to be named `avg_Fourier_nnls_ntc_0.005_1000_0.0_300.0_600.0_{threshold}_{Ntarget}_after500000.nc`.

The parameters of the simulations (size of cells, etc.) are set at the top of the script after the imports;
then plotting parameters (font sizes, font families) are set. If `savefigs` is set to true, the produced figures will be saved as PDFs.

The file produces 5 plots.