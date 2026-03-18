# Scripts for post-processing of data for the paper "Moment-preserving particle merging via non-negative least squares" (2026)
The scripts in this directory are
for post-processing the numerical results for the paper ["Moment-preserving particle merging via non-negative least squares"](TODO: link)
by G. Oblapenko and M. Torrilhon.

The results were computed using `Merzbild.jl` version `0.7.8`, Julia version `1.12`. The Python
scripts require `numpy`, `scipy`, `matplotlib`, `netCDF4`.

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
are multiplied by 1e15.
* `process_ionization.py`: takes converted files and produces plot

## Post-processing of 1D Fourier results
The numerical results produced by `fourier_fw.jl`, `fourier_varweight_nnls.jl` and `fourier_varweight_octree.jl`
can be processed with the `process_fourier.py` script. The script loads the files automatically: output files are assumed to be located in
`scratch/data`, octree results are assumed to be named `avg_Fourier_octree_ntc_0.005_1000_0.0_300.0_600.0_{threshold}_{Ntarget}_after500000.nc`,
NNLS results are assumed to be named `avg_Fourier_nnls_ntc_0.005_1000_0.0_300.0_600.0_{threshold}_{Ntarget}_after500000.nc`.

The parameters of the simulations (size of cells, etc.) are set at the top of the script after the imports;
then plotting parameters (font sizes, font families) are set. If `savefigs` is set to true, the produced figures will be saved as PDFs.
