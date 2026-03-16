# Scripts for post-processing of data for the paper "Moment-preserving particle merging via non-negative least squares" (2026)
The scripts in this directory are
for post-processing the numerical results for the paper ["Moment-preserving particle merging via non-negative least squares"](TODO: link)
by G. Oblapenko and M. Torrilhon.

The results were computed using `Merzbild.jl` version `0.7.8`, Julia version `1.12`. The scripts require `numpy`, `scipy`, `matplotlib`, netCDF4`.

## Post-processing of BKW results
The numerical results produced by `bkw_varweight_octree.jl` and `bkw_varweight_nnls.jl` can be post-processed
with the `process_bkw.py` script. The script loads the files automatically: output files are assumed to be located in
`scratch/data`, octree results are assumed to be named `octree_mid_{threshold}_{Ntarget}_{seed}.nc`, 
NNLS results are assumed to be named `bkw_nnls_$(n_full_up_to_total)full_$(threshold)_$(seed).nc`.

The parameters of the simulations (number of random seeds, timestep, etc.) are set at the top of the script after the imports;
then plotting parameters (font sizes, font families) are set. If `savefigs` is set to true, the produced figures will be saved as PDFs.

The file produces 3 plots.