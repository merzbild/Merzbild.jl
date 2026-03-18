import numpy as np
from netCDF4 import Dataset
from matplotlib import pyplot as plt
from scipy import constants
from scipy.special import gamma

# convert ionization results to a set of files containing instantaneous values of
# 1) ionization rate coefficients 2) electron temperature 3) number of electron particles
# that can then be used in post-processing/plotting scripts
# files are assumed to be named as {path_to_dir}/{fname}_seed{seed}.nc (if seed > 1)
# or {path_to_dir}/{fname}.nc (if seed == 0)
# first, a run_names list is constructed (one run_name per parameter set), and then iterated over
# a second list (n_seeds_for_run) contains the number of different starting random seeds for each corresponding
# run in run_names (if seeds = 0,...,N were used, n_seeds_for_run should be [N-1])
# output is written to {path_to_dir}/{fname}_seed{seed}_rate_data_only.nc

field_Tn = 400
dt = 5e-14

path_to_dir = "scratch/data/"


# octree simulation parameters, uncomment these and comment out the NNLS run_names to use the former
# octree_runs = [[41, 38], [62, 58], [95, 88], [131, 122], [178, 166], [236, 220]]
# n_seeds_for_run = [63, 63, 63, 15, 15, 15]
# run_names = [f"{path_to_dir}ionization_Ar_{field_Tn}Tn_octree_mid_{run[0]}_to_{run[1]}_es" for run in octree_runs]

# NNLS simulation parameters, uncomment these and comment out the octree run_names to use the former
# nnlstypes refers to whether rate preservation is exact, approximate, or turned off
nnls_runs = [[4, 41, 38], [5, 62, 58], [6, 95, 88], [7, 131, 122], [8, 178, 166], [9, 236, 220]] 
n_seeds_for_run = [63, 63, 63, 15, 15, 15]
nnlstypes = ["", "rate_exact", "rate_approx"]
run_names = [f"{path_to_dir}ionization_Ar_{field_Tn}Tn_NNLS{rp}_{run[0]}_to_{run[1]}_es" for run in nnls_runs for rp in nnlstypes]

print(run_names)

# for ns, run in zip(n_seeds_for_run, octree_runs):
#     post_process(f"/home/georgii/Data/Sciebo/PIC_DSMC/NNLS_Paper/ionization/Ar_{field_Tn}Tn_octree_mid_{run[0]}_to_{run[1]}", nseeds=ns)

# compute rate given an array of densities of shape [n_time, n_species, 1], the indices of the species, and array
# of time values
def get_rate(n, species_n, species_i, species_e, time):
    dn_ion = n[1:, species_i, 0] - n[:-1, species_i, 0]
    dt = time[1:] - time[:-1]
    
    dn_ion = dn_ion / (dt * n[1:, species_n, 0] * n[1:, species_e, 0])
    return dn_ion


# convert a single netCDF file with densities, etc. to a smaller file
# containing only a single dimension ("time")
def post_process_single_run(fname, species_e=2):
    ds = Dataset(fname + ".nc")
    
    ndens = np.asarray(ds.variables["ndens"])
    ts = np.asarray(ds.variables["timestep"])
    np_e = np.asarray(ds.variables["np"])[1:, species_e, 0]
    T = np.asarray(ds.variables["T"])[1:]
    
    ds.close()

    start_t = 0
    end_t = len(ts)

    k_ion = get_rate(ndens, 0, 1, 2, ts * dt) * 1e15
    
    rootgrp = Dataset(fname + "_rate_data_only.nc", "w")
    timedim = rootgrp.createDimension("time", end_t-1)
    
    v_rates = rootgrp.createVariable("k_ion", "f8", ("time",))
    v_np_e = rootgrp.createVariable("np_e", "f8", ("time",))
    v_T_e = rootgrp.createVariable("T_e", "f8", ("time",))
    
    v_rates[:] = k_ion
    v_np_e[:] = np_e
    v_T_e[:] = T[:,species_e,0]
    
    rootgrp.close()

# process a set of files that have a naming scheme: {fname}_seed{seed}.nc (or {fname}.nc in case seed=0)
def post_process(fname, nseeds):
    print(f"Processing {fname} with {nseeds} seeds")
    post_process_single_run(fname)
    
    for adds in range(nseeds):
        post_process_single_run(fname + f"_seed{adds+1}")

# for ns, run in zip(n_seeds_for_run, run_names):
#     post_process(run, nseeds=ns)