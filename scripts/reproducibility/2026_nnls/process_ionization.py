import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from netCDF4 import Dataset
from math import floor
import json

dt = 5e-14
field_Tn = 100


path_to_dir = "scratch/data/"

# assumed to be equal across all runs
n_seeds_for_run = [63, 63, 63, 15, 15, 15]

# parameters of runs with octree merging
octree_runs = [[41, 38], [62, 58], [95, 88], [131, 122], [178, 166], [236, 220]]

# parameters of runs with nnls merging without rate preservation
nnls_runs = [[4, 41, 38], [5, 62, 58], [6, 95, 88], [7, 131, 122], [8, 178, 166], [9, 236, 220]]

# parameters of runs with nnls merging with approximate rate preservation
nnls_rp_runs = [[4, 41, 38], [5, 62, 58], [6, 95, 88], [7, 131, 122], [8, 178, 166], [9, 236, 220]]

# parameters of runs with nnls merging with exact rate preservation
nnls_erp_runs = [[4, 41, 38], [5, 62, 58], [6, 95, 88], [7, 131, 122], [8, 178, 166], [9, 236, 220]]


# in which time window is averaging performed for the two different field strengths
times = {400: (0.75e-8, 2.5e-8), 100: (0.75e-7, 2.55e-7)}

# reference values for rate and temperature (eV) from octree simulations for the two different field strengths
ref_vals = {400: 4.461519882565042e-15, 100: 3.7675363518883525e-16}
ref_T_vals = {400: 6.7947153861840635, 100: 4.956735}

ref_val_mean = ref_vals[field_Tn]
ref_T_val = ref_T_vals[field_Tn]
t_min = times[field_Tn][0]
t_max = times[field_Tn][1]
ts_min = round(t_min / dt)
ts_max = round(t_max / dt)
print(f"Averaging over {ts_max-ts_min} timesteps")

octree_data = {}
nnls_data = {}
nnls_rp_data = {}
nnls_erp_data = {}

# process files with ionization rate data produced by convert_ionization_data.py
def get_bias_mean_noise_np_from_rate_file(ref_rate, start_t, end_t, fname, nseeds):
    ds = Dataset(fname + "_rate_data_only.nc")
    
    ts = ds.dimensions["time"].size

    start_t = max(0, start_t)
    end_t = min(end_t, ts)
    
    k_ion = np.asarray(ds.variables["k_ion"])[start_t:end_t]  # already * 1e15    
    T_e = np.asarray(ds.variables["T_e"])[start_t:end_t]
    npmean = np.mean(np.asarray(ds.variables["np_e"])[start_t:end_t])
    ds.close()

    noise_k_tmp = np.std(k_ion)
    avg_of_bias = abs(np.mean(k_ion) - ref_rate * 1e15)
    
    for adds in range(nseeds):
        ds = Dataset(fname + f"_seed{adds+1}_rate_data_only.nc")
        k_ion_tmp = np.asarray(ds.variables["k_ion"])[start_t:end_t]  # already * 1e15    
        npmean += np.mean(np.asarray(ds.variables["np_e"])[start_t:end_t])
        T_e += np.asarray(ds.variables["T_e"])[start_t:end_t]
        ds.close()

        k_ion += k_ion_tmp
        noise_k_tmp += np.std(k_ion_tmp)
        avg_of_bias += abs(np.mean(k_ion_tmp) - ref_rate * 1e15)
    
    k_ion /= (nseeds+1)
    noise_k_tmp /= (nseeds+1)
    avg_of_bias /= (nseeds+1)
    npmean /= (nseeds+1)
    T_e /= (nseeds+1)
    
    mean_k_tmp = np.mean(k_ion)
    bias_of_avg = abs(mean_k_tmp - ref_rate * 1e15)
    
    return avg_of_bias / 1e15, bias_of_avg / 1e15, mean_k_tmp / 1e15, noise_k_tmp / 1e15, npmean, np.mean(T_e), np.std(T_e)

for field_Tn in [100, 400]:
    ref_val_mean = ref_vals[field_Tn]
    ref_T_val = ref_T_vals[field_Tn]
    t_min = times[field_Tn][0]
    t_max = times[field_Tn][1]
    
    ts_min = round(t_min / dt)
    ts_max = round(t_max / dt)

    octree_data[field_Tn] = {x[0]: {} for x in octree_runs}
    for ns, run in zip(n_seeds_for_run, octree_runs):
        print(run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_octree_mid_{run[0]}_to_{run[1]}"
        avg_of_bias, bias_of_avg, mk, nk, npm, Te, nTe = get_bias_mean_noise_np_from_rate_file(ref_val_mean,
                                                                                               ts_min, ts_max,
                                                                                               filename,
                                                                                               nseeds=ns)
        octree_data[field_Tn][run[0]] = {"avg_of_bias": avg_of_bias, "bias_of_avg": bias_of_avg, "mean": mk,
                                         "noise": nk, "np": npm, "T_e": Te, "std_T_e": nTe}

    nnls_data[field_Tn] = {x[0]: {} for x in nnls_runs}
    for ns, run in zip(n_seeds_for_run, nnls_runs):
        print(run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_NNLS_{run[0]}full_{run[1]}"
        avg_of_bias, bias_of_avg, mk, nk, npm, Te, nTe = get_bias_mean_noise_np_from_rate_file(ref_val_mean,
                                                                                               ts_min, ts_max, 
                                                                                               filename,
                                                                                               nseeds=ns)
        nnls_data[field_Tn][run[0]] = {"avg_of_bias": avg_of_bias, "bias_of_avg": bias_of_avg, "mean": mk,
                                       "noise": nk, "np": npm, "T_e": Te, "std_T_e": nTe}

    nnls_rp_data[field_Tn] = {x[0]: {} for x in nnls_rp_runs}
    for ns, run in zip(n_seeds_for_run, nnls_rp_runs):
        print(run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_NNLS_approx_{run[0]}full_{run[1]}"
        avg_of_bias, bias_of_avg, mk, nk, npm, Te, nTe = get_bias_mean_noise_np_from_rate_file(ref_val_mean,
                                                                                               ts_min, ts_max, 
                                                                                               filename,
                                                                                               nseeds=ns)
        nnls_rp_data[field_Tn][run[0]] = {"avg_of_bias": avg_of_bias, "bias_of_avg": bias_of_avg, "mean": mk,
                                          "noise": nk, "np": npm, "T_e": Te, "std_T_e": nTe}

    nnls_erp_data[field_Tn] = {x[0]: {} for x in nnls_erp_runs}
    for ns, run in zip(n_seeds_for_run, nnls_erp_runs):
        print(run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_NNLS_exact_{run[0]}full_{run[1]}"
        avg_of_bias, bias_of_avg, mk, nk, npm, Te, nTe = get_bias_mean_noise_np_from_rate_file(ref_val_mean,
                                                                                               ts_min, ts_max, 
                                                                                               filename,
                                                                                               nseeds=ns)
        nnls_erp_data[field_Tn][run[0]] = {"avg_of_bias": avg_of_bias, "bias_of_avg": bias_of_avg, "mean": mk,
                                           "noise": nk, "np": npm, "T_e": Te, "std_T_e": nTe}