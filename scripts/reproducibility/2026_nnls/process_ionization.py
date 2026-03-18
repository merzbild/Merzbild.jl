import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from netCDF4 import Dataset

savefigs = False
# savefigs = True  # - uncomment to save plots

dt = 5e-14
field_Tn = 400

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

# plotting parameters
label_size = 24
tick_size = 20
legend_size = 20

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["axes.linewidth"] = 0.8

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

    print(f"Processing data for E = {field_Tn}Tn; averaging over {ts_max-ts_min} timesteps")

    octree_data[field_Tn] = {x[0]: {} for x in octree_runs}
    for ns, run in zip(n_seeds_for_run, octree_runs):
        print("Octree: ", run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_octree_mid_{run[0]}_to_{run[1]}_es"
        avg_of_bias, bias_of_avg, mk, nk, npm, Te, nTe = get_bias_mean_noise_np_from_rate_file(ref_val_mean,
                                                                                               ts_min, ts_max,
                                                                                               filename,
                                                                                               nseeds=ns)
        octree_data[field_Tn][run[0]] = {"avg_of_bias": avg_of_bias, "bias_of_avg": bias_of_avg, "mean": mk,
                                         "noise": nk, "np": npm, "T_e": Te, "std_T_e": nTe}

    nnls_data[field_Tn] = {x[0]: {} for x in nnls_runs}
    for ns, run in zip(n_seeds_for_run, nnls_runs):
        print("NNLS :", run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_NNLS_{run[0]}full_{run[1]}_es"
        avg_of_bias, bias_of_avg, mk, nk, npm, Te, nTe = get_bias_mean_noise_np_from_rate_file(ref_val_mean,
                                                                                               ts_min, ts_max, 
                                                                                               filename,
                                                                                               nseeds=ns)
        nnls_data[field_Tn][run[0]] = {"avg_of_bias": avg_of_bias, "bias_of_avg": bias_of_avg, "mean": mk,
                                       "noise": nk, "np": npm, "T_e": Te, "std_T_e": nTe}

    nnls_rp_data[field_Tn] = {x[0]: {} for x in nnls_rp_runs}
    for ns, run in zip(n_seeds_for_run, nnls_rp_runs):
        print("NNLS (ARP): ", run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_NNLSrate_approx_{run[0]}full_{run[1]}_es"
        avg_of_bias, bias_of_avg, mk, nk, npm, Te, nTe = get_bias_mean_noise_np_from_rate_file(ref_val_mean,
                                                                                               ts_min, ts_max, 
                                                                                               filename,
                                                                                               nseeds=ns)
        nnls_rp_data[field_Tn][run[0]] = {"avg_of_bias": avg_of_bias, "bias_of_avg": bias_of_avg, "mean": mk,
                                          "noise": nk, "np": npm, "T_e": Te, "std_T_e": nTe}

    nnls_erp_data[field_Tn] = {x[0]: {} for x in nnls_erp_runs}
    for ns, run in zip(n_seeds_for_run, nnls_erp_runs):
        print("NNLS (RP): ", run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_NNLSrate_exact_{run[0]}full_{run[1]}_es"
        avg_of_bias, bias_of_avg, mk, nk, npm, Te, nTe = get_bias_mean_noise_np_from_rate_file(ref_val_mean,
                                                                                               ts_min, ts_max, 
                                                                                               filename,
                                                                                               nseeds=ns)
        nnls_erp_data[field_Tn][run[0]] = {"avg_of_bias": avg_of_bias, "bias_of_avg": bias_of_avg, "mean": mk,
                                           "noise": nk, "np": npm, "T_e": Te, "std_T_e": nTe}
        
def get_bias_of_avg_data(rundata, runs):
    xmean_vals_ = [rundata[run[0]]["np"] for run in runs]
    y_vals_ = np.asarray([rundata[run[0]][f"bias_of_avg"] / ref_val_mean * 100 for run in runs])
    
    return xmean_vals_, y_vals_

def get_avg_of_bias_data(rundata, runs):
    xmean_vals_ = [rundata[run[0]]["np"] for run in runs]
    y_vals_ = np.asarray([rundata[run[0]][f"avg_of_bias"] / ref_val_mean * 100 for run in runs])
    
    return xmean_vals_, y_vals_

def get_mean_rate(rundata, runs):
    xmean_vals_ = [rundata[run[0]]["np"] for run in runs]
    y_vals_ = np.asarray([rundata[run[0]][f"mean"] for run in runs])
    
    return xmean_vals_, y_vals_

def get_noise_data(rundata, runs):
    xmean_vals_ = [rundata[run[0]]["np"] for run in runs]
    y_vals_ = [rundata[run[0]][f"noise"] / rundata[run[0]][f"mean"] * 100 for run in runs]
    
    return xmean_vals_, y_vals_

def get_temperature_data(rundata, runs):
    xmean_vals_ = [rundata[run[0]]["np"] for run in runs]
    y_vals_ = np.asarray([rundata[run[0]][f"T_e"] / 11605 for run in runs])
    
    return xmean_vals_, y_vals_

# bias in rate vs number of particles
fig = plt.figure(figsize=(18,6))

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

for ax, field_Tn in zip([ax1, ax2], [100, 400]):
    xmean_vals, y_vals = get_avg_of_bias_data(octree_data[field_Tn], octree_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"Octree")

    xmean_vals, y_vals = get_avg_of_bias_data(nnls_data[field_Tn], nnls_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"NNLS")

    xmean_vals, y_vals = get_avg_of_bias_data(nnls_erp_data[field_Tn], nnls_erp_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"NNLS, RP")

    xmean_vals, y_vals = get_avg_of_bias_data(nnls_rp_data[field_Tn], nnls_rp_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"NNLS, ARP")


ax1.legend(fontsize=legend_size, framealpha=1.0, title="E = 100 Tn",
    title_fontsize=legend_size)

ax2.legend([],
    [],
    title="E = 400 Tn",
    framealpha=1.0,
    title_fontsize=legend_size
)

for ax in [ax1, ax2]:
    ax.grid()
    ax.tick_params(axis='both', labelsize=tick_size,)
    ax.set_xlabel(r"$\overline{N_p}$", fontsize=label_size)
ax1.set_ylabel(r"$\overline{\mathcal{B}}(k_{ion})$, \%", fontsize=label_size)

if savefigs:
    fig.savefig(f"ionization_bias_k_ion.pdf", bbox_inches="tight")


# noise in rate vs number of particles
fig = plt.figure(figsize=(18,6))

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

for ax, field_Tn in zip([ax1, ax2], [100, 400]):
    xmean_vals, y_vals = get_noise_data(octree_data[field_Tn], octree_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"Octree")

    xmean_vals, y_vals = get_noise_data(nnls_data[field_Tn], nnls_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"NNLS")

    xmean_vals, y_vals = get_noise_data(nnls_erp_data[field_Tn], nnls_erp_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"NNLS, RP")

    xmean_vals, y_vals = get_noise_data(nnls_rp_data[field_Tn], nnls_rp_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"NNLS, ARP")

ax1.legend(fontsize=legend_size, framealpha=1.0, title="E = 100 Tn",
           title_fontsize=legend_size)

ax2.legend([],
    [],
    title="E = 400 Tn",
    framealpha=1.0,
    title_fontsize=legend_size
)

for ax in [ax1, ax2]:
    ax.grid()
    ax.tick_params(axis='both', labelsize=tick_size,)
    ax.set_xlabel(r"$\overline{N_p}$", fontsize=label_size)
ax1.set_ylabel(r"$\overline{\mathcal{N}}(k_{ion})$, \%", fontsize=label_size)

if savefigs:
    fig.savefig(f"ionization_noise_k_ion.pdf", bbox_inches="tight")


# electron temperature vs number of particles
fig = plt.figure(figsize=(18,6))

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

for ax, field_Tn in zip([ax1, ax2], [100, 400]):
    ax.plot([20, 250], [ref_T_vals[field_Tn], ref_T_vals[field_Tn]], 'k', linewidth=2, label="Reference")
    
    xmean_vals, y_vals = get_temperature_data(octree_data[field_Tn], octree_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"Octree")

    xmean_vals, y_vals = get_temperature_data(nnls_data[field_Tn], nnls_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"NNLS")

    xmean_vals, y_vals = get_temperature_data(nnls_erp_data[field_Tn], nnls_erp_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"NNLS, RP")

    xmean_vals, y_vals = get_temperature_data(nnls_rp_data[field_Tn], nnls_rp_runs)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"NNLS, ARP")

    # print(ax.get_xlim()
    ax.set_xlim([20, 250])
    
ax1.legend(fontsize=legend_size, framealpha=1.0, title="E = 100 Tn",
    title_fontsize=legend_size)

ax2.legend([],
    [],
    title="E = 400 Tn",
    framealpha=1.0,
    title_fontsize=legend_size
)

for ax in [ax1, ax2]:
    # ax.set_xlim(ax.get_xlim())
    print(ax.get_xlim())
    ax.grid()
    ax.tick_params(axis='both', labelsize=tick_size,)
    ax.set_xlabel(r"$\overline{N_p}$", fontsize=label_size)
ax1.set_ylabel(r"$T_e$, eV", fontsize=label_size)

if savefigs:
    fig.savefig(f"ionization_t_electron.pdf", bbox_inches="tight")


# now process bias in first timesteps after a merging event
# process single file
def get_bias_post_merge_single(ref_rate_mean, ref_T_mean, start_t, end_t, fname, window_size):
    ds = Dataset(fname)
    
    k_ion = np.asarray(ds.variables["k_ion"])
    ts = np.asarray(range(1,len(k_ion)+1))
    npart = np.asarray(ds.variables["np_e"])
    Te = np.asarray(ds.variables["T_e"])
    
    
    merge_timesteps_mask = npart[1:] < npart[0:-1]
    
    
    merge_timesteps = ts[1:][merge_timesteps_mask]
    
    if len(merge_timesteps) == 0:
        print(f"0 for {fname}")
    
    ds.close()

    start_t = max(0, start_t)
    end_t = min(end_t, len(ts))

    merge_timesteps = np.array([x for x in merge_timesteps if (x >= ts_min) and (x <= ts_max)])
    
    if len(merge_timesteps) == 0:
        print(f"0 for {fname} after ts_min/max, {ts_min}, {ts_max}")
    
    n_merges = len(merge_timesteps)
    n_between = sum(merge_timesteps[1:] - merge_timesteps[:-1]) / (n_merges - 1)
    
    k_ion_post = np.array([np.mean(k_ion[t0+1:t0+1+window_size])/(1e15 * ref_rate_mean) for t0 in merge_timesteps])
    Te_post = np.array([np.mean(Te[t0+1:t0+1+window_size])/ref_T_mean for t0 in merge_timesteps])
    
    return np.mean(np.abs(k_ion_post - 1.0)), np.mean(np.abs(Te_post - 1.0)), n_merges, n_between

# process a group of files with different random seeds
def get_bias_post_merge(ref_rate_mean, ref_T_mean, start_t, end_t, fname, nseeds, window_size):
    bias_k_tmp, bias_T_tmp, n_merge_tmp, n_between_tmp = get_bias_post_merge_single(ref_rate_mean, ref_T_mean,
                                                                                    start_t, end_t,
                                                                                    fname + "_rate_data_only.nc",
                                                                                    window_size)
    
    for adds in range(nseeds):
        bias_k_tmp2, bias_T_tmp2, n_merge_tmp2, n_between_tmp2 = get_bias_post_merge_single(ref_rate_mean, ref_T_mean,
                                                                                            start_t, end_t,
                                                                                            fname + f"_seed{adds+1}_rate_data_only.nc",
                                                                                            window_size)
        bias_k_tmp += bias_k_tmp2
        bias_T_tmp += bias_T_tmp2
        n_merge_tmp += n_merge_tmp2
        n_between_tmp += n_between_tmp2
    
    bias_k_tmp /= (nseeds+1)
    bias_T_tmp /= (nseeds+1)
    n_merge_tmp /= (nseeds+1)
    n_between_tmp /= (nseeds+1)
    
    return bias_k_tmp, bias_T_tmp, n_merge_tmp, n_between_tmp


octree_data_window = {}
nnls_data_window = {}
nnls_rp_data_window = {}
nnls_erp_data_window = {}
ws = 50 # window size

for field_Tn in [100, 400]:
    print(f"Processing windowed data for E = {field_Tn}Tn")
    ref_val_mean = ref_vals[field_Tn]
    ref_T_val = ref_T_vals[field_Tn]
    t_min = times[field_Tn][0]
    t_max = times[field_Tn][1]
    
    ts_min = round(t_min / dt)
    ts_max = round(t_max / dt)

    octree_data_window[field_Tn] = {x[0]: {} for x in octree_runs}

    for ns, run in zip(n_seeds_for_run, octree_runs):
        print("Octree: ", run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_octree_mid_{run[0]}_to_{run[1]}_es"
        bias_k_w, bias_T_w, n_merge_avg, n_between = get_bias_post_merge(ref_val_mean, ref_T_val * 11605.0,
                                                                         ts_min, ts_max, 
                                                                         filename,
                                                                         ns, ws)
        octree_data_window[field_Tn][run[0]] = {"bias_k_w": bias_k_w, "bias_T_w": bias_T_w, "n_merge_avg": n_merge_avg, "n_between": n_between}

    nnls_data_window[field_Tn] = {x[0]: {} for x in nnls_runs}

    for ns, run in zip(n_seeds_for_run, nnls_runs):
        print("NNLS: ", run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_NNLS_{run[0]}full_{run[1]}_es"
        bias_k_w, bias_T_w, n_merge_avg, n_between = get_bias_post_merge(ref_val_mean, ref_T_val * 11605.0,
                                                                         ts_min, ts_max, 
                                                                         filename,
                                                                         ns, ws)
        nnls_data_window[field_Tn][run[0]] = {"bias_k_w": bias_k_w, "bias_T_w": bias_T_w, "n_merge_avg": n_merge_avg, "n_between": n_between}

    nnls_rp_data_window[field_Tn] = {x[0]: {} for x in nnls_rp_runs}

    for ns, run in zip(n_seeds_for_run, nnls_rp_runs):
        print("NNLS ARP: ", run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_NNLSrate_approx_{run[0]}full_{run[1]}_es"
        bias_k_w, bias_T_w, n_merge_avg, n_between = get_bias_post_merge(ref_val_mean, ref_T_val * 11605.0,
                                                                         ts_min, ts_max, 
                                                                         filename,
                                                                         ns, ws)
        nnls_rp_data_window[field_Tn][run[0]] = {"bias_k_w": bias_k_w, "bias_T_w": bias_T_w, "n_merge_avg": n_merge_avg, "n_between": n_between}
        
    nnls_erp_data_window[field_Tn] = {x[0]: {} for x in nnls_erp_runs}

    for ns, run in zip(n_seeds_for_run, nnls_erp_runs):
        print("NNLS RP: ", run)
        filename = f"{path_to_dir}ionization_Ar_{field_Tn}Tn_NNLSrate_exact_{run[0]}full_{run[1]}_es"
        bias_k_w, bias_T_w, n_merge_avg, n_between = get_bias_post_merge(ref_val_mean, ref_T_val * 11605.0,
                                                                         ts_min, ts_max, 
                                                                         filename,
                                                                         ns, ws)
        nnls_erp_data_window[field_Tn][run[0]] = {"bias_k_w": bias_k_w, "bias_T_w": bias_T_w, "n_merge_avg": n_merge_avg, "n_between": n_between}


def get_window_bias_k(rundata, runs):
    y_vals_ = np.asarray([rundata[run[0]][f"bias_k_w"] * 100 for run in runs])
    return y_vals_

def get_window_bias_T(rundata, runs):
    y_vals_ = np.asarray([rundata[run[0]][f"bias_T_w"] * 100 for run in runs])
    return y_vals_

# plot bias in rate in 50 steps after merging event
fig = plt.figure(figsize=(18,6))

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

for ax, field_Tn in zip([ax1, ax2], [100, 400]):
    xmean_vals, y_vals = get_avg_of_bias_data(octree_data[field_Tn], octree_runs)
    y_vals = get_window_bias_k(octree_data_window[field_Tn], octree_runs)
    
    xmean_vals = np.asarray(xmean_vals)
    y_vals = np.asarray(y_vals)
    mask = ~np.isnan(y_vals)
    
    ax.plot(xmean_vals[mask], y_vals[mask], '-o', linewidth=2, label=f"Octree")

    xmean_vals, y_vals = get_avg_of_bias_data(nnls_data[field_Tn], nnls_runs)
    y_vals = get_window_bias_k(nnls_data_window[field_Tn], nnls_runs)
    xmean_vals = np.asarray(xmean_vals)
    y_vals = np.asarray(y_vals)
    mask = ~np.isnan(y_vals)
    ax.plot(xmean_vals[mask], y_vals[mask], '-o', linewidth=2, label=f"NNLS")

    xmean_vals, y_vals = get_avg_of_bias_data(nnls_erp_data[field_Tn], nnls_erp_runs)
    y_vals = get_window_bias_k(nnls_erp_data_window[field_Tn], nnls_erp_runs)
    xmean_vals = np.asarray(xmean_vals)
    y_vals = np.asarray(y_vals)
    mask = ~np.isnan(y_vals)
    ax.plot(xmean_vals[mask], y_vals[mask], '-o', linewidth=2, label=f"NNLS, RP")

    xmean_vals, y_vals = get_avg_of_bias_data(nnls_rp_data[field_Tn], nnls_rp_runs)
    y_vals = get_window_bias_k(nnls_rp_data_window[field_Tn], nnls_rp_runs)
    xmean_vals = np.asarray(xmean_vals)
    y_vals = np.asarray(y_vals)
    mask = ~np.isnan(y_vals)
    ax.plot(xmean_vals[mask], y_vals[mask], '-o', linewidth=2, label=f"NNLS, ARP")


ax1.legend(fontsize=legend_size, framealpha=1.0, title="E = 100 Tn",
    title_fontsize=legend_size)

ax2.legend([],
    [],
    title="E = 400 Tn",
    framealpha=1.0,
    title_fontsize=legend_size
)

for ax in [ax1, ax2]:
    ax.grid()
    ax.tick_params(axis='both', labelsize=tick_size,)
    ax.set_xlabel(r"$\overline{N_p}$", fontsize=label_size)
ax1.set_ylabel(r"$\overline{\mathcal{B}}_{50}(k_{ion})$, \%", fontsize=label_size)

# ax1.text(x=150
if savefigs:
    fig.savefig(f"ionization_bias_k_ion_50.pdf", bbox_inches="tight")


# plot bias in temperature in 50 steps after merging event
fig = plt.figure(figsize=(18,6))

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

for ax, field_Tn in zip([ax1, ax2], [100, 400]):
    xmean_vals, y_vals = get_avg_of_bias_data(octree_data[field_Tn], octree_runs)
    y_vals = get_window_bias_T(octree_data_window[field_Tn], octree_runs)
    
    xmean_vals = np.asarray(xmean_vals)
    y_vals = np.asarray(y_vals)
    mask = ~np.isnan(y_vals)
    
    ax.plot(xmean_vals[mask], y_vals[mask], '-o', linewidth=2, label=f"Octree")

    xmean_vals, y_vals = get_avg_of_bias_data(nnls_data[field_Tn], nnls_runs)
    y_vals = get_window_bias_T(nnls_data_window[field_Tn], nnls_runs)
    xmean_vals = np.asarray(xmean_vals)
    y_vals = np.asarray(y_vals)
    mask = ~np.isnan(y_vals)
    ax.plot(xmean_vals[mask], y_vals[mask], '-o', linewidth=2, label=f"NNLS")

    xmean_vals, y_vals = get_avg_of_bias_data(nnls_erp_data[field_Tn], nnls_erp_runs)
    y_vals = get_window_bias_T(nnls_erp_data_window[field_Tn], nnls_erp_runs)
    xmean_vals = np.asarray(xmean_vals)
    y_vals = np.asarray(y_vals)
    mask = ~np.isnan(y_vals)
    ax.plot(xmean_vals[mask], y_vals[mask], '-o', linewidth=2, label=f"NNLS, RP")

    xmean_vals, y_vals = get_avg_of_bias_data(nnls_rp_data[field_Tn], nnls_rp_runs)
    y_vals = get_window_bias_T(nnls_rp_data_window[field_Tn], nnls_rp_runs)
    xmean_vals = np.asarray(xmean_vals)
    y_vals = np.asarray(y_vals)
    mask = ~np.isnan(y_vals)
    ax.plot(xmean_vals[mask], y_vals[mask], '-o', linewidth=2, label=f"NNLS, ARP")


ax1.legend(fontsize=legend_size, framealpha=1.0, title="E = 100 Tn",
    title_fontsize=legend_size)

ax2.legend([],
    [],
    title="E = 400 Tn",
    framealpha=1.0,
    title_fontsize=legend_size
)

for ax in [ax1, ax2]:
    ax.grid()
    ax.tick_params(axis='both', labelsize=tick_size,)
    ax.set_xlabel(r"$\overline{N_p}$", fontsize=label_size)
ax1.set_ylabel(r"$\overline{\mathcal{B}}_{50}(T_e)$, \%", fontsize=label_size)

if savefigs:
    fig.savefig(f"ionization_bias_T_ion_50.pdf", bbox_inches="tight")