import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from netCDF4 import Dataset

# uncomment savefigs = False to turn off saving of figures
savefigs = True
# savefigs = False

dt = 5e-14

pref = "scratch/data/"

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
#
# the bias is the signed deviation of the ENSEMBLE-MEAN rate from the reference.
# The absolute value has to come after the ensemble average, or not at all -- the signed
# value is plotted here so that the sign of the bias stays visible.
#
# bias_se is the standard error of the ensemble mean.
def get_bias_mean_noise_np_from_rate_file(ref_rate, start_t, end_t, fname, nseeds):
    seed_means = []  # window-mean rate, one entry per seed
    seed_stds = []   # temporal std of the rate within the window, one entry per seed
    npmean = 0.0
    T_e = None

    for adds in range(nseeds + 1):
        ds = Dataset(fname + (f"_seed{adds}" if adds > 0 else "") + "_rate_data_only.nc")

        if adds == 0:
            start_t = max(0, start_t)
            end_t = min(end_t, ds.dimensions["time"].size)

        k_ion = np.asarray(ds.variables["k_ion"])[start_t:end_t]  # already * 1e15
        T_seed = np.asarray(ds.variables["T_e"])[start_t:end_t]
        npmean += np.mean(np.asarray(ds.variables["np_e"])[start_t:end_t])
        ds.close()

        seed_means.append(np.mean(k_ion))
        seed_stds.append(np.std(k_ion))
        T_e = T_seed.copy() if T_e is None else T_e + T_seed

    nfiles = nseeds + 1
    seed_means = np.asarray(seed_means) / 1e15
    mean_k = np.mean(seed_means)
    T_e /= nfiles

    return {"bias_of_avg": mean_k - ref_rate,
            "bias_se": np.std(seed_means, ddof=1) / np.sqrt(nfiles) if nfiles > 1 else np.nan,
            "mean": mean_k,
            "noise": np.mean(seed_stds) / 1e15,
            "np": npmean / nfiles,
            "T_e": np.mean(T_e),
            "std_T_e": np.std(T_e)}

for field_Tn in [100, 400]:
    ref_val_mean = ref_vals[field_Tn]
    ref_T_val = ref_T_vals[field_Tn]
    t_min = times[field_Tn][0]
    t_max = times[field_Tn][1]
    
    ts_min = round(t_min / dt)
    ts_max = round(t_max / dt)

    print(f"Processing data for E = {field_Tn}Tn; averaging over {ts_max-ts_min} timesteps")

    for label, data, runs, tag in [("Octree",     octree_data,   octree_runs,   "octree_mid_{0}_to_{1}"),
                                   ("NNLS",       nnls_data,     nnls_runs,     "NNLS_{0}full_{1}"),
                                   ("NNLS (ARP)", nnls_rp_data,  nnls_rp_runs,  "NNLSrate_approx_{0}full_{1}"),
                                   ("NNLS (RP)",  nnls_erp_data, nnls_erp_runs, "NNLSrate_exact_{0}full_{1}")]:
        data[field_Tn] = {}
        for ns, run in zip(n_seeds_for_run, runs):
            print(f"{label}: ", run)
            filename = f"{pref}ionization_Ar_{field_Tn}Tn_" + tag.format(*run[:2]) + "_es"
            data[field_Tn][run[0]] = get_bias_mean_noise_np_from_rate_file(ref_val_mean,
                                                                          ts_min, ts_max,
                                                                          filename,
                                                                          nseeds=ns)


def get_np_data(rundata, runs):
    return np.asarray([rundata[run[0]]["np"] for run in runs])

# ref is passed in rather than read off the module-level ref_val_mean: that global is left holding
# the value of the last field strength processed, so both panels would be normalised by it
def get_bias_of_avg_data(rundata, runs, ref):
    xmean_vals_ = get_np_data(rundata, runs)
    y_vals_ = np.asarray([rundata[run[0]][f"bias_of_avg"] / ref * 100 for run in runs])
    y_err_ = np.asarray([rundata[run[0]][f"bias_se"] / ref * 100 for run in runs])

    return xmean_vals_, y_vals_, y_err_

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
    for data, runs, label in [(octree_data, octree_runs, "Octree"), (nnls_data, nnls_runs, "NNLS"),
                              (nnls_erp_data, nnls_erp_runs, "NNLS, RP"),
                              (nnls_rp_data, nnls_rp_runs, "NNLS, ARP")]:
        xmean_vals, y_vals, y_err = get_bias_of_avg_data(data[field_Tn], runs, ref_vals[field_Tn])
        ax.plot(xmean_vals, y_vals, marker='o', linewidth=2, label=label)

    ax.axhline(0.0, color='k', linewidth=0.8)


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
ax1.set_ylabel(r"$\mathcal{B}(k_{ion})$, \%", fontsize=label_size)

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


# now process the effect of a merging event on the following window_size timesteps
#
# indexing: a merge at t0 means the electron count dropped between array indices t0-2 and t0-1.
# np_e and T_e are written after merging, so t0-1 is the first post-merge sample of those, whereas
# k_ion[j] spans the interval (j -> j+1), so k_ion[t0-1] is the interval *containing* the merge and
# k_ion[t0] is the first interval sampled from the merged population.
#
# the estimator is a difference of differences,
#   (post-window - pre-window) at merge steps  -  (post-window - pre-window) at control steps
# and both halves of that are needed:
#
#  - signed rather than mean(|x - 1|): only ~1% of timesteps contain an ionization event at all, so
#    a short window mean is dominated by whether any event happened to land in it. mean(|x - 1|)
#    reports that shot noise and cannot average it away -- evaluated at randomly chosen timesteps it
#    comes out the same as at merge steps, i.e. it carries no merge signal whatsoever.
#
#  - control steps rather than arbitrary timesteps: merging triggers on the particle count crossing
#    the threshold, which only happens on a step that just ionized, so k_ion[t0-1] is selected for
#    being non-zero (it averages ~1/P(ionization) ~ 75x the mean rate) and it is excluded from both
#    windows. merges also cluster in locally hot stretches of the run, which lifts T_e symmetrically
#    on *both* sides of the merge. the controls -- steps that ionized but did not trigger a merge --
#    carry the same selection, so differencing them out isolates what the merge itself did.
#
# events whose windows would overlap another merge are dropped, so that one merge is not measuring
# the tail of the previous one and the controls stay genuinely merge-free.
def get_post_merge_effect_single(ref_rate_mean, ref_T_mean, start_t, end_t, fname, window_size):
    ds = Dataset(fname)
    k_ion = np.asarray(ds.variables["k_ion"]) / (1e15 * ref_rate_mean)
    npart = np.asarray(ds.variables["np_e"])
    Te = np.asarray(ds.variables["T_e"]) / ref_T_mean
    ds.close()

    n_ts = len(k_ion)
    ts = np.arange(1, n_ts + 1)

    start_t = max(0, start_t)
    end_t = min(end_t, n_ts)

    merges = ts[1:][npart[1:] < npart[:-1]]
    # k_ion[j] > 0 means ionization happened during the interval ending at step j+1
    ionized = np.flatnonzero(k_ion > 0.0) + 1
    controls = np.setdiff1d(ionized, merges)

    in_window = merges[(merges >= start_t) & (merges <= end_t)]
    n_merges = len(in_window)
    n_between = np.mean(np.diff(in_window)) if n_merges > 1 else np.nan

    if n_merges == 0:
        print(f"0 merges for {fname} in [{start_t}, {end_t}]")
        return np.nan, np.nan, n_merges, n_between

    def select(t0):
        # both windows must fit in the array, and t0 must lie in the averaging window
        t0 = t0[(t0 >= window_size + 1) & (t0 + window_size <= n_ts)
                & (t0 >= start_t) & (t0 <= end_t)]
        # no *other* merge within +-window_size (a merge point finds itself, hence the allowance)
        n_near = (np.searchsorted(merges, t0 + window_size, side="right")
                  - np.searchsorted(merges, t0 - window_size, side="left"))
        return t0[n_near <= np.where(np.isin(t0, merges), 1, 0)]

    # prefix sums so every window mean is O(1); the arrays are ~5e5 long and there are ~1e4 events
    cs_k = np.concatenate(([0.0], np.cumsum(k_ion)))
    cs_T = np.concatenate(([0.0], np.cumsum(Te)))

    def wmean(cs, a, b):
        return (cs[b] - cs[a]) / (b - a)

    def deltas(t0):
        # the pre-window ends just before the merge interval, the post-window starts at it
        d_k = wmean(cs_k, t0, t0 + window_size) - wmean(cs_k, t0 - 1 - window_size, t0 - 1)
        d_T = wmean(cs_T, t0 - 1, t0 - 1 + window_size) - wmean(cs_T, t0 - 1 - window_size, t0 - 1)
        return d_k, d_T

    merge_sel, control_sel = select(merges), select(controls)

    if len(merge_sel) == 0 or len(control_sel) == 0:
        print(f"no isolated merges/controls for {fname}")
        return np.nan, np.nan, n_merges, n_between

    d_k_merge, d_T_merge = deltas(merge_sel)
    d_k_control, d_T_control = deltas(control_sel)

    return (np.mean(d_k_merge) - np.mean(d_k_control),
            np.mean(d_T_merge) - np.mean(d_T_control),
            n_merges, n_between)

# process a group of files with different random seeds
def get_post_merge_effect(ref_rate_mean, ref_T_mean, start_t, end_t, fname, nseeds, window_size):
    per_seed = np.asarray([get_post_merge_effect_single(ref_rate_mean, ref_T_mean, start_t, end_t,
                                                        fname + (f"_seed{adds}" if adds > 0 else "")
                                                        + "_rate_data_only.nc",
                                                        window_size)
                           for adds in range(nseeds + 1)])

    nfiles = nseeds + 1
    sem = np.std(per_seed, axis=0, ddof=1) / np.sqrt(nfiles) if nfiles > 1 else np.full(4, np.nan)

    return {"effect_k": np.mean(per_seed[:, 0]), "effect_k_se": sem[0],
            "effect_T": np.mean(per_seed[:, 1]), "effect_T_se": sem[1],
            "n_merge_avg": np.mean(per_seed[:, 2]), "n_between": np.nanmean(per_seed[:, 3])}


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

    for label, data, runs, tag in [("Octree",   octree_data_window,   octree_runs,   "octree_mid_{0}_to_{1}"),
                                   ("NNLS",     nnls_data_window,     nnls_runs,     "NNLS_{0}full_{1}"),
                                   ("NNLS ARP", nnls_rp_data_window,  nnls_rp_runs,  "NNLSrate_approx_{0}full_{1}"),
                                   ("NNLS RP",  nnls_erp_data_window, nnls_erp_runs, "NNLSrate_exact_{0}full_{1}")]:
        data[field_Tn] = {}
        for ns, run in zip(n_seeds_for_run, runs):
            print(f"{label}: ", run)
            filename = f"{pref}ionization_Ar_{field_Tn}Tn_" + tag.format(*run[:2]) + "_es"
            data[field_Tn][run[0]] = get_post_merge_effect(ref_val_mean, ref_T_val * 11605.0,
                                                           ts_min, ts_max,
                                                           filename,
                                                           ns, ws)


def get_window_effect_k(rundata, runs):
    y_vals_ = np.asarray([rundata[run[0]][f"effect_k"] * 100 for run in runs])
    y_err_ = np.asarray([rundata[run[0]][f"effect_k_se"] * 100 for run in runs])
    return y_vals_, y_err_

def get_window_effect_T(rundata, runs):
    y_vals_ = np.asarray([rundata[run[0]][f"effect_T"] * 100 for run in runs])
    y_err_ = np.asarray([rundata[run[0]][f"effect_T_se"] * 100 for run in runs])
    return y_vals_, y_err_

# plot bias in rate in 50 steps after merging event
fig = plt.figure(figsize=(18,6))

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

for ax, field_Tn in zip([ax1, ax2], [100, 400]):
    for data, data_w, runs, label in [(octree_data, octree_data_window, octree_runs, "Octree"),
                                      (nnls_data, nnls_data_window, nnls_runs, "NNLS"),
                                      (nnls_erp_data, nnls_erp_data_window, nnls_erp_runs, "NNLS, RP"),
                                      (nnls_rp_data, nnls_rp_data_window, nnls_rp_runs, "NNLS, ARP")]:
        xmean_vals = get_np_data(data[field_Tn], runs)
        y_vals, y_err = get_window_effect_k(data_w[field_Tn], runs)
        mask = ~np.isnan(y_vals)

        ax.plot(xmean_vals[mask], y_vals[mask], marker='o', linewidth=2, label=label)

    ax.axhline(0.0, color='k', linewidth=0.8)


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
ax1.set_ylabel(r"$\Delta_{50}(k_{ion})$, \%", fontsize=label_size)

# ax1.text(x=150
if savefigs:
    fig.savefig(f"ionization_bias_k_ion_50.pdf", bbox_inches="tight")


# plot bias in temperature in 50 steps after merging event
fig = plt.figure(figsize=(18,6))

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

for ax, field_Tn in zip([ax1, ax2], [100, 400]):
    for data, data_w, runs, label in [(octree_data, octree_data_window, octree_runs, "Octree"),
                                      (nnls_data, nnls_data_window, nnls_runs, "NNLS"),
                                      (nnls_erp_data, nnls_erp_data_window, nnls_erp_runs, "NNLS, RP"),
                                      (nnls_rp_data, nnls_rp_data_window, nnls_rp_runs, "NNLS, ARP")]:
        xmean_vals = get_np_data(data[field_Tn], runs)
        y_vals, y_err = get_window_effect_T(data_w[field_Tn], runs)
        mask = ~np.isnan(y_vals)

        ax.plot(xmean_vals[mask], y_vals[mask], marker='o', linewidth=2, label=label)

    ax.axhline(0.0, color='k', linewidth=0.8)


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
ax1.set_ylabel(r"$\Delta_{50}(T_e)$, \%", fontsize=label_size)

if savefigs:
    fig.savefig(f"ionization_bias_T_ion_50.pdf", bbox_inches="tight")