import numpy as np
from netCDF4 import Dataset

# convert the tail diagnostics written by the _tail.jl simulations in simulations/0D/ionization
# into one small per-run summary file that process_ionization_tails.py can plot directly.
#
# the raw {fname}_tail.nc files hold ~19 values per timestep, i.e. ~76 MB per 500k-step run, and a
# full sweep is ~1000 runs; reducing each one to a few dozen numbers here means the plotting script
# no longer has to read ~90 GB every time a figure is tweaked.
#
# input, per run:
#   {pref}{fname}_tail.nc  -- tail_pre/tail_post, w_total_pre/post, w_ratio/sigma_w/sigma_logw
#                             pre and post, sigma_g_w_max, n_coll, n_eq_w_coll, merge_kind
#   {pref}{fname}.nc       -- only the np variable, for the mean electron particle count
# output:
#   {pref}{fname}_tail_summary.nc
#
# files are assumed to be named as {pref}{fname}_seed{seed}_tail.nc (if seed > 0)
# or {pref}{fname}_tail.nc (if seed == 0), as in convert_ionization_data.py: a run_names list is
# built (one run_name per parameter set) along with n_seeds_for_run, and the two are zipped into
# run_names_and_seeds, each element a (file prefix, number of extra seeds) tuple
#
# everything is reduced over a single averaging window, chosen so that the run has reached a
# quasi-steady state. the window is stored in the output as an attribute; change it here and re-run
# rather than trying to re-window downstream

field_Tn = 400
dt = 5e-14

# path to directory with simulation results
pref = "scratch/data2/"

# in which time window is averaging performed for the two different field strengths
times = {400: (0.75e-8, 2.5e-8), 100: (0.75e-7, 2.55e-7)}

# threshold = 1.075 * target
octree_runs = [[41, 38], [62, 58], [95, 88], [131, 122], [178, 166], [236, 220]]
nnls_runs = [[4, 41, 38], [5, 62, 58], [6, 95, 88], [7, 131, 122], [8, 178, 166], [9, 236, 220]]
n_seeds_for_run = [15, 15, 15, 15, 15, 15]

# threshold = 1.5 * target
# octree_runs = [[57, 38], [87, 58], [132, 88], [183, 122], [250, 166], [330, 220]]
# nnls_runs = [[4, 57, 38], [5, 87, 58], [6, 132, 88], [7, 183, 122], [8, 250, 166], [9, 330, 220]]
# n_seeds_for_run = [31, 31, 31, 7, 7, 7]

run_names_and_seeds = [(f"{pref}ionization_Ar_{field_Tn}Tn_octree_mid_{run[0]}_to_{run[1]}_es", ns)
                       for run, ns in zip(octree_runs, n_seeds_for_run)]

for rp in ["", "rate_exact", "rate_approx"]:
    run_names_and_seeds += [(f"{pref}ionization_Ar_{field_Tn}Tn_NNLS{rp}_{run[0]}full_{run[1]}_es", ns)
                            for run, ns in zip(nnls_runs, n_seeds_for_run)]

# the reference run: enough particles that merging barely acts. it merges rarely or not at all, so
# the merge-conditioned entries of its summary may come out as NaN, which is expected
run_names_and_seeds += [(f"{pref}ionization_Ar_{field_Tn}Tn_octree_mid_24000_to_12000_es", 7)]

print(run_names_and_seeds)

# quantities recorded pre- and post-merge alongside the tail functions
weight_names = ["w_ratio", "sigma_w", "sigma_logw"]

# quantities recorded once per timestep
step_names = ["sigma_g_w_max", "n_coll", "n_eq_w_coll"]


# reduce a single run to its summary values over [start_t, end_t)
def summarize_single_run(fname, start_t, end_t):
    ds = Dataset(fname + "_tail.nc")

    energies = np.asarray(ds.variables["cutoff_energy_eV"])
    start_t = max(0, start_t)
    end_t = min(end_t, ds.dimensions["time"].size)
    sl = slice(start_t, end_t)

    # stored as (cutoff, time) in the simulation, which reads back as (time, cutoff) here
    tail_pre = np.asarray(ds.variables["tail_pre"])[sl]
    tail_post = np.asarray(ds.variables["tail_post"])[sl]
    w_total_pre = np.asarray(ds.variables["w_total_pre"])[sl]
    w_total_post = np.asarray(ds.variables["w_total_post"])[sl]
    merge_kind = np.asarray(ds.variables["merge_kind"])[sl]

    wstats = {n: (np.asarray(ds.variables[n + "_pre"])[sl], np.asarray(ds.variables[n + "_post"])[sl])
              for n in weight_names}
    steps = {n: np.asarray(ds.variables[n])[sl] for n in step_names}
    ds.close()

    # the tail file does not carry the particle count, so take it from the main output; only the np
    # variable is read, not the whole file
    ds = Dataset(fname + ".nc")
    npmean = np.mean(np.asarray(ds.variables["np"])[1:, 2, 0][sl])
    ds.close()

    # the tail functions are only meaningful as fractions of the total weight: the electron density
    # grows by three orders of magnitude over a run, so the raw weights cannot be averaged over time
    frac_pre = tail_pre / w_total_pre[:, None]
    frac_post = tail_post / w_total_post[:, None]

    merged = merge_kind > 0
    n_merges = int(merged.sum())

    if n_merges == 0:
        print(f"  no merges in [{start_t}, {end_t}] for {fname}")

    def at_merge(a):
        return a[merged].mean(axis=0) if n_merges > 0 else np.full(a.shape[1:], np.nan)

    out = {"np": npmean,
           "n_merges": n_merges,
           # fraction of merges that fell back to a lower-order NNLS or to octree
           "fallback_frac": np.count_nonzero(merge_kind >= 2) / n_merges if n_merges > 0 else np.nan,
           "F": frac_pre.mean(axis=0),
           "F_merge_pre": at_merge(frac_pre),
           "F_merge_post": at_merge(frac_post)}

    for n in weight_names:
        pre, post = wstats[n]
        out[n] = pre.mean()
        out[n + "_merge_pre"] = at_merge(pre)
        out[n + "_merge_post"] = at_merge(post)

    for n in step_names:
        out[n] = steps[n].mean()

    return energies, out


# write the summary of a single run. the derived per-merge changes are stored alongside the raw
# merge-step means so that the file can be plotted as-is, and renormalised if needed.
#
# the change is normalised by the pre-merge value *at the merge steps*, not by the window mean:
# merges fire when the weight spread has grown most, so the two differ a lot for w_ratio and
# normalising by the window mean would put the change past -100%
def write_summary(fname, energies, summary, start_t, end_t):
    rootgrp = Dataset(fname + "_tail_summary.nc", "w")
    rootgrp.createDimension("cutoff", len(energies))

    rootgrp.start_t = start_t
    rootgrp.end_t = end_t
    rootgrp.COMMENT = ("per-run summary of {fname}_tail.nc over [start_t, end_t). F is the mean "
                       "tail weight fraction at each cutoff energy; d<name> is the mean signed "
                       "fractional change of <name> across a single merging event, normalised by "
                       "its pre-merge value at the merge steps")

    rootgrp.createVariable("cutoff_energy_eV", "f8", ("cutoff",))[:] = energies

    for name in ["F"] + weight_names:
        dims = ("cutoff",) if name == "F" else ()
        pre, post = summary[name + "_merge_pre"], summary[name + "_merge_post"]

        rootgrp.createVariable(name, "f8", dims)[...] = summary[name]
        rootgrp.createVariable(name + "_merge_pre", "f8", dims)[...] = pre
        rootgrp.createVariable(name + "_merge_post", "f8", dims)[...] = post
        rootgrp.createVariable("d" + name, "f8", dims)[...] = (post - pre) / pre

    for name in step_names + ["np", "fallback_frac"]:
        rootgrp.createVariable(name, "f8", ())[...] = summary[name]

    rootgrp.createVariable("n_merges", "i8", ())[...] = summary["n_merges"]

    rootgrp.close()


# process a set of files that have a naming scheme: {fname}_seed{seed}_tail.nc
# (or {fname}_tail.nc in case seed=0)
def summarize(fname, nseeds, start_t, end_t):
    print(f"Processing {fname} with {nseeds} seeds ({nseeds+1} files in total)")

    for adds in range(nseeds + 1):
        name = fname + (f"_seed{adds}" if adds > 0 else "")
        energies, summary = summarize_single_run(name, start_t, end_t)
        write_summary(name, energies, summary, start_t, end_t)


t_min, t_max = times[field_Tn]
ts_min = round(t_min / dt)
ts_max = round(t_max / dt)

print(f"Summarizing tail data for E = {field_Tn}Tn over timesteps [{ts_min}, {ts_max})")

for (run, ns) in run_names_and_seeds:
    summarize(run, ns, ts_min, ts_max)
