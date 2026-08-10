import numpy as np
import os
from matplotlib import pyplot as plt
from netCDF4 import Dataset

# plot the per-run tail summaries produced by convert_ionization_tails.py, which has already
# reduced each run over the averaging window -- nothing here reads the raw per-timestep files.
#
# for every quantity two things are plotted against the mean number of electron particles:
#  1) its mean over the averaging window (F for the tail weight fractions, w_ratio / sigma_w /
#     sigma_logw for the electron weight spread, n_coll / sigma_g_w_max for the collision work)
#  2) dF, dw_ratio, ...: its mean signed fractional change across a single merging event. this is
#     measured at the merge itself, so it needs no correction for the fact that merges only fire on
#     steps that just ionized
#
# both are averaged over the ensemble here, with the spread over seeds giving the error bars

# uncomment savefigs = False to turn off saving of figures
savefigs = True
# savefigs = False

pref = "scratch/data2/"

# assumed to be equal across all runs; number of extra seeds, i.e. nseeds+1 files per run
n_seeds_for_run = [15, 15, 15, 15, 15, 15]

octree_runs = [[41, 38], [62, 58], [95, 88], [131, 122], [178, 166], [236, 220]]
nnls_runs = [[4, 41, 38], [5, 62, 58], [6, 95, 88], [7, 131, 122], [8, 178, 166], [9, 236, 220]]
nnls_rp_runs = [[4, 41, 38], [5, 62, 58], [6, 95, 88], [7, 131, 122], [8, 178, 166], [9, 236, 220]]
nnls_erp_runs = [[4, 41, 38], [5, 62, 58], [6, 95, 88], [7, 131, 122], [8, 178, 166], [9, 236, 220]]

# run used as the reference: enough particles that merging barely acts. set to None to drop the
# reference line. it needs a _tail.nc of its own, i.e. it has to be run from a _tail.jl
reference_run = "octree_mid_24000_to_12000"
reference_nseeds = 7

# the averaging window is baked into the summaries by convert_ionization_tails.py and recorded
# there as the start_t / end_t attributes; change it in that script and re-run it

field_Tn = 400

# which two of the recorded cutoff energies to plot. the simulations record four, placed at
# [0.33, 0.75, 1.25, 2.0] x the ionization threshold, i.e. bulk / just below / just above / far
# tail. indices 1 and 2 straddle the ionization threshold, which is the pair the rate responds to
cutoff_indices = [1, 2]

# plotting parameters
label_size = 24
tick_size = 20
legend_size = 20

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["axes.linewidth"] = 0.8

# the weight-spread quantities, which are recorded pre- and post-merge like the tail functions
weight_names = ["w_ratio", "sigma_w", "sigma_logw"]


# read one run's summary, returning a flat dict of scalars: F{i} / dF{i} are the mean tail fraction
# and its change per merge at cutoff i, and likewise {name} / d{name} for the other quantities
def read_single_run(fname):
    ds = Dataset(fname + "_tail_summary.nc")

    energies = np.asarray(ds.variables["cutoff_energy_eV"])
    out = {"energies": energies}

    for i in range(len(energies)):
        out[f"F{i}"] = float(ds.variables["F"][i])
        out[f"dF{i}"] = float(ds.variables["dF"][i])

    for name in weight_names + ["n_coll", "n_eq_w_coll", "sigma_g_w_max",
                                "np", "n_merges", "fallback_frac"]:
        out[name] = float(ds.variables[name][...])
        if "d" + name in ds.variables:
            out["d" + name] = float(ds.variables["d" + name][...])

    ds.close()
    return out


# average a run over its seeds; the spread over seeds gives the error bars. returns None if any
# member of the ensemble is missing, so a sweep that is still running can be plotted as it goes
def read_run(fname, nseeds):
    names = [fname + (f"_seed{adds}" if adds > 0 else "") for adds in range(nseeds + 1)]

    missing = [n for n in names if not os.path.exists(n + "_tail_summary.nc")]
    if missing:
        print(f"  skipping {os.path.basename(fname)}: {len(missing)}/{len(names)} summaries missing")
        return None

    per_seed = [read_single_run(n) for n in names]
    nfiles = len(per_seed)

    out = {"energies": per_seed[0]["energies"]}
    for key in per_seed[0]:
        if key == "energies":
            continue
        vals = np.array([r[key] for r in per_seed])
        out[key] = np.mean(vals)
        out[key + "_se"] = np.std(vals, ddof=1) / np.sqrt(nfiles) if nfiles > 1 else np.nan

    return out


print(f"Plotting tail summaries for E = {field_Tn}Tn")

schemes = [("Octree", octree_runs, "octree_mid_{0}_to_{1}", "o"),
           ("NNLS", nnls_runs, "NNLS_{0}full_{1}", "d"),
           ("NNLS, RP", nnls_erp_runs, "NNLSrate_exact_{0}full_{1}", "s"),
           ("NNLS, ARP", nnls_rp_runs, "NNLSrate_approx_{0}full_{1}", "^")]

data = {}
for label, runs, tag, _ in schemes:
    print(f"{label}:")
    data[label] = []
    for ns, run in zip(n_seeds_for_run, runs):
        filename = f"{pref}ionization_Ar_{field_Tn}Tn_" + tag.format(*run[:2]) + "_es"
        data[label].append(read_run(filename, ns))

reference = None
if reference_run is not None:
    print("Reference:")
    reference = read_run(f"{pref}ionization_Ar_{field_Tn}Tn_{reference_run}_es", reference_nseeds)

energies = next(d["energies"] for runs in data.values() for d in runs if d is not None)

print()
hdr = f"{'run':<24}{'Np':>8}{'merges':>9}{'fallb':>7}{'n_coll':>9}"
for i in cutoff_indices:
    hdr += f"{f'F({energies[i]:.1f})':>11}{'dF %':>9}"
hdr += f"{'w_max/w_min':>13}{'dw %':>9}{'sig_lnw':>10}{'dsig %':>9}"
print(hdr)
for label, _, _, _ in schemes:
    for d in data[label]:
        if d is None:
            continue
        line = f"{label:<24}{d['np']:8.1f}{d['n_merges']:9.0f}{d['fallback_frac']*100:6.1f}%{d['n_coll']:9.2f}"
        for i in cutoff_indices:
            line += f"{d[f'F{i}']:11.5f}{d[f'dF{i}']*100:9.3f}"
        line += f"{d['w_ratio']:13.4g}{d['dw_ratio']*100:9.2f}{d['sigma_logw']:10.4f}{d['dsigma_logw']*100:9.2f}"
        print(line)
    print()


# plot one quantity per subplot against the mean number of particles, one curve per scheme
def plot_vs_np(keys, ylabels, fname, ref=False, logy=False, zero_line=False, scale=1.0):
    fig = plt.figure(figsize=(18, 6))
    axes = [fig.add_subplot(1, 2, 1), fig.add_subplot(1, 2, 2)]

    for ax, key, ylabel in zip(axes, keys, ylabels):
        if ref and reference is not None:
            ax.plot([xl1, xl2], [reference[key] * scale] * 2, color="k", linewidth=2, label="Reference")

        for label, _, _, marker in schemes:
            present = [d for d in data[label] if d is not None]
            if not present:
                continue
            x = np.asarray([d["np"] for d in present])
            y = np.asarray([d[key] for d in present]) * scale
            e = np.asarray([d[key + "_se"] for d in present]) * scale
            ax.errorbar(x, y, yerr=e, marker=marker, capsize=4, linewidth=2, label=label)

        if zero_line:
            ax.axhline(0.0, color='k', linewidth=0.8)
        if logy:
            ax.set_yscale("log")

        ax.margins(y=0.08)
        ax.set_xlim([xl1, xl2])
        ax.grid()
        ax.tick_params(axis='both', labelsize=tick_size,)
        ax.set_xlabel(r"$\overline{N_p}$", fontsize=label_size)
        ax.set_ylabel(ylabel, fontsize=label_size)

    axes[1].legend(fontsize=legend_size, framealpha=1.0, title=f"E = {field_Tn} Tn",
                   title_fontsize=legend_size)

    if savefigs:
        fig.savefig(fname, bbox_inches="tight")


x_all = np.concatenate([[d["np"] for d in data[label] if d is not None] for label, _, _, _ in schemes])
xl1 = x_all.min() - 2
xl2 = x_all.max() + 2

i1, i2 = cutoff_indices

# mean tail functions, and their change per merging event
plot_vs_np([f"F{i1}", f"F{i2}"],
           [rf"$F({energies[i]:.1f}\,\mathrm{{eV}})$" for i in cutoff_indices],
           f"ionization_tail_functions_{field_Tn}Tn.pdf", ref=True)

plot_vs_np([f"dF{i1}", f"dF{i2}"],
           [rf"$\Delta F({energies[i]:.1f}\,\mathrm{{eV}})$, \%" for i in cutoff_indices],
           f"ionization_tail_change_per_merge_{field_Tn}Tn.pdf", zero_line=True, scale=100.0)

# weight spread, and its change per merging event
plot_vs_np(["w_ratio", "sigma_logw"],
           [r"$w_{\mathrm{max}}/w_{\mathrm{min}}$", r"$\sigma_{\ln w}$"],
           f"ionization_weight_spread_{field_Tn}Tn.pdf", ref=True, logy=True)

plot_vs_np(["dw_ratio", "dsigma_logw"],
           [r"$\Delta (w_{\mathrm{max}}/w_{\mathrm{min}})$, \%", r"$\Delta \sigma_{\ln w}$, \%"],
           f"ionization_weight_change_per_merge_{field_Tn}Tn.pdf", zero_line=True, scale=100.0)

# collision work: candidate pairs tested per timestep, and the NTC majorant that sets it
plot_vs_np(["n_coll", "sigma_g_w_max"],
           [r"$\overline{N_{\mathrm{coll}}}$", r"$\overline{(\sigma g w)_{\mathrm{max}}}$"],
           f"ionization_collision_work_{field_Tn}Tn.pdf", ref=False)
