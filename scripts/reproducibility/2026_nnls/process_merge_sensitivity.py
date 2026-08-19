import numpy as np
from matplotlib import pyplot as plt

# plot the amplification factor (ratio of RMS of chamfer distance to sqrt(3)*perturbation sigma)
# measured by simulations/0D/basic/merge_sensitivity.jl: how far a
# merging scheme moves its output particles in response to a small perturbation of its input,
# relative to the size of that perturbation.

# uncomment savefigs = False to turn off saving of figures
savefigs = True
# savefigs = False

# directory where data is located
# files assumed to be named {pref}merge_sensitivity_{stype}.log
# where stype is either equalweight or weighted
pref = "scratch/data/"

# set plotting parameters
label_size = 24
tick_size = 20
legend_size = 16

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["axes.linewidth"] = 0.8


# each configuration is written as a block of "key: value" lines headed by a "method:" line
def parse_log(path):
    runs = []
    cur = None

    with open(path) as f:
        for line in f:
            stripped = line.strip()

            if stripped.startswith("method:"):
                # method: nnls, merge parameter: 4, sampling: equal_weight
                parts = stripped.split(",")
                cur = {"method": parts[0].split(":")[-1].strip(),
                       "merge_parameter": int(parts[1].split(":")[-1]),
                       "sampling": parts[2].split(":")[-1].strip()}
                runs.append(cur)
            elif cur is None:
                continue
            elif stripped.startswith("Npost:"):
                cur["npost"] = float(stripped.split(":")[-1])
            elif stripped.startswith("noise sigma"):
                cur["noise"] = float(stripped.split(":")[-1])
            elif stripped.startswith("amplification:"):
                cur["amplification"] = float(stripped.split(":")[-1])

    return runs


# the zero-noise configurations are a check that the two merges are handed identical inputs, so
# their amplification is 0/0; they carry no information here
def select(runs, method, noise):
    sel = [r for r in runs
           if r["method"] == method and abs(r["noise"] - noise) < 1e-12
           and np.isfinite(r["amplification"])]

    sel.sort(key=lambda r: r["npost"])

    return np.asarray([r["npost"] for r in sel]), np.asarray([r["amplification"] for r in sel])


colors = {"octree": "tab:blue", "nnls": "tab:orange"}
markers = {"octree": "o", "nnls": "d"}
labels = {"octree": "Octree", "nnls": "NNLS"}
linestyles = ["-", "--", ":", "-."]

fig = plt.figure(figsize=(18, 6))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

for ax, stype, title in zip([ax1, ax2], ["equalweight", "weighted"],
                            ["Fixed-weight samples", "Variable-weight samples"]):
    runs = parse_log(f"{pref}merge_sensitivity_{stype}.log")

    # the noise levels actually present in the log, smallest first
    noise_levels = sorted({r["noise"] for r in runs if np.isfinite(r["amplification"])})

    for method in ["octree", "nnls"]:
        for noise, ls in zip(noise_levels, linestyles):
            npost, amplification = select(runs, method, noise)

            if len(npost) == 0:
                continue

            if title == "Fixed-weight samples":
                ax.plot(npost, amplification, marker=markers[method], color=colors[method],
                        linestyle=ls, linewidth=2)
            else:
                ax.plot(npost, amplification, marker=markers[method], color=colors[method],
                        linestyle=ls, linewidth=2,
                        label=rf"{labels[method]}, $\sigma/v_{{ref}} = {noise:g}$")

    ax.set_yscale("log")
    ax.grid()
    ax.tick_params(axis='both', labelsize=tick_size,)
    ax.set_xlabel(r"$N_{post}$", fontsize=label_size)

    if title == "Fixed-weight samples":
        ax.legend([], [], fontsize=legend_size, framealpha=1.0, title=title, title_fontsize=legend_size)
    else:
        ax.legend(fontsize=legend_size, framealpha=1.0, title=title, title_fontsize=legend_size)

ax1.set_ylabel(r"$\overline{d}_c/(\sqrt{3}\sigma_v)$", fontsize=label_size)

if savefigs:
    fig.savefig("merge_sensitivity_amplification.pdf", bbox_inches="tight")
