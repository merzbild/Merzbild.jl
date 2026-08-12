import re

import numpy as np
from matplotlib import pyplot as plt

# plot the cost of merging in the 1D Fourier simulations
# (simulations/1D/fourier_varweight_octree.jl and simulations/1D/fourier_varweight_nnls.jl)
# as a function of the average number of particles in the simulation.
# the input is the stdout of those simulations: for each run it contains a TimerOutputs.jl
# table, followed by the average number of particles in that run.
# the octree and the NNLS runs are of different length, so all per-timestep costs are
# normalized by the number of timesteps of the run they come from

# uncomment savefigs = False to turn off saving of figures
savefigs = True
# savefigs = False

# directory where the logs are located
pref = "scratch/"

# logs to parse; all runs found in all of them are pooled together
logfiles = [f"{pref}fourier_octree.log", f"{pref}fourier_nnls.log"]

# annotate the first and last NNLS points with the order of the moment system used
annotate_L = True

# TimerOutputs sections holding the merging cost, per merging method; the first section of a
# method is the one whose per-merge cost is plotted, the remaining ones are the fall-back
# merges performed when an NNLS merge fails, which are included in the full timestep cost
merge_sections = {"octree": ["merge"],
                  "nnls": ["merge NNLS", "merge NNLS backup", "merge octree"]}

# sections which, together with the merging, make up the cost of a timestep
extra_sections = ["collide", "convect", "sort", "restore ordering"]

# section called once per timestep, used to count the number of timesteps of a run
timestep_section = "convect"

# set plotting parameters
label_size = 24
tick_size = 20
legend_size = 20

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["axes.linewidth"] = 0.8

# a row of a TimerOutputs table: "restore ordering    30.0k    10.2s   13.8%   339μs  ..."
# the section name may contain single spaces, the columns are separated by at least two
section_re = re.compile(r"^\s*(?P<name>\S(?:.*?\S)?)\s{2,}"
                        r"(?P<ncalls>\d[\d.]*[kMGT]?)\s+"
                        r"(?P<time>\d[\d.]*\s*[a-zμµ]+)\s+"
                        r"(?P<pct>\d[\d.]*)%\s+"
                        r"(?P<avg>\d[\d.]*\s*[a-zμµ]+)\s")

# the average number of particles, printed on its own line after the table
navg_re = re.compile(r"^\s*n_p_avg\s*=\s*(?P<navg>\d+\.\d+(?:[eE][-+]?\d+)?)\s*$")

# the parameters of an NNLS run, printed before the table as
# "{n_species}, [{order of the moment system}, {target number of particles}, {n_moments}]"
config_re = re.compile(r"^\s*\d+,\s*\[(?P<order>\d+),\s*(?P<target>\d+),\s*(?P<nmoments>\d+)\]\s*$")

time_units = {"ns": 1e-9, "μs": 1e-6, "µs": 1e-6, "us": 1e-6, "ms": 1e-3, "s": 1.0,
              "min": 60.0, "h": 3600.0}
count_units = {"": 1, "k": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}


def parse_time(s):
    value, unit = re.match(r"(\d[\d.]*)\s*([a-zμµ]+)", s).groups()

    return float(value) * time_units[unit]


def parse_count(s):
    value, unit = re.match(r"(\d[\d.]*?)([kMGT]?)$", s).groups()

    return float(value) * count_units[unit]


# a run is a TimerOutputs table terminated by the line holding the average number of
# particles; the sections timed at t=0 are the warm-up merges performed on the initially
# sampled particles and are dropped
def parse_log(path):
    runs = []
    sections = {}
    order = None

    with open(path) as f:
        for line in f:
            config_match = config_re.match(line)
            navg_match = navg_re.match(line)

            if config_match is not None:
                order = int(config_match.group("order"))
                sections = {}
                continue

            if navg_match is not None:
                if len(sections) > 0:
                    runs.append({"navg": float(navg_match.group("navg")),
                                 "order": order,
                                 "sections": sections})

                sections = {}
                order = None
                continue

            section_match = section_re.match(line)

            if section_match is None:
                continue

            name = section_match.group("name")

            if name.endswith("(t=0)") or name.endswith("t=0"):
                continue

            sections[name] = {"ncalls": parse_count(section_match.group("ncalls")),
                              "time": parse_time(section_match.group("time")),
                              "avg": parse_time(section_match.group("avg"))}

    return runs


# total time spent in the sections that are present in the run, in seconds
def total_time(run, names):
    return sum(run["sections"][name]["time"] for name in names if name in run["sections"])


# a run belongs to a merging method if its primary merging section was timed; the NNLS runs
# also contain an octree section, but only as a fall-back for the failed NNLS merges
def select(runs, method):
    sections = merge_sections[method]
    sel = [r for r in runs if sections[0] in r["sections"]]

    sel.sort(key=lambda r: r["navg"])

    navg = np.asarray([r["navg"] for r in sel])
    nsteps = np.asarray([r["sections"][timestep_section]["ncalls"] for r in sel])

    # the cost of a single merge event, i.e. the "avg" column of the merging section
    t_merge = np.asarray([r["sections"][sections[0]]["avg"] for r in sel])

    # the cost of a full timestep is normalized per timestep instead, so that the merging
    # cost, which is incurred once per cell per merging interval, is additive with the cost
    # of the collisions, of the convection and of the sorting of the particles
    t_full = np.asarray([total_time(r, sections + extra_sections) for r in sel]) / nsteps

    return navg, 1e6 * t_merge, 1e6 * t_full, [r["order"] for r in sel]


runs = []

for logfile in logfiles:
    runs += parse_log(logfile)

colors = {"octree": "tab:blue", "nnls": "tab:orange"}
markers = {"octree": "o", "nnls": "d"}
labels = {"octree": "Octree", "nnls": "NNLS"}

# cost of a single merge and cost of a full timestep vs the average number of particles
fig = plt.figure(figsize=(18, 6))

ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

for method in ["octree", "nnls"]:
    navg, t_merge, t_full, orders = select(runs, method)

    if len(navg) == 0:
        print(f"no runs found for {labels[method]}")
        continue

    for ax, y in zip([ax1, ax2], [t_merge, t_full]):
        ax.plot(navg, y, marker=markers[method], color=colors[method], linewidth=2,
                label=labels[method])

        if annotate_L and method == "nnls":
            # the first label is placed to the right of its point, the last one to the left,
            # so that neither of them is pushed outside of the axes
            for i, offset, align in [(0, 1.03, "left"), (-1, 0.97, "right")]:
                if orders[i] is not None:
                    ax.text(navg[i] * offset, y[i], f"L={orders[i]}",
                            fontsize=legend_size - 2, ha=align, va="center")

for ax in [ax1, ax2]:
    ax.set_yscale("log")
    ax.grid()
    ax.grid(which="minor", linewidth=0.4)
    ax.tick_params(axis='both', labelsize=tick_size,)
    ax.set_xlabel(r"$\overline{N}_p$", fontsize=label_size)

ax1.legend(fontsize=legend_size, framealpha=1.0, loc="upper left")

ax1.set_ylabel(r"$t_{merge}$ per merge, $\mu$s", fontsize=label_size)
ax2.set_ylabel(r"$t_{full}$ per timestep, $\mu$s", fontsize=label_size)

if savefigs:
    fig.savefig("fourier_merge_performance.pdf", bbox_inches="tight")
