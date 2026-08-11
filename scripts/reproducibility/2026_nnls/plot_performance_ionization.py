import re

import numpy as np
from matplotlib import pyplot as plt

# plot the cost of electron merging in the 0D ionization simulations
# (simulations/0D/ionization/0D_ionization_1neutralspecies_es.jl and
# simulations/0D/ionization/0D_ionization_1neutralspecies_nnls_es.jl) as a function of the
# average number of electron particles in the simulation.
# the input is the stdout of those simulations: for each run it contains a TimerOutputs.jl
# table, directly followed by the average number of electron particles in that run.

# uncomment savefigs = False to turn off saving of figures
savefigs = True
# savefigs = False

# directory where the logs are located
pref = "scratch/"

# logs to parse; all runs found in all of them are pooled together
logfiles = [f"{pref}ion_octree.log", f"{pref}ion_nnls.log"]

# only keep runs performed with the electrostatic (ES) versions of the scripts,
# i.e. runs in which the electron-neutral collisions are timed as "coll n-e ES"
es_only = True

# target numbers of particles of the run configurations to plot, taken from the "_{L}full_{Np}"
# part of the name of the output file of a run; set to None to plot all runs found in the logs
# (this only filters the NNLS runs, the octree runs are not identified by an output file)
allowed_targets = [41, 62, 95, 131, 178, 236]

# annotate the first and last NNLS points with the order of the moment system used
annotate_L = True

# TimerOutputs sections holding the electron merging cost, per merging method
merge_sections = {"octree": "merge e",
                  "nnls": "NNLSmerge e",
                  "nnls_arp": "NNLSmergeARP e",
                  "nnls_erp": "NNLSmergeERP e"}

# sections added to the merging cost in the second subplot
extra_sections = ["coll n-e", "coll n-e ES", "acc e"]

# section used to count the number of timesteps of a run
timestep_section = "props"

# set plotting parameters
label_size = 24
tick_size = 20
legend_size = 20

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["axes.linewidth"] = 0.8

# a row of a TimerOutputs table: "coll n-e ES       500k   33.0ms    3.6%  66.0ns  ..."
# the section name may contain single spaces, the columns are separated by at least two
section_re = re.compile(r"^\s*(?P<name>\S(?:.*?\S)?)\s{2,}"
                        r"(?P<ncalls>\d[\d.]*[kMGT]?)\s+"
                        r"(?P<time>\d[\d.]*\s*[a-zμµ]+)\s+"
                        r"(?P<pct>\d[\d.]*)%\s+"
                        r"(?P<avg>\d[\d.]*\s*[a-zμµ]+)\s")

# the average number of electron particles, printed on its own line after the table
navg_re = re.compile(r"^\s*(?P<navg>\d+\.\d+(?:[eE][-+]?\d+)?)\s*$")

# the netCDF file the run wrote to, printed before the table; the "{L}full" part of the name
# is the order of the moment system used by the NNLS merges
outfile_re = re.compile(r"^\s*\S*ionization_\S+\.nc\s*$")
order_re = re.compile(r"_(\d+)full_(\d+)")

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
# electron particles; sections timed at t=0 ("merge e (t=0)" or "merge e t=0") are the
# warm-up merges performed on the initially sampled particles and are dropped
def parse_log(path):
    runs = []
    sections = {}
    outfile = None

    with open(path) as f:
        for line in f:
            navg_match = navg_re.match(line)

            if outfile_re.match(line):
                outfile = line.strip()
                sections = {}
                continue

            if navg_match is not None:
                if len(sections) > 0:
                    runs.append({"navg": float(navg_match.group("navg")),
                                 "outfile": outfile,
                                 "sections": sections})

                sections = {}
                outfile = None
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


# the order of the moment system and the target number of particles of a run, taken from the
# name of its output file; the octree runs do not write one, so both are unknown for them
def run_config(run):
    if run["outfile"] is None:
        return None, None

    match = order_re.search(run["outfile"])

    return (None, None) if match is None else (int(match.group(1)), int(match.group(2)))


# a run belongs to a merging method if the corresponding electron merging section was timed;
# the NNLS runs also contain an octree "merge e t=0" section, but that one is dropped above
def select(runs, method):
    section = merge_sections[method]
    sel = [r for r in runs if section in r["sections"]]

    if es_only:
        sel = [r for r in sel if "coll n-e ES" in r["sections"]]

    if allowed_targets is not None:
        sel = [r for r in sel
               if run_config(r)[1] is None or run_config(r)[1] in allowed_targets]

    sel.sort(key=lambda r: r["navg"])

    navg = np.asarray([r["navg"] for r in sel])
    nsteps = np.asarray([r["sections"][timestep_section]["ncalls"] for r in sel])

    # the cost of a single merge event, i.e. the "avg" column of the merging section
    t_merge = np.asarray([r["sections"][section]["avg"] for r in sel])

    # the combined cost is normalized per timestep instead, so that the merging cost, which
    # is incurred only every few hundred timesteps, is additive with the per-timestep cost of
    # the electron collisions and of the acceleration of the electrons by the field
    t_total = np.asarray([total_time(r, [section] + extra_sections) for r in sel]) / nsteps

    return navg, 1e6 * t_merge, 1e6 * t_total, [run_config(r)[0] for r in sel]


runs = []

for logfile in logfiles:
    runs += parse_log(logfile)

colors = {"octree": "tab:blue", "nnls": "tab:orange",
          "nnls_arp": "tab:green", "nnls_erp": "tab:red"}
markers = {"octree": "o", "nnls": "d", "nnls_arp": "s", "nnls_erp": "^"}
labels = {"octree": "Octree", "nnls": "NNLS",
          "nnls_arp": "NNLS, ARP", "nnls_erp": "NNLS, RP"}

# cost of a single electron merge event, and cost of electron merging plus electron collisions
# and acceleration per timestep, vs the average number of electron particles
fig = plt.figure(figsize=(18, 6))

ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

for method in ["octree", "nnls", "nnls_arp", "nnls_erp"]:
    navg, t_merge, t_total, orders = select(runs, method)

    if len(navg) == 0:
        print(f"no runs found for {labels[method]}")
        continue

    for ax, y in zip([ax1, ax2], [t_merge, t_total]):
        ax.plot(navg, y, marker=markers[method], color=colors[method], linewidth=2,
                label=labels[method])

        if annotate_L and method == "nnls":
            # the first label is placed to the right of its point, the last one to the left,
            # so that neither of them is pushed outside of the axes
            for i, offset, align in [(0, 1.05, "left"), (-1, 0.95, "right")]:
                if orders[i] is not None:
                    ax.text(navg[i] * offset, y[i], f"L={orders[i]}",
                            fontsize=legend_size - 2, ha=align, va="center")

for ax in [ax1, ax2]:
    # ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid()
    ax.grid(which="minor", linewidth=0.4)
    ax.tick_params(axis='both', labelsize=tick_size,)
    ax.set_xlabel(r"$\overline{N_e}$", fontsize=label_size)

# the data fills the whole axes, so room for the legend is made above the curves
ax1.set_ylim(top=8.0 * ax1.get_ylim()[1])
ax1.legend(fontsize=legend_size, framealpha=1.0, title="Merging", loc="upper left", ncol=2,
           title_fontsize=legend_size)

ax2.legend([],
           [],
           title="Merging + collisions + acceleration",
           framealpha=1.0,
           loc="upper left",
           title_fontsize=legend_size
           )

ax1.set_ylabel(r"$t_{merge}$ per merge, $\mu$s", fontsize=label_size)
ax2.set_ylabel(r"$t_{merge}+t_{coll}+t_{acc}$ per timestep, $\mu$s", fontsize=label_size)

if savefigs:
    fig.savefig("ionization_merge_performance.pdf", bbox_inches="tight")
