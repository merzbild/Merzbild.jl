import re

import numpy as np
from netCDF4 import Dataset
from matplotlib import pyplot as plt

# plot the error in the wall pressure of the 1D Fourier simulations against the cost of a
# timestep of those simulations, so that the accuracy of the octree and of the NNLS merging
# can be compared at equal cost rather than at an equal number of particles.
# the error is computed exactly as in process_fourier.py, the cost is taken from the
# TimerOutputs.jl tables written by the simulations (see plot_performance_fourier.py)

# uncomment savefigs = False to turn off saving of figures
savefigs = True
# savefigs = False

# directory with the simulation results
pref = "scratch/data/"

# directory where the logs of those simulations are located
log_pref = "scratch/"

# logs to parse; all runs found in all of them are pooled together
logfiles = [f"{log_pref}fourier_octree.log", f"{log_pref}fourier_nnls.log"]

# this is the timestep after which the averaging of results was turned on in the simulation
avg_start = 500000

# the runs of each merging method, named by the "{target Np}_{n_moments}" part of the file name
runs = {"octree": ["45_38", "69_58", "105_88", "146_122", "200_166", "266_220"],
        "nnls": ["45_38", "69_58", "105_88", "146_122", "200_166"]}

# how each method is named in the file names; "fw" is the fixed-weight reference simulation
file_tags = {"octree": "octree", "nnls": "NNLS", "fw": "FW"}

# the wall the pressure error is computed at
wall = 1

# annotate the first and last NNLS points with the order of the moment system used
annotate_L = True

# TimerOutputs sections holding the merging cost, per merging method; the first section of a
# method is the primary one, the remaining ones are the fall-back merges performed when an
# NNLS merge fails
merge_sections = {"octree": ["merge"],
                  "nnls": ["merge NNLS", "merge NNLS backup", "merge octree"]}

# sections which, together with the merging, make up the cost of a timestep
extra_sections = ["collide", "convect", "sort", "restore ordering"]

# section called once per timestep, used to count the number of timesteps of a run
timestep_section = "convect"

# a run of the log is matched to a run of the simulation results by the average number of
# particles; the log averages over the whole simulation and the results only over the part
# of it that was averaged over, so the two differ slightly
navg_rtol = 0.02

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
    parsed = []
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
                    parsed.append({"navg": float(navg_match.group("navg")),
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

    return parsed


# total time spent in the sections that are present in the run, in seconds
def total_time(run, names):
    return sum(run["sections"][name]["time"] for name in names if name in run["sections"])


# cost of a timestep of a run, in microseconds; the octree and the NNLS simulations are of
# different length, so every run is normalized by its own number of timesteps
def timestep_cost(run, method):
    nsteps = run["sections"][timestep_section]["ncalls"]

    return 1e6 * total_time(run, merge_sections[method] + extra_sections) / nsteps


# the log holds no name of the run it timed, so the run of the log and the run of the
# simulation results are matched by their average number of particles
def find_run(parsed, method, navg):
    candidates = [r for r in parsed if merge_sections[method][0] in r["sections"]]
    closest = min(candidates, key=lambda r: abs(r["navg"] - navg))

    if abs(closest["navg"] - navg) > navg_rtol * navg:
        print(f"no timings found for the {method} run with Np = {navg:.1f}")
        return None

    return closest


# load the results of the merging simulations and of the fixed-weight reference simulation
def load(method, run):
    base = f"{pref}avg_Fourier_{file_tags[method]}_ntc_0.005_1000_0.0_300.0_600.0"
    base = base if run is None else f"{base}_{run}"

    ds = Dataset(f"{base}_after{avg_start}.nc")
    n_p = np.asarray(ds["np"][:][0, 0, :])
    ds.close()

    ds = Dataset(f"{base}_surf_after{avg_start}.nc")
    s_p = np.asarray(ds["normal_pressure"][0, 0, :])
    ds.close()

    return np.mean(n_p), s_p


ref_s_p = load("fw", None)[1]

parsed = []

for logfile in logfiles:
    parsed += parse_log(logfile)

colors = {"octree": "tab:blue", "nnls": "tab:orange"}
markers = {"octree": "o", "nnls": "d"}
labels = {"octree": "Octree", "nnls": "NNLS"}

# error in the wall pressure vs cost of a timestep
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(1, 1, 1)

for method in ["octree", "nnls"]:
    x_t = []
    y_p = []
    orders = []

    for run in runs[method]:
        navg, s_p = load(method, run)
        timings = find_run(parsed, method, navg)

        if timings is None:
            continue

        x_t.append(timestep_cost(timings, method))
        # the error in the wall pressure, computed as in process_fourier.py
        y_p.append(np.sqrt(np.sum((s_p[wall] / ref_s_p[wall] - 1)**2)) * 100)
        orders.append(timings["order"])

    ax.plot(x_t, y_p, marker=markers[method], color=colors[method], linewidth=2,
            label=labels[method])

    if annotate_L and method == "nnls":
        # the first label is placed to the right of its point, the last one to the left,
        # so that neither of them is pushed outside of the axes
        for i, offset, align in [(0, 1.03, "left"), (-1, 0.97, "right")]:
            if orders[i] is not None:
                ax.text(x_t[i] * offset, y_p[i], f"$L={orders[i]}$",
                        fontsize=legend_size - 2, ha=align, va="center")

ax.set_xscale("log")
ax.set_yscale("log")
ax.grid()
ax.grid(which="minor", linewidth=0.4)
ax.tick_params(axis='both', labelsize=tick_size)

ax.legend(fontsize=legend_size, framealpha=1.0)

ax.set_xlabel(r"$t_{full}$ per timestep, $\mu$s", fontsize=label_size)
ax.set_ylabel(r"$\overline{\mathcal{B}}(p_s)$, \%", fontsize=label_size)

if savefigs:
    fig.savefig("fourier_error_vs_cost.pdf", bbox_inches="tight")
