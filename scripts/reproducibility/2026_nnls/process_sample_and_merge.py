import numpy as np
from matplotlib import pyplot as plt

# uncomment savefigs = False to turn off saving of figures
savefigs = True
# savefigs = False

# directory where data is located
# files assumed to be named {pref}octree_{stype} and {pref}nnls_{stype}
# where stype is either equalweight or weighted
pref = "scratch/data/"

# set plotting parameters
label_size = 24
tick_size = 20
legend_size = 20

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["axes.linewidth"] = 0.8

def parse_log(path):
    w_std = []
    w_log_std = []
    w_r = []
    w_std_orig = []
    w_log_std_orig = []
    w_r_orig = []
    np_mean = []
    tailf_500 = []
    tailf_750 = []
    tailf_500_orig = []
    tailf_750_orig = []
    
    with open(path) as f:
        for line in f:
            if line.startswith("Npost"):
                np_mean.append(float(line.split("=")[-1]))
            elif line.startswith("pre-merge weight ratio"):
                w_r_orig.append(float(line.split(":")[-1]))
            elif line.startswith("pre-merge weight std"):
                w_std_orig.append(float(line.split(":")[-1]))
            elif line.startswith("pre-merge weight log std"):
                w_log_std_orig.append(float(line.split(":")[-1]))
            elif line.startswith("weight ratio"):
                w_r.append(float(line.split(":")[-1]))
            elif line.startswith("weight std"):
                w_std.append(float(line.split(":")[-1]))
            elif line.startswith("weight log std"):
                w_log_std.append(float(line.split(":")[-1]))
            elif line.startswith("f_tail(500.0)"):
                lsp = line.split()
                tailf_500.append(float(lsp[-1]))
                tailf_500_orig.append(float(lsp[-3]))
            elif line.startswith("f_tail(750.0)"):
                lsp = line.split()
                tailf_750.append(float(lsp[-1]))
                tailf_750_orig.append(float(lsp[-3]))
    
    return {"np_mean": np_mean,
            "w_r": w_r, "w_std": w_std, "w_log_std": w_log_std,
            "w_r_orig": w_r_orig, "w_std_orig": w_std_orig, "w_log_std_orig": w_log_std_orig,
            "tailf_500": tailf_500, "tailf_750": tailf_750,
            "tailf_500_orig": tailf_500_orig, "tailf_750_orig": tailf_750_orig}

for stype in ["equalweight", "weighted"]:
    octree_data = parse_log(f"{pref}octree_{stype}.log")
    nnls_data = parse_log(f"{pref}nnls_{stype}.log")

    xl1 = np.min(octree_data["np_mean"]) - 2
    xl2 = np.max(nnls_data["np_mean"]) + 2

    # plot tail functions
    fig = plt.figure(figsize=(18, 6))
    ax1 = fig.add_subplot(1,2,1)

    ax1.plot(octree_data["np_mean"], octree_data["tailf_500"], marker="o", linewidth=2, label=f"Octree")
    ax1.plot(nnls_data["np_mean"], nnls_data["tailf_500"], marker="d", linewidth=2, label=f"NNLS")
    t500_mean = np.mean(octree_data["tailf_500_orig"])
    ax1.plot([xl1, xl2], [t500_mean, t500_mean], color="k", linewidth=2)
    ax1.set_xlim([xl1, xl2])

    ax2 = fig.add_subplot(1,2,2)
    t750_mean = np.mean(octree_data["tailf_750_orig"])
    ax2.plot([xl1, xl2], [t750_mean, t750_mean], color="k", linewidth=2, label="Sample value")
    ax2.plot(octree_data["np_mean"], octree_data["tailf_750"], marker="o", linewidth=2, label=f"Octree")
    ax2.plot(nnls_data["np_mean"], nnls_data["tailf_750"], marker="d", linewidth=2, label=f"NNLS")
    ax2.set_xlim([xl1, xl2])

    for ax in [ax1, ax2]:
        ax.grid()
        ax.tick_params(axis='both', labelsize=tick_size,)
        ax.set_xlabel(r"$N_{post}$", fontsize=label_size)
        
    ax2.legend(fontsize=legend_size, framealpha=1.0)   
    ax1.set_ylabel(r"$F(500)$", fontsize=label_size)
    ax2.set_ylabel(r"$F(750)$", fontsize=label_size)

    if savefigs:
        fig.savefig(f"tail_functions_{stype}.pdf", bbox_inches="tight")

    # plot weight deviations
    fig = plt.figure(figsize=(18, 6))
    ax1 = fig.add_subplot(1,2,1)

    if stype == "weighted":
        ax1.plot([xl1, xl2], [octree_data["w_std_orig"][0], octree_data["w_std_orig"][0]], linewidth=2, color="k", label="Sample value")
    ax1.plot(octree_data["np_mean"], octree_data["w_std"], marker="o", linewidth=2, color="tab:blue", label="Octree")
    ax1.plot(nnls_data["np_mean"], nnls_data["w_std"], marker="d", linewidth=2, color="tab:orange", label="NNLS")
    ax1.set_xlim([xl1, xl2])

    ax2 = fig.add_subplot(1,2,2)

    if stype == "weighted":
        ax2.plot([xl1, xl2], [octree_data["w_log_std_orig"][0], octree_data["w_log_std_orig"][0]], linewidth=2, color="k", label="Sample value")
    ax2.plot(octree_data["np_mean"], octree_data["w_log_std"], marker="o", linewidth=2, color="tab:blue", label="Octree")
    ax2.plot(nnls_data["np_mean"], nnls_data["w_log_std"], marker="d", linewidth=2, color="tab:orange", label="NNLS")
    ax2.set_xlim([xl1, xl2])

    for ax in [ax1, ax2]:
        ax.grid()
        ax.tick_params(axis='both', labelsize=tick_size,)
        ax.set_xlabel(r"$N_{post}$", fontsize=label_size)
        
    ax2.legend(fontsize=legend_size, framealpha=1.0)
    ax1.set_ylabel(r"$\sigma_{w}$", fontsize=label_size)
    ax2.set_ylabel(r"$\sigma_{\ln w}$", fontsize=label_size)

    if savefigs:
        fig.savefig(f"sigma_weights_{stype}.pdf", bbox_inches="tight")

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1,1,1)

    if stype == "weighted":
        ax.plot([xl1, xl2], [octree_data["w_r_orig"][0], octree_data["w_r_orig"][0]], linewidth=2, color="k", label="Sample value")
    ax.plot(octree_data["np_mean"], octree_data["w_r"], marker="o", linewidth=2, label=f"Octree")
    ax.plot(nnls_data["np_mean"], nnls_data["w_r"], marker="o", linewidth=2, label=f"NNLS")
    ax.set_xlim([xl1, xl2])

    ax.grid()
    ax.tick_params(axis='both', labelsize=tick_size,)
    ax.set_xlabel(r"$N_{post}$", fontsize=label_size)

    ax.legend(fontsize=legend_size, framealpha=1.0)
    ax.set_yscale("log")
    ax.set_ylabel(r"$w_{\mathrm{max}}/w_{\mathrm{min}}$", fontsize=label_size)

    if savefigs:
        fig.savefig(f"ratio_weights_{stype}.pdf", bbox_inches="tight")