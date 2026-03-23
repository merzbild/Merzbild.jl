import numpy as np
from netCDF4 import Dataset
from matplotlib import pyplot as plt
from scipy import constants
from scipy.special import gamma

savefigs = False
# savefigs = True  # - uncomment to save plots

# number of timesteps (including the 0th timestep) and timestep size in each simulation
n_t = 601
dt = 0.025

# path to directory with simulation results
pref = "scratch/data/"

# the octree merging simulations' parameters ([threshold, Ntarget]) and number of seeds for each set of parameters
octree_mid_runs = [[42, 36], [66, 55], [100, 85], [142, 120], [200, 164], [264, 220]]
octree_mid_seeds = [10800,    4800,     2160,      1040,       560,        400]
# code assumes files have names {pref}bkw_octree_mid_{threshold}_{Ntarget}_{seed}.nc

# the NNLS merging simulations' parameters ([n_full_up_to_total, threshold, ntarget_octree]) and number of seeds for each set of parameters
nnls_runs = [[4, 42, 36], [5, 66, 55], [6, 100, 85], [7, 142, 120], [8, 200, 164], [9, 264, 220]] #,
nnls_seeds = [10800,       4800,        2160,         1040,          560,           400]
# code assumes files have names {pref}bkw_nnls_$(n_full_up_to_total)full_$(threshold)_$(seed).nc

# plotting parameters
label_size = 24
tick_size = 20
legend_size = 20

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["axes.linewidth"] = 0.8

def analytic(time, magic_factor, N):
    # The time scaling in the analytical solution is different due to differences in definition of the
    # reference cross-section; hence the magic_factor which accounts for that
    C = 1.0 - 0.4 * np.exp(-time * magic_factor / 6)
    kk = N / 2
    return C**(kk - 1) * (kk - (kk - 1) * C)

Tref = 273.0
mref = 66.3e-27 
mcd = mref / 2.0
dref = 4.11e-10
nref = 1e23
Lref = 1.0 / (nref * constants.pi * dref**2)
vref = ((2 * constants.k * Tref) / mref)**0.5
time_ref = Lref / vref

kappa_mult = constants.pi * dref**2 * (mcd / (2 * constants.k * Tref))**(-0.5) / gamma(5/2 - 1.0)
ttt_bkw = 1 / (4 * constants.pi * nref * kappa_mult)

# compute the time-scaling factor that accounts for different definitions of reference cross-sections
mf = time_ref / ttt_bkw / (4 * constants.pi)

# load, ensemble-average octree data
octree_mid_data = {x[0]: {"M4_avg": np.zeros(n_t), "M4_std": np.zeros(n_t),
                          "M6_avg": np.zeros(n_t), "M6_std": np.zeros(n_t),
                          "M8_avg": np.zeros(n_t), "M8_std": np.zeros(n_t),
                          "np_avg": np.zeros(n_t), "np_min": np.zeros(n_t), "np_max": np.zeros(n_t), "t": np.zeros(n_t)} for x in octree_mid_runs}

for octree_mid_seed_n, run in zip(octree_mid_seeds, octree_mid_runs):
    print(run)
    
    moment_arrs = [np.zeros((octree_mid_seed_n, n_t)), np.zeros((octree_mid_seed_n, n_t)), np.zeros((octree_mid_seed_n, n_t))]
    np_arr = np.zeros((octree_mid_seed_n, n_t-1))
    
    for seed in range(1,octree_mid_seed_n+1):
        threshold = run[0]
        Ntarget = run[1]
        
        ds = Dataset(f"{pref}bkw_octree_mid_{threshold}_{Ntarget}_{seed}.nc")
        
        octree_mid_data[threshold]["t"] = np.asarray(ds.variables["timestep"][:].data)
        np_arr[seed-1, :] = np.asarray(ds.variables["np"][:].data[1:, 0, 0])
        
        for (mom_id, mom) in enumerate([4, 6, 8]):
            moment_arrs[mom_id][seed-1, :] = np.asarray(ds.variables["moments"][:].data[:, 0, 0, mom_id])

        ds.close()
    
    octree_mid_data[threshold]["t"] *= dt
    octree_mid_data[threshold][f"np_avg"] = np.mean(np_arr)
    octree_mid_data[threshold][f"np_min"] = np.mean(np.min(np_arr, axis=1))
    octree_mid_data[threshold][f"np_max"] = np.mean(np.max(np_arr, axis=1))
    
    for (mom_id, mom) in enumerate([4, 6, 8]):
        octree_mid_data[threshold][f"M{mom}_avg"] = np.mean(moment_arrs[mom_id], axis=0)
        octree_mid_data[threshold][f"M{mom}_std"] = np.std(moment_arrs[mom_id], axis=0)
            
for run in octree_mid_runs:
    threshold = run[0]
    for M in [4, 6, 8]:
        octree_mid_data[threshold][f"M{M}_bias"] = np.mean(np.abs(analytic(octree_mid_data[threshold]["t"], mf, M) - octree_mid_data[threshold][f"M{M}_avg"]))

# load, ensemble-average NNLS data
nnls_data = {x[0]: {"M4_avg": np.zeros(n_t), "M4_std": np.zeros(n_t),
                    "M6_avg": np.zeros(n_t), "M6_std": np.zeros(n_t),
                    "M8_avg": np.zeros(n_t), "M8_std": np.zeros(n_t),
                    "np_avg": np.zeros(n_t), "np_min": np.zeros(n_t), "np_max": np.zeros(n_t), "t": np.zeros(n_t)} for x in nnls_runs}

for (nnls_seed_n, run) in zip(nnls_seeds, nnls_runs):
    print(run)
    
    moment_arrs = [np.zeros((nnls_seed_n, n_t)), np.zeros((nnls_seed_n, n_t)), np.zeros((nnls_seed_n, n_t))]
    np_arr = np.zeros((nnls_seed_n, n_t-1))
    
    for seed in range(1,nnls_seed_n+1):
        threshold = run[0]
        
        ds = Dataset(f"{pref}bkw_nnls_{run[0]}full_{run[1]}_{seed}.nc")
        
        nnls_data[threshold]["t"] = np.asarray(ds.variables["timestep"][:].data)
        np_arr[seed-1, :] = np.asarray(ds.variables["np"][:].data[1:, 0, 0])
        
        for (mom_id, mom) in enumerate([4, 6, 8]):
            moment_arrs[mom_id][seed-1, :] = np.asarray(ds.variables["moments"][:].data[:, 0, 0, mom_id])

        ds.close()
    
    nnls_data[threshold]["t"] *= dt
    nnls_data[threshold][f"np_avg"] = np.mean(np_arr)
    nnls_data[threshold][f"np_min"] = np.mean(np.min(np_arr, axis=1))
    nnls_data[threshold][f"np_max"] = np.mean(np.max(np_arr, axis=1))
    
    for (mom_id, mom) in enumerate([4, 6, 8]):
        nnls_data[threshold][f"M{mom}_avg"] = np.mean(moment_arrs[mom_id], axis=0)
        nnls_data[threshold][f"M{mom}_std"] = np.std(moment_arrs[mom_id], axis=0)
            
for run in nnls_runs:
    threshold = run[0]
    for M in [4, 6, 8]:
        nnls_data[threshold][f"M{M}_bias"] = np.mean(np.abs(analytic(nnls_data[threshold]["t"], mf, M) - nnls_data[threshold][f"M{M}_avg"]))

# plot evolution of 4th and 8th moments for 2 cases: NNLS with up to 4th-order moments conserved; NNLS with up to 9th-order moments conserved
# octree runs chosen accordingly
for nnls_upto in [4,9]:
    fig = plt.figure(figsize=(18,6))

    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    thresh_post = {p[0]: [p[1], p[2]] for p in nnls_runs}

    for M, ax in zip([4,8], [ax1,ax2]):
        np_thr = thresh_post[nnls_upto][0]
        np_post = thresh_post[nnls_upto][1]

        ax.plot(octree_mid_data[np_thr]["t"], analytic(octree_mid_data[np_thr]["t"], mf, M), label="Analytical", linewidth=2, color="k") 
        ax.plot(octree_mid_data[np_thr]["t"], octree_mid_data[np_thr][f"M{M}_avg"], label=f"Octree, {np_thr}:{np_post}", linewidth=2) 
        
        ax.plot(nnls_data[nnls_upto]["t"], nnls_data[nnls_upto][f"M{M}_avg"], label=f"NNLS, {nnls_upto} mixed", linewidth=2) 

        if M == 4:
            ax.legend(fontsize=legend_size, framealpha=1.0)
        ax.grid()
        ax.tick_params(axis='both', labelsize=tick_size,)

        ax.set_xlabel(r"$\hat{t}$", fontsize=label_size)
        ax.set_ylabel(r"$\hat{M}_"+ f"{M}" +"$", fontsize=label_size)

        ax.set_xlim([-0.05, 15.0])

    if savefigs:
        fig.savefig(f"bkw_M_both_nnlsupto{nnls_upto}.pdf", bbox_inches="tight")

def get_bias_data(rundata, runs, MM):
    xmin_vals_ = [rundata[run[0]]["np_min"] for run in runs]
    xmean_vals_ = [rundata[run[0]]["np_avg"] for run in runs]
    xmax_vals_ = [rundata[run[0]]["np_max"] for run in runs]
    y_vals_ = [rundata[run[0]][f"M{MM}_bias"] for run in runs]
    
    return xmin_vals_, xmean_vals_, xmax_vals_, y_vals_

# plot bias in moments as function of average number of particles
fig = plt.figure(figsize=(18,6))

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

for M, ax in zip([4,8], [ax1,ax2]):
    xmin_vals, xmean_vals, xmax_vals, y_vals = get_bias_data(octree_mid_data, octree_mid_runs, M)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"Octree")

    xmin_vals, xmean_vals, xmax_vals, y_vals = get_bias_data(nnls_data, nnls_runs, M)
    ax.plot(xmean_vals, y_vals, '-o', linewidth=2, label=f"NNLS")

    if M == 4:
        ax.legend(fontsize=legend_size, framealpha=1.0)
    ax.grid()
    ax.tick_params(axis='both', labelsize=tick_size,)

    ax.set_xlabel(r"$\overline{N_p}$", fontsize=label_size)
    ax.set_ylabel(r"$\mathcal{B}(\hat{M}_"+ f"{M}" +")$", fontsize=label_size)

if savefigs:
    fig.savefig(f"bkw_bias_M_both.pdf", bbox_inches="tight")
