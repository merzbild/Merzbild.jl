import numpy as np
from netCDF4 import Dataset
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

savefigs = False
savefigs = True  # - uncomment to save plots

# this is the timestep after which the averaging of results was turned on in the simulation
# assumes filenames end with "_after{avg_start}.nc"
avg_start = 1000

# x-axis (mm)
x_ax = np.linspace(0, 5, 1000)
# output is # of physical particles in cell, we need number density, hence the scaling factor
nd_c = 1.0 / (5e-3/1000)

# this is the list of all octree runs for error-vs-average number of particles plots
oc_runs = ["45_38", "69_58", "105_88", "146_122", "200_166"]
# this is the list of all NNLS runs for error-vs-average number of particles plots
nnls_runs = ["45_38", "69_58", "105_88", "146_122", "200_166"]

# prefix to directory with data
pref = "scratch/data/"

# 3 datasets used to plot example results, 69:58 (approximately) merging
data_nnls = f"{pref}avg_Fourier_NNLS_ntc_0.005_1000_0.0_300.0_600.0_69_58_after{avg_start}.nc"
data_oc = f"{pref}avg_Fourier_octree_ntc_0.005_1000_0.0_300.0_600.0_69_58_after{avg_start}.nc"
data_fw = f"{pref}avg_Fourier_FW_ntc_0.005_1000_0.0_300.0_600.0_after{avg_start}.nc"

# plotting parameters
label_size = 24
tick_size = 20
legend_size = 20

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["axes.linewidth"] = 0.8

# load octree data for bulk and surface properties
oc_T = {}
oc_np = {}
oc_n = {}
oc_s_p = {}
oc_s_ke = {}
oc_s_fi = {}
for oc in oc_runs:    
    ds = Dataset(f"{pref}avg_Fourier_octree_ntc_0.005_1000_0.0_300.0_600.0_{oc}_after{avg_start}.nc")
    oc_T[oc] = np.asarray(ds["T"][:][0,0,:])
    oc_n[oc] = np.asarray(ds["ndens"][:][0,0,:])
    oc_np[oc] = np.asarray(ds["np"][:][0,0,:])
    ds.close()

    ds = Dataset(f"{pref}avg_Fourier_octree_ntc_0.005_1000_0.0_300.0_600.0_{oc}_surf_after{avg_start}.nc")
    oc_s_p[oc] = np.asarray(ds["normal_pressure"][0,0,:])
    oc_s_ke[oc] = np.asarray(ds["kinetic_energy_flux"][0,0,:])
    oc_s_fi[oc] = np.asarray(ds["flux_incident"][0,0,:])
    ds.close()

# load NNLS data for bulk and surface properties
nnls_T = {}
nnls_n = {}
nnls_np = {}
nnls_s_p = {}
nnls_s_ke = {}
nnls_s_fi = {}

for nnls in nnls_runs:
    ds = Dataset(f"{pref}avg_Fourier_NNLS_ntc_0.005_1000_0.0_300.0_600.0_{nnls}_after{avg_start}.nc")
    nnls_T[nnls] = np.asarray(ds["T"][:][0,0,:])
    nnls_n[nnls] = np.asarray(ds["ndens"][:][0,0,:])
    nnls_np[nnls] = np.asarray(ds["np"][:][0,0,:])
    ds.close()

    ds = Dataset(f"{pref}avg_Fourier_NNLS_ntc_0.005_1000_0.0_300.0_600.0_{nnls}_surf_after{avg_start}.nc")
    nnls_s_p[nnls] = np.asarray(ds["normal_pressure"][0,0,:])
    nnls_s_ke[nnls] = np.asarray(ds["kinetic_energy_flux"][0,0,:])
    nnls_s_fi[nnls] = np.asarray(ds["flux_incident"][0,0,:])
    ds.close()

# load reference values for bulk and surface properties
ds = Dataset(f"{pref}avg_Fourier_FW_ntc_0.005_1000_0.0_300.0_600.0_after{avg_start}.nc")
ref_T = np.asarray(ds["T"][:][0,0,:])
ref_n = np.asarray(ds["ndens"][:][0,0,:])
ds.close()

ds = Dataset(f"{pref}avg_Fourier_FW_ntc_0.005_1000_0.0_300.0_600.0_surf_after{avg_start}.nc")
ref_s_p = np.asarray(ds["normal_pressure"][0,0,:])
ref_s_ke = np.asarray(ds["kinetic_energy_flux"][0,0,:])
ref_s_fi = np.asarray(ds["flux_incident"][0,0,:])
ds.close()

# plot example profiles of T and number density
fig = plt.figure(figsize=(20,6))

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

# plot results
for (ax, v, cf) in zip([ax1,ax2], ["T", "ndens"], [1.0, nd_c]):
    ds = Dataset(data_fw)
    ax.plot(x_ax, ds[v][:][0,0,:]*cf, label="Fixed weight", linewidth=2, color="k") 
    ds.close()

    ds = Dataset(data_oc)
    ax.plot(x_ax, ds[v][:][0,0,:]*cf, label="Octree", linewidth=2, color="tab:blue") 
    ds.close()


    ds = Dataset(data_nnls)
    ax.plot(x_ax, ds[v][:][0,0,:]*cf, label="NNLS", linewidth=2, color="tab:orange") 
    ds.close()

    # ax.legend(fontsize=legend_size, framealpha=1.0)
    ax.set_xlabel(r"$x$, mm", fontsize=label_size)
    ax.grid()
    ax.tick_params(axis='both', labelsize=tick_size,)
    ax.set_xlim([-0.05, 5.05])

ax1.legend(fontsize=legend_size, framealpha=1.0)
ax1.set_ylabel(r"$T$, K", fontsize=label_size)
ax2.set_ylabel(r"$n$, m$^{-3}$", fontsize=label_size)

# plot zoomed-in insets
# Create a zoomed inset axes
axins = zoomed_inset_axes(ax1, zoom=2, loc='lower right', borderpad=3.5)  # zoom factor and location
ds = Dataset(data_fw)
axins.plot(x_ax, ds["T"][:][0,0,:], label="Fixed weight", linewidth=2, color="k") 
ds.close()

ds = Dataset(data_oc)
axins.plot(x_ax, ds["T"][:][0,0,:], label="Octree", linewidth=2, color="tab:blue") 
ds.close()

ds = Dataset(data_nnls)
axins.plot(x_ax, ds["T"][:][0,0,:], label="NNLS", linewidth=2, color="tab:orange") 
ds.close()

# Set the limits for the zoomed region
x1, x2, y1, y2 = 900*5e-3, 1000*5e-3, 560, 600  # Define the region to zoom in
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# # Mark the zoomed region on the main plot
mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")

# Create a zoomed inset axes
axins2 = zoomed_inset_axes(ax2, zoom=2, loc='upper right', borderpad=3.5)  # zoom factor and location

ds = Dataset(data_fw)
axins2.plot(x_ax, ds["ndens"][:][0,0,:]*nd_c, label="Fixed weight", linewidth=2, color="k") 
ds.close()

ds = Dataset(data_oc)
axins2.plot(x_ax, ds["ndens"][:][0,0,:]*nd_c, label="Octree", linewidth=2, color="tab:blue") 
ds.close()

ds = Dataset(data_nnls)
axins2.plot(x_ax, ds["ndens"][:][0,0,:]*nd_c, label="NNLS", linewidth=2, color="tab:orange") 
ds.close()

# # Set the limits for the zoomed region
x1, x2, y1, y2 = 900*5e-3, 1000*5e-3, 1.15e17*nd_c, 1.3e17*nd_c  # Define the region to zoom in
axins2.set_xlim(x1, x2)
axins2.set_ylim(y1, y2)

# # Mark the zoomed region on the main plot
mark_inset(ax2, axins2, loc1=1, loc2=3, fc="none", ec="0.5")

if savefigs:
    fig.savefig(f"fourier_T_ndens_profile.pdf", bbox_inches="tight")

# plot bias in temperature and density as function of number of particles
fig = plt.figure(figsize=(20,6))

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

avg_T = np.mean(ref_T)
avg_n = np.mean(ref_n)

x_np = []
y_n = []
for ocrun in oc_runs:
    x_np.append(np.mean(oc_np[ocrun]))
    y_n.append(np.sqrt(np.sum((oc_T[ocrun]/ref_T - 1)**2)/1000) * 100)
ax1.plot(x_np, y_n, label="Octree", color="tab:blue", linewidth=2, marker="o")

x_np = []
y_n = []
for nnls in nnls_runs:
    x_np.append(np.mean(nnls_np[nnls]))
    y_n.append(np.sqrt(np.sum((nnls_T[nnls]/ref_T - 1)**2)/1000) * 100)
ax1.plot(x_np, y_n, label="NNLS", color="tab:orange", linewidth=2, marker="o")

x_np = []
y_n = []
for ocrun in oc_runs:
    x_np.append(np.mean(oc_np[ocrun]))
    y_n.append(np.sqrt(np.sum((oc_n[ocrun]/ref_n - 1)**2)/1000) * 100)
ax2.plot(x_np, y_n, label="Octree", color="tab:blue", linewidth=2, marker="o")

x_np = []
y_n = []
for nnls in nnls_runs:
    x_np.append(np.mean(nnls_np[nnls]))
    y_n.append(np.sqrt(np.sum((nnls_n[nnls]/ref_n - 1)**2)/1000) * 100)
ax2.plot(x_np, y_n, label="NNLS", color="tab:orange", linewidth=2, marker="o")
    

ax1.legend(fontsize=legend_size, framealpha=1.0)

for ax in [ax1,ax2]:
    ax.grid()
    ax.tick_params(axis='both', labelsize=tick_size)

    ax.set_xlabel(r"$\overline{N}_p$", fontsize=label_size)
ax1.set_ylabel(r"$\overline\mathcal{B}(T)$, \%", fontsize=label_size)
ax2.set_ylabel(r"$\overline\mathcal{B}(n)$, \%", fontsize=label_size)

if savefigs:
    fig.savefig(f"fourier_bias_T_N.pdf", bbox_inches="tight")

# plot bias in incident number density flux
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(1,1,1)

wall = 1

x_np = []
y_n = []
for ocrun in oc_runs:
    x_np.append(np.mean(oc_np[ocrun]))
    y_n.append(np.sqrt(np.sum((oc_s_fi[ocrun][wall]/ref_s_fi[wall] - 1)**2)) * 100)
ax.plot(x_np, y_n, label="Octree", color="tab:blue", linewidth=2, marker="o")

x_np = []
y_n = []
for nnls in nnls_runs:
    x_np.append(np.mean(nnls_np[nnls]))
    y_n.append(np.sqrt(np.sum((nnls_s_fi[nnls][wall]/ref_s_fi[wall] - 1)**2)) * 100)
ax.plot(x_np, y_n, label="NNLS", color="tab:orange", linewidth=2, marker="o")

ax.legend(fontsize=legend_size, framealpha=1.0)

ax.grid()
ax.tick_params(axis='both', labelsize=tick_size)

ax.set_xlabel(r"$\overline{N}_p$", fontsize=label_size)
ax.set_ylabel(r"$\overline\mathcal{B}(\dot{n}_{in})$, \%", fontsize=label_size)

if savefigs:
    fig.savefig(f"fourier_bias_incidentflux.pdf", bbox_inches="tight")


# plot bias in pressure
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(1,1,1)

wall = 1

x_np = []
y_n = []
for ocrun in oc_runs:
    x_np.append(np.mean(oc_np[ocrun]))
    y_n.append(np.sqrt(np.sum((oc_s_p[ocrun][wall]/ref_s_p[wall] - 1)**2)) * 100)
ax.plot(x_np, y_n, label="Octree", color="tab:blue", linewidth=2, marker="o")

x_np = []
y_n = []
for nnls in nnls_runs:
    x_np.append(np.mean(nnls_np[nnls]))
    y_n.append(np.sqrt(np.sum((nnls_s_p[nnls][wall]/ref_s_p[wall] - 1)**2)) * 100)

ax.plot(x_np, y_n, label="NNLS", color="tab:orange", linewidth=2, marker="o")

ax.legend(fontsize=legend_size, framealpha=1.0)

ax.grid()
ax.tick_params(axis='both', labelsize=tick_size)

ax.set_xlabel(r"$\overline{N}_p$", fontsize=label_size)
ax.set_ylabel(r"$\overline\mathcal{B}(p_s)$, \%", fontsize=label_size)

if savefigs:
    fig.savefig(f"fourier_bias_pressure.pdf", bbox_inches="tight")

# plot bias in kinetic energy flux
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(1,1,1)

wall = 1

x_np = []
y_n = []
for ocrun in oc_runs:
    x_np.append(np.mean(oc_np[ocrun]))
    y_n.append(np.sqrt(np.sum((oc_s_ke[ocrun][wall]/ref_s_ke[wall] - 1)**2)) * 100)
ax.plot(x_np, y_n, label="Octree", color="tab:blue", linewidth=2, marker="o")

x_np = []
y_n = []
for nnls in nnls_runs:
    x_np.append(np.mean(nnls_np[nnls]))
    y_n.append(np.sqrt(np.sum((nnls_s_ke[nnls][wall]/ref_s_ke[wall] - 1)**2)) * 100)

ax.plot(x_np, y_n, label="NNLS", color="tab:orange", linewidth=2, marker="o")

ax.legend(fontsize=legend_size, framealpha=1.0)

ax.grid()
ax.tick_params(axis='both', labelsize=tick_size)

ax.set_xlabel(r"$\overline{N}_p$", fontsize=label_size)
ax.set_ylabel(r"$\overline\mathcal{B}(q_s)$, \%", fontsize=label_size)

if savefigs:
    fig.savefig(f"fourier_bias_heat.pdf", bbox_inches="tight")