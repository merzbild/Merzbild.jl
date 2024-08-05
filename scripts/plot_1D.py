# Usage example:
# python scripts/plot_1D.py --file scratch/data/couette_0.0005_50_500.0_300.0_1000.nc --propname T --startt 1
# --plotname plots/couette_T.png

from matplotlib import pyplot as plt
import argparse
from netCDF4 import Dataset
import numpy as np

parser = argparse.ArgumentParser(description='Plot data from a 1D single-species simulation')
parser.add_argument("--file", required=True)
parser.add_argument("--propname", required=True)
parser.add_argument("--plotname", required=True)
parser.add_argument("--startt", required=True)
args = parser.parse_args()


propname = args.propname

label_size = 22
tick_size = 18
legend_size = 18

fig = plt.figure(figsize=(12,10))

ax = fig.add_subplot(1,1,1)

v_map = {
    "vx": 0,
    "vy": 1,
    "vz": 2,
    "v_x": 0,
    "v_y": 1,
    "v_z": 2,
}

ds = Dataset(args.file)
startt = int(args.startt)

if propname in ["ndens", "np", "T"]:
    nx = np.shape(ds.variables[propname][:].data)[-1]
    nt = np.shape(ds.variables[propname][:].data)[0]
    x_arr = np.linspace(1, nx, nx)

    data_arr = np.asarray(ds.variables[propname][:].data[:, 0, :])
    for i in range(startt, nt):
        ax.plot(x_arr, data_arr[i, :], label=f"nt={i}", linewidth=2)
elif propname.lower() in v_map:
    nx = np.shape(ds.variables["v"][:].data)[-2]
    nt = np.shape(ds.variables["v"][:].data)[0]
    x_arr = np.linspace(1, nx, nx)

    data_arr = np.asarray(ds.variables["v"][:].data[:, 0, :, v_map[propname]])
    voffset = np.linspace(-500, 500, nx)
    for i in range(startt,nt):
        ax.plot(x_arr, data_arr[i, :] - voffset, label=f"nt={i}", linewidth=2)
    # ax.plot(ds.variables["timestep"][:].data, ds.variables["v"][:].data[:, 0, 0, v_map[propname]], label=label, linewidth=2)

ds.close()
    
ax.legend(fontsize=legend_size, ncol=2, framealpha=1.0)

ax.grid()
ax.tick_params(axis='both', labelsize=tick_size)

ax.set_xlabel(r"Timestep", fontsize=label_size)
ax.set_ylabel(f"{propname}", fontsize=label_size)

fig.savefig(f"{args.plotname}", bbox_inches="tight")
