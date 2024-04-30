# Usage example:
# python scripts/plot_0D_1species.py --files test/data/bkw_20k_seed1234.nc test/data/bkw_vw_grid_seed1234.nc --propname=T --plotname plots/T_BKW_and_vw.png --labels "Fixed weight" "VW+merging"

from matplotlib import pyplot as plt
import argparse
from netCDF4 import Dataset
import numpy as np

parser = argparse.ArgumentParser(description='Plot data from multiple 0D 1-species files over time')
parser.add_argument("--files", nargs='+', required=True)
parser.add_argument("--propname", required=True)
parser.add_argument("--plotname", required=True)
parser.add_argument("--labels", nargs='+')
args = parser.parse_args()

labels = args.labels
if labels == None:
    labels = args.files
elif len(labels) != len(args.files):
    raise ValueError("Length of labels list should be the same as of the files' list!")

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

for label, file in zip(labels, args.files):

    ds = Dataset(file)

    if propname in ["ndens", "np", "T"]:
        ax.plot(ds.variables["timestep"][:].data, ds.variables[propname][:].data[:, 0, 0], label=label, linewidth=2)
    elif propname.lower() in v_map:
        ax.plot(ds.variables["timestep"][:].data, ds.variables["v"][:].data[:, 0, 0, v_map[propname]], label=label, linewidth=2)
    elif propname.startswith("M"):
        momno = int(propname[1:])
        momno_list = np.asarray(ds.variables["moment_powers"])
        mom_id = np.where(momno_list == momno)

        ax.plot(ds.variables["timestep"][:].data, ds.variables["moments"][:].data[:, 0, 0, mom_id[0]], label=label, linewidth=2)

    ds.close()
    
ax.legend(fontsize=legend_size, ncol=2, framealpha=1.0)

ax.grid()
ax.tick_params(axis='both', labelsize=tick_size)

ax.set_xlabel(r"Timestep", fontsize=label_size)
ax.set_ylabel(f"{propname}", fontsize=label_size)

fig.savefig(f"{args.plotname}", bbox_inches="tight")
