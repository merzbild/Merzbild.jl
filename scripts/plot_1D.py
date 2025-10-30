# Usage example:
# python scripts/plot_1D.py --file scratch/data/couette_0.0005_50_500.0_300.0_1000.nc --propname T --startt 1
# --plotname plots/couette_T.png --labels "Couette"
# velocity is plotted with an offset of a linear profile that goes from -500 to 500 m/s
# spartafile: optional path to SPARTA output file
# spartavarid: number of the variable in the file to plot
# it is assumed that the sparta cells are sorted, are the same ones as used in Merzbild.jl
# the output should begin at line 10, line 9 is the header line
# ITEM: CELLS id f_1[1] f_1[2] f_1[3] f_1[4] f_1[5] f_1[6] - to plot id set varid to 0, f_1[1] - set to 1, etc.

from matplotlib import pyplot as plt
import argparse
from netCDF4 import Dataset
import numpy as np

parser = argparse.ArgumentParser(description='Plot data from a 1D single-species simulation')
parser.add_argument("--files", nargs='+', required=True)
parser.add_argument("--propname", required=True)
parser.add_argument("--plotname", required=True)
parser.add_argument("--startt", required=True)
parser.add_argument("--labels", nargs='+')
parser.add_argument("--spartafile", required=False)
parser.add_argument("--spartavarid", required=False)
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

labels = args.labels
if labels == None:
    labels = args.files
elif len(labels) != len(args.files):
    raise ValueError("Length of labels list should be the same as of the files' list!")

startt = int(args.startt)

for label, file in zip(labels, args.files):
    ds = Dataset(file)

    if propname in ["ndens", "np", "T"]:
        nx = np.shape(ds.variables[propname][:].data)[-1]
        nt = np.shape(ds.variables[propname][:].data)[0]
        x_arr = np.linspace(1, nx, nx)

        data_arr = np.asarray(ds.variables[propname][:].data[:, 0, :])
        for i in range(startt, nt):
            ax.plot(x_arr, data_arr[i, :], label=f"{label}, nt={i}", linewidth=2)
    elif propname.lower() in v_map:
        nx = np.shape(ds.variables["v"][:].data)[-2]
        nt = np.shape(ds.variables["v"][:].data)[0]
        x_arr = np.linspace(1, nx, nx)

        data_arr = np.asarray(ds.variables["v"][:].data[:, 0, :, v_map[propname]])
        voffset = np.linspace(-500, 500, nx)
        for i in range(startt,nt):
            ax.plot(x_arr, data_arr[i, :] - voffset, label=f"{label}, nt={i}", linewidth=2)
        # ax.plot(ds.variables["timestep"][:].data, ds.variables["v"][:].data[:, 0, 0, v_map[propname]], label=label, linewidth=2)

    ds.close()
        
if args.spartafile:
    if not args.spartavarid:
        print("No SPARTA variable number given!")
    else:
        varid = int(args.spartavarid)
        sp_data = []
        with open(args.spartafile) as f:
            for i, line in enumerate(f):
                if i>=9:
                    lsp = line.split()
                    sp_data.append(float(lsp[varid]))
        if propname in ["ndens", "np", "T"]:
            ax.plot(x_arr, sp_data, label=f"SPARTA, nt={i}", linewidth=2)
        else:
            voffset = np.linspace(-500, 500, nx)
            ax.plot(x_arr, np.array(sp_data) - voffset, label=f"SPARTA, nt={i}", linewidth=2)

ax.legend(fontsize=legend_size, ncol=2, framealpha=1.0)

ax.grid()
ax.tick_params(axis='both', labelsize=tick_size)

ax.set_xlabel(r"x", fontsize=label_size)
ax.set_ylabel(f"{propname}", fontsize=label_size)

fig.savefig(f"{args.plotname}", bbox_inches="tight")
