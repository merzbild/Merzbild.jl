# Plot 1D PIC data (potential/charge density/E-field) from a NetCDF file
# Usage example:
# python scripts/plot_1D_PIC.py --files scratch/data/picdata.nc --propname E --startt 1
# --plotname plots/Efield.png --labels "100 ppc" --scalex 1
# possible propnames: E or electric_field (E field), rho or charge_density (charge density), Phi or potential (potential)
# scalex is an optional list of values to scale the x axis of each plotted result by (since plotting is done via cell indices and not actual grid size)

from matplotlib import pyplot as plt
import argparse
from netCDF4 import Dataset
import numpy as np

parser = argparse.ArgumentParser(description='Plot electric field data from a 1D PIC simulation')
parser.add_argument("--files", nargs='+', required=True)
parser.add_argument("--propname", required=True)
parser.add_argument("--plotname", required=True)
parser.add_argument("--startt", required=True)
parser.add_argument("--labels", nargs='+')
parser.add_argument("--scalex", nargs='+', required=False)
args = parser.parse_args()


propname = args.propname

label_size = 22
tick_size = 18
legend_size = 18

fig = plt.figure(figsize=(12,10))

ax = fig.add_subplot(1,1,1)

labels = args.labels
if labels == None:
    labels = args.files
elif len(labels) != len(args.files):
    raise ValueError("Length of labels list should be the same as of the files' list!")

scales = [1 for _ in args.files]

propname_map = {
    "electric_field": "electric_field",
    "e": "electric_field",
    "potential": "potential",
    "phi": "potential",
    "charge_density": "charge_density",
    "rho": "charge_density"
}

if args.scalex is not None:
    scales = [float(s) for s in args.scalex]

startt = int(args.startt)

for label, file, sc in zip(labels, args.files, scales):
    ds = Dataset(file)

    try:
        pn = propname_map[propname.lower()]
        nx = np.shape(ds.variables[pn][:].data)[-1]
        nt = np.shape(ds.variables[pn][:].data)[0]
        x_arr = np.linspace(0, nx - 1, nx)

        data_arr = np.asarray(ds.variables[pn][:].data[:, :])
        for i in range(startt, nt):
            ax.plot(x_arr * sc, data_arr[i, :], label=f"{label}, nt={i}", linewidth=2)
    except KeyError:
        print(f"Variable {propname} not found in file {file}")
    ds.close()

ax.legend(fontsize=legend_size, ncol=2, framealpha=1.0)

ax.grid()
ax.tick_params(axis='both', labelsize=tick_size)

ax.set_xlabel(r"x", fontsize=label_size)
ax.set_ylabel(f"{propname}", fontsize=label_size)

fig.savefig(f"{args.plotname}", bbox_inches="tight")
