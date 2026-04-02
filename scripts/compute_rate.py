import argparse
from netCDF4 import Dataset
import numpy as np

parser = argparse.ArgumentParser(description='Compute ionization rate and electron temperature' \
                                             ' of species nid colliding with electrons with index eid' \
                                             ' in a 0D simulation; averaging over timesteps from tstart to tend')
parser.add_argument("--filename", required=True)
parser.add_argument("--tstart", required=True)
parser.add_argument("--tend", required=True)
parser.add_argument("--nid", required=True)
parser.add_argument("--ionid", required=True)
parser.add_argument("--eid", required=True)
parser.add_argument("--dt", required=True)

args = parser.parse_args()

ds = Dataset(args.filename)

tstart = int(args.tstart)
tend = int(args.tend)

nid = int(args.nid)
ionid = int(args.ionid)
eid = int(args.eid)

T_e = np.mean(np.asarray(ds["T"])[tstart:tend, eid, 0]) / 11605.0

densities = np.asarray(ds["ndens"])

ds.close()

print(densities.shape)

# dn_ion / dt = k_ion * n_neutral * n_e
dn_ion = densities[tstart:tend-1, ionid, 0] - densities[tstart-1:tend-2, ionid, 0]
k_ion = dn_ion / (float(args.dt) * densities[tstart-1:tend-2, nid, 0] * densities[tstart-1:tend-2, eid, 0])

k_ion_avg = np.mean(k_ion)

print(k_ion_avg, T_e)