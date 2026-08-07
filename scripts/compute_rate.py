# python scripts/compute_rate.py --filename test/data/tmp_ion_no_es.nc --tstart 100000 --tend 290000 
# --nid 0 --ionid 1 --eid 2 --dt 5e-14
import argparse
from netCDF4 import Dataset
import numpy as np

parser = argparse.ArgumentParser(description='Compute ionization rate and electron temperature' \
                                             ' of species nid colliding with electrons with index eid' \
                                             ' in a 0D simulation; averaging over timesteps from tstart to tend; if nseed is provided' \
                                             ' additional files with names filename_seed<seed>.nc are read and averaged, where 1 <= seed < nseed-1')
parser.add_argument("--filename", required=True)
parser.add_argument("--tstart", required=True)
parser.add_argument("--tend", required=True)
parser.add_argument("--nid", required=True)
parser.add_argument("--ionid", required=True)
parser.add_argument("--eid", required=True)
parser.add_argument("--dt", required=True)
parser.add_argument("-nseed", default=0)

args = parser.parse_args()


tstart = int(args.tstart)
tend = int(args.tend)

nid = int(args.nid)
ionid = int(args.ionid)
eid = int(args.eid)
nseed = int(args.nseed)
print(nseed)

if nseed == 0:
    ds = Dataset(args.filename)
    T_e = np.mean(np.asarray(ds["T"])[tstart:tend, eid, 0]) / 11605.0
    T_e_std = np.std(np.asarray(ds["T"])[tstart:tend, eid, 0]) / 11605.0

    densities = np.asarray(ds["ndens"])

    ds.close()

    # dn_ion / dt = k_ion * n_neutral * n_e
    dn_ion = densities[tstart:tend-1, ionid, 0] - densities[tstart-1:tend-2, ionid, 0]
    k_ion = dn_ion / (float(args.dt) * densities[tstart-1:tend-2, nid, 0] * densities[tstart-1:tend-2, eid, 0])

    k_ion_avg = np.mean(k_ion)
    k_ion_std = np.std(k_ion)

    print(k_ion_avg, "+-", k_ion_std, T_e, "+-", T_e_std)
else:
    ds = Dataset(args.filename)
    T_e = np.asarray(ds["T"])[tstart:tend, eid, 0] / 11605.0

    densities = np.asarray(ds["ndens"])

    ds.close()

    for i in range(1, nseed):
        ds = Dataset(args.filename.replace(".nc", f"_seed{i}.nc"))
        T_e += np.asarray(ds["T"])[tstart:tend, eid, 0] / 11605.0
        # T_e_std += np.std(np.asarray(ds["T"])[tstart:tend, eid, 0]) / 11605.0
    
        densities += np.asarray(ds["ndens"])
        ds.close()

    T_e_std = np.std(T_e)
    T_e = np.mean(T_e)
    T_e /= nseed
    T_e_std /= nseed
    densities /= nseed

    # dn_ion / dt = k_ion * n_neutral * n_e
    dn_ion = densities[tstart:tend-1, ionid, 0] - densities[tstart-1:tend-2, ionid, 0]
    k_ion = dn_ion / (float(args.dt) * densities[tstart-1:tend-2, nid, 0] * densities[tstart-1:tend-2, eid, 0])

    k_ion_avg = np.mean(k_ion)
    k_ion_std = np.std(k_ion)

    print(k_ion_avg, "+-", k_ion_std, T_e, "+-", T_e_std)
    