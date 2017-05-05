"""
Split the halo catalogs into mass bins.
"""
import numpy as np
import matplotlib.pyplot as plt
import os, sys

indices   = [6, 7, 8, 9]
redshifts = [1.0, 0.5, 0.25, 0.0]
lM_edges  = [13.0, 13.1, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.5, 15.0, 16.0]

#First read in the smallest catalog and look at some shit in there
#datapath = "/calvin1/tmcclintock/fox_data/richness_halos/rich_snapdir_%03d/reduced_richness_halos_%03d" #Only valid on CALVIN
datapath = "./reduced_catalogs/reduced_richness_halos_%03d"

#The path to the output
outpath  = "./mass_catalogs/SV_halos_m%d_z%.2f.txt"

header = "X[Mpc/h] Y[Mpc/h] Z[Mpc/h] Np M200b[Msun/h] Richness"
fmt = "%.5f %.5f %.5f %d %.4e %.2f"

for i in range(len(indices)):
    ind = indices[i]
    red = redshifts[i]
    #data = np.genfromtxt(datapath%(ind,ind)) #For on CALVIN
    data = np.genfromtxt(datapath%ind)
    M = data[:,4]
    lM = np.log10(M)
    for j in range(len(lM_edges)-1):
        lo = lM_edges[j]
        hi = lM_edges[j+1]
        sample = (lM>=lo)&(lM < hi)
        outdata = data[sample]
        np.savetxt(outpath%(j,red), outdata, fmt=fmt, header=header)
