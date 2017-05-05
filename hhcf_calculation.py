"""
Calculate the halo-halo correlation functions.

This is in advance of calculating the HMCFs, since that takes forever.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc("text", usetex=True, fontsize=24)
from Corrfunc.theory.xi import xi

halopath = "./mass_catalogs/SV_halos_m%d_z%.2f.txt"

boxsize = 1050.0
nthreads = 8
nbins = 50
#NOTE - these are much finer than what we will end up using for the deltasigmas
#this is because we can calculate the bin-averaged deltasigmas later
bins = np.logspace(np.log10(0.01), np.log10(150.0), nbins+1)
R = (bins[:-1]+bins[1:])/2.

do_hhcf = False
do_calc = True
do_jk   = True
do_jkcalc = True

indices   = [6, 7, 8, 9]
redshifts = [1.0, 0.5, 0.25, 0.0]
lM_edges  = [13.0, 13.1, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.5, 15.0, 16.0]
zstrings = ["1.0", "0.5", "0.25", "0.0"]

c = np.linspace(1.0, 0.0, len(lM_edges)-1)
cmaps = ['Reds','Oranges','Greens', 'Blues']

if do_hhcf:
    for i in range(len(indices)):
        ind = indices[i]
        red = redshifts[i]
        cmap = plt.get_cmap(cmaps[i])
        for j in range(len(lM_edges)-1):
            if do_calc:
                data = np.genfromtxt(halopath%(j,red))
                if not data.size: continue
                x,y,z,N,M,lam = data.T
                result = xi(boxsize, nthreads, bins, x, y, z)
                np.savetxt("results/hhcf_m%d_z%.2f.txt"%(j,red), result)
            result = np.loadtxt("results/hhcf_m%d_z%.2f.txt"%(j,red))
            hhcf = result[:,3]
            plt.loglog(R, hhcf, c=cmap(c[j]), label=r"$z=%.2f\ m%d$"%(red,j))
        plt.legend(loc=0, fontsize=10)
        plt.xlabel(r"$R\ [{\rm Mpc/h}]$", fontsize=24)
        plt.ylabel(r"$\xi_{\rm hh}$", fontsize=24)
        plt.subplots_adjust(bottom=0.17, left=0.2)
        #plt.gcf().savefig("figures/hhcf_richnesses_z%0.2f.png"%z)
        plt.show()
        plt.clf()

if do_jk:
    ndivs = 8
    step  = boxsize/ndivs
    for i in range(len(indices)):
        ind = indices[i]
        red = redshifts[i]
        cmap = plt.get_cmap(cmaps[i])
        for j in range(len(lM_edges)-1):
            if do_jkcalc:
                data = np.genfromtxt(halopath%(j,red))
                if not data.size: continue
                x,y,z,N,M,lam = data.T
                #Now need split the data into little sub boxes
                #Then do auto correlations
                #Then do cross correlations
                #Then resum to get the totals
                #Then remove boxes one at a time to get the LOO samples
                #Then we can find the covariance matrix
