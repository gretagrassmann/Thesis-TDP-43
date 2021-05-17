import os, sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as mpl
import math
import pandas as pd
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt

import ZernikeFunc as ZF
import SurfaceFunc as SF
import my_functions

Npixel = 25  # the plane edge in pixels...

Daa = 3.  # the maximum distance between antigen and antibody binding site atoms..
Dpp = 1.  # the distance between points of the same patch (needed to remove islands)

verso = float(1)
if (verso == 1):
    ver = "up"
elif (verso == -1):
    ver = "down"
ZOrder = 20  # the Zernike expansion order..

pdb_file = "..\RRM2\cluster2.dms"
respath = "..\RRM2\cluster2\R_s_8"
name_ = "Cluster 2 of RRM2 EXPERIMENT"
Npoint = int(1)

surf_ = pd.read_csv(pdb_file)
lag = len(surf_["x"])
surf = np.zeros((lag, 6))
surf[:, :] = surf_[["x", "y", "z", "Nx", "Ny", "Nz"]]
surf_obj = SF.Surface(surf[:, :], patch_num=0, r0=.6, theta_max=45)
ltmp = np.shape(surf)[0]

            # SCREENING E ZERNIKE PER OGNI SINGOLO PUNTO
print("Vuoi fare Zernike per ogni singolo punto? y o n?")
o = input()
if o == "y":
    index_all_possible_area = np.arange(ltmp)[::Npoint]

    res_inv_ = my_functions.ZernikeCoefficients(index_all_possible_area, surf_obj,verso, Npixel,ZOrder, respath, "total", 1)

   # np.savetxt("%s/zernike_total.dat" % (respath), res_inv_,fmt="%.4e")
   # np.savetxt("%s/Selected_points_1out%d_indexes_%s_verso_%s.dat" % (respath, Npoint, name_, ver), index_all_possible_area)


            # SCREENING PER DIVERSI VALORI DEI PARAMETRI

print("Vuoi fare lo screening per diversi valori di alpha? y o n?")
oo = input()
if oo == "y":

    Rs_select = 8.
    #alpha = [.01, .1, 1]
    #alpha = [.02, .03, .04, .05, .06, .08, .09]
    #alpha = [.065, .07, .071, .072, .073, .074, .075, .076, .077, .078, .079, .085]
    alpha = [1.1]
    step =1

    for i in alpha:
        index_possible_area = my_functions.CosScanning(surf,surf_obj,Rs_select,i,step,respath)



#zernike = my_functions.ZernikeCoefficients(index_possible_area,surf_obj,verso,Npixel,ZOrder)
#np.savetxt("{}\zernike_alpha_{}.txt".format(respath,alpha),zernike)


        ## PCA PER CIASCUNA LISTA DI COEFFICIENTI DI ZERNIKE

#pca = my_functions.PCA(zernike)
#np.savetxt("{}\pca_alpha_{}.txt".format(respath, Rs_select, bin), pca)

#x_grid = np.linspace(min(pca[0]),max(pca[0]), num=int(len(pca[0])/20))
#y_grid = np.linspace(min(pca[1]),max(pca[1]), num=int(len(pca[1])/20))

#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.set_title("PCA: 2D projection of the Zernike coefficients")
#ax1.set_xlabel('Projection on eigenvector 1 (nm)')
#ax1.set_ylabel('Projection on eigenvector 2 (nm)')
#ax1.xaxis.set_ticks(x_grid)
#ax1.yaxis.set_ticks(y_grid)
#ax1.grid(True)
    #ax1.set_yticklabels([])
    #ax1.set_xticklabels([])
#sc = plt.scatter(pca[0], pca[1])
#leg = ax1.legend()
#plt.show()