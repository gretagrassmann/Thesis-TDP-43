#!/usr/bin/env python
# coding: utf-8


import os, sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as mpl

import pandas as pd
from scipy.spatial import distance_matrix

import ZernikeFunc as ZF
import SurfaceFunc as SF

import time

SAVE_FIG = 0
PLOT = 0

Npixel = 25  # the plane edge in pixels...

Daa = 3.  # the maximum distance between antigen and antibody binding site atoms..
Rs = 1.  # the radius of the sphere that includes the patch..
Dpp = 1.  # the distance between points of the same patch (needed to remove islands)

ZOrder = 5  # the Zernike expansion order..

pdb_file = "..\RRM2\cluster2.dms"
verso = float(1)
#respath = "..\RRM2\cluster2"
name_ = "Cluster 2 of RRM2 EXPERIMENT"
Npoint = int(1)


if (verso == 1):
    ver = "up"
elif (verso == -1):
    ver = "down"

# loading surface file..

surf_name = pdb_file

print("Processing protein %s with verso: (%s)\n" % (name_, verso))

try:
    surf_ = pd.read_csv(surf_name)
except:
    print("Given protein is not present!")
    exit()

lag = len(surf_["x"])


surf = np.zeros((lag, 6))
surf[:, :] = surf_[["x", "y", "z", "Nx", "Ny", "Nz"]]

surf_obj = SF.Surface(surf[:, :], patch_num=0, r0=Rs, theta_max=45)

ltmp = np.shape(surf)[0]

index_possible_area = np.arange(ltmp)[::Npoint]  #MM]
ltmp = len(index_possible_area)
###     QUA DEVO FARE LA COSA DEL PRODOTTO SCALARE, E POI DIRE IN ltmp
###     QUANTI PUNTI STO GUARDANDO


for i in index_possible_area:
    sys.stderr.write(("\r Processing point %i" % (i)))
    sys.stderr.flush()
    patch, mask = surf_obj.BuildPatch(point_pos=index_possible_area[i], Dmin=.5)
    surf_obj.real_br = mask
    prod = patch[0,3:6]
    for j in range(1,len(patch)):
        prod = prod*patch[j,3:6]
    prod_scal_patch = sum(prod) #QUESTO MI DICE QUANTO LA PATCH E' PIATTA

    if abs(prod_scal_patch) > 0.01: # SALVO GLI INDICI DEI PUNTI DOVE SERVE CAMPIONAMENTO PIU' FITTO
        slope_index = np.zeros((len(patch)))
        for k in range(len(patch)):
            slope_index[k] = np.where((surf[:,0] == patch[k,0])&(surf[:,1] == patch[k,1])&(surf[:,2] == patch[k,2]))[0]

        index_possible_area = [i for i in index_possible_area if i not in slope_index]


###     index_possible_area DEVO DEFINIRLO IN MANIERA DIVERSA,
###     METTENDOCI GLI INDICI DEI PUNTI CHE NON HO ELIMINATO CON IL PRODOTTO SCALARE

#index_possible_area =

###     DA RIGA 84 DI SurfaceScreening_general.py