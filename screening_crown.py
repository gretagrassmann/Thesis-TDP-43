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

Npixel = 25  # the plane edge in pixels...

Daa = 3.  # the maximum distance between antigen and antibody binding site atoms..
Dpp = 1.  # the distance between points of the same patch (needed to remove islands)

ZOrder = 20  # the Zernike expansion order..

pdb_file = "..\RRM2\cluster2.dms"
respath = "..\RRM2\cluster2\positive_crown_sampling"
name_ = "Cluster 2 of RRM2 EXPERIMENT"
Npoint = int(1)

Rs_select = 4.  # the radius of the sphere that includes the patch..
Max_scalar  = 1.e-60
R_crown = 1.e60

# loading surface file..

surf_name = pdb_file
surf_ = pd.read_csv(surf_name)

lag = len(surf_["x"])

surf = np.zeros((lag, 6))
surf[:, :] = surf_[["x", "y", "z", "Nx", "Ny", "Nz"]]

surf_obj = SF.Surface(surf[:, :], patch_num=0, r0=Rs_select, theta_max=45)

ltmp = np.shape(surf)[0]

index_possible_area = np.arange(1,ltmp,3)
#index_possible_area = np.arange(ltmp)[:Npoint]  #MM]

for i in index_possible_area:
    ltmp = len(index_possible_area)
    sys.stderr.write(("\r Processing point %i" % (i)))
    sys.stderr.flush()
    patch, mask = surf_obj.BuildPatch(point_pos=i, Dmin=.5)
    surf_obj.real_br = mask
    prod = patch[0,3:6]
    for j in range(1,len(patch)):
        prod = prod*patch[j,3:6]
    prod_scal_patch = sum(prod) #QUESTO MI DICE QUANTO LA PATCH E' PIATTA
    if abs(prod_scal_patch) > Max_scalar: # SALVO GLI INDICI DEI PUNTI DOVE SERVE CAMPIONAMENTO PIU' FITTO
       # A SECONDA DI QUANTO LA PATCH E' PIATTA, SALVO GLI INDICI SOLO DI UNA CORONA STRETTA SUL BORDO DELLA PATCH
        point_index = np.zeros((len(patch)))
        d2 = np.zeros(len(patch))
        center_index = []
        for k in range(len(patch)):
            point_index[k] = np.where((surf[:, 0] == patch[k, 0]) & (surf[:, 1] == patch[k, 1]) & (surf[:, 2] == patch[k, 2]))[0]
            d2[k] = (surf[int(point_index[k]),0] - surf[i,0])**2 + (surf[int(point_index[k]),1] - surf[i,1])**2 + (surf[int(point_index[k]),2] - surf[i,2])**2
            #print("d2[{}]=".format(k), d2[k])
###CAMBIA DA > A <
            if d2[k] < abs(prod_scal_patch)*R_crown: #SALVO I PUNTI NEL CENTRO
                center_index.append(point_index[k])


        # STO TENENDO COME PUNTI DI ZERNIKE ANCHE QUELLI NELLA CORONA DELLE PARTI PIATTE
        index_possible_area = [int(i) for i in index_possible_area if i not in center_index]


np.savetxt("{}\index_possible_area.txt".format(respath), index_possible_area)
