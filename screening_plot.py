import os, sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as mpl

import pandas as pd
from scipy.spatial import distance_matrix

import ZernikeFunc as ZF
import SurfaceFunc as SF

Rs = 2.  # the radius of the sphere that includes the patch..


pdb_file = "..\RRM2\cluster2.dms"
surf_name = pdb_file

#TO VISUALIZE THE WHOLE MOLECULE
surf_total = pd.read_csv(surf_name)
tot_len = len(surf_total["x"])
surf_tot = np.zeros((tot_len, 6))
surf_tot[:, :] = surf_total[["x", "y", "z", "Nx", "Ny", "Nz"]]

with open("..\RRM2\cluster2\positive_new_sampling\index_possible_area.txt") as f:
    index_possible_area = [int(float(x)) for x in f.read().split()]

###PLOT TO VISUALIZE WHICH POINT WHERE SELECTED AS CENTER OF THE PATCHES
l_a = len(index_possible_area)-1
print("l_a=", l_a)
r = [i for i in range(total_len) if i not in index_possible_area]
surf_select = pd.read_csv(surf_name, skiprows=r)
surf_a = np.zeros((l_a, 6))
surf_a[:, :] = surf_select[["x", "y", "z", "Nx", "Ny", "Nz"]]




### To inizialize the Surface class:
surf_a_obj = SF.Surface(surf_a[:, :], patch_num=0, r0=Rs, theta_max=45)
surf_tot_obj = SF.Surface(surf_tot[:, :], patch_num=0, r0=Rs, theta_max=45)
### To plot surface in 3D:
res1, c = SF.ConcatenateFigPlots([surf_tot_obj.surface[:, :3], surf_a_obj.surface[:, :3]])

SF.Plot3DPoints(res1[:, 0], res1[:, 1], res1[:, 2], c, 0.3)
