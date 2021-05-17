import os, sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as mpl
import matplotlib.pyplot as plt

import pandas as pd
from scipy.spatial import distance_matrix

import ZernikeFunc as ZF
import SurfaceFunc as SF
import my_functions

Rs_select = 1.  # the radius of the sphere that includes the patch..


pdb_file = "..\RRM2\cluster2.dms"
surf_name = pdb_file

with open("..\RRM2\cluster2\R_s_8\index\index_possible_area_R_s_8.0_alpha_0.078_step_1_points_2161.txt") as f:
    index_possible_area = [int(float(x)) for x in f.read().split()]

#TO VISUALIZE THE WHOLE MOLECULE
surf_total1 = pd.read_csv(surf_name)
tot_len = len(surf_total1["x"])
#r_t = [i for i in range(tot_len) if (i in index_possible_area) & (i != 0)]
r_t = [i+1 for i in range(tot_len) if (i in index_possible_area)] #la prima riga me la vede vuota


surf_total = pd.read_csv(surf_name, skiprows=r_t)
print(surf_total)
len_tot_partial = len(surf_total["x"])

###PLOT TO VISUALIZE WHICH POINT WHERE SELECTED AS CENTER OF THE PATCHES
l_a = len(index_possible_area)
print("number of considered points=", l_a)
print("number of total points in the surface=", tot_len)
r = [i for i in range(tot_len) if (i not in index_possible_area) & (i!=0)]

#colnames=["x", "y", "z", "Nx", "Ny", "Nz"]
#surf_select = pd.read_csv(surf_name, names=colnames, header=None, skiprows=r)
surf_select = pd.read_csv(surf_name, skiprows=r)
l = len(surf_select["x"])
surf_a = np.zeros((l, 6))
surf_a[:, :] = surf_select[["x", "y", "z", "Nx", "Ny", "Nz"]]

surf_tot = np.zeros((len_tot_partial, 6))
surf_tot[:, :] = surf_total[["x", "y", "z", "Nx", "Ny", "Nz"]]



### To inizialize the Surface class:
surf_a_obj = SF.Surface(surf_a[:, :], patch_num=0, r0=Rs_select, theta_max=45)
surf_tot_obj = SF.Surface(surf_tot[:, :], patch_num=0, r0=Rs_select, theta_max=45)
### To plot surface in 3D:
res1, c = SF.ConcatenateFigPlots([surf_tot_obj.surface[:, :3], surf_a_obj.surface[:, :3]])

SF.Plot3DPoints(res1[:, 0], res1[:, 1], res1[:, 2], c, 0.3)
