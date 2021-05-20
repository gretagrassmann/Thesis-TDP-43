import numpy as np
import pandas as pd
import SurfaceFunc as SF

    ##########  PARAMETERS  ##########

with open('configuration.txt') as f:
          for line in f:
              exec(line)
alpha = 0.075

#Rs_select = 6  # the radius of the sphere that includes the patch..
#R_zernike = 6
#alpha = 0.075
#fragment = 208
#luster = 1
#tep = 1



    ############### END PARAMETERS  ##############
pdb_file = "..\\{}\cluster{}.dms".format(fragment, cluster)
respath = "..\\{}\cluster{}\R_zernike_{}\R_s_{}".format(fragment, cluster, R_zernike, Rs_select)

with open("{}\index\index_possible_area_R_s_{}_alpha_{}_step_{}.txt".format(respath,Rs_select, alpha, step)) as f:
    index_possible_area = [int(float(x)) for x in f.read().split()]

#TO VISUALIZE THE WHOLE MOLECULE
surf_total1 = pd.read_csv(pdb_file)
tot_len = len(surf_total1["x"])
#r_t = [i for i in range(tot_len) if (i in index_possible_area) & (i != 0)]
r_t = [i+1 for i in range(tot_len) if (i in index_possible_area)] #la prima riga me la vede vuota


surf_total = pd.read_csv(pdb_file, skiprows=r_t)
len_tot_partial = len(surf_total["x"])

###PLOT TO VISUALIZE WHICH POINT WHERE SELECTED AS CENTER OF THE PATCHES
l_a = len(index_possible_area)
print("number of considered points=", l_a)
print("number of total points in the surface=", tot_len)
r = [i for i in range(tot_len) if (i not in index_possible_area) & (i!=0)]

#colnames=["x", "y", "z", "Nx", "Ny", "Nz"]
#surf_select = pd.read_csv(surf_name, names=colnames, header=None, skiprows=r)
surf_select = pd.read_csv(pdb_file, skiprows=r)
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
