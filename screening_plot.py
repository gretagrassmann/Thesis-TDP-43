import numpy as np
import pandas as pd
import SurfaceFunc as SF
import random
    ##########  PARAMETERS  ##########

with open('configuration.txt') as f:
          for line in f:
              exec(line)

#screening = "index_continuous_distribution_exp"
#alpha = 11.0
#Rs_select = 4  # the radius of the sphere that includes the patch..
##R_zernike = 6
fragment = 208
cluster = 1
##step = 1
a = 1.
b = 1.
g = 10.0
d = 0.


    ############### END PARAMETERS  ##############
pdb_file = "..\\{}\cluster{}.dms".format(fragment, cluster)
respath = "..\\{}\cluster{}\R_zernike_{}\R_s_{}".format(fragment, cluster, R_zernike, Rs_select)

with open("{}\EDOTTIA\\a{}_b{}_g{}_d{}.txt".format(respath, a, b, g, d)) as f:
    index_possible_area = list(set([int(float(x)) for x in f.read().split()]))

#with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_{}.txt".format(respath, screening, Rs_select, alpha, step)) as f:
#    index_possible_area = [int(float(x)) for x in f.read().split()]

#TO VISUALIZE THE WHOLE MOLECULE
surf_total1 = pd.read_csv(pdb_file)
tot_len = len(surf_total1["x"])
r_t = [i+1 for i in index_possible_area]


TRE = 1
if TRE == 1:
    ## TOTAL
    surf_a_ = pd.read_csv(pdb_file)
    l_a = len(surf_a_["x"])
    print("Npoints", l_a)
    surf_a = np.zeros((l_a, 6))
    surf_a[:, :] = surf_a_[["x", "y", "z", "Nx", "Ny", "Nz"]]
    surf_a_obj = SF.Surface(surf_a[:, :], patch_num=0, r0=Rs_select, theta_max=45)
    res1, c = SF.ConcatenateFigPlots([surf_a_obj.surface[:, :3]])
    SF.Plot3DPoints(res1[:, 0], res1[:, 1], res1[:, 2], c, 0.3)

    ## SAMPLING
    r = [i + 1 for i in range(tot_len) if (i not in index_possible_area)]
    surf_select = pd.read_csv(pdb_file, skiprows=r)
    l = len(surf_select["x"])
    surf_a = np.zeros((l, 6))
    surf_a[:, :] = surf_select[["x", "y", "z", "Nx", "Ny", "Nz"]]
    surf_a_obj = SF.Surface(surf_a[:, :], patch_num=0, r0=Rs_select, theta_max=45)
    res1, c = SF.ConcatenateFigPlots([surf_a_obj.surface[:, :3]])
    SF.Plot3DPoints(res1[:, 0], res1[:, 1], res1[:, 2], c, 0.3)

    ### RANDOM
    random_index = sorted(random.sample(list(np.arange(0, l_a)), l))
    r_r = [i + 1 for i in range(l_a) if (i not in random_index)]
    surf_random = pd.read_csv(pdb_file, skiprows=r_r)
    l_r = len(surf_random["x"])
    surf_r = np.zeros((l_r, 6))
    surf_r[:, :] = surf_random[["x", "y", "z", "Nx", "Ny", "Nz"]]
    surf_r_obj = SF.Surface(surf_r[:, :], patch_num=0, r0=Rs_select, theta_max=45)
    res1, c = SF.ConcatenateFigPlots([surf_r_obj.surface[:, :3]])
    SF.Plot3DPoints(res1[:, 0], res1[:, 1], res1[:, 2], c, 0.3)







surf_total = pd.read_csv(pdb_file, skiprows=r_t)
len_tot_partial = len(surf_total["x"])

###PLOT TO VISUALIZE WHICH POINT WHERE SELECTED AS CENTER OF THE PATCHES
l_a = len(index_possible_area)
print("number of considered points=", l_a)
print("number of total points in the surface=", tot_len)
#r = [i for i in range(tot_len) if (i not in index_possible_area) & (i!=0)]
r = [i+1 for i in range(tot_len) if (i not in index_possible_area)]

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
