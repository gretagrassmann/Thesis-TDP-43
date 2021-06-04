import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import ZernikeFunc as ZF
import SurfaceFunc as SF
import my_functions
import shutil
        ######################## PARAMETERS #########
with open('configuration.txt') as f:
    for line in f:
        exec(line)
#alpha = [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1.]

#alpha = [ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.]
#alpha = [21., 22., 23., 24., 25., 26., 27., 28., 29., 30.]
#alpha = [40., 50., 60., 70., 80., 90., 100.]
#alpha = [71., 72., 73., 74., 75., 76., 77., 78., 79.,
#         81., 82., 83., 84., 85., 86., 87., 88., 89.,
#         91., 92., 93., 94., 95., 96.]
alpha = [103]


screening = "index_continuous_distribution_exp"

#number_crowns = 20
#R_zernike = 6
Rs_select = 4
#fragment = 208
cluster = 1
#step = 1
#verso = float(1)
#Npoint = int(1)

    ###################################   END PARAMETERS    ###########

pdb_file = "..\\{}\cluster{}.dms".format(fragment, cluster)
respath = "..\\{}\cluster{}\R_zernike_{}\R_s_{}".format(fragment, cluster, R_zernike, Rs_select)
if (verso == 1):
    ver = "up"
elif (verso == -1):
    ver = "down"
ZOrder = 20  # the Zernike expansion order..
Npixel = 25  # the plane edge in pixels...
Daa = 3.  # the maximum distance between antigen and antibody binding site atoms..
Dpp = 1.  # the distance between points of the same patch (needed to remove islands)

surf_ = pd.read_csv(pdb_file)
lag = len(surf_["x"])
surf = np.zeros((lag, 6))
surf[:, :] = surf_[["x", "y", "z", "Nx", "Ny", "Nz"]]
surf_obj = SF.Surface(surf[:, :], patch_num=0, r0=R_zernike, theta_max=45)
surf_obj_scan = SF.Surface(surf[:, :], patch_num=0, r0=Rs_select, theta_max=45)
ltmp = np.shape(surf)[0]



            ################## ZERNIKE PER OGNI SINGOLO PUNTO   #####################
print("Vuoi fare ZERNIKE per OGNI {} PUNTI? y o n?".format(Npoint))
o = input()
if o == "y":
    index_all_possible_area = np.arange(ltmp)[::Npoint]
    ltmp = len(index_all_possible_area)
    zernike_sampling_inv_a = np.zeros((121, ltmp))

    all_patches = []
    for i in range(ltmp):
        sys.stderr.write("\r Processing %i out of %i point" % (i, ltmp))
        sys.stderr.flush()
        # finding antigen patch, plane and zernike descriptors..

        patch, mask = surf_obj.BuildPatch(point_pos=index_all_possible_area[i], Dmin=.5)
        surf_obj.real_br = mask

        rot_patch, rot_ag_patch_nv = surf_obj.PatchReorientNew(patch, verso)

        z = surf_obj.FindOrigin(rot_patch)
        plane, weigths, dist_plane, thetas = surf_obj.CreatePlane(patch=rot_patch, z_c=z, Np=Npixel)
        new_plane = plane.copy()
        new_plane___ = plane.copy()

        if (np.shape(rot_patch)[1] == 4):

            new_plane_re = surf_obj.FillTheGap_everywhere(plane_=np.real(plane))
            new_plane_im = surf_obj.FillTheGap_everywhere(plane_=np.imag(plane))

            new_plane_re_ = surf_obj.EnlargePixels(new_plane_re)
            new_plane_im_ = surf_obj.EnlargePixels(new_plane_im)

            new_plane_ = new_plane_re_ + 1j * new_plane_im_ / np.max(np.abs(new_plane_im_))
        else:
            new_plane = surf_obj.FillTheGap_everywhere(plane_=plane)
            ## enlarging plane..
            new_plane_ = surf_obj.EnlargePixels(new_plane)

        try:
            zernike_env.img = new_plane_
        except:
            zernike_env = ZF.Zernike2d(new_plane_)
        # br_recon, br_coeff = zernike_env.ZernikeReconstruction(order=ZOrder, PLOT=0)
        br_coeff = zernike_env.ZernikeDecomposition(order=ZOrder)

        zernike_sampling_inv_a[:, i] = np.absolute(br_coeff)


    res_inv_ = np.row_stack([index_all_possible_area, zernike_sampling_inv_a])
    if verso == 1:
        np.savetxt("{}/zernike/zernike_positive/zernike_total.dat".format(respath), res_inv_, fmt="%.4e")

    else:
        np.savetxt("{}/zernike/zernike_negative/zernike_total.dat".format(respath), res_inv_, fmt="%.4e")

    #############################   ROUGHNESS DELLE PATCH  ###########
print("Vuoi calcolare la MEDIA DEL COSENO di ciascuna possibile patch? y o n?")
ooo = input()
if ooo == "y":
    cosine = my_functions.CosWithoutScaling(surf, surf_obj_scan, respath, Rs_select, step)
   # with open("{}\COSINE_step_{}_Rs_{}.txt".format(respath, step, Rs_select)) as f:
   #     mean_cos = [2 * (float(x)) - 1 for x in f.read().split()]
   # plt.hist(mean_cos, bins='auto')  # arguments are passed to np.histogram
   # plt.title("Mean cosine value of the surface's {} patches".format(len(mean_cos)))
   # plt.show()

print("Vuoi vedere la media dei coseni per ciascuna patch? y o n?")
o = input()
if o == "y":
    with open("{}\COSINE_step_{}_Rs_{}.txt".format(respath, step, Rs_select)) as f:
#        mean_cos = [2 * (float(x)) - 1 for x in f.read().split()]
        mean_cos = [float(x) for x in f.read().split()]

    plt.hist(mean_cos, bins='auto')  # arguments are passed to np.histogram
    plt.title("Mean cosine value of the surface's {} patches".format(len(mean_cos)))
    plt.show()

    ############### SCREENING PARTENDO DA COSINE.txt    #################
print("Vuoi fare lo SCREENING partendo da cosine.txt per diversi valori di alpha? y o n?")
o = input()
if o == "y":
    with open("{}\COSINE_step_{}_Rs_{}.txt".format(respath, step, Rs_select)) as f:
        #mean_cos = [2 * (float(x)) - 1 for x in f.read().split()]
        mean_cos = [float(x) for x in f.read().split()]
    for i in alpha:
#        index_possible_area = my_functions.NewCosScanning(mean_cos, surf, surf_obj_scan, Rs_select, i, step, respath)
#        index_possible_area = my_functions.PercentageScreening(mean_cos, surf, surf_obj_scan, Rs_select, i, step, respath)
#        index_possible_area = my_functions.CosDistributionScreening(mean_cos, surf, surf_obj_scan, Rs_select, i, step, respath, number_crowns)
        index_possible_area = my_functions.ContinuousDistribution(mean_cos, surf, surf_obj_scan, Rs_select, i, step, respath, screening)



print("Vuoi salvare i coefficenti di Zernike dello screening migliore? Se si con che alpha? Altrimenti digita n")
alpha = input()
if alpha != "n":
    # Apro il file con gli indici delle patch di cui voglio salvare Zernike
    shutil.copy("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, alpha),
                "{}\zernike\index_alpha.txt".format(respath))


    with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, alpha)) as f:
        index_possible_points = [int(float(x)) for x in f.read().split()]

    #Apro il file con tutti gli indici di Zernike, compresa la riga 1 che contiene gli indici
    if verso == 1:
        zernike_total = np.loadtxt("{}\zernike\zernike_positive\zernike_total.dat".format(respath), delimiter=" ")
    else:
        zernike_total = np.loadtxt("{}\zernike\zernike_negative\zernike_total.dat".format(respath), delimiter=" ")



    zernike = zernike_total[:,index_possible_points]



    if verso == 1:
        np.savetxt("{}/zernike/zernike_positive/zernike_alpha.dat".format(respath, alpha), zernike, fmt="%.4e")
    else:
        np.savetxt("{}/zernike/zernike_negative/zernike_alpha.dat".format(respath, alpha), zernike, fmt="%.4e")

