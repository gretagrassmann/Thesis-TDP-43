import ZernikeFunc as ZF

import os, sys
import numpy as np
import math
import statistics

def find_nearest_vector(array,value):
    idx = np.array([np.linalg.norm(x+y+z) for (x,y,z) in array-value]).argmin()
    return idx

def ZernikeCoefficients(index_possible_area, surf_obj, verso, Npixel, ZOrder, respath, kind, alpha):
    """
    Questa viene da Zernike.py.
    La uso per trovare i coefficienti di Zernike sia per tutti i punti (kind=total),
     sia una volta noti gli indici  (index_possible_area) dei punti che mi servono (kind=partial).

    """
    ltmp = len(index_possible_area)
    zernike_sampling_inv_a = np.zeros((121, ltmp))
    for i in range(ltmp):
        sys.stderr.write("\r Processing %i out of %i point" % (i, ltmp))
        sys.stderr.flush()
        patch, mask = surf_obj.BuildPatch(point_pos=index_possible_area[i], Dmin=.5)
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
        br_coeff = zernike_env.ZernikeDecomposition(order=ZOrder)

        zernike_sampling_inv_a[:, i] = np.absolute(br_coeff)

    res_inv_ = np.row_stack([index_possible_area, zernike_sampling_inv_a])
    if kind == "total":
        if verso ==1:
            np.savetxt("%s/zernike/zernike_positive/zernike_total.dat" % (respath), res_inv_,fmt="%.4e")
        else:
            np.savetxt("%s/zernike/zernike_negative/zernike_total.dat" % (respath), res_inv_,fmt="%.4e")

    else:
        if verso == 1:
            np.savetxt("{}/zernike/zernike_positive/zernike_alpha_{}_point_{}.dat".format(respath, alpha, ltmp), res_inv_,fmt="%.4e")
        else:
            np.savetxt("{}/zernike/zernike_negative/zernike_alpha_{}_point_{}.dat".format(respath, alpha, ltmp), res_inv_,fmt="%.4e")


    return(res_inv_)


def CosScanning(surf, surf_obj, Rs_select, alpha, step, respath):
    ltmp = np.shape(surf)[0]
    index_possible_area_total = np.arange(0, ltmp, step)
    # index_possible_area = np.arange(ltmp)
    index_possible = []

    mean_cos = []
    for i in index_possible_area_total:
        sys.stderr.write(("\r Processing point {} out of {}".format(i, ltmp)))
        sys.stderr.flush()
        patch, mask = surf_obj.BuildPatch(point_pos=i, Dmin=.5)

        if len(patch) <= 1:
            mean_cos = 1
        else:
            # PER CIASCUN PATCH TROVO LA MEDIA DEI COSENI
            normal_v = patch[:, 3:6]
            cos = []
            for ii in range(len(patch) - 1):
                norm_i = np.sqrt(np.sum(np.power(normal_v[ii], 2)))
                for j in range(ii + 1, len(patch)):
                    norm_j = np.sqrt(np.sum(np.power(normal_v[j], 2)))

                    dotproduct = 0
                    for a, b in zip(normal_v[ii], normal_v[j]):
                        dotproduct = dotproduct + a * b

                    cos.append(dotproduct / (norm_i * norm_j))

            mean_cos = (statistics.mean(cos) + 1) / 2.  # PER TENERE CONTO DEI COSENI NEGATIVI, SCALO TUTTO

        R_c = mean_cos * Rs_select * alpha
        # SE MEAN_COS->1+1 IL PATCH E' PIANO.
        # DEFINISCO CHE PUNTI SALVARE
        #point_index = np.zeros((len(patch)))
        d2 = np.zeros(len(patch))

        patch_center = \
        np.where((patch[:, 0] == surf[i, 0]) & (patch[:, 1] == surf[i, 1]) & (patch[:, 2] == surf[i, 2]))[0]
        if patch_center.size == 0:
            patch_center = find_nearest_vector(patch[:, 0:3], surf[i, 0:3])

        center_index = []
        for k in range(len(patch)):
            # point_index[k] = \
            # np.where((surf[:, 0] == patch[k, 0]) & (surf[:, 1] == patch[k, 1]) & (surf[:, 2] == patch[k, 2]))[0]
            # d2[k] = (surf[int(point_index[k]), 0] - surf[i, 0]) ** 2 + (surf[int(point_index[k]), 1] - surf[i, 1]) ** 2 + (
            #            surf[int(point_index[k]), 2] - surf[i, 2]) ** 2
            d2[k] = (patch[k, 0] - patch[patch_center, 0]) ** 2 + (patch[k, 1] - patch[patch_center, 1]) ** 2 + (
                    patch[k, 2] - patch[patch_center, 2]) ** 2
            if d2[k] > R_c ** 2:  # SALVO I PUNTI NEL CENTRO
                # center_index.append(point_index[k])
                indx = \
                np.where((surf[:, 0] == patch[k, 0]) & (surf[:, 1] == patch[k, 1]) & (surf[:, 2] == patch[k, 2]))[0]
                if indx not in index_possible:
                    index_possible.append(indx[0])

        # index_possible_area = [int(i) for i in index_possible_area if i not in center_index]
    index_possible_area = sorted(index_possible)

    np.savetxt("{}\index\index_possible_area_R_s_{}_alpha_{}_step_{}_points_{}.txt".format(respath,Rs_select, alpha, step, len(index_possible_area)), index_possible_area)
    print("\r number of points for R_s={},alpha={},step={} =".format(Rs_select,alpha,step), len(index_possible_area))

    return(index_possible_area)


def PCAgrid(x):
    X = x.transpose()

    cov_mat = np.cov(X, rowvar=False)
    # Calculating Eigenvalues and Eigenvectors of the covariance matrix
    eigen_values, eigen_vectors = np.linalg.eigh(cov_mat)
    # sort the eigenvalues in descending order
    sorted_index = np.argsort(eigen_values)[::-1]
    sorted_eigenvalue = eigen_values[sorted_index]
    # similarly sort the eigenvectors
    sorted_eigenvectors = eigen_vectors[:, sorted_index]

    # select the first n eigenvectors, n is desired dimension
    # of our final reduced data.
    n_components = 2  # you can select any number of components.
    eigenvector_subset = sorted_eigenvectors[:, 0:n_components]

    # Transform the data
    # X_reduced = np.dot(eigenvector_subset.transpose(), X_meaned.transpose()).transpose()
    X_reduced = np.dot(eigenvector_subset.transpose(), X.transpose()).transpose()
    x_ax = X_reduced[:, 0]
    y_ax = X_reduced[:, 1]

    number_division = int(math.log10(len(x[0, :]))) * 5

    x_grid = np.linspace(min(x_ax), max(x_ax), num=number_division)
    y_grid = np.linspace(min(y_ax), max(y_ax), num=number_division)
    cell_w = x_grid[1] - x_grid[0]
    cell_h = y_grid[1] - y_grid[0]

    count = np.zeros((len(x_grid), len(y_grid)))
    for i in range(len(x_grid)):
        for j in range(len(y_grid)):
            index_cell = np.where(
                (X_reduced[:, 0] >= x_grid[0] + i * cell_w) & (X_reduced[:, 0] < x_grid[0] + (i + 1) * cell_w) & \
                (X_reduced[:, 1] >= y_grid[0] + j * cell_h) & (X_reduced[:, 1] < y_grid[0] + (j + 1) * cell_h))
            count_cell = len(index_cell[0])
            count[i, j] = count_cell

    return(count)
