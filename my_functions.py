import ZernikeFunc as ZF
import random
import os, sys
import numpy as np
import math
import statistics
import operator
import random
import pandas as pd
import SurfaceFunc as SF


def find_nearest_vector2D(array,value):
    dist_2 = np.sum((array - value)**2, axis=1)
    return np.argmin(dist_2)


#def find_nearest_vector(array,value):
#    idx = np.array([np.linalg.norm(x+y+z) for (x,y,z) in array-value]).argmin()
#    return idx

def find_nearest_vector(array,value):
    array = np.asarray(array)
    dist_3 = np.sum((array - value)**2, axis=1)
    return np.argmin(dist_3)





################ POINTS SCREENING ########################


def CosWithoutScaling(surf, surf_obj_scan, respath, Rs_select, step ):
    ltmp = np.shape(surf)[0]
    index_possible_area_total = np.arange(0, ltmp, step)

    mean_cos_surface = []
    for i in index_possible_area_total:
        sys.stderr.write(("\r Processing point {} out of {} for the cosine".format(i, ltmp)))
        sys.stderr.flush()
        patch, mask = surf_obj_scan.BuildPatch(point_pos=i, Dmin=.5)
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

            if statistics.mean(cos) >= 0 :
                mean_cos = statistics.mean(cos)
            else:
                mean_cos = 0.
        mean_cos_surface.append(mean_cos)
    np.savetxt("{}\COSINE_step_{}_Rs_{}.txt".format(respath, step, Rs_select), mean_cos_surface)
    print("Number of considered patches={}".format(len(mean_cos_surface)))
    return()


def ContinuousDistribution(mean_cos, surf, surf_obj_scan, Rs_select, alpha, step, respath, screening):
    ltmp = np.shape(surf)[0]
    index_possible_area_total = np.arange(0, ltmp, step)

    index_possible = []
    delete_index = []
    for i in index_possible_area_total:
        if i not in np.array(delete_index):
            sys.stderr.write(("\r Processing point {} out of {} for the screening with alpha={}".format(i, ltmp, alpha)))
            sys.stderr.flush()
            patch, mask = surf_obj_scan.BuildPatch(point_pos=i, Dmin=.5)

            # Ogni volta prendo il centro: qua trovo il centro
            patch_center_indx_inpatch = np.where((patch[:, 0] == surf[i, 0]) & (patch[:, 1] == surf[i, 1]) & (patch[:, 2] == surf[i, 2]))[0]
            # Se quando ho costruito la patch ho cancellato il punto di partenza per pulirla, prendo come centro il punto piu' vicino
            if patch_center_indx_inpatch.size == 0:
                patch_center_indx_inpatch = find_nearest_vector(patch[:, 0:3], surf[i, 0:3])
                patch_center = np.where((surf[:, 0] == patch[patch_center_indx_inpatch,0]) & (surf[:, 1] == patch[patch_center_indx_inpatch,1]) & (surf[:, 2] == patch[patch_center_indx_inpatch,2]))[0]
                index_possible.append(patch_center)
                delete_index.append(patch_center)
            else:
                index_possible.append(i)
                delete_index.append(i)


            # Ora calcolo per ciascun punto la distanza dal centro e la normalizzo. Poi estraggo un numero casuale: se e' minore
            # della prob considero quel punto

            d2 = np.sqrt((patch[:, 0] - patch[patch_center_indx_inpatch, 0]) ** 2 + (patch[:, 1] - patch[patch_center_indx_inpatch, 1]) ** 2 + (
                        patch[:, 2] - patch[patch_center_indx_inpatch, 2]) ** 2)/Rs_select
            prob = (d2**alpha)*(1-mean_cos[i])
            rand = np.random.random_sample(size=len(patch))
            mask2 = rand < prob

            mask_reject = rand >= prob
            patch_reject = patch[mask_reject,:]
            for k in range(len(patch_reject)):
                index_in_surface = \
                    np.where((surf[:, 0] == patch_reject[k, 0]) & (surf[:, 1] == patch_reject[k, 1]) & (
                                surf[:, 2] == patch_reject[k, 2]))[0][0]
                delete_index.append(index_in_surface)


            patch_selected = patch[mask2, :]
            for k in range(len(patch_selected)):
                index_in_surface = \
                np.where((surf[:, 0] == patch_selected[k, 0]) & (surf[:, 1] == patch_selected[k, 1]) & (surf[:, 2] == patch_selected[k, 2]))[0][0]
                if index_in_surface not in np.array(index_possible):
                    index_possible.append(index_in_surface)
                delete_index.append(index_in_surface)

            #if (alpha == -4.) & (p % 10 == 0):
            #    fragment = 208
            #    cluster = 1
            #    r_t = [i + 1 for i in list(index_possible)]
            #    pdb_file = "..\\{}\cluster{}.dms".format(fragment, cluster)
            #    surf_total = pd.read_csv(pdb_file, skiprows=r_t)
            #    surf_tot = np.zeros((len(surf_total['x']), 6))
            #    surf_tot[:, :] = surf_total[["x", "y", "z", "Nx", "Ny", "Nz"]]
            #    r = [i + 1 for i in range(ltmp) if (i not in list(index_possible))]
            #    surf_select = pd.read_csv(pdb_file, skiprows=r)
            #    surf_a = np.zeros((len(surf_select['x']), 6))
            #    surf_a[:, :] = surf_select[["x", "y", "z", "Nx", "Ny", "Nz"]]
                ### To inizialize the Surface class:
            #    surf_a_obj = SF.Surface(surf_a[:, :], patch_num=0, r0=Rs_select, theta_max=45)
            #    surf_tot_obj = SF.Surface(surf_tot[:, :], patch_num=0, r0=Rs_select, theta_max=45)
                ### To plot surface in 3D:
            #    res1, c = SF.ConcatenateFigPlots([surf_tot_obj.surface[:, :3], surf_a_obj.surface[:, :3]])
            #    SF.Plot3DPoints(res1[:, 0], res1[:, 1], res1[:, 2], c, 0.3)

            #for k in range(len(patch)):
            #    distance = (math.sqrt((patch[k, 0] - patch[patch_center_indx_inpatch, 0]) ** 2 + (patch[k, 1] - patch[patch_center_indx_inpatch, 1]) ** 2 + (
            #            patch[k, 2] - patch[patch_center_indx_inpatch, 2]) ** 2))/Rs_select
            #    if screening == "index_continuous_distribution":
            #        prob = (alpha*(distance**2)+(1-alpha)*distance)*(1-mean_cos[i])
            #    else:
            #        prob = (distance**alpha)*(1-mean_cos[i])
            #    rand = random.random()
            #    index_in_surface = \
            #    np.where((surf[:, 0] == patch[k, 0]) & (surf[:, 1] == patch[k, 1]) & (surf[:, 2] == patch[k, 2]))[0]
            #    if rand < prob:
            #        if index_in_surface not in np.array(index_possible):
            #            index_possible.append(index_in_surface)
            #    delete_index.append(index_in_surface)

    index_possible_area = sorted(index_possible)

    np.savetxt("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_{}.txt".format(respath, screening, Rs_select, alpha, step),
               index_possible_area)
    print("\r number of points for R_s={},alpha={},step={} =".format(Rs_select, alpha, step), len(index_possible_area))

    return ()






def CosDistributionScreening(mean_cos, surf, surf_obj_scan, Rs_select, alpha, step, respath, number_crowns):
    ltmp = np.shape(surf)[0]
    index_possible_area_total = np.arange(0, ltmp, step)

    crown_limits = np.linspace(.0, Rs_select, num=number_crowns+1)

    index_possible = []
    delete_index = []
    for i in index_possible_area_total:
        if i not in np.array(delete_index):
            sys.stderr.write(("\r Processing point {} out of {} for the screening with alpha=".format(i, ltmp, alpha)))
            sys.stderr.flush()
            patch, mask = surf_obj_scan.BuildPatch(point_pos=i, Dmin=.5)

            # Ogni volta per il disco interno prendo solo in centro: qua trovo il centro
            patch_center_indx_inpatch = np.where((patch[:, 0] == surf[i, 0]) & (patch[:, 1] == surf[i, 1]) & (patch[:, 2] == surf[i, 2]))[0]
            patch_center = np.where((surf[:, 0] == patch[patch_center_indx_inpatch, 0]) & (surf[:, 1] == patch[patch_center_indx_inpatch, 1]) & (surf[:, 2] == patch[patch_center_indx_inpatch, 2]))[0]
            # Se quando ho costruito la patch ho cancellato il punto di partenza per pulirla, prendo come centro il punto piu' vicino
            if patch_center.size == 0:
                patch_center_indx_inpatch = find_nearest_vector(patch[:, 0:3], surf[i, 0:3])
                patch_center = np.where((surf[:, 0] == patch[patch_center_indx_inpatch,0]) & (surf[:, 1] == patch[patch_center_indx_inpatch,1]) & (surf[:, 2] == patch[patch_center_indx_inpatch,2]))[0]

            index_possible.append(patch_center)

            # Ora invece riempio per ciascuna corona le liste con i punti che le appartengono
            distance = []
            for k in range(len(patch)):
                distance.append(math.sqrt((patch[k, 0] - patch[patch_center_indx_inpatch, 0]) ** 2 + (patch[k, 1] - patch[patch_center_indx_inpatch, 1]) ** 2 + (
                        patch[k, 2] - patch[patch_center_indx_inpatch, 2]) ** 2))

            sorted_distance = sorted(enumerate(distance), key=operator.itemgetter(1))

            first_point = 0
            for ii in range(1, number_crowns):
                points_in_crown = []
                for i in range(first_point, len(patch)):
                    if (sorted_distance[i][1] > crown_limits[ii]) & (sorted_distance[i][1] <= crown_limits[ii+1]):
                        points_in_crown.append(sorted_distance[i][0])
                    if sorted_distance[i][1] > crown_limits[ii+1]:
                        first_point = i
                        break

                tot_points = len(points_in_crown)
                n_points = int(tot_points*(1-mean_cos[i])*(crown_limits[ii+1]/Rs_select)*alpha)
                patch_points = random.sample(points_in_crown, n_points)
                for n in range(n_points):
                    index_possible.append(np.where((surf[:, 0] == patch[patch_points[n], 0]) & (surf[:, 1] == patch[patch_points[n], 1]) & (surf[:, 2] == patch[patch_points[n], 2]))[0])



            #Ora faccio in modo di non andare a considerare piu' nessuno dei punti di questa patch
            #delete_index_patch = np.zeros(len(patch))
            for k in range(len(patch)):
                delete = np.where((surf[:, 0] == patch[k,0]) & (surf[:, 1] == patch[k,1]) & (surf[:, 2] == patch[k,2]))[0]
                delete_index.append(delete)
              #  delete_index_patch[k] = (np.where((surf[:, 0] == patch[k,0]) & (surf[:, 1] == patch[k,1]) & (surf[:, 2] == patch[k,2]))[0])

            #delete_index.append(delete_index_patch)


    index_possible_area = sorted(index_possible)
    np.savetxt("{}\index_distribution_{}_crowns\index_possible_area_R_s_{}_alpha_{}_step_{}.txt".format(respath, number_crowns, Rs_select, alpha, step),
               index_possible_area)
    print("\r number of points for R_s={},alpha={},step={} =".format(Rs_select, alpha, step), len(index_possible_area))

    return ()



def PercentageScreening(mean_cos, surf, surf_obj_scan, Rs_select, alpha, step, respath):
    ltmp = np.shape(surf)[0]
    index_possible_area_total = np.arange(0, ltmp, step)
    index_possible = []
    delete_index = []
    iter = 0
    for i in index_possible_area_total:
        if i not in np.array(delete_index):
            iter += 1
            sys.stderr.write(("\r Processing point {} out of {} for the screening with alpha=".format(i, ltmp, alpha)))
            sys.stderr.flush()
            patch, mask = surf_obj_scan.BuildPatch(point_pos=i, Dmin=.5)
            number_points = int((1.-mean_cos[i])*len(patch)*alpha)
            # Se e' completamente piatta (cos=1) prendo solo il punto al centro
            if number_points == 0:
                patch_center_indx_inpatch = np.where((patch[:, 0] == surf[i, 0]) & (patch[:, 1] == surf[i, 1]) & (patch[:, 2] == surf[i, 2]))[0]
                patch_center = np.where((surf[:, 0] == patch[patch_center_indx_inpatch, 0]) & (surf[:, 1] == patch[patch_center_indx_inpatch, 1]) & (surf[:, 2] == patch[patch_center_indx_inpatch, 2]))[0]
                # Se quando ho costruito la patch ho cancellato il punto di partenza per pulirla, prendo come centro il punto piu' vicino
                if patch_center.size == 0:
                    patch_center_indx_inpatch = find_nearest_vector(patch[:, 0:3], surf[i, 0:3])
                    patch_center = np.where((surf[:, 0] == patch[patch_center_indx_inpatch,0]) & (surf[:, 1] == patch[patch_center_indx_inpatch,1]) & (surf[:, 2] == patch[patch_center_indx_inpatch,2]))[0]

                if patch_center not in index_possible:
                    index_possible.append(patch_center)

            else:
                patch_points = random.sample(list(patch), number_points)

                for k in range(len(patch_points)):
                    indx = np.where((surf[:, 0] == patch_points[k][0]) & (surf[:, 1] == patch_points[k][1]) & (surf[:, 2] == patch_points[k][2]))[0]
                    if indx not in index_possible:
                        index_possible.append(indx[0])


            #Ora faccio in modo di non andare a considerare piu' nessuno dei punti di questa patch
            #delete_index_patch = np.zeros(len(patch))
            for k in range(len(patch)):
                delete = np.where((surf[:, 0] == patch[k,0]) & (surf[:, 1] == patch[k,1]) & (surf[:, 2] == patch[k,2]))[0]
                if delete not in delete_index:
                    delete_index.append(delete)
              #  delete_index_patch[k] = (np.where((surf[:, 0] == patch[k,0]) & (surf[:, 1] == patch[k,1]) & (surf[:, 2] == patch[k,2]))[0])

            #delete_index.append(delete_index_patch)


    index_possible_area = sorted(index_possible)
    np.savetxt("{}\index_percentage\index_possible_area_R_s_{}_alpha_{}_step_{}.txt".format(respath, Rs_select, alpha, step),
               index_possible_area)
    print("\r number of points for R_s={},alpha={},step={} =".format(Rs_select, alpha, step), len(index_possible_area))

    return ()


def NewCosScanning(mean_cos, surf, surf_obj_scan, Rs_select, alpha, step, respath):
    ltmp = np.shape(surf)[0]
    index_possible_area_total = np.arange(0, ltmp, step)
    index_possible = []

    for i in index_possible_area_total:
        sys.stderr.write(("\r Processing point {} out of {} for the screening with alpha=".format(i, ltmp, alpha)))
        sys.stderr.flush()
        patch, mask = surf_obj_scan.BuildPatch(point_pos=i, Dmin=.5)

        R_c = mean_cos[i] * Rs_select * alpha
        # SE MEAN_COS->1 IL PATCH E' PIANO.
        # DEFINISCO CHE PUNTI SALVARE
        d2 = np.zeros(len(patch))

        patch_center = \
        np.where((patch[:, 0] == surf[i, 0]) & (patch[:, 1] == surf[i, 1]) & (patch[:, 2] == surf[i, 2]))[0]
        if patch_center.size == 0:
            patch_center = find_nearest_vector(patch[:, 0:3], surf[i, 0:3])

        center_index = []
        for k in range(len(patch)):
            d2[k] = (patch[k, 0] - patch[patch_center, 0]) ** 2 + (patch[k, 1] - patch[patch_center, 1]) ** 2 + (
                    patch[k, 2] - patch[patch_center, 2]) ** 2
            if d2[k] > R_c ** 2:  # SALVO I PUNTI NON NEL CENTRO
                indx = \
                np.where((surf[:, 0] == patch[k, 0]) & (surf[:, 1] == patch[k, 1]) & (surf[:, 2] == patch[k, 2]))[0]
                if indx not in index_possible:
                    index_possible.append(indx[0])

    index_possible_area = sorted(index_possible)

    np.savetxt("{}\index\index_possible_area_R_s_{}_alpha_{}_step_{}.txt".format(respath,Rs_select, alpha, step), index_possible_area)
    print("\r number of points for R_s={},alpha={},step={} =".format(Rs_select,alpha,step), len(index_possible_area))

    return()

        ##################### LOSS FUNCTION MINIMIZATION    #########################

def cart2polar(pca, centroid):

    x_centered = pca[:,0] - centroid[0, 0]
    y_centered = pca[:,1] - centroid[0, 1]

    radius= []
    theta = []
    for i in range(len(x_centered)):
        radius.append(math.sqrt(x_centered[i] * x_centered[i] + y_centered[i] * y_centered[i]))
        #    theta.append(math.atan(y_tot_centered[i]/x_tot_centered[i])*180/math.pi)
        theta.append(math.atan2(x_centered[i], y_centered[i]) * 180 / math.pi)

    return(theta, radius)


def PCA(x, n_components):
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
  # you can select any number of components.
    eigenvector_subset = sorted_eigenvectors[:, 0:n_components]
    eigen_values = sorted_eigenvalue[0:n_components]

    # Transform the data
    # X_reduced = np.dot(eigenvector_subset.transpose(), X_meaned.transpose()).transpose()
    X_reduced = np.dot(eigenvector_subset.transpose(), X.transpose()).transpose()

    return(X_reduced, eigenvector_subset, eigen_values, sorted_eigenvalue)






##########################  INUTILI ###################################################



def Cos(surf, surf_obj_scan, respath, Rs_select, step ):
    ltmp = np.shape(surf)[0]
    index_possible_area_total = np.arange(0, ltmp, step)

    mean_cos_surface = []
    for i in index_possible_area_total:
        sys.stderr.write(("\r Processing point {} out of {} for the cosine".format(i, ltmp)))
        sys.stderr.flush()
        patch, mask = surf_obj_scan.BuildPatch(point_pos=i, Dmin=.5)
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
        mean_cos_surface.append(mean_cos)
    np.savetxt("{}\COSINE_step_{}_Rs_{}.txt".format(respath, step, Rs_select), mean_cos_surface)
    print("Number of considered patches={}".format(len(mean_cos_surface)))
    return()



def CosScanning(surf, surf_obj_scan, Rs_select, alpha, step, respath):
    ltmp = np.shape(surf)[0]
    index_possible_area_total = np.arange(0, ltmp, step)
    # index_possible_area = np.arange(ltmp)
    index_possible = []

    if alpha == 1.:
        mean_cos_surface = []
    for i in index_possible_area_total:
        sys.stderr.write(("\r Processing point {} out of {}".format(i, ltmp)))
        sys.stderr.flush()
        patch, mask = surf_obj_scan.BuildPatch(point_pos=i, Dmin=.5)
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
        if alpha == 1.:
            mean_cos_surface.append(mean_cos)

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

    np.savetxt("{}\index\index_possible_area_R_s_{}_alpha_{}_step_{}.txt".format(respath,Rs_select, alpha, step), index_possible_area)
    print("\r number of points for R_s={},alpha={},step={} =".format(Rs_select,alpha,step), len(index_possible_area))
    if alpha == 1:
        np.savetxt("{}\cos\COS.txt".format(respath), mean_cos_surface)

    return()


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



def ZernikeCoefficients(index_possible_area, surf_obj, verso, Npixel, ZOrder, respath, kind, alpha):
    """
        NON FUNZIONA DA QUA IN CERTI CASI, LA METTO DIRETTAMENTE IN COMPLETE.PY
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
            np.savetxt("{}/zernike/zernike_positive/zernike_alpha_{}.dat".format(respath, alpha), res_inv_,fmt="%.4e")
        else:
            np.savetxt("{}/zernike/zernike_negative/zernike_alpha_{}.dat".format(respath, alpha), res_inv_,fmt="%.4e")


    return()

