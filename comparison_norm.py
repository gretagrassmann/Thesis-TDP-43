import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import shutil
import my_functions
import sys
from scipy.signal import savgol_filter
import random
import statistics
from scipy.spatial.distance import cdist
from kneed import KneeLocator

screening = "index_continuous_distribution_exp"

    ################### PARAMETERS  ########
with open('configuration.txt') as f:
    for line in f:
        exec(line)
fragment = 220
cluster = 2
R_zernike = 6
Rs_select = 4
ver = 1
alpha = [ -3., -2., -1.,
           1., 2., 3., 4., 5., 6., 7., 8., 9.
         , 10., 11., 12., 13., 14., 15., 16., 17., 18., 19.
         , 20., 21., 22., 23., 24., 25., 26., 27., 28., 29.
         , 30., 31., 32., 33., 34., 35., 36., 37., 38., 39.
         , 40., 41., 42., 43., 44., 45., 46., 47., 48., 49.
         , 50., 51., 52., 53., 54., 55., 56., 57., 58., 59.
         , 60., 61., 62., 63., 64., 65., 66., 67., 68., 69.
         , 70., 71., 72., 73., 74., 75., 76., 77., 78., 79.
         , 80., 81., 82., 83., 84., 85., 86., 87., 88., 89.
         , 90., 91., 92., 93., 94., 95., 96.
         ]
    ############### END PARAMETERS  ################
respath = "..\\{}\cluster{}\R_zernike_{}\R_s_{}".format(fragment, cluster, R_zernike, Rs_select)
pdb_file = "..\\{}\cluster{}.dms".format(fragment, cluster)

if ver == 1:
    total = np.loadtxt("{}\zernike\zernike_positive\zernike_total.dat".format(respath), delimiter=" ", skiprows=1)
else:
    total = np.loadtxt("{}\zernike\zernike_negative\zernike_total.dat".format(respath), delimiter=" ", skiprows=1)

surf_ = pd.read_csv(pdb_file)
lag = len(surf_["x"])
surf = np.zeros((lag, 6))
surf[:, :] = surf_[["x", "y", "z", "Nx", "Ny", "Nz"]]

   ####################################################3
print('Vuoi vedere la media dei punti selezionati nel caso random e per il sampling? y o n?')
o = input()
if o == 'y':
    print("len total=", total.shape)
    iteration = 100
    n_test_point = lag//50
    iter = []
    iter_random = []
    for j in np.arange(iteration):
        zernike_alpha = []
        random_zernike_alpha = []
        for a in alpha:
            with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
                index_possible_points = list(set([int(float(x)) for x in f.read().split()]))
                number_points = len(index_possible_points)
            random_index_selected = sorted(random.sample(list(np.arange(0, lag)), number_points))

            surface_random = surf[[i for i in random_index_selected], :]
            zernike_random = total[:, [i for i in random_index_selected]]

            surface = surf[[i for i in index_possible_points], :]
            zernike = total[:, [i for i in index_possible_points]]

            ##### !!!RANDOM!!! DISTANCE BETWEEN THE MEAN SELECTED AND EXCLUDED ZERNIKE IN A PATCH
            random_zernike_point_mean = []
            zernike_point_mean = []
            random_test_point = sorted(random.sample(list(np.arange(0,lag)), n_test_point))
            for i in random_test_point:
            #for i in np.arange(lag)[::50]:
                sys.stderr.write("\r Processing %i out of %i point for alpha = %i, iteration %i" % (i, lag, a, j))
                sys.stderr.flush()
                d2_random = (surface_random[:, 0] - surf[i, 0]) ** 2 + (
                        surface_random[:, 1] - surf[i, 1]) ** 2 + (
                             surface_random[:, 2] - surf[i, 2]) ** 2
                mask_random = d2_random <= Rs_select ** 2
                zernike_patch_random = zernike_random[:, mask_random]

                d2 = (surface[:, 0] - surf[i, 0]) ** 2 + (
                        surface[:, 1] - surf[i, 1]) ** 2 + (
                             surface[:, 2] - surf[i, 2]) ** 2
                mask = d2 <= Rs_select ** 2
                zernike_patch = zernike[:, mask]

                if len(zernike_patch_random[1]) != 0:
                    zernike_patch_mean_random = np.mean(np.array(zernike_patch_random), axis=1)
                else:
                    zernike_patch_index_random = my_functions.find_nearest_vector(surface_random[:, 0:3], surf[i, 0:3])
                    zernike_patch_mean_random = np.mean(np.array(zernike_random[:,[zernike_patch_index_random]]), axis=1)

                random_zernike_point_mean.append(zernike_patch_mean_random.mean())

                if len(zernike_patch[1]) != 0:
                    zernike_patch_mean = np.mean(np.array(zernike_patch), axis=1)
                else:
                    zernike_patch_index = my_functions.find_nearest_vector(surface[:,0:3],surf[i,0:3])
                    zernike_patch_mean = np.mean(np.array(zernike[:,[zernike_patch_index]]),axis=1)

                zernike_point_mean.append(zernike_patch_mean.mean())

            zernike_alpha.append(np.mean(zernike_point_mean))
            random_zernike_alpha.append(np.mean(random_zernike_point_mean))

        iter.append(zernike_alpha)
        iter_random.append(random_zernike_alpha)

    final = np.mean(iter,axis=0)
    final_random = np.mean(iter_random,axis=0)
    np.savetxt("{}\mean_zernike_mean_coeff_iter_{}_points_{}.txt".format(respath, iteration, n_test_point), final)
    np.savetxt("{}\mean_zernike_mean_coeff_rand_iter_{}_points_{}.txt".format(respath, iteration, n_test_point), final_random)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Mean value of the mean selected Zernike coefficients")
    plt.plot(alpha, final, label='Points selected between the sampled, averaged over {} iteration'.format(iteration))
    plt.plot(alpha, final_random, label='Points selected randomly, averaged over {} iteration'.format(iteration))
    leg = ax1.legend()
    plt.show()
    leg = ax1.legend()
    plt.show()



        ###############################
print("Vuoi vedere i grafici gia' ottenuti? y o n?")
o = input()
if o == 'y':
    iteration = 100
    n_test_point = lag//50

    sample = np.loadtxt("{}\mean_zernike_mean_coeff_iter_{}_points_{}.txt".format(respath, iteration, n_test_point))
    randomm = np.loadtxt("{}\mean_zernike_mean_coeff_rand_iter_{}_points_{}.txt".format(respath, iteration, n_test_point))


    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Mean value of the mean selected Zernike coefficients")
    plt.plot(alpha, sample, label='Points selected between the sampled, averaged over {} iteration'.format(iteration))
    plt.plot(alpha, randomm, label='Points selected randomly, averaged over {} iteration'.format(iteration))
    leg = ax1.legend()
    plt.show()
    leg = ax1.legend()
    plt.show()

    sample_smooth = savgol_filter(sample,51, 5)  # window size 51, polynomial order 3

    kn = KneeLocator(alpha, sample_smooth, curve='concave', direction='increasing')

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Mean value of the mean selected Zernike coefficients")
    plt.plot(alpha, sample, label='Points selected between the sampled, averaged over {} iteration'.format(iteration))
    plt.plot(alpha, sample_smooth, label='Smoothed function'.format(iteration))
    plt.plot(kn.knee, kn.knee_y, '*', label='Elbow point at $\\alpha$={}'.format(kn.knee))
    leg = ax1.legend()
    plt.vlines(kn.knee, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
    plt.show()


        ############################################
print("Qual'e' il miglior valore per alpha?n se non vuoi vedere")
best_alpha = input()
if best_alpha != 'n':
    with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select,
                                                                            best_alpha)) as f:
        index_possible_points = [int(float(x)) for x in f.read().split()]
    x = total[:, index_possible_points]
    z_total = np.loadtxt("{}\COSINE_step_{}_Rs_{}.txt".format(respath, step, Rs_select), delimiter=" ")
    z = z_total[index_possible_points]

    TOT_reduced, total_eigenvector_subset, total_eigenvalues_subset, total_eigenvalues = my_functions.PCA(total, 2)
    tot_x_ax = TOT_reduced[:, 0]
    tot_y_ax = TOT_reduced[:, 1]

    X = x.transpose()
    number_points = x.shape[1]
    reduced = np.dot(total_eigenvector_subset.transpose(), X.transpose()).transpose()
    x_ax = reduced[:, 0]
    y_ax = reduced[:, 1]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("2D projection of the Zernike coefficients")
    ax1.set_xlabel('Projection on eigenvector 1')
    ax1.set_ylabel('Projection on eigenvector 2')
    sc = plt.scatter(tot_x_ax, tot_y_ax, c='black')
    cm = plt.cm.get_cmap('RdYlBu')
    sc = plt.scatter(x_ax, y_ax, c=z, cmap=cm)
    plt.colorbar(sc, format='%.1f',
                 label="Mean cosine value of the patches selected with $\\alpha$={}".format(best_alpha))

    plt.legend(["All {} points of the surface".format(lag),
                "{} points sampled with alpha={}".format(number_points, best_alpha)])

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("2D projection of the Zernike coefficients")
    ax1.set_xlabel('Projection on eigenvector 1')
    ax1.set_ylabel('Projection on eigenvector 2')
    cm = plt.cm.get_cmap('RdYlBu')
    sc = plt.scatter(tot_x_ax, tot_y_ax, c=z_total, cmap=cm)
    plt.colorbar(sc, format='%.1f', label="Mean cosine value of the patch")

    plt.show()



        #######################################################
print('Vuoi vedere la proiezione dei punti esclusi su quelli selezionati?')
o = input()
if o == 'y':
    alpha = [1., 10., 20., 30., 40., 50., 60., 70., 80., 90.]

    maximum_min_sampling_alpha = []
    maximum_min_random_alpha = []
    for a in alpha:

        with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
            index_possible_points = list(set([int(float(x)) for x in f.read().split()]))
            number_points = len(index_possible_points)
        random_index_selected = sorted(random.sample(list(np.arange(0, lag)), number_points))
        r_selected_surface = surf[[i for i in random_index_selected], :]
        r_selected_zernike = total[:, [i for i in random_index_selected]]
        r_excluded_surface = surf[[i for i in range(lag) if i not in random_index_selected], :]
        r_excluded_zernike = total[:, [i for i in range(lag) if i not in random_index_selected]]

        selected_surface = surf[[i for i in index_possible_points], :]
        selected_zernike = total[:, [i for i in index_possible_points]]
        excluded_surface = surf[[i for i in range(lag) if i not in index_possible_points], :]
        excluded_zernike = total[:,[i for i in range(lag) if i not in index_possible_points]]

            ##### MAXIMUM BETWEEN THE MINIMUM DISTANCE BETWEEN THE NON SELECTED AND SELECTED ZERNIKE VECTORS
        maximum_min_sampling = []
        maximum_min_random = []
        for i in np.arange(lag)[::50]:
            sys.stderr.write("\r Processing %i out of %i point for alpha = %i" % (i, lag, a))
            sys.stderr.flush()
            selected_index = i

            d_sampling = (selected_surface[:, 0] - surf[selected_index, 0]) ** 2 + (
                    selected_surface[:, 1] - surf[selected_index, 1]) ** 2 + (
                                 selected_surface[:, 2] - surf[selected_index, 2]) ** 2
            mask_sampling = d_sampling <= Rs_select ** 2
            zernike_patch_sampling = selected_zernike[:, mask_sampling]

            d_sampling_excluded = (excluded_surface[:, 0] - surf[selected_index, 0]) ** 2 + (
                    excluded_surface[:, 1] - surf[selected_index, 1]) ** 2 + (
                                 excluded_surface[:, 2] - surf[selected_index, 2]) ** 2
            mask_sampling_excluded = d_sampling_excluded <= Rs_select ** 2
            zernike_patch_sampling_excluded = excluded_zernike[:, mask_sampling_excluded]

            d_random = (r_selected_surface[:, 0] - surf[selected_index, 0]) ** 2 + (
                    r_selected_surface[:, 1] - surf[selected_index, 1]) ** 2 + (
                               r_selected_surface[:, 2] - surf[selected_index, 2]) ** 2
            mask_random = d_random <= Rs_select ** 2
            zernike_patch_random = r_selected_zernike[:, mask_random]

            d_random_excluded = (r_excluded_surface[:, 0] - surf[selected_index, 0]) ** 2 + (
                    r_excluded_surface[:, 1] - surf[selected_index, 1]) ** 2 + (
                               r_excluded_surface[:, 2] - surf[selected_index, 2]) ** 2
            mask_random_excluded = d_random_excluded <= Rs_select ** 2
            zernike_patch_random_excluded = r_excluded_zernike[:, mask_random_excluded]

            zernike_patch_sampling = np.array(zernike_patch_sampling).transpose()
            zernike_patch_random = np.array(zernike_patch_random).transpose()
            zernike_patch_sampling_excluded = np.array(zernike_patch_sampling_excluded).transpose()
            zernike_patch_random_excluded = np.array(zernike_patch_random_excluded).transpose()

            sampling = cdist(zernike_patch_sampling_excluded, zernike_patch_sampling)
            randomm = cdist(zernike_patch_random_excluded, zernike_patch_random)

            sampling = sampling.min(axis=1)
            randomm = randomm.min(axis=1)
            print('sampl', sampling.mean())
            print('rand', randomm.mean())


            maximum_min_sampling.append(sampling.mean())
            maximum_min_random.append(randomm.mean())

        maximum_min_sampling_alpha.append(np.mean(maximum_min_sampling))
        maximum_min_random_alpha.append(np.mean(maximum_min_random))
    # np.savetxt("{}\zernike_scalar_dist_total_sampl.txt".format(respath),mean_nearest_total_sampling_alpha)
    # np.savetxt("{}\zernike_scalar_dist_total_rand.txt".format(respath), mean_nearest_total_random_alpha)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Distance between Zernike coefficients")
    plt.plot(alpha, maximum_min_sampling_alpha, label='Sampling')
    plt.plot(alpha, maximum_min_random_alpha, label='Random')
    leg = ax1.legend()
    plt.show()

        ####################################################

print('Vuoi vedere la differenza dello Zernike variance per patch total-sampling e total-random?')
o = input()
if o == 'y':
    alpha = [1., 10., 20., 30., 40., 50., 60., 70., 80., 90.]

    mean_total_sampling_alpha = []
    mean_total_random_alpha = []
    for a in alpha:
        with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
            index_possible_points = list(set([int(float(x)) for x in f.read().split()]))
            number_points = len(index_possible_points)
        random_index_selected = sorted(random.sample(list(np.arange(0, lag)), number_points))
        r_selected_surface = surf[[i for i in random_index_selected], :]
        r_selected_zernike = total[:, [i for i in random_index_selected]]

        selected_surface = surf[[i for i in index_possible_points], :]
        selected_zernike = total[:, [i for i in index_possible_points]]

        ##### DISTANCE BETWEEN EACH ZERNIKE VECTOR OF THE TOTAL PATCH (NO SAMPLING) AND THE NEAREST VECTOR IN THE SAMPLING AND RANDOM
        mean_total_sampling = []
        mean_total_random = []
        for i in np.arange(lag)[::50]:
            sys.stderr.write("\r Processing %i out of %i point for alpha = %i" % (i, lag, a))
            sys.stderr.flush()
            #selected_index = index_possible_points[i]
            selected_index = i
            d_total = (surf[:, 0] - surf[selected_index, 0]) ** 2 + (
                    surf[:, 1] - surf[selected_index, 1]) ** 2 + (
                         surf[:, 2] - surf[selected_index, 2]) ** 2
            mask_total = d_total <= Rs_select ** 2
            zernike_patch_total = total[:, mask_total]

            d_sampling = (selected_surface[:, 0] - surf[selected_index, 0]) ** 2 + (
                    selected_surface[:, 1] - surf[selected_index, 1]) ** 2 + (
                              selected_surface[:, 2] - surf[selected_index, 2]) ** 2
            mask_sampling = d_sampling <= Rs_select ** 2
            zernike_patch_sampling = selected_zernike[:, mask_sampling]

            d_random = (r_selected_surface[:, 0] - surf[selected_index, 0]) ** 2 + (
                    r_selected_surface[:, 1] - surf[selected_index, 1]) ** 2 + (
                              r_selected_surface[:, 2] - surf[selected_index, 2]) ** 2
            mask_random = d_random <= Rs_select ** 2
            zernike_patch_random = r_selected_zernike[:, mask_random]

            zernike_patch_total = np.array(zernike_patch_total)
            zernike_patch_sampling = np.array(zernike_patch_sampling)
            zernike_patch_random = np.array(zernike_patch_random)

            zernike_patch_total_mean = zernike_patch_total.mean(axis=1)
            zernike_patch_sampling_mean = zernike_patch_sampling.mean(axis=1)
            zernike_patch_random_mean = zernike_patch_random.mean(axis=1)

            zernike_patch_total_var = zernike_patch_total.var(axis=1)
            zernike_patch_sampling_var = zernike_patch_sampling.var(axis=1)
            zernike_patch_random_var = zernike_patch_random.var(axis=1)

            zernike_patch_total_div = np.true_divide(zernike_patch_total_var, zernike_patch_total_mean).reshape(-1, 1)
            zernike_patch_sampling_div = np.true_divide(zernike_patch_sampling_var, zernike_patch_sampling_mean).reshape(-1 ,1)
            zernike_patch_random_div = np.true_divide(zernike_patch_random_var, zernike_patch_random_mean).reshape(-1, 1)

            total_sampling = sum(abs(zernike_patch_total_div-zernike_patch_sampling_div))
            total_random = sum(abs(zernike_patch_total_div-zernike_patch_random_div))

            mean_total_sampling.append(total_sampling)
            mean_total_random.append(total_random)

        mean_total_sampling_alpha.append(np.mean(mean_total_sampling))
        mean_total_random_alpha.append(np.mean(mean_total_random))
    #np.savetxt("{}\zernike_scalar_dist_total_sampl.txt".format(respath),mean_nearest_total_sampling_alpha)
    #np.savetxt("{}\zernike_scalar_dist_total_rand.txt".format(respath), mean_nearest_total_random_alpha)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Distance between Zernike coefficients")
    plt.plot(alpha, mean_total_sampling_alpha, label='Selected by sampling')
    plt.plot(alpha, mean_total_random_alpha, label='Selected randomly')
    leg = ax1.legend()
    plt.show()
            #####################################

print('Vuoi vedere la differenza sampling-random dellla distanza minima rispetto a ciascuno dei Zernike totali?')
o = input()
if o == 'y':
    alpha = [-3., 1., 10., 20., 30., 40., 50., 60., 70., 80., 90.]

    mean_nearest_total_sampling_alpha = []
    mean_nearest_total_random_alpha = []
    for a in alpha:
        with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
            index_possible_points = list(set([int(float(x)) for x in f.read().split()]))
            number_points = len(index_possible_points)
        random_index_selected = sorted(random.sample(list(np.arange(0, lag)), number_points))
        r_selected_surface = surf[[i for i in random_index_selected], :]
        r_selected_zernike = total[:, [i for i in random_index_selected]]

        selected_surface = surf[[i for i in index_possible_points], :]
        selected_zernike = total[:, [i for i in index_possible_points]]

        ##### DISTANCE BETWEEN EACH ZERNIKE VECTOR OF THE TOTAL PATCH (NO SAMPLING) AND THE NEAREST VECTOR IN THE SAMPLING AND RANDOM
        mean_nearest_total_sampling = []
        mean_nearest_total_random = []
        #for i in range(number_points):
        for i in np.arange(lag)[::50]:
            sys.stderr.write("\r Processing %i out of %i point for alpha = %i" % (i, number_points, a))
            sys.stderr.flush()
            #selected_index = index_possible_points[i]
            selected_index = i
            zernike_point = total[:, selected_index]

            d_total = (surf[:, 0] - surf[selected_index, 0]) ** 2 + (
                    surf[:, 1] - surf[selected_index, 1]) ** 2 + (
                         surf[:, 2] - surf[selected_index, 2]) ** 2
            mask_total = d_total <= Rs_select ** 2
            zernike_patch_total = total[:, mask_total]

            d_sampling = (selected_surface[:, 0] - surf[selected_index, 0]) ** 2 + (
                    selected_surface[:, 1] - surf[selected_index, 1]) ** 2 + (
                              selected_surface[:, 2] - surf[selected_index, 2]) ** 2
            mask_sampling = d_sampling <= Rs_select ** 2
            zernike_patch_sampling = selected_zernike[:, mask_sampling]

            d_random = (r_selected_surface[:, 0] - surf[selected_index, 0]) ** 2 + (
                    r_selected_surface[:, 1] - surf[selected_index, 1]) ** 2 + (
                              r_selected_surface[:, 2] - surf[selected_index, 2]) ** 2
            mask_random = d_random <= Rs_select ** 2
            zernike_patch_random = r_selected_zernike[:, mask_random]

            zernike_patch_total = np.array(zernike_patch_total).transpose()
            zernike_patch_sampling = np.array(zernike_patch_sampling).transpose()
            zernike_patch_random = np.array(zernike_patch_random).transpose()


            total_sampling = cdist(zernike_patch_total, zernike_patch_sampling)
            total_random = cdist(zernike_patch_total, zernike_patch_random)

            total_sampling = total_sampling.min(axis=1)
            total_random = total_random.min(axis=1)

            mean_nearest_total_sampling.append(total_sampling.mean())
            mean_nearest_total_random.append(total_random.mean())

        mean_nearest_total_sampling_alpha.append(np.average(mean_nearest_total_sampling))
        mean_nearest_total_random_alpha.append(np.average(mean_nearest_total_random))
    #np.savetxt("{}\zernike_comparison_total_sampl.txt".format(respath),mean_nearest_total_sampling_alpha)
    #np.savetxt("{}\zernike_comparison_total_rand.txt".format(respath), mean_nearest_total_random_alpha)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Distance between Zernike coefficients")
    plt.plot(alpha, mean_nearest_total_sampling_alpha, label='Selected by sampling')
    plt.plot(alpha, mean_nearest_total_random_alpha, label='Selected randomly')
    leg = ax1.legend()
    plt.show()

print('Vuoi studiare il caso con tante iterazioni?')
o = input()
iteration = 10
if o == 'y':
    zernike_euclidean = []
    for iter in np.arange(iteration):
        zernike_euclidean.append(np.loadtxt("{}\zernike_euclidean_iter\zernike_euclidean_iter_{}.txt".format(respath, iter)))

    print('Vuoi vedere il grafico per ogni iter?')
    o = input()
    if o == 'y':
        for iter in np.arange(iteration):
            min_norm = min(zernike_euclidean[iter])
            min_alpha = alpha[np.where(min_norm == zernike_euclidean[iter])[0][0]]

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.set_xlabel('$\\alpha$ value')
            ax1.set_ylabel("Distance between Zernike coefficients")
            plt.plot(alpha, zernike_euclidean[iter], label='Selected by sampling, iteration {}'.format(iter))
            plt.plot(min_alpha, min_norm, '*', label='Minimum at $\\alpha$={}'.format(min_alpha))
            leg = ax1.legend()
        plt.show()

    zernike_euclidean_mean = np.mean(np.array(zernike_euclidean), axis=0)
    np.savetxt("{}\zernike_euclidean_averaged.txt".format(respath),zernike_euclidean_mean)

    min_norm = min(zernike_euclidean_mean)
    min_alpha = alpha[np.where(min_norm == zernike_euclidean_mean)[0][0]]
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Distance between Zernike coefficients")
    plt.plot(alpha, zernike_euclidean_mean, label='Selected by sampling, averaged over {} iterations'.format(iteration))
    plt.plot(min_alpha, min_norm, '*', label='Minimum at $\\alpha$={}'.format(min_alpha))
    leg = ax1.legend()
    plt.show()


print("Vuoi vedere il numero di punti in funzione di alpha? y o n?")
o = input()
if o == "y":
    alpha = [
        #-25., -24., -23., -22., -21., -20.,
        #-19., -18., -17., -16., -15., -14., -13., -11., -10.,
        -9., -8., -7., -6., -5., -4., -3., -2., -1.,
         1., 2., 3., 4., 5., 6., 7., 8., 9.
        , 10., 11., 12., 13., 14., 15., 16., 17., 18., 19.
        , 20., 21., 22., 23., 24., 25., 26., 27., 28., 29.
        , 30., 31., 32., 33., 34., 35., 36., 37., 38., 39.
        , 40., 41., 42., 43., 44., 45., 46., 47., 48., 49.
        , 50., 51., 52., 53., 54., 55., 56., 57., 58., 59.
        , 60., 61., 62., 63., 64., 65., 66., 67., 68., 69.
        , 70., 71., 72., 73., 74., 75., 76., 77., 78., 79.
        , 80., 81., 82., 83., 84., 85., 86., 87., 88., 89.
        , 90., 91., 92., 93., 94., 95., 96.
    ]
    alpha_points_plot = []
    for a in alpha:
        with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
            index_possible_points = list(set([int(float(x)) for x in f.read().split()]))
            alpha_points_plot.append(len(index_possible_points))

    print(alpha_points_plot)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Number of sampled points over {} points".format(lag))
    plt.plot(alpha, alpha_points_plot)
    plt.show()


print("Vuoi vedere i grafici gia' ottenuti? y o n?")
o = input()
if o == 'y':
    alpha_random = [ 1., 10., 20., 30., 40., 50., 60., 70., 80., 90.]

    zernike_euclidean = np.loadtxt("{}\zernike_euclidean.txt".format(respath))
    z_euclidean_alpha = np.loadtxt("{}\z_euclidean_alpha.txt".format(respath))

    min_norm = min(zernike_euclidean)
    min_alpha = z_euclidean_alpha[np.where(min_norm==zernike_euclidean)[0][0]]

    random_norm_alpha = np.loadtxt("{}\\random_norm_alpha.txt".format(respath))

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Distance between Zernike coefficients")
    plt.plot(z_euclidean_alpha, zernike_euclidean, label='Selected by sampling')
    #plt.plot(min_alpha,min_norm, '*', label='Minimum at $\\alpha$={}'.format(min_alpha))
    leg = ax1.legend()


    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Distance between Zernike coefficients")
    plt.plot(z_euclidean_alpha, zernike_euclidean, label='Selected by sampling')
    plt.plot(alpha_random, random_norm_alpha, label='Randomly selected')
    leg = ax1.legend()

    plt.show()

print('Vuoi trovare l\'alpha migliore? y o n?')
o = input()
if o == 'y':
    print("len total=", total.shape)
    zernike_alpha = []
    random_zernike_alpha = []
    #alpha = [-3., 10., 20., 30., 40., 50., 60., 70., 80., 90.]
    #zernike_alpha_excluded = []
    for a in alpha:
        with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
            index_possible_points = list(set([int(float(x)) for x in f.read().split()]))
            number_points = len(index_possible_points)
        random_index_selected = sorted(random.sample(list(np.arange(0, lag)), number_points))


        removed_surface_random = surf[[i for i in range(lag) if i not in random_index_selected], :]
        removed_zernike_random = total[:, [i for i in range(lag) if i not in random_index_selected]]
        surface_random = surf[[i for i in random_index_selected], :]
        zernike_random = total[:, [i for i in random_index_selected]]

        removed_surface = surf[[i for i in range(lag) if i not in index_possible_points], :]
        removed_zernike = total[:, [i for i in range(lag) if i not in index_possible_points]]
        surface = surf[[i for i in index_possible_points], :]
        zernike = total[:, [i for i in index_possible_points]]

        ##### !!!RANDOM!!! DISTANCE BETWEEN THE MEAN SELECTED AND EXCLUDED ZERNIKE IN A PATCH
        random_zernike_point_mean = []
        #for i in random_index_selected:
        for i in np.arange(lag)[::50]:
            sys.stderr.write("\r Processing %i out of %i point for alpha = %i" % (i, len(index_possible_points), a))
            sys.stderr.flush()
            #zernike_point = total[:, [i]]
            d2_removed = (removed_surface_random[:, 0] - surf[i, 0]) ** 2 + (
                    removed_surface_random[:, 1] - surf[i, 1]) ** 2 + (
                         removed_surface_random[:, 2] - surf[i, 2]) ** 2
            mask_removed = d2_removed <= Rs_select ** 2
            zernike_patch_removed = removed_zernike_random[:, mask_removed]

            d2 = (surface_random[:, 0] - surf[i, 0]) ** 2 + (
                    surface_random[:, 1] - surf[i, 1]) ** 2 + (
                         surface_random[:, 2] - surf[i, 2]) ** 2
            mask = d2 <= Rs_select ** 2
            zernike_patch = zernike_random[:, mask]

            if len(zernike_patch_removed[1]) != 0:
                zernike_patch_removed_mean = np.mean(np.array(zernike_patch_removed), axis=1)
            else:
                zernike_patch_removed_index = my_functions.find_nearest_vector(removed_surface_random[:, 0:3], surf[i, 0:3])
                zernike_patch_removed_mean = np.mean(np.array(removed_zernike_random[:, [zernike_patch_removed_index]]), axis=1)


            if len(zernike_patch[1]) != 0:
                zernike_patch_mean = np.mean(np.array(zernike_patch), axis=1)
            else:
                zernike_patch_index = my_functions.find_nearest_vector(surface_random[:, 0:3], surf[i, 0:3])
                zernike_patch_mean = np.mean(np.array(zernike_random[:,[zernike_patch_index]]),axis=1)


            euc_dist_point = zernike_patch_mean.mean()
            #euc_dist_point = sum(abs(zernike_patch_mean-zernike_patch_removed_mean))
            #euc_dist_point = np.linalg.norm(zernike_patch_mean - zernike_patch_removed_mean)
            random_zernike_point_mean.append(euc_dist_point)

        #if len(zernike_patch[1]) != 0:
            #    euc_dist_point = np.linalg.norm(np.array(zernike_patch) - np.array(zernike_point), axis=0)
            #    random_zernike_point_mean.append(np.min(euc_dist_point))
            #else:
            #    zernike_patch_index = my_functions.find_nearest_vector(removed_surface_random[:, 0:3], surf[i, 0:3])
            #    zernike_patch = removed_zernike_random[:, [zernike_patch_index]]
            #    euc_dist_point = np.linalg.norm(np.array(zernike_patch) - np.array(zernike_point), axis=0)
            #    random_zernike_point_mean.append(np.min(euc_dist_point))

        #random_zernike_alpha.append(np.mean(random_zernike_point_mean))

            ##### DISTANCE BETWEEN THE ZERNIKE VECTOR OF EACH SELECTED POINT AND THE SURROUNDING EXCLUDED POINTS
        zernike_point_mean = []
        #for i in index_possible_points:
        for i in np.arange(lag)[::50]:

            sys.stderr.write("\r Processing %i out of %i point for alpha = %i" % (i, len(index_possible_points), a))
            sys.stderr.flush()
            #zernike_point = total[:, [i]]
            d2_removed = (removed_surface[:, 0] - surf[i, 0]) ** 2 + (
                        removed_surface[:, 1] - surf[i, 1]) ** 2 + (
                             removed_surface[:, 2] - surf[i, 2]) ** 2
            mask_removed = d2_removed <= Rs_select ** 2
            zernike_patch_removed = removed_zernike[:, mask_removed]
            #print('patch removed', np.array(zernike_patch_removed).shape)

            d2 = (surface[:, 0] - surf[i, 0]) ** 2 + (
                        surface[:, 1] - surf[i, 1]) ** 2 + (
                             surface[:, 2] - surf[i, 2]) ** 2
            mask = d2 <= Rs_select ** 2
            zernike_patch = zernike[:, mask]

            if len(zernike_patch_removed[1]) != 0:
                zernike_patch_removed_mean = np.mean(np.array(zernike_patch_removed), axis=1)
            else:
                zernike_patch_removed_index = my_functions.find_nearest_vector(removed_surface[:, 0:3], surf[i, 0:3])
                #print('index',zernike_patch_removed_index)
                #print('sele',np.array(removed_zernike[:, [zernike_patch_removed_index]]).shape )
                zernike_patch_removed_mean = np.mean(np.array(removed_zernike[:, [zernike_patch_removed_index]]), axis=1)

            if len(zernike_patch[1]) != 0:
                zernike_patch_mean = np.mean(np.array(zernike_patch), axis=1)
            else:
                zernike_patch_index = my_functions.find_nearest_vector(surface[:,0:3],surf[i,0:3])
                zernike_patch_mean = np.mean(np.array(zernike[:,[zernike_patch_index]]),axis=1)

            euc_dist_point = zernike_patch_mean.mean()
            #euc_dist_point = sum(abs(zernike_patch_mean-zernike_patch_removed_mean))
            #euc_dist_point = np.linalg.norm(zernike_patch_mean - zernike_patch_removed_mean)
            zernike_point_mean.append(euc_dist_point)



            #if len(zernike_patch[1]) != 0:
            #    euc_dist_point = np.linalg.norm(np.array(zernike_patch) - np.array(zernike_point), axis=0)
            #    zernike_point_mean.append(np.min(euc_dist_point))
            #else:
            #    zernike_patch_index = my_functions.find_nearest_vector(removed_surface[:, 0:3], surf[i, 0:3])
            #    zernike_patch = removed_zernike[:,[zernike_patch_index]]
            #    euc_dist_point = np.linalg.norm(np.array(zernike_patch) - np.array(zernike_point), axis=0)
            #    zernike_point_mean.append(np.min(euc_dist_point))

            #        euc_dist_point = np.linalg.norm(zernike_patch_mean-zernike_patch_excluded_mean)
            #        zernike_excluded_point_mean.append(euc_dist_point)

        zernike_alpha.append(np.mean(zernike_point_mean))
        random_zernike_alpha.append(np.mean(random_zernike_point_mean))
    #np.savetxt("{}\zernike_euclidean_prova.txt".format(respath),zernike_alpha)

                ###### DIFFERENCE BETWEEN MEAN ZERNIKE VECTOR OF EXCLUDED AND SELECTED POINTS IN A PATCH
    #    selected_surface = surf[[i for i in index_possible_points], :]
    #    selected_zernike = total[:, [i for i in index_possible_points]]

    #    zernike_excluded_point_mean = []
    #    for i in range(len(selected_zernike[1])):
    #        sys.stderr.write("\r Processing %i out of %i point for alpha = %i" % (i, len(selected_zernike[1]), a))
    #        sys.stderr.flush()
    #        d_excluded = (removed_surface[:, 0]-selected_surface[i, 0]) ** 2 + (
    #                removed_surface[:, 1] - selected_surface[i, 1]) ** 2 + (
    #                     removed_surface[:, 2] - selected_surface[i, 2]) ** 2
    #        mask_excluded = d_excluded <= Rs_select ** 2
    #        zernike_patch_excluded = removed_zernike[:, mask_excluded]
    #        if len(zernike_patch_excluded[1]) != 0:
    #           zernike_patch_excluded_mean = np.mean(np.array(zernike_patch_excluded), axis=1)
    #        else:
    #            zernike_patch_excluded_index = my_functions.find_nearest_vector(removed_surface[:, 0:3], selected_surface[i, 0:3])
    #            zernike_patch_excluded_mean = np.mean(np.array(removed_zernike[:, [zernike_patch_excluded_index]]),axis=1)

    #        d = (selected_surface[:, 0] - selected_surface[i, 0]) ** 2 + (
    #                selected_surface[:, 1] - selected_surface[i, 1]) ** 2 + (
    #                             selected_surface[:, 2] - selected_surface[i, 2]) ** 2
    #        mask = d <= Rs_select ** 2
    #        zernike_patch = selected_zernike[:, mask]
    #        zernike_patch_mean = np.mean(np.array(zernike_patch), axis=1)
    #        euc_dist_point = np.linalg.norm(zernike_patch_mean-zernike_patch_excluded_mean)
    #        zernike_excluded_point_mean.append(euc_dist_point)


     #   zernike_alpha_excluded.append(np.mean(zernike_excluded_point_mean))

#    np.savetxt("{}\zernike_euclidean_diff.txt".format(respath),zernike_alpha)
#    np.savetxt("{}\z_euclidean_alpha.txt".format(respath),alpha)
    print('MIN=', alpha[np.where(np.array(zernike_alpha)==min(zernike_alpha))[0][0]])

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Zernike coefficients' mean euclidean distance")
    plt.plot(alpha, zernike_alpha, label = 'Sampling')
    plt.plot(alpha, random_zernike_alpha, label='Random')
    leg = ax1.legend()
    plt.show()


print('Vuoi studiare il caso random? y o n?')
o = input()
if o == 'y':
    alpha_random = [ 1., 10., 20., 30., 40., 50., 60., 70., 80., 90.]
    random_norm_alpha = []
    zernike_alpha = []
    for a in alpha_random:
        with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
            index_possible_points = [int(float(x)) for x in f.read().split()]
        number_points = len(index_possible_points)

        random_index_selected = random.sample(list(np.arange(0, lag)), number_points)
        removed_surface = surf[[i for i in range(lag) if i not in random_index_selected], :]
        removed_zernike = total[:, [i for i in range(lag) if i not in random_index_selected]]

        random_norm_point = []
        for i in random_index_selected:
            sys.stderr.write("\r Processing %i out of %i point for alpha = %i" % (i, len(index_possible_points), a))
            sys.stderr.flush()
            zernike_random_selected = total[:, [i]]

            d2 = (removed_surface[:, 0] - surf[i, 0]) ** 2 + (
                    removed_surface[:, 1] - surf[i, 1]) ** 2 + (
                         removed_surface[:, 2] - surf[i, 2]) ** 2
            mask = d2 <= Rs_select ** 2
            zernike_patch = removed_zernike[:, mask]

            euc_dist_point = np.linalg.norm(np.array(zernike_patch) - np.array(zernike_random_selected), axis=0)
            random_norm_point.append(np.mean(euc_dist_point))
        random_norm_alpha.append(np.mean(random_norm_point))

            ######################################################
       # removed_surface = surf[[i for i in range(lag) if i not in random_index_selected], :]
       # removed_zernike = total[:, [i for i in range(lag) if i not in random_index_selected]]

       # zernike_point_mean = []
       # for i in random_index_selected:
       #     sys.stderr.write("\r Processing %i out of %i point for alpha = %i" % (i, len(random_index_selected), a))
       #     sys.stderr.flush()
       #     zernike_point = total[:, [i]]
       #     d2 = (removed_surface[:, 0] - surf[i, 0]) ** 2 + (
       #             removed_surface[:, 1] - surf[i, 1]) ** 2 + (
       #                  removed_surface[:, 2] - surf[i, 2]) ** 2
       #     mask = d2 <= Rs_select ** 2
       #     zernike_patch = removed_zernike[:, mask]
       #     if len(zernike_patch[1]) != 0:
       #         euc_dist_point = np.linalg.norm(np.array(zernike_patch) - np.array(zernike_point), axis=0)
       #         zernike_point_mean.append(np.mean(euc_dist_point))
       #     else:
       #         zernike_patch_index = my_functions.find_nearest_vector(removed_surface[:, 0:3], surf[i, 0:3])
       #         zernike_patch = removed_zernike[:, [zernike_patch_index]]
       #         euc_dist_point = np.linalg.norm(np.array(zernike_patch) - np.array(zernike_point), axis=0)
       #         zernike_point_mean.append(np.mean(euc_dist_point))

       # zernike_alpha.append(np.mean(zernike_point_mean))
    np.savetxt("{}\\random_norm_alpha.txt".format(respath), random_norm_alpha)
    np.savetxt("{}\\alpha_random_norm_alpha.txt".format(respath), alpha_random)


    zernike_euclidean = np.loadtxt("{}\zernike_euclidean_prova.txt".format(respath))
    z_euclidean_alpha = np.loadtxt("{}\z_euclidean_alpha.txt".format(respath))

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Distance between Zernike coefficients")
    plt.plot(z_euclidean_alpha, zernike_euclidean, label='Selected by sampling')
    plt.plot(alpha_random, random_norm_alpha, label='Randomly selected')
    leg = ax1.legend()

    plt.show()



print('Vuoi fare lo smoothing? y o n?')
o = input()
if o == 'y':
    zernike_euclidean = np.loadtxt("{}\zernike_euclidean.txt".format(respath))
    z_euclidean_alpha = np.loadtxt("{}\z_euclidean_alpha.txt".format(respath))

    z_euclidean_smooth = savgol_filter(zernike_euclidean, 9, 3) # window size 51, polynomial order 3

    min = min(z_euclidean_smooth)
    min_alpha = z_euclidean_alpha[np.where(min == z_euclidean_smooth)[0]]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Zernike coefficients' mean euclidean distance")
    plt.plot(z_euclidean_alpha, zernike_euclidean, label='Original function')
    plt.plot(z_euclidean_alpha, z_euclidean_smooth, label='Smoothed function' )
    plt.plot(min_alpha, min, '*', label='Minimum at $\\alpha$={}'.format(min_alpha[0]))
    leg = ax1.legend()
    plt.show()

print("Qual'e' il miglior valore per alpha?")
best_alpha = input()

with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, best_alpha)) as f:
    index_possible_points = [int(float(x)) for x in f.read().split()]
x = total[:, index_possible_points]
z_total = np.loadtxt("{}\COSINE_step_{}_Rs_{}.txt".format(respath, step, Rs_select), delimiter=" ")
z = z_total[index_possible_points]

TOT_reduced, total_eigenvector_subset, total_eigenvalues_subset, total_eigenvalues = my_functions.PCA(total,2)
tot_x_ax = TOT_reduced[:, 0]
tot_y_ax = TOT_reduced[:, 1]

X = x.transpose()
number_points = x.shape[1]
reduced = np.dot(total_eigenvector_subset.transpose(), X.transpose()).transpose()
x_ax = reduced[:, 0]
y_ax = reduced[:, 1]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("2D projection of the Zernike coefficients")
ax1.set_xlabel('Projection on eigenvector 1')
ax1.set_ylabel('Projection on eigenvector 2')
sc = plt.scatter(tot_x_ax, tot_y_ax, c='black')
cm = plt.cm.get_cmap('RdYlBu')
sc = plt.scatter(x_ax, y_ax, c=z, cmap=cm)
plt.colorbar(sc, format='%.1f', label="Mean cosine value of the patches selected with $\\alpha$={}".format(best_alpha))

plt.legend(["All {} points of the surface".format(lag),
            "{} points sampled with alpha={}".format(number_points, best_alpha)])



fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("2D projection of the Zernike coefficients")
ax1.set_xlabel('Projection on eigenvector 1')
ax1.set_ylabel( 'Projection on eigenvector 2')
cm = plt.cm.get_cmap('RdYlBu')
sc = plt.scatter(tot_x_ax, tot_y_ax, c=z_total, cmap=cm)
plt.colorbar(sc, format='%.1f', label= "Mean cosine value of the patch")


plt.show()



#print("Vuoi fare la media di tutti i casi? y o n?")
#o = input()
#if o == 'y':
#    z_euclidean_alpha = np.loadtxt("{}\z_euclidean_alpha.txt".format(respath))

#    zernike_euclidean_all = np.array(5,len(z_euclidean_alpha))
#    for i in [1,2,3,4,5]:
#        respath = "..\\{}\cluster{}\R_zernike_{}\R_s_{}".format(fragment, i, R_zernike, Rs_select)
#        zernike_euclidean_all(i) = np.loadtxt("{}\zernike_euclidean.txt".format(respath))

#    print(np.array(zernike_euclidean_all).shape)
#    zernike_euclidean = statistics.mean(zernike_euclidean_all)
#    print(np.array(zernike_euclidean).shape)
#    z_euclidean_alpha = np.loadtxt("{}\z_euclidean_alpha.txt".format(respath))

#    min_norm = min(zernike_euclidean)
#    min_alpha = z_euclidean_alpha[np.where(min_norm==zernike_euclidean)[0][0]]


#    fig = plt.figure()
#    ax1 = fig.add_subplot(111)
#    ax1.set_xlabel('$\\alpha$ value')
#    ax1.set_ylabel("Mean distance between Zernike coefficients")
#    plt.plot(z_euclidean_alpha, zernike_euclidean, label='Averaged over all the configurations of fragment A')
#    plt.plot(min_alpha,min_norm, '*', label='Minimum at $\\alpha$={}'.format(min_alpha))
#    leg = ax1.legend()
#    plt.show()




print("Vuoi salvare i coefficenti di Zernike dello screening migliore? Se si con che alpha? Altrimenti digita n")
alpha = input()
if alpha != "n":
    # Apro il file con gli indici delle patch di cui voglio salvare Zernike
    shutil.copy("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, alpha),
                "{}\zernike\index_alpha.txt".format(respath))


    with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, alpha)) as f:
        index_possible_points = [int(float(x)) for x in f.read().split()]

    #Apro il file con tutti gli indici di Zernike, compresa la riga 1 che contiene gli indici
    zernike_total_positive = np.loadtxt("{}\zernike\zernike_positive\zernike_total.dat".format(respath), delimiter=" ")
    zernike_total_negative = np.loadtxt("{}\zernike\zernike_negative\zernike_total.dat".format(respath), delimiter=" ")



    zernike_positive_alpha = zernike_total_positive[:,index_possible_points]
    zernike_negative_alpha = zernike_total_negative[:,index_possible_points]



    np.savetxt("{}/zernike/zernike_positive/zernike_alpha.dat".format(respath, alpha), zernike_positive_alpha, fmt="%.4e")
    np.savetxt("{}/zernike/zernike_negative/zernike_alpha.dat".format(respath, alpha), zernike_negative_alpha, fmt="%.4e")



