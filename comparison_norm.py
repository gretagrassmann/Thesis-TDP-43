import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import shutil
import my_functions
import sys
from scipy.signal import savgol_filter
import random
import statistics

screening = "index_continuous_distribution_exp"

    ################### PARAMETERS  ########
with open('configuration.txt') as f:
    for line in f:
        exec(line)
fragment = 208
cluster = 1
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
    plt.plot(min_alpha,min_norm, '*', label='Minimum at $\\alpha$={}'.format(min_alpha))
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
    for a in alpha:
        with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
            index_possible_points = list(set([int(float(x)) for x in f.read().split()]))

        removed_surface = surf[[i for i in range(lag) if i not in index_possible_points], :]
        removed_zernike = total[:, [i for i in range(lag) if i not in index_possible_points]]

        zernike_point_mean = []
#        for i in range(len(index_possible_points)):
        for i in index_possible_points:
            sys.stderr.write("\r Processing %i out of %i point for alpha = %i" % (i, len(index_possible_points), a))
            sys.stderr.flush()
            zernike_point = total[:, [i]]
            d2 = (removed_surface[:, 0] - surf[i, 0]) ** 2 + (
                        removed_surface[:, 1] - surf[i, 1]) ** 2 + (
                             removed_surface[:, 2] - surf[i, 2]) ** 2
            mask = d2 <= Rs_select ** 2
            zernike_patch = removed_zernike[:, mask]
            if len(zernike_patch[1]) != 0:
                euc_dist_point = np.linalg.norm(np.array(zernike_patch) - np.array(zernike_point), axis=0)
                zernike_point_mean.append(np.mean(euc_dist_point))

        zernike_alpha.append(np.mean(zernike_point_mean))


    np.savetxt("{}\zernike_euclidean.txt".format(respath),zernike_alpha)
    np.savetxt("{}\z_euclidean_alpha.txt".format(respath),alpha)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Zernike coefficients' mean euclidean distance")
    plt.plot(alpha, zernike_alpha)
    leg = ax1.legend()
    plt.show()


print('Vuoi studiare il caso random? y o n?')
o = input()
if o == 'y':
    alpha_random = [ 1., 10., 20., 30., 40., 50., 60., 70., 80., 90.]
    random_norm_alpha = []
    for a in alpha_random:
        with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
            index_possible_points = [int(float(x)) for x in f.read().split()]
        number_points = len(index_possible_points)

        random_index_selected = random.sample(list(np.arange(0, lag)), number_points)
        random_norm_point = []
        for i in random_index_selected:
            sys.stderr.write("\r Processing %i out of %i point for alpha = %i" % (i, len(index_possible_points), a))
            sys.stderr.flush()
            zernike_random_selected = total[:, [i]]

            d2 = (surf[:, 0] - surf[i, 0]) ** 2 + (
                    surf[:, 1] - surf[i, 1]) ** 2 + (
                         surf[:, 2] - surf[i, 2]) ** 2
            mask = d2 <= Rs_select ** 2
            mask = d2 != 0
            zernike_patch = total[:, mask]

            euc_dist_point = np.linalg.norm(np.array(zernike_patch) - np.array(zernike_random_selected), axis=0)
            random_norm_point.append(np.mean(euc_dist_point))
        random_norm_alpha.append(np.mean(random_norm_point))

    np.savetxt("{}\\random_norm_alpha.txt".format(respath), random_norm_alpha)
    np.savetxt("{}\\alpha_random_norm_alpha.txt".format(respath), alpha_random)


    zernike_euclidean = np.loadtxt("{}\zernike_euclidean.txt".format(respath))
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



