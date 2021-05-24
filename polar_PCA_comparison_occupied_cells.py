import numpy as np
import matplotlib.pyplot as plt
import my_functions
from scipy.stats import pearsonr
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
from sklearn.cluster import KMeans


    ################### PARAMETERS  ########
with open('configuration.txt') as f:
    for line in f:
        exec(line)

#number_division = 60
number_division = 30

#alpha = [1.5, 1.6, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
#         3., 4., 5., 10.,
#         3.1, 3.2, 3.3, 3.4, 3.6, 3.7, 3.8, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.5]
alpha = [.01, .1, .2, 0.3, .4, .5, .6, .7, .8, .9, 1.]

alpha_plot = [.1]
total_plot = 0 #1 se vuoi vederlo
#fragment = 208
#cluster = 1
#R_zernike = 6
#Rs_select = 2

    ############### END PARAMETERS  ################
respath = "..\\{}\cluster{}\R_zernike_{}\R_s_{}".format(fragment, cluster, R_zernike, Rs_select)

total = np.loadtxt("{}\zernike\zernike_positive\zernike_total.dat".format(respath), delimiter=" ", skiprows=1)
#total = np.loadtxt("{}\zernike\zernike_positive\zernike_total.dat".format(respath), delimiter=" ")

print("len total=", total.shape)
tot_points = total.shape[1]

### PCA of the Zernike coefficients for all patches
TOT_reduced, total_eigenvector_subset = my_functions.PCA(total)
tot_x_ax = TOT_reduced[:, 0]
tot_y_ax = TOT_reduced[:, 1]

### PCA of the Zernike coefficients for all patches in polar coordinates
clusterer = KMeans(n_clusters=1, random_state=10)
cluster_labels_total = clusterer.fit_predict(TOT_reduced)
centroid_total = clusterer.cluster_centers_
theta_total, radius_total = my_functions.cart2polar(TOT_reduced,centroid_total)
                        ### GRID ###
x_grid = np.linspace(min(theta_total), max(theta_total), num=number_division)
y_grid = np.linspace(min(radius_total), max(radius_total), num=number_division)
cell_w = x_grid[1] - x_grid[0]
cell_h = y_grid[1] - y_grid[0]
                        ### END GRID###

if total_plot == 1:
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("2D projection of the Zernike coefficients in polar coordinates")
    ax1.set_xlabel('$\\theta$ (degree)')
    ax1.set_ylabel('r')
    sc = plt.scatter(theta_total, radius_total)
    ax1.xaxis.set_ticks(x_grid)
    ax1.yaxis.set_ticks(y_grid)
    #ax1.set_yticklabels([])
   # ax1.set_xticklabels(x_grid)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%1.f'))
    ax1.grid(True)
    plt.tight_layout()
    plt.show()

reward = []
alpha_points = []
for a in alpha:

    with open("{}\index_percentage\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, Rs_select, a)) as f:
        index_possible_points = [int(float(x)) for x in f.read().split()]
    x = total[:,index_possible_points]
    #x = np.loadtxt("{}\zernike\zernike_positive\zernike_alpha_{}.dat".format(respath,a), delimiter=" ", skiprows=1)

    print("len alpha {}=".format(a), x.shape)

     ### PROIEZIONE SULLE DUE PC DEI COEFFICIENTI DI ZERNIKE PER OGNI PATCH
    X = x.transpose()
    number_points = x.shape[1]
    alpha_points.append(number_points)
    reduced = np.dot(total_eigenvector_subset.transpose(), X.transpose()).transpose()
    x_ax = reduced[:,0]
    y_ax = reduced[:,1]

    #### PROIEZIONE SULLE DUE PC IN COORDINATE POLARI
    theta_alpha, radius_alpha = my_functions.cart2polar(reduced,centroid_total)

    if a in alpha_plot:
    #print("Vuoi fare il plot delle PCA per lo screening totale e con alpha={}? y o n?".format(a))
    #o = input()
    #if o == "y":
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title("PCA: 2D projection of the Zernike coefficients")
        ax1.set_xlabel('Projection on eigenvector 1')
        ax1.set_ylabel('Projection on eigenvector 2)')
        sc = plt.scatter(tot_x_ax, tot_y_ax)
        sc = plt.scatter(x_ax, y_ax, alpha=.3)
        plt.legend(["All {} points of the surface".format(tot_points), "{} points found with alpha={}".format(number_points,a)])
        plt.show()

    #print("Vuoi fare il plot delle PCA per lo screening totale e con alpha={} in COORDINATE POLARI? y o n?".format(a))
    #o = input()
    #if o == "y":
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title("2D projection of the Zernike coefficients in polar coordinates")
        ax1.set_xlabel('$\\theta$')
        ax1.set_ylabel('r')
        ax1.xaxis.set_ticks(x_grid)
        ax1.yaxis.set_ticks(y_grid)
        # ax1.set_yticklabels([])
        ax1.set_xticklabels([])
        ax1.grid(True)
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%1.f'))

        sc = plt.scatter(theta_total, radius_total)
        sc = plt.scatter(theta_alpha, radius_alpha)
        plt.legend(["All {} points of the surface".format(tot_points), "{} points found with alpha={}".format(number_points,a)])
        plt.tight_layout()
        plt.show()


    ##### PEARSON CORRELATION PER IL NUMERO DI PUNTI NELLE GRID #######
    count = np.zeros((len(x_grid), len(y_grid)))
    count_cell_total=np.zeros(len(x_grid)*len(y_grid))
    count_cell_alpha = np.zeros(len(x_grid)*len(y_grid))
    common_points = []

    k = 0
    occupied_cells = 0
    for i in range(len(x_grid)):
        for j in range(len(y_grid)):
            index_cell_total = np.where(
                (theta_total >= x_grid[0] + i * cell_w) & (theta_total < x_grid[0] + (i + 1) * cell_w) & \
                (radius_total>= y_grid[0] + j * cell_h) & (radius_total < y_grid[0] + (j + 1) * cell_h))

            index_cell_alpha = np.where(
                (theta_alpha >= x_grid[0] + i * cell_w) & (theta_alpha < x_grid[0] + (i + 1) * cell_w) & \
                (radius_alpha>= y_grid[0] + j * cell_h) & (radius_alpha < y_grid[0] + (j + 1) * cell_h))

            #punti blu nella cella su tutti i punti blu
            count_cell_total[k] = len(index_cell_total[0])/tot_points

            #punti arancioni nella cella su tutti i punti blu
            count_cell_alpha[k] = len(index_cell_alpha[0])/number_points
            if count_cell_alpha[k] != 0:
                occupied_cells +=1
#                common_points.append((count_cell_total[k]/tot_points)-(count_cell_alpha[k]/number_points))
#                common_points.append((count_cell_total[k]/tot_points)-(count_cell_alpha[k]/tot_points))
                common_points.append(count_cell_total[k]-count_cell_alpha[k])

            k += 1

  #  print("len common points=", np.shape(common_points))
  #  print("number occupied cells=", occupied_cells)
  #  print("sum *=", sum(common_points)*occupied_cells)
    reward.append(sum(common_points)*occupied_cells/number_points)


#np.savetxt("{}\\reward_function.txt".format(respath),reward)


fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax1.set_title("Cost function")
ax1.set_xlabel('$\\alpha$ value')
ax1.set_ylabel("$Common points \\times \\frac{Occupied~cells}{Number~of~points_for~\\alpha}$")

for i in range(len(alpha)):
    plt.plot(alpha[i], reward[i], 'bo')
    plt.text(alpha[i] * (1 + 0.01), reward[i] * (1 + 0.01), "{} points".format(alpha_points[i]))
leg = ax1.legend()
plt.show()