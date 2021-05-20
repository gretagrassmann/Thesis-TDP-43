import numpy as np
import matplotlib.pyplot as plt
import math
from sklearn.cluster import KMeans

    ################### PARAMETERS  ########
with open('configuration.txt') as f:
    for line in f:
        exec(line)

alpha = [.01]
#alpha = [.01, .1, .5, 1.]
#fragment = 208
#cluster = 1
#R_zernike = 6
#Rs_select = 2

    ############### END PARAMETERS  ################
respath = "..\\{}\cluster{}\R_zernike_{}\R_s_{}".format(fragment, cluster, R_zernike, Rs_select)

total = np.loadtxt("{}\zernike\zernike_positive\zernike_total.dat".format(respath), delimiter=" ", skiprows=1)

print("len total=", total.shape)
TOT = total.transpose()
tot_points = total.shape[1]

# TOTAL
cov_mat = np.cov(TOT, rowvar=False)
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
#TOT_reduced_with0 = np.dot(eigenvector_subset.transpose(), TOT.transpose()).transpose()
#TOT_reduced = TOT_reduced_with0[np.all(TOT_reduced_with0 != 0.0, axis=1)]
TOT_reduced = np.dot(eigenvector_subset.transpose(), TOT.transpose()).transpose()
tot_x_ax = TOT_reduced[:, 0]
tot_y_ax = TOT_reduced[:, 1]

    ####### find centroid and inertia of the total  #####
clusterer_total = KMeans(n_clusters=1, random_state=10)
cluster_labels_total = clusterer_total.fit_predict(TOT_reduced)
centroid_total = clusterer_total.cluster_centers_
inertia_total = clusterer_total.inertia_
print("center=",centroid_total)
print("inertia=", inertia_total)

alpha_points = []
centroid_dist = []
inertia_dist = []
for a in alpha:
    x = np.loadtxt("{}\zernike\zernike_positive\zernike_alpha_{}.dat".format(respath,a), delimiter=" ", skiprows=1)


    print("len alpha {}=".format(a), x.shape)
    X = x.transpose()
    number_points = x.shape[1]
    alpha_points.append(number_points)

    reduced = np.dot(eigenvector_subset.transpose(), X.transpose()).transpose()
    x_ax = reduced[:,0]
    y_ax = reduced[:,1]

    #########   trovo centroide e inerzia ################
    clusterer_alpha = KMeans(n_clusters=1, random_state=10)
    cluster_labels_alpha = clusterer_alpha.fit_predict(reduced)
    centroid_alpha = clusterer_alpha.cluster_centers_
    inertia_alpha = clusterer_alpha.inertia_
    print("centroid for alpha={}".format(a), centroid_alpha)
    print("inertia for alpha={}".format(a), inertia_alpha)

    centroid_dist.append(math.sqrt((centroid_total[0,0]-centroid_alpha[0,0])**2)+(centroid_total[0,1]-centroid_alpha[0,1])**2)
    inertia_dist.append(inertia_total-inertia_alpha)



    # TROVO LA MINIMA DISTANZA
#    dist = []
#    for i in range(number_points):
#        sys.stderr.write(("\r Processing point {} out of {} for alpha = {}".format(i, number_points, a)))
#        sys.stderr.flush()
#        # punto della PCA totale piu' vicino
#        index = my_functions.find_nearest_vector2D(TOT_reduced,reduced[i])
#        # distanza tra il punto piu' vicino e quello che sto considerando
#        dist.append(math.sqrt((TOT_reduced[index,0]-reduced[i,0])**2)+(TOT_reduced[index, 1]-reduced[i,1]**2))
#    distance.append(sum(dist)/alpha_points)

    print("Vuoi fare il plot delle PCA per lo screening totale e con alpha={}? y o n?".format(a))
    o = input()
    if o == "y":
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title("PCA: 2D projection of the Zernike coefficients")
        ax1.set_xlabel('Projection on eigenvector 1 (nm)')
        ax1.set_ylabel('Projection on eigenvector 2 (nm)')
        sc = plt.scatter(tot_x_ax, tot_y_ax)
        sc = plt.scatter(x_ax, y_ax, alpha=.3)
        plt.legend(["All {} points of the surface".format(tot_points), "{} points found with alpha={}".format(number_points,a)])
        plt.show()

y = []
x = []
for i in range(len(alpha)):
    y.append(centroid_dist[i]*inertia_dist[i]*alpha_points[i])
    x.append(alpha[i])

cost_func = np.column_stack([x,y])
np.savetxt("{}\cost_function.txt".format(respath),cost_func)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Cost function")
ax1.set_xlabel('$\\alpha$ value')
ax1.set_ylabel("Centroid distance $\\times$ inertia difference $\\times$ number of point")
for i in range(len(alpha)):
    plt.plot(x[i], y[i], 'bo')
    plt.text(x[i] * (1 + 0.01), y[i] * (1 + 0.01), "{} points".format(alpha_points[i]))
leg = ax1.legend()
plt.show()