from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
from scipy.cluster.vq import vq

range_n_clusters = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]


print("Di che parte di RRM2 si tratta: W per intera, 208 per taglio a 208, 219 per taglio a 219")
proteina=input()

if proteina == "W":
    loc = "../RRM2/"
elif proteina == "208":
    loc = "../cutfrom208_RRM2/"
elif proteina == "219":
    loc = "../cut220_RRM2/"

file = "%s2dproj.xvg" % loc
x1 = []
y = []
with open(file, 'r+') as f:
    for line in f:
        cols = line.split()

        if len(cols) == 2:
            x1.append(float(cols[0]))
            y.append(float(cols[1]))
d = {'x': x1,'y': y}
X = pd.DataFrame(d)



silhouette_avg1 = []
clusterer1 = []
cluster_labels1 = []
for k in range_n_clusters:

    clusterer1.append(KMeans(n_clusters=k, random_state=10))
    cluster_labels1.append(clusterer1[k-2].fit_predict(X))

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg1.append(silhouette_score(X, cluster_labels1[k-2]))
    print("For n_clusters =", k,
          "The average silhouette_score is :", silhouette_avg1[k-2])


print("Quanti cluster vuoi usare?")
n_clusters = int(input())

fig = plt.figure()
ax1 = fig.add_subplot(111)

clusterer = clusterer1[n_clusters-2]
cluster_labels = cluster_labels1[n_clusters-2]
# Compute the silhouette scores for each sample
sample_silhouette_values = silhouette_samples(X, cluster_labels)

silhouette_avg = silhouette_avg1[k-2]

y_lower = 10
for i in range(n_clusters):
    # Aggregate the silhouette scores for samples belonging to
    # cluster i, and sort them


    ith_cluster_silhouette_values = \
        sample_silhouette_values[cluster_labels == i]

    ith_cluster_silhouette_values.sort()

    size_cluster_i = ith_cluster_silhouette_values.shape[0]
    y_upper = y_lower + size_cluster_i

    color = cm.nipy_spectral(float(i) / n_clusters)
    ax1.fill_betweenx(np.arange(y_lower, y_upper),
                      0, ith_cluster_silhouette_values,
                      facecolor=color, edgecolor=color, alpha=0.7)

    # Label the silhouette plots with their cluster numbers at the middle
    ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

    # Compute the new y_lower for next plot
    y_lower = y_upper + 10  # 10 for the 0 samples

ax1.set_title("Silhouette plot for the various clusters.")
ax1.set_xlabel("Silhouette coefficient values")
ax1.set_ylabel("Cluster label")

# The vertical line for average silhouette score of all the values
ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

ax1.set_yticks([])  # Clear the yaxis labels / ticks
ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
plt.show()

# 2nd Plot showing the actual clusters formed

fig = plt.figure()
ax2 = fig.add_subplot(111)

colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
ax2.scatter(X['x'], X['y'], marker='.', s=30, lw=0, alpha=0.7,
        c=colors, edgecolor='k')

# Labeling the clusters
centers = clusterer.cluster_centers_
print(centers)
np.savetxt("{}centers.txt".format(loc), centers)



# Draw white circles at cluster centers
ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
        c="white", alpha=1, s=200, edgecolor='k')

for i, c in enumerate(centers):
    ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
                s=50, edgecolor='k')

ax2.set_title("Visualization of the clustered data.")
ax2.set_xlabel("Projection on eigenvector 1 (nm)")
ax2.set_ylabel("Projection on eigenvector 2 (nm)")


plt.show()

    ################    PRESO DA CONFORMATION.PY #######################3
c1 = []
c2 = []
with open('%scenters.txt' % loc, 'r+') as f:
    for line in f:
        cols = line.split()

        if len(cols) == 2:
            c1.append(float(cols[0]))
            c2.append(float(cols[1]))
c1_c2 = {'c1': c1,'c2': c2}
pca_coor = pd.DataFrame(c1_c2)

closest, distances = vq(pca_coor, X)
index = X.index[closest].tolist()

print(index)

i = 1
for j in index:
    with open('{}frame{}_index.ndx'.format(loc, i), 'w+') as f:
        f.write("[frames]"+"\n")
        f.write(str(j))
    i = i+1
