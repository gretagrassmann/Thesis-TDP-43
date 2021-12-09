import numpy as np
import matplotlib.pyplot as plt
import statistics
import math
import matplotlib.ticker as mtick
from matplotlib.ticker import MaxNLocator
import scipy.stats as st

x, y = [], []

print("Di che parte di RRM2 si tratta: W per intera, 208 per taglio a 208, 219 per taglio a 219, 219 rmsd per la sezione")
proteina=input()

print("Di cosa fare il grafico: E per Energy minimization, T per Temperature, P per Pressure, R per RMSD, D per Density, G for gyratio,"
      "PCA per PCA, EVR per EVR, ClUSTER per cluster.")
grafico=input()

if proteina == "W":
    loc = "../RRM2/"
elif proteina == "208":
    loc = "../cutfrom208_RRM2/"
elif proteina == "219":
    loc = "../cut220_RRM2/"
elif proteina == "219 rmsd":
    loc = "../cut220_RRM2/RMSD/"


if grafico == "E":
    file = "%spotential_11.xvg" % loc
elif grafico == "T":
    file = "%stemperature.xvg" % loc
elif grafico == "P":
    if proteina == 'W':
        file = "%spressure2.xvg" % loc
    else:
        file = "%spressure.xvg" % loc
elif grafico == "D":
    file = "%sdensity.xvg" % loc
elif grafico == "R":
    file = "%srmsd.xvg" % loc #equilibrated
    file2 = "%srmsd_xtal.xvg" % loc #crystal
elif grafico == "Rh":
    file = "%srmsd.xvg" % loc  # equilibrated
elif grafico == "G":
    file = "%sgyrate.xvg" % loc
elif grafico == "PCA":
    file = "%s2dproj.xvg" % loc
elif grafico == "EVR":
    file = "%seigenval.xvg" % loc
elif grafico == 'CLUSTER':
    print("Quanti cluster vuoi usare?")
    n_clusters = int(input())

    sample_silhouette_values = np.loadtxt('%ssample_silhouette_values.txt' % loc)
    cluster_labels = np.loadtxt('%scluster_labels.txt' % loc)
    clusterer=np.loadtxt("%sclusterer.txt" % loc)
    clusterer_labels = np.loadtxt('%sclusterer_labels.txt' % loc)
    silhouette_avg = np.loadtxt('%ssilhouette_avg.txt' % loc)
    centers = np.loadtxt('%scenters.txt' % loc)
    X = np.loadtxt('%sX.txt' % loc)

elif grafico == 'R e G':
    file = "%srmsd.xvg" % loc  # equilibrated
    file2 = "%sgyrate.xvg" % loc

else:
    print("Error")

if grafico == "G":
    with open(file, "r+") as f:
        for line in f:
            cols = line.split()

            if len(cols) == 5:
                x.append(float(cols[0])*0.001) #from pico to nano
                y.append(float(cols[1]))
else:
    with open(file, "r+") as f:
        for line in f:
            cols = line.split()

            if len(cols) == 2:
                x.append(float(cols[0]))
                y.append(float(cols[1]))
    if grafico == "R":
        x2, y2 = [], []
        with open(file2, "r+") as f:
            for line in f:
                cols = line.split()
                x2.append(float(cols[0]))
                y2.append(float(cols[1]))

    if grafico == 'R e G':
        x2, y2 = [], []
        with open(file2, "r+") as f:
            for line in f:
                cols = line.split()
                x2.append(float(cols[0])*0.001)  # from pico to nano
                y2.append(float(cols[1]))


std = math.sqrt(statistics.variance(y))
mean = statistics.mean(y)
print("sigma=",std)
print("mean=", mean)

fig = plt.figure()
ax1 = fig.add_subplot(111)
if grafico == "E":
    ax1.set_title("Energy minimization")
    ax1.set_xlabel('Energy minimization step (ps)')
    ax1.set_ylabel('Potential energy (kj/mol)')
    ax1.plot(x, y, c='r', label='Potential')
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))


elif grafico == "T":
    ax1.set_title("Temperature")
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Temperature (K)')
    ax1.plot(x, y, c='r', label='')
elif grafico == "P":
    ax1.set_title("Pressure")
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Pressure (bar)')
    ax1.plot(x, y, c='r', label='')
elif grafico == "D":
    ax1.set_title("Density")
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Density (kg/$m^3$)')
    ax1.plot(x, y, c='r', label='')
elif grafico == "R":
    print('Vuoi vedere solo rispetto alla struttura equilibrata o anche quella iniziale? 1 o 2')
    o =input()
    print('Vuoi vedere i limiti di unfolding di Fragment B? y o n')
    oo = input()
    ax1.set_title("RMSD",fontsize =30)
    ax1.set_xlabel('Time (ns)',fontsize =20)
    ax1.set_ylabel('RMSD (nm)',fontsize =20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    if o == '1':
        ax1.plot(x, y, c='grey') #alpha=0.3, label='Equilibrated'
    else:
        ax1.plot(x, y, c='grey', alpha=0.3, label='Equilibrated')
        ax1.plot(x2,y2, c='r', label='Crystal')
    print('Vuoi vedere la posizione dei centroidi? y o n')
    centroid = input()
    if centroid == 'y':
        if proteina == '208':
            centroid_frame = [570494, 588955, 178412, 524193, 119199]
        else:
            centroid_frame = [514432, 672666]
        for i in centroid_frame:
            print(i,x[i])
            plt.vlines(x[i], plt.ylim()[0], plt.ylim()[1], linestyles='dashed', linewidth=2)

    if oo == 'y':
        plt.vlines(5500, plt.ylim()[0], plt.ylim()[1], linestyles='dashed', linewidth=2)
        plt.vlines(6500, plt.ylim()[0], plt.ylim()[1], linestyles='dashed', linewidth=2)
        print('MAX RMSD=', max(y))
        print(np.where(np.array(y)==max(y))[0][0])
        print('at time=', x[np.where(np.array(y)==max(y))[0][0]])
    X_detail = list(range(400,421))
    Y_detail = y[400:421]
    Y2_detail = y2[400:421]
    ax1.yaxis.set_major_locator(MaxNLocator(6))

    print('Vuoi vedere la porzione zoomata? y o n')
    o = input()
    if o == 'y':
        if proteina == "W":
            sub_axes = plt.axes([0.58,.16,.25,.2]) #rect = [left, bottom, width, height]
        if proteina == "219":
            sub_axes = plt.axes([0.2,.55,.3,.3]) #rect = [left, bottom, width, height]
        if proteina == "208":
            sub_axes = plt.axes([0.25, .17, .25, .25])  # rect = [left, bottom, width, height]
        sub_axes.plot(X_detail, Y_detail, c='b')
        sub_axes.plot(X_detail, Y2_detail, c='r')

elif grafico == "R e G":
    ax1.set_title("RMSD and $R_g$",fontsize=22)
    if proteina == '219' or proteina == '219 rmsd':
        ax1.set_ylim(.0,2.3)
    elif proteina == '208':
        ax1.set_ylim(.0,1.6)
    else:
        ax1.set_ylim(0.,1.3)
    ax1.set_xlabel('Time (ns)',fontsize=22)
    ax1.set_ylabel('RMSD (nm)',fontsize=22)
    plt.yticks(fontsize=22)
    plt.xticks(fontsize=22)
    ax1.plot(x, y, c='r') #alpha=0.3,
    if proteina == '219 rmsd':
        plt.vlines(5500, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
        plt.vlines(6500, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
    # location for the gyratio
    if proteina == "W":
       sub_axes = plt.axes([0.22, .58, .3, .25])  # rect = [left, bottom, width, height]
    if proteina == "219" or proteina == '219 rmsd':
       sub_axes = plt.axes([0.22, .58, .3, .25])  # rect = [left, bottom, width, height]
    if proteina == "208":
       sub_axes = plt.axes([0.22, .58, .3, .25])  # rect = [left, bottom, width, height]
    sub_axes.set_xlabel('Time (ns)',fontsize=22)
    sub_axes.set_ylabel('$R_g$ (nm)',fontsize=22)
    sub_axes.set_xticklabels([])
    sub_axes.set_yticklabels([])
    sub_axes.plot(x2, y2, c='r')
    if proteina == '219 rmsd':
        plt.vlines(5500, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
        plt.vlines(6500, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')


elif grafico == "G":
    ax1.set_title("Radius of total gyration",fontsize=30)
    ax1.set_xlabel('Time (ns)',fontsize=20)
    ax1.set_ylabel('$R_g$ (nm)',fontsize=20)
    ax1.plot(x, y, c='grey', label='')
    plt.vlines(5500, plt.ylim()[0], plt.ylim()[1], linestyles='dashed', linewidth=2)
    plt.vlines(6500, plt.ylim()[0], plt.ylim()[1], linestyles='dashed', linewidth=2)
    ax1.yaxis.set_major_locator(MaxNLocator(5))

    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)

elif grafico == "PCA":
    print('Vuoi vedere la densita? y o n?')
    o = input()
    if o =='y':
        # Creating bins
        x_min = np.min(x)
        x_max = np.max(x)
        y_min = np.min(y)
        y_max = np.max(y)
        x_bins = np.linspace(x_min, x_max, 100)
        y_bins = np.linspace(y_min, y_max, 100)
        fig, ax = plt.subplots()
        # Creating plot
        counts, binx, biny = np.histogram2d(x,y, bins=[x_bins,y_bins])
        counts = counts.T
        counts = np.log(counts)
        plt.imshow(counts, interpolation='nearest', origin='lower',
                   extent=[binx[0], binx[-1], biny[0], biny[-1]])

        with open('%scenters.txt' % loc, 'r+') as f:
            c1,c2 = [],[]
            for line in f:
                cols = line.split()
                if len(cols) == 2:
                    c1.append(float(cols[0]))
                    c2.append(float(cols[1]))
        ax.scatter(c1, c2, marker='o',
                    c="white", alpha=1, s=20, edgecolor='k')
        #for i, c1_i in enumerate(c1):
        #    ax.scatter(c1[i], c2[i], marker='$A%d$' % (i+1), alpha=1,
        #                s=20, edgecolor='k')
        #cfset = ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax], aspect='auto')
        ax.xaxis.set_tick_params(labelsize=20)
        ax.yaxis.set_tick_params(labelsize=20)
        ax.set_xlabel('PC 1', fontsize=30, linespacing=3.1)
        ax.set_ylabel('PC 2', fontsize=30, linespacing=3.1)
        ax.set_title('PCA density', fontsize=30, linespacing=3.1)
        #cbar = plt.colorbar(cfset, orientation="vertical")
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax1.xaxis.set_major_locator(MaxNLocator(4))
        #cbar.ax.tick_params(labelsize=20)
        plt.tight_layout()
        plt.show()

    cm = plt.cm.get_cmap('RdYlBu')
    ax1.set_title("PCA", fontsize=22)
    ax1.set_xlabel('Projection on PC 1', fontsize=25)
    ax1.set_ylabel('Projection on PC 2', fontsize=25)
    plt.yticks(fontsize=25)
    plt.xticks(fontsize=25)
    z = list(range(len(x)))
    sc = plt.scatter(x, y, c=z, vmin=0, vmax=len(x), s=.3, cmap=cm)
    print('Vuoi vedere i centroid?y o n?')
    o = input()
    if o == 'y':
        with open('%scenters.txt' % loc, 'r+') as f:
            c1,c2 = [],[]
            for line in f:
                cols = line.split()
                if len(cols) == 2:
                    c1.append(float(cols[0]))
                    c2.append(float(cols[1]))
        ax1.scatter(c1, c2, marker='o',
                    c="white", alpha=1, s=700, edgecolor='k')
        for i, c1_i in enumerate(c1):
            ax1.scatter(c1[i], c2[i], marker='B$%d$' % (i+1), alpha=1,
                        s=250, edgecolor='k')

    cbar = plt.colorbar(sc)
    ax1.yaxis.set_major_locator(MaxNLocator(5))
    cbar.ax.tick_params(labelsize=22)
    cbar.set_ticks([0,250000,500000,750000,1000000])
    cbar.set_ticklabels(['{}e+00'.format(0.0),'{}e+03'.format(2.5),'{}e+03'.format(5.0),'{}e+03'.format(7.5),'{}e+04'.format(1.0)])
    cbar.set_label('Time (ns)',fontsize=25)


elif grafico == "EVR":
    ax1.set_title("Explained Variance Ratio")
    ax1.set_xlabel('Eigenvector index')
    ax1.set_ylabel('Explained Variance Ratio')

    evr = np.divide(y, sum(y))
    ax1.scatter(x[2:], evr[2:], s=.5, c='r', label='')
    ax1.scatter(x[0:2], evr[0:2], s=50, c='y', marker='*', label='First two eigenvalues\' EVRs')



leg = ax1.legend()
plt.tight_layout()
plt.show()

if grafico == "Rh":
    fig = plt.figure()
    ax = fig.add_subplot(111)
    print('Vuoi vedere la posizione dei centroidi? y o n')
    centroid = input()
    plt.hist(y, 1000,density=True,color='lightgray')
    if centroid == 'y':
        if proteina == '208':
            centroid_frame = [570494, 588955, 178412, 524193, 119199]
            centroid_color = ['dimgrey','b','cyan','lawngreen','orange' ]
        else:
            centroid_frame = [514432, 672666]
            centroid_color = ['firebrick','lawngreen']
        j=0
        for i in centroid_frame:
            print(y[i])
            plt.vlines(y[i],ymin=0, ymax=6, linestyles='dashed', linewidth=2,colors=centroid_color[j],label="B{}".format(j+1))
            j+=1

    plt.xlabel("RMSD (nm)",fontsize=25)
    plt.title("RMSD distribution",fontsize=25)
    plt.yticks(fontsize=25)
    plt.xticks(fontsize=25)
    plt.legend(fontsize=25)
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.xaxis.set_major_locator(MaxNLocator(5))
    plt.tight_layout()
    plt.show()



if grafico == 'CLUSTER':
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
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

