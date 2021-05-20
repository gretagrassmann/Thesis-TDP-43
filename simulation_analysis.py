import numpy as np
import matplotlib.pyplot as plt
import statistics
import math
import matplotlib.ticker as mtick

x, y = [], []

print("Di che parte di RRM2 si tratta: W per intera, 208 per taglio a 208, 219 per taglio a 219")
proteina=input()

print("Di cosa fare il grafico: E per Energy minimization, T per Temperature, P per Pressure, R per RMSD, D per Density, G for gyratio,"
      "PCA per PCA, EVR per EVR.")
grafico=input()

if proteina == "W":
    loc = "../RRM2/"
elif proteina == "208":
    loc = "../cutfrom208_RRM2/"
elif proteina == "219":
    loc = "../cut220_RRM2/"


if grafico == "E":
    file = "%spotential_11.xvg" % loc
elif grafico == "T":
    file = "%stemperature.xvg" % loc
elif grafico == "P":
    file = "%spressure.xvg" % loc #se stai studiando W ricorda di cambiare a pressure2.xvg!!!!
elif grafico == "D":
    file = "%sdensity.xvg" % loc
elif grafico == "R":
    file = "%srmsd.xvg" % loc #equilibrated
    file2 = "%srmsd_xtal.xvg" % loc #crystal
elif grafico == "G":
    file = "%sgyrate.xvg" % loc
elif grafico == "PCA":
    file = "%s2dproj.xvg" % loc
elif grafico == "EVR":
    file = "%seigenval.xvg" % loc
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
    ax1.set_title("RMSD")
    ax1.set_xlabel('Time (ns)')
    ax1.set_ylabel('RMSD (nm)')
    ax1.plot(x, y, c='b', label='Equilibrated') #alpha=0.3,
    ax1.plot(x2,y2, c='r', label='Crystal')
    X_detail = list(range(400,421))
    Y_detail = y[400:421]
    Y2_detail = y2[400:421]
    # location for the zoomed portion
    if proteina == "W":
        sub_axes = plt.axes([0.55,.2,.3,.3]) #rect = [left, bottom, width, height]
    else :
        sub_axes = plt.axes([0.25, .17, .25, .25])  # rect = [left, bottom, width, height]



    # plot the zoomed portion
    sub_axes.plot(X_detail, Y_detail, c='b')
    sub_axes.plot(X_detail, Y2_detail, c='r')

elif grafico == "G":
    ax1.set_title("Radius of total gyration")
    ax1.set_xlabel('Time (ns)')
    ax1.set_ylabel('$R_g$ (nm)')
    ax1.plot(x, y, c='r', label='')

elif grafico == "PCA":
    cm = plt.cm.get_cmap('RdYlBu')

    ax1.set_title("PCA: 2D projection of the trajectory")
    ax1.set_xlabel('Projection on eigenvector 1 (nm)')
    ax1.set_ylabel('Projection on eigenvector 2 (nm)')

    z = list(range(len(x))) #len(x)=155606 for W
    sc = plt.scatter(x, y, c=z, vmin=0, vmax=len(x), s=.3, cmap=cm)

    plt.colorbar(sc, format='%.1e')

elif grafico == "EVR":
    ax1.set_title("Explained Variance Ratio")
    ax1.set_xlabel('Eigenvector index')
    ax1.set_ylabel('Explained Variance Ratio')

    evr = np.divide(y, sum(y))
    ax1.scatter(x[2:], evr[2:], s=.5, c='r', label='')
    ax1.scatter(x[0:2], evr[0:2], s=50, c='y', marker='*', label='First two eigenvalues\' EVRs')



leg = ax1.legend()
plt.show()

if grafico == "R":
    plt.hist(y,100)
    plt.xlabel("RMSD (nm)")
    plt.title("RMSD distribution")
    plt.show()