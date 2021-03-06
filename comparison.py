import numpy as np
import matplotlib.pyplot as plt
import my_functions
from matplotlib.ticker import FormatStrFormatter
import shutil
import statistics

    ################### PARAMETERS  ########
with open('configuration.txt') as f:
    for line in f:
        exec(line)

screening = "index_continuous_distribution_exp"

#number_crowns = 5

number_division = 20

alpha = [30., 60., 90., 120., 150., 180., 210., 240., 270., 300., 330., 360., 390.,
         420., 450., 480., 510., 540., 570., 600., 630., 660.,
         690., 720., 750., 780., 810.]




alpha_plot = []
total_plot = 1 #1 se vuoi vederlo
#fragment = 208
#cluster = 1
#R_zernike = 6
Rs_select = 4

    ############### END PARAMETERS  ################
respath = "..\\{}\cluster{}\R_zernike_{}\R_s_{}".format(fragment, cluster, R_zernike, Rs_select)

total = np.loadtxt("{}\zernike\zernike_positive\zernike_total.dat".format(respath), delimiter=" ", skiprows=1)

print("len total=", total.shape)
tot_points = total.shape[1]

### PCA of the Zernike coefficients for all patches
TOT_reduced, total_eigenvector_subset, total_eigenvalues_subset, total_eigenvalues = my_functions.PCA(total)
tot_elips = total_eigenvalues_subset[0]/total_eigenvalues_subset[1]
tot_x_ax = TOT_reduced[:, 0]
tot_y_ax = TOT_reduced[:, 1]

######  EVR #######
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("Explained Variance Ratio")
ax1.set_xlabel('Eigenvector index')
ax1.set_ylabel('Explained Variance Ratio')
evr = np.divide(total_eigenvalues, sum(total_eigenvalues))
x = np.arange(1,len(total_eigenvalues)+1)
ax1.scatter(x[2:], evr[2:], s=5, c='r', label='')
ax1.scatter(x[0:2], evr[0:2], s=50, c='y', marker='*', label='First two eigenvalues\' EVRs')
leg = ax1.legend()
plt.show()

z_total = np.loadtxt("{}\COSINE_step_{}_Rs_{}.txt".format(respath, step, Rs_select), delimiter=" ")
print(np.shape(z_total))
### PCA of the Zernike coefficients for all patches in polar coordinates
                        ### GRID ###
x_grid = np.linspace(min(tot_x_ax), max(tot_x_ax), num=number_division)
y_grid = np.linspace(min(tot_y_ax), max(tot_y_ax), num=number_division)
cell_w = x_grid[1] - x_grid[0]
cell_h = y_grid[1] - y_grid[0]
                        ### END GRID###

if total_plot == 1:
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_title("2D projection of the Zernike coefficients for each patch")
    ax1.set_xlabel('Projection on eigenvector 1')
    ax1.set_ylabel( 'Projection on eigenvector 1')
    cm = plt.cm.get_cmap('RdYlBu')
    sc = plt.scatter(tot_x_ax, tot_y_ax, c=z_total, cmap=cm)
    plt.colorbar(sc, format='%.1f', label= "Mean cosine value of the patch")
    #plt.show()


elipse_tot = []
variance_tot = []
free_cells_tot = []

reward = []
alpha_points = []
print("INIZIALE+",sum(z_total)/tot_points)
points_cosine = []
for a in alpha:

    with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
        index_possible_points = [int(float(x)) for x in f.read().split()]
    x = total[:,index_possible_points]
    z = z_total[index_possible_points]
    print("len alpha {}=".format(a), x.shape)

     ### PROIEZIONE SULLE DUE PC DEI COEFFICIENTI DI ZERNIKE PER OGNI PATCH
    X = x.transpose()
    number_points = x.shape[1]
    alpha_points.append(number_points)
    reduced = np.dot(total_eigenvector_subset.transpose(), X.transpose()).transpose()
    x_ax = reduced[:,0]
    y_ax = reduced[:,1]

    if a in alpha_plot:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title("2D projection of the Zernike coefficients for each patch")
        ax1.set_xlabel('Projection on eigenvector 1')
        ax1.set_ylabel('Projection on eigenvector 1')

        sc = plt.scatter(tot_x_ax, tot_y_ax, c='black')

        cm = plt.cm.get_cmap('RdYlBu')
        sc = plt.scatter(x_ax, y_ax,  c=z, cmap=cm)
        plt.colorbar(sc, format='%.1f', label="Mean cosine value of the patches selected with $\\alpha$={}".format(a))

        plt.legend(["All {} points of the surface".format(tot_points),
                    "{} points found with alpha={}".format(number_points, a)])
        ax1.xaxis.set_ticks(x_grid)
        ax1.yaxis.set_ticks(y_grid)
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%1.f'))
        ax1.grid(True)
        plt.tight_layout()
        #plt.show()

    _, _, eigenvalues_subset, _ = my_functions.PCA(x)
    alpha_elips = eigenvalues_subset[0] / eigenvalues_subset[1]



    count = np.zeros((len(x_grid), len(y_grid)))
    count_cell_total=np.zeros(len(x_grid)*len(y_grid))
    count_cell_alpha = np.zeros(len(x_grid)*len(y_grid))
    common_points = []
    alpha_points_for_variance = []



    k = 0
    occupied_cells = 0
    for i in range(len(x_grid)):
        for j in range(len(y_grid)):
            index_cell_total = np.where(
                (tot_x_ax >= x_grid[0] + i * cell_w) & (tot_x_ax < x_grid[0] + (i + 1) * cell_w) & \
                (tot_y_ax>= y_grid[0] + j * cell_h) & (tot_y_ax < y_grid[0] + (j + 1) * cell_h))

            index_cell_alpha = np.where(
                (x_ax >= x_grid[0] + i * cell_w) & (x_ax < x_grid[0] + (i + 1) * cell_w) & \
                (y_ax>= y_grid[0] + j * cell_h) & (y_ax < y_grid[0] + (j + 1) * cell_h))

            #punti blu nella cella su tutti i punti blu
            count_cell_total[k] = len(index_cell_total[0])/tot_points

            #punti arancioni nella cella su tutti i punti arancioni
            count_cell_alpha[k] = len(index_cell_alpha[0])/number_points
            if count_cell_alpha[k] != 0:
                occupied_cells +=1
                common_points.append(abs(count_cell_total[k]-count_cell_alpha[k]))
                alpha_points_for_variance.append(len(index_cell_alpha[0]))
            k += 1


    #reward.append(abs(tot_elips-alpha_elips)*(sum(z)/number_points)*sum(common_points)*((len(x_grid)*(len(y_grid)))-occupied_cells)/number_points)

#    reward.append(abs(tot_elips-alpha_elips)*(sum(z)/number_points)*(statistics.variance(alpha_points_for_variance)/number_points)*((len(x_grid)*len(y_grid))-occupied_cells)/(len(x_grid)*(len(y_grid))))
    reward.append((sum(z)/number_points)*(statistics.variance(alpha_points_for_variance)/number_points)*((len(x_grid)*len(y_grid))-occupied_cells)/(len(x_grid)*(len(y_grid))))
    points_cosine.append(sum(z)/number_points)

    elipse_tot.append(abs(tot_elips-alpha_elips))
    variance_tot.append(statistics.variance(alpha_points_for_variance)/number_points)
    free_cells_tot.append(((len(x_grid)*len(y_grid))-occupied_cells)/(len(x_grid)*(len(y_grid))))

#np.savetxt("{}\\reward_function_continuous_distribution_exp.txt".format(respath),reward)

difference = []
alpha_index = []
scaled_reward = []
for i in range(1,len(alpha_points)):
    difference.append(abs(alpha_points[i]-alpha_points[i-1]))
    alpha_index.append(alpha[i])
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_ylabel("Number of discarded points")
ax1.set_xlabel('$\\alpha$ value')
ax1.set_xticks(alpha_index)
plt.plot(alpha_index, difference)
leg = ax1.legend()


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel('$\\alpha$ value')
ax1.set_ylabel("Components")
plt.plot(alpha, variance_tot, label = 'variance')
#plt.plot(alpha, elipse_tot, label = 'elipse')
plt.plot(alpha, free_cells_tot, label = 'not occupied')
plt.plot(alpha, points_cosine, label = 'mean cosine')
leg = ax1.legend()




######### Qua trovo i minimi locali ##############
#alpha_min_loc = []
#loss_min_loc = []
#for i in range(1,len(alpha)-1):
#    if (reward[i-1] > reward[i]) & (reward[i]<reward[i+1]):
#        alpha_min_loc.append(alpha[i])
#        loss_min_loc.append(reward[i])

loss_min = np.min(reward)
alpha_min = np.where(loss_min == reward)[0][0]
mask = np.array(reward) == loss_min
color = np.where(mask, 'red', 'blue')

fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax1.set_title("Cost function")
ax1.set_xlabel('$\\alpha$ value')
ax1.set_ylabel("Loss function ")
ax1.set_xticks(alpha)
#ax1.scatter(alpha, reward, color=color)
ax1.plot(alpha, reward)
#ax1.plot(alpha[alpha_min], reward[alpha_min], '*')
#ax1.legend(["The minimum of the Loss function is at $\\alpha$={}".format(alpha[alpha_min])])


fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax1.set_title("Cost function")
ax1.set_xlabel('$\\alpha$ value')
ax1.set_ylabel("Mean cosine value of the selected patches")
#ax1.set_xticks(np.arange(0, 1.05, .05))
for i in range(len(alpha)):
    plt.plot(alpha[i], points_cosine[i], "*")
leg = ax1.legend()
#plt.show()

plt.show()




print("Qual'e' il miglior valore per alpha?")
best_alpha = input()

with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, best_alpha)) as f:
    index_possible_points = [int(float(x)) for x in f.read().split()]
x = total[:, index_possible_points]
z = z_total[index_possible_points]

X = x.transpose()
number_points = x.shape[1]
alpha_points.append(number_points)
reduced = np.dot(total_eigenvector_subset.transpose(), X.transpose()).transpose()
x_ax = reduced[:, 0]
y_ax = reduced[:, 1]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("2D projection of the Zernike coefficients for each patch")
ax1.set_xlabel('Projection on eigenvector 1')
ax1.set_ylabel('Projection on eigenvector 1')

sc = plt.scatter(tot_x_ax, tot_y_ax, c='black')

cm = plt.cm.get_cmap('RdYlBu')
sc = plt.scatter(x_ax, y_ax, c=z, cmap=cm)
plt.colorbar(sc, format='%.1f', label="Mean cosine value of the patches selected with $\\alpha$={}".format(best_alpha))

plt.legend(["All {} points of the surface".format(tot_points),
            "{} points found with alpha={}".format(number_points, best_alpha)])
ax1.xaxis.set_ticks(x_grid)
ax1.yaxis.set_ticks(y_grid)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%1.f'))
ax1.grid(True)
plt.tight_layout()

plt.show()


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
