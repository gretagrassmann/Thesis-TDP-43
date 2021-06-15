import numpy as np
import matplotlib.pyplot as plt
import my_functions
from matplotlib.ticker import FormatStrFormatter
import shutil
import statistics
import random
from random import randrange

    ################### PARAMETERS  ########
with open('configuration.txt') as f:
    for line in f:
        exec(line)

screening = "index_continuous_distribution_exp"

nr = 10

alpha = [1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19.
         , 20., 21., 22., 23., 24., 25., 26., 27., 28., 29.
         , 30., 31., 32., 33., 34., 35., 36., 37., 38., 39.
         ,40., 41., 42., 43., 44., 45., 46., 47., 48., 49.
         ,50., 51., 52., 53., 54., 55., 56., 57., 58., 59.
         ,60., 61., 62., 63., 64., 65., 66., 67., 68., 69.
         ,70., 71., 72., 73., 74., 75., 76., 77., 78., 79.
         ,80., 81., 82., 83., 84., 85., 86., 87., 88., 89.
         ,90., 91., 92., 93., 94., 95., 96.
     #  , 100., 120., 150., 180., 210.
     #   , 240., 270., 300., 330., 360., 390.,
     #    420., 450., 480., 510.
     #         , 540., 570., 600., 630., 660.,
     #    690., 720., 750., 780., 810.
     #    ,840., 870., 900., 930., 960., 990., 1000., 2000., 3000.
         ]

alpha_plot = []
#fragment = 208
#cluster = 1
#R_zernike = 6
Rs_select = 4

    ############### END PARAMETERS  ################
respath = "..\\{}\cluster{}\R_zernike_{}\R_s_{}".format(fragment, cluster, R_zernike, Rs_select)



print("Vuoi vedere i grafici gia' ottenuti? y o n?")
o = input()
if o == 'y':
    zernike_missed_var =  np.loadtxt("{}\zernike_missed_var.txt".format(respath))
    zernike_missed_var_random = np.loadtxt("{}\zernike_missed_var_random.txt".format(respath))
    zernike_var = np.loadtxt("{}\zernike_var.txt".format(respath))
    zernike_var_random = np.loadtxt("{}\zernike_var_random.txt".format(respath))
    alpha_plot = np.loadtxt("{}\zm_alpha.txt".format(respath))

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Zernike coefficients' variance")
    plt.plot(alpha_plot, zernike_var, label='Selected by sampling')
    plt.plot(alpha_plot, zernike_missed_var, label='Missed by sampling')
    plt.plot(alpha_plot, zernike_var_random, label='Randomly selected')
    plt.plot(alpha_plot, zernike_missed_var_random, label='Randomly missed')
    leg = ax1.legend()

    plt.show()



total = np.loadtxt("{}\zernike\zernike_positive\zernike_total.dat".format(respath), delimiter=" ", skiprows=1)
print("len total=", total.shape)
tot_points = total.shape[1]
z_total = np.loadtxt("{}\COSINE_step_{}_Rs_{}.txt".format(respath, step, Rs_select), delimiter=" ")

zernike_missed_var = []
zernike_missed_var_random = []
zernike_var = []
zernike_var_random = []

print("Vuoi vedere il numero di punti in funzione di alpha? y o n?")
o = input()
if o == "y":
    alpha_points_plot = []
    for a in alpha:
        with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
            index_possible_points = list(set([int(float(x)) for x in f.read().split()]))
            alpha_points_plot.append(len(index_possible_points))

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('$\\alpha$ value')
    ax1.set_ylabel("Number of points")
    plt.plot(alpha, alpha_points_plot)
    plt.show()

alpha_points = []
for a in alpha:
    with open("{}\{}\index_possible_area_R_s_{}_alpha_{}_step_1.txt".format(respath, screening, Rs_select, a)) as f:
        index_possible_points = list(set([int(float(x)) for x in f.read().split()]))

    x = total[:, index_possible_points]
    x_missed = total[:, [i for i in range(tot_points) if (i not in index_possible_points)]]

    print("len alpha {}=".format(a), x.shape)
    number_points = x.shape[1]
    alpha_points.append(number_points)

    zernike_var_sum = 0
    for i in range(len(x)):
        zernike_var_sum += np.var(x[i])
    zv = zernike_var_sum
    zernike_var.append(zv)

    zernike_var_miss_sum = 0
    for i in range(len(x_missed)):
        zernike_var_miss_sum += np.var(x_missed[i])
    zvm = zernike_var_miss_sum
    zernike_missed_var.append(zvm)


        ############    RANDOM  #####################
    zernike_var_random_nr = 0
    zernike_missed_var_random_nr = 0
    for i in range(nr):
        random_index = random.sample(list(np.arange(0, tot_points)), number_points)
        x_random = total[:, random_index]
        x_missed_random = total[:, [i for i in range(tot_points) if (i not in random_index)]]


        zernike_var_rand_sum = 0
        for i in range(len(x_random)):
            zernike_var_rand_sum += np.var(x_random[i])
        zernike_var_random_nr += zernike_var_rand_sum

        zernike_var_rand_miss_sum = 0
        for i in range(len(x_missed_random)):
            zernike_var_rand_miss_sum += np.var(x_missed_random[i])
        zernike_missed_var_random_nr += zernike_var_rand_miss_sum

    zvr = zernike_var_random_nr/nr
    zernike_var_random.append(zvr)

    zmvr = zernike_missed_var_random_nr/nr
    zernike_missed_var_random.append(zmvr)


np.savetxt("{}\zernike_missed_var.txt".format(respath),zernike_missed_var)
np.savetxt("{}\zernike_missed_var_random.txt".format(respath),zernike_missed_var_random)
np.savetxt("{}\zernike_var.txt".format(respath),zernike_var)
np.savetxt("{}\zernike_var_random.txt".format(respath),zernike_var_random)
np.savetxt("{}\zm_alpha.txt".format(respath),alpha)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel('$\\alpha$ value')
ax1.set_ylabel("Zernike coefficients' variance")
plt.plot(alpha, zernike_var, label='Selected by sampling')
plt.plot(alpha, zernike_missed_var, label='Missed by sampling')
plt.plot(alpha, zernike_var_random, label='Randomly selected')
plt.plot(alpha, zernike_missed_var_random, label='Randomly missed')
leg = ax1.legend()

plt.show()

