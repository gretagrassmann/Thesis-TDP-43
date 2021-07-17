import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import DataFrame
import random
import matplotlib.ticker as mtick

from scipy import interpolate
from matplotlib.ticker import FuncFormatter
import ZernikeFunc as ZF
import SurfaceFunc as SF
import seaborn as sns
import my_functions
import shutil
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
        ######################## PARAMETERS #########
with open('configuration.txt') as f:
    for line in f:
        exec(line)
#alpha = [-9., -8., -7., -6., -5., -4., -3., -2., -1., 1., 2., 3., 4., 5., 6., 7., 8., 9.
#         , 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29.
#         , 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49.
#         , 50., 51., 52., 53., 54., 55., 56., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69.
#         ,  70., 71., 72., 73., 74., 75., 76., 77., 78., 79., 80., 81., 82., 83., 84., 85., 86., 87., 88., 89.
#         , 90., 91., 92., 93., 94., 95., 96.]
screening = "index_continuous_distribution_exp_FAST"

#R_zernike = 6
Rs_select = 4
fragment = 208
cluster = 1
#step = 1
verso = float(1)
#Npoint = int(1)

    ###################################   END PARAMETERS    ###########

pdb_file = "..\\{}\cluster{}.dms".format(fragment, cluster)
respath = "..\\{}\cluster{}\R_zernike_{}\R_s_{}".format(fragment, cluster, R_zernike, Rs_select)
if (verso == 1):
    ver = "up"
elif (verso == -1):
    ver = "down"
ZOrder = 20  # the Zernike expansion order..
Npixel = 25  # the plane edge in pixels...
Daa = 3.  # the maximum distance between antigen and antibody binding site atoms..
Dpp = 1.  # the distance between points of the same patch (needed to remove islands)

surf_ = pd.read_csv(pdb_file)
lag = len(surf_["x"])
surf = np.zeros((lag, 6))
surf[:, :] = surf_[["x", "y", "z", "Nx", "Ny", "Nz"]]
surf_obj = SF.Surface(surf[:, :], patch_num=0, r0=R_zernike, theta_max=45)
surf_obj_scan = SF.Surface(surf[:, :], patch_num=0, r0=Rs_select, theta_max=45)
ltmp = np.shape(surf)[0]

print('Vuoi fare il nuovo sampling EDOTTIA?')
o = input()
if o == 'y':
    with open("{}\COSINE_step_{}_Rs_{}_FAST.txt".format(respath, step, Rs_select)) as f:
        mean_cos = [float(x) for x in f.read().split()]
    alpha1 = [.1, .2, .4, .6, .8, 1.]
    beta = [0., .2, .4, .6, .8, 1.]
    gamma = [0., 2., 4., 6., 8., 10.]
    delta = [0., 2., 4., 6., 8., 10.]

    for a in alpha1:
        for b in beta:
            for g in gamma:
                for d in delta:
                    my_functions.NewDistribution(a,b,g,d,ltmp, mean_cos, surf, surf_obj_scan, Rs_select, step, respath)

print('Vuoi vedere i grafici per il sampling EDOTTIA?')
o = input()
if o == 'y':
    alpha1 = [.1, .4, 1.]
    beta = [0.,.4, 1.]
    gamma = [0., 2., 4., 6., 8., 10.]
    delta = [0., 2., 4., 6., 8., 10.]

    for a in alpha1:
        for b in beta:
            i = 0
            fig = plt.figure(figsize=(12,18))
            ax1 = fig.add_subplot()
            ax1.set_xlabel('$\\delta$ value',fontsize=40)
            ax1.set_ylabel("Fraction of sampled points",fontsize=40)
            ax1.set_ylim([0,1.1])
            plt.title('$\\alpha$={}, $\\beta$={}'.format(a, b),fontsize=40)
            number_points_gamma = []
            for g in gamma:
                number_points_delta = []
                for d in delta:
                    with open("{}\EDOTTIA\\a{}_b{}_g{}_d{}.txt".format(respath,a,b,g,d)) as f:
                        index_possible_points = list(set([int(float(x)) for x in f.read().split()]))
                        number_points_delta.append(round(len(index_possible_points)/lag,4))
                number_points_gamma.append(number_points_delta)

                plt.plot(delta, number_points_gamma[i],linewidth=5, label='$\gamma$={}'.format(g))
                plt.yticks(fontsize=40)
                plt.xticks(fontsize=40)
                plt.legend(fontsize=40)
                plt.tight_layout
                i += 1
            plt.show()

######################## PARAMETRI PER I PROSSIMI DUE STEP
frac_test = 200
PLOT = 1
alpha1 = [.1, .4, 1.]
beta = [0., .2, .4, .6, .8, 1.]
gamma = [0., 2., 4., 6., 8., 10.]
delta = [0., 2., 4., 6., 8., 10.]
#total_sampling = []
#total_random = []
###################    GRAPHS OF THE MEAN(TOT VS R)-MEAN(TOT VS S) AS A FUNCTION OF DELTA, AND OF THE LOWEST ONES LABELD WITH THE MEAN(TOT VS S)
print('Vuoi vedere tutti i grafici del mondo?')
o = input()
if o == 'y':
    file = "{}\EDOTTIA\\difference_{}_{}_over_{}points.csv".format(respath,fragment,cluster,int(lag/frac_test))
    results = pd.DataFrame()
    results['alpha'] = pd.read_csv(file,usecols=['alpha'],squeeze=True)
    results['beta'] = pd.read_csv(file,usecols=['beta'],squeeze=True)
    results['gamma'] = pd.read_csv(file,usecols=['gamma'],squeeze=True)
    results['delta'] = pd.read_csv(file,usecols=['delta'],squeeze=True)
    results['mean(tot vs s)'] = pd.read_csv(file,usecols=['m(tot vs s)'],squeeze=True)
    results['mean(tot vs r)'] = pd.read_csv(file,usecols=['m(tot vs r)'],squeeze=True)
    results['selected points'] = pd.read_csv(file,usecols=['selected points'],squeeze=True)
    results['difference'] = results['mean(tot vs r)']-results['mean(tot vs s)']
    results['loss']=results['mean(tot vs s)']*results['difference']/results['selected points']


    # fissato alpha, provo per un paio di beta a fare un grafico 3D della Loss=diff/selected points rispetto a gamma e delta

    a_fix = 1.
    a_results = pd.DataFrame()
    a_results = results[results.alpha.isin([a_fix])]
    max_d = []
    for b in beta:
        a_b_results = pd.DataFrame()
        a_b_results = a_results[a_results.beta.isin([b])]
        g_diff = []
        for g in gamma:
            a_b_g_results = pd.DataFrame()
            a_b_g_results = a_b_results[a_b_results.gamma.isin([g])]
            d_diff = []
            for d in delta:
                a_b_g_d_results = pd.DataFrame
                a_b_g_d_results = a_b_g_results[a_b_g_results.delta.isin([d])]
                d_diff.append(a_b_g_d_results['difference'])
            g_diff.append(d_diff)

        X_rough = delta
        Y_rough = gamma
        X_rough, Y_rough = np.meshgrid(X_rough,Y_rough)
        Z_rough = np.array(g_diff).reshape(6,6)

        X, Y = np.mgrid[0:10:.5, 0:10:.5]
        tck = interpolate.bisplrep(X_rough, Y_rough, Z_rough, s=0)
        Z = interpolate.bisplev(X[:, 0], Y[0, :], tck)
        fig = plt.figure(figsize=(12, 12))
        ax = fig.gca(projection='3d')
        ax.plot_surface(X, Y, Z, cmap='summer', rstride=1, cstride=1, alpha=None, antialiased=True)
        # Creating color map
        #my_cmap = plt.get_cmap('hot')
        #fig = plt.figure(figsize=(14, 9))
        #ax = plt.axes(projection='3d')
        #surf = ax.plot_surface(X, Y, Z, cmap=my_cmap, edgecolor='none')
        cset = ax.contour(X, Y, Z, zdir='z', offset=-0.1, cmap='summer')
        cset = ax.contour(X, Y, Z, zdir='x', offset=-2, cmap='summer')
        cset = ax.contour(X, Y, Z, zdir='y', offset=12, cmap='summer')
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.zaxis.set_tick_params(labelsize=20)
        ax.xaxis.set_tick_params(labelsize=20)
        ax.yaxis.set_tick_params(labelsize=20)

        ax.set_xlabel('$\\delta$', fontsize=30, linespacing=3.1)
        ax.set_ylabel('$\\gamma$', fontsize=30, linespacing=3.1)
        ax.set_title('$d$ for $\\alpha$={} and $\\beta$={}'.format(a_fix,b), fontsize=30, linespacing=3.1)

        # calc index of max Z value
        xmax, ymax = np.unravel_index(np.argmax(Z_rough), Z_rough.shape)
        # min max points in 3D space (x,y,z)
        max_d.append([X_rough[xmax, ymax], Y_rough[xmax, ymax], Z_rough.max(),b])
    plt.show()

    print(np.array(max_d))

    #d_max = list(np.array(max_d)[:,0])
    #g_max = list(np.array(max_d)[:,1])
    #b_max = list(np.array(max_d)[:,3])
    #d_max, g_max = np.meshgrid(d_max, g_max)
    #dist_max = np.array(max_d)[:,2].reshape(6, 6)



        ############## MIGLIORI 20 SAMPLING x=Z(tot-r)-Z(tot-s) y=#points, label=parameters ##############
            # GUARDA I PARAMETRI MIGLIORI E FACCI IL BOX PLOT
    best = results.nlargest(20, columns=['difference'])
    best.reset_index(drop=True, inplace=True)
    x = best['difference']
    y = best['selected points']
    color = best['gamma']
    scale = best['delta']*10
    annotations = []
    for i in range(best.shape[0]):
        annotations.append(r'Dist= %.2f, with $\alpha$=%.1f,$\beta$=%.1f,$\gamma$=%.1f,$\delta$=%.1f' % (best.loc[i,'mean(tot vs s)'], best.loc[i,'alpha'], best.loc[i,'beta'], best.loc[i,'gamma'], best.loc[i,'delta']))

    plt.figure(figsize=(8, 6))
    plt.xlabel("Difference between the random and sampling mean distance from the total surface")
    plt.ylabel("Number of sampled points")
    plt.title("Best performing samplings", fontsize=15)
    cm = plt.cm.get_cmap('RdYlBu')
    #for i, label in enumerate(annotations):
        #plt.scatter(x[i], y[i], s=100, label=annotations[i]))

    sc = plt.scatter(x,y,c=color, s=scale,cmap=cm)


    #sc = plt.scatter(x, y, c=color, label=annotations, s=scale, cmap=cm)
    plt.colorbar(sc, format='%.1f',
                 label="$\\gamma$")
    sc.set_sizes(scale)

    plt.legend()
    plt.show()

        ########## Z(tot-r)-Z(tot-s) AS A FUNCTION OF D, FOR FIXED A, B AND ALL THE G IN A GRAPH
    for a in alpha1:
        for b in beta:
            i = 0
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.set_xlabel('$\\delta$ value')
            ax1.set_ylabel("Difference between the random and sampling mean distance from the total surface")
            plt.title('$\\alpha$={}, $\\beta$={}'.format(a, b))
            number_points_gamma = []
            for g in gamma:
                diff = []
                for d in delta:
                    diff_p = results.loc[(results['alpha'] == a)&(results['beta'] == b)&(results['gamma']==g),'difference']
                    print('diff p', np.array(diff_p).shape)
                    print('d', np.array(delta).shape)
                    print('alpha={},beta={},gamma={}'.format(a,b,g))
                    diff.append(diff_p)
                plt.plot(delta, diff[i], label='$\gamma$={}'.format(g))
                i += 1
            leg = ax1.legend()
            plt.show()

        ################## MIGLIOR SAMPLING PER CIASCUN RANGE DI PUNTI #######
    min_points = results['selected points'].min()
    max_points = results['selected points'].max()
    div = 10
    number_limits = np.arange(min_points, max_points + (max_points - min_points) / div, (max_points - min_points) / div)
    for i in np.arange(div):
        diff = results[
            (results['selected points'] >= number_limits[i]) & (results['selected points'] < number_limits[i + 1])]
        diff_max = diff.nlargest(1, columns=['difference'])
        print('For the range {} - {}'.format(number_limits[i], number_limits[i + 1]))
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', -1)
        print(diff_max)



###################  COMPARISON BETWEEN ZERNIKE TOTAL AND THE NEWLY DEFINED SAMPLING/RANDOM ZERNIKE #########
print('Vuoi vedere la differenza Zernike totale vs sampling/random?')
o = input()
if o == 'y':
    results = pd.DataFrame(columns=['alpha', 'beta', 'gamma', 'delta', 'selected points', 'm(tot vs s)', 'm(tot vs r)'])
    #results.to_csv(
    #    "{}\EDOTTIA\\difference_{}_{}_over_{}points.csv".format(respath, fragment, cluster, int(lag / frac_test)),
    #    mode='a')

    zernike = np.loadtxt("{}\zernike\zernike_positive\zernike_total.dat".format(respath), delimiter=" ", skiprows=1)
    alpha1 = [1.]
    for a in [1.0]:
        for b in [0.0]:
            for g in [4.]:
                for d in [4.]:
                    with open("{}\EDOTTIA\\a{}_b{}_g{}_d{}.txt".format(respath, a,b,g,d)) as f:
                        sampling_index = list(set([int(float(x)) for x in f.read().split()]))
                        number_points = len(sampling_index)
                    print(number_points)
                    random_index = sorted(random.sample(list(np.arange(0, lag)), number_points))
                    surface_sampling = surf[[i for i in sampling_index], :]
                    surface_random = surf[[i for i in random_index], :]
                    if PLOT == 1:
                        res, c = SF.ConcatenateFigPlots(list_=[surface_sampling[:, :3]])
                        SF.Plot3DPoints(res[:, 0], res[:, 1], res[:, 2], c)
                        res_r, c_r = SF.ConcatenateFigPlots(list_=[surface_random[:, :3]])
                        SF.Plot3DPoints(res_r[:, 0], res_r[:, 1], res_r[:, 2], c_r)

                    test_index = sorted(random.sample(list(np.arange(0, lag)), int(lag/frac_test)))
                    zernike_total = zernike[:, [i for i in test_index]]

                    surf_sampling = SF.Surface(surface_sampling[:, :], patch_num=0, r0=R_zernike, theta_max=45)
                    surf_random = SF.Surface(surface_random[:,:], patch_num=0, r0=R_zernike, theta_max=45)

                    zernike_sampling = np.zeros((121, int(lag/frac_test)))
                    zernike_random = np.zeros((121, int(lag/frac_test)))
                    j = 0
                    for i in test_index:
                        sys.stderr.write("\r Processing {} out of {} point for alpha={}, beta={}, gamma={}, delta={}".format(j,int(lag/frac_test),a,b,g,d))
                        sys.stderr.flush()
                        # finding antigen patch, plane and zernike descriptors..
                        patch, mask = surf_sampling.BuildMixedPatch(point_pos=surf[i,:], Dmin=.5)
                        surf_sampling.real_br = mask
                        rot_patch, rot_ag_patch_nv = surf_sampling.PatchReorientNew(patch, verso)
                        z = surf_sampling.FindOrigin(rot_patch)
                        plane, weights, dist_plane, thetas = surf_sampling.CreatePlane(patch=rot_patch, z_c=z, Np=Npixel)
                        new_plane = plane.copy()
                        new_plane___ = plane.copy()
                        if (np.shape(rot_patch)[1] == 4):
                            new_plane_re = surf_sampling.FillTheGap_everywhere(plane_=np.real(plane))
                            new_plane_im = surf_sampling.FillTheGap_everywhere(plane_=np.imag(plane))
                            new_plane_re_ = surf_sampling.EnlargePixels(new_plane_re)
                            new_plane_im_ = surf_sampling.EnlargePixels(new_plane_im)
                            new_plane_ = new_plane_re_ + 1j * new_plane_im_ / np.max(np.abs(new_plane_im_))
                        else:
                            new_plane = surf_sampling.FillTheGap_everywhere(plane_=plane)
                            ## enlarging plane..
                            new_plane_ = surf_sampling.EnlargePixels(new_plane)
                        try:
                            zernike_env.img = new_plane_
                        except:
                            zernike_env = ZF.Zernike2d(new_plane_)
                        # br_recon, br_coeff = zernike_env.ZernikeReconstruction(order=ZOrder, PLOT=0)
                        br_coeff = zernike_env.ZernikeDecomposition(order=ZOrder)
                        zernike_sampling[:, j] = np.absolute(br_coeff)
                        #############   RANDOM #############
                        patch_r, mask_r = surf_random.BuildMixedPatch(point_pos=surf[i,:], Dmin=.5)
                        surf_random.real_br = mask_r
                        rot_patch_r, rot_ag_patch_nv_r = surf_random.PatchReorientNew(patch_r, verso)
                        z_r = surf_random.FindOrigin(rot_patch_r)
                        plane_r, weights_r, dist_plane_r, thetas_r = surf_random.CreatePlane(patch=rot_patch_r, z_c=z_r, Np=Npixel)
                        new_plane_r = plane_r.copy()
                        new_plane___r = plane_r.copy()
                        if (np.shape(rot_patch_r)[1] == 4):
                            new_plane_re_r = surf_random.FillTheGap_everywhere(plane_=np.real(plane_r))
                            new_plane_im_r = surf_random.FillTheGap_everywhere(plane_=np.imag(plane_r))
                            new_plane_re_r = surf_random.EnlargePixels(new_plane_re_r)
                            new_plane_im_r = surf_random.EnlargePixels(new_plane_im_r)
                            new_plane_r = new_plane_re_r + 1j * new_plane_im_r / np.max(np.abs(new_plane_im_r))
                        else:
                            new_plane_r = surf_random.FillTheGap_everywhere(plane_=plane_r)
                            ## enlarging plane..
                            new_plane_r = surf_random.EnlargePixels(new_plane_r)
                        try:
                            zernike_env_r.img = new_plane_r
                        except:
                            zernike_env_r = ZF.Zernike2d(new_plane_r)
                        # br_recon, br_coeff = zernike_env.ZernikeReconstruction(order=ZOrder, PLOT=0)
                        br_coeff_r = zernike_env_r.ZernikeDecomposition(order=ZOrder)
                        zernike_random[:, j] = np.absolute(br_coeff_r)
                        j += 1

                    #total_sampling.append(np.mean(np.linalg.norm(np.array(zernike_total) - np.array(zernike_sampling), axis=0)))
                    #total_random.append(np.mean(np.linalg.norm(np.array(zernike_total) - np.array(zernike_random), axis=0)))
                    results.loc[0,'alpha'] = a
                    results.loc[0,'beta'] = b
                    results.loc[0,'gamma'] = g
                    results.loc[0,'delta'] = d
                    results.loc[0,'selected points'] = number_points
                    results.loc[0,'m(tot vs s)'] = np.mean(np.linalg.norm(np.array(zernike_total) - np.array(zernike_sampling), axis=0))
                    results.loc[0,'m(tot vs r)'] = np.mean(np.linalg.norm(np.array(zernike_total) - np.array(zernike_random), axis=0))
                    #results.to_csv("{}\EDOTTIA\\difference_{}_{}_over_{}points.csv".format(respath,fragment,cluster,int(lag/frac_test)),mode='a', header=False)

################### !!BOX PLOT!!!  COMPARISON BETWEEN ZERNIKE TOTAL AND THE NEWLY DEFINED SAMPLING/RANDOM ZERNIKE #########
print('Vuoi vedere il BOX PLOT per la differenza Zernike totale vs sampling/random?')
o = input()
if o == 'y':
    n_bins = 10
    PLOT = 0
    a = 1.
    b = .0
    g = 6.0
    d = 0.

    with open("{}\EDOTTIA\\a{}_b{}_g{}_d{}.txt".format(respath, a, b, g, d)) as f:
        sampling_index = list(set([int(float(x)) for x in f.read().split()]))
        number_points = len(sampling_index)
    random_index = sorted(random.sample(list(np.arange(0, lag)), number_points))
    surface_sampling = surf[[i for i in sampling_index], :]
    surface_random = surf[[i for i in random_index], :]
    if PLOT == 1:
        res, c = SF.ConcatenateFigPlots(list_=[surface_sampling[:, :3]])
        SF.Plot3DPoints(res[:, 0], res[:, 1], res[:, 2], c)
        res_r, c_r = SF.ConcatenateFigPlots(list_=[surface_random[:, :3]])
        SF.Plot3DPoints(res_r[:, 0], res_r[:, 1], res_r[:, 2], c_r)

    zernike = np.loadtxt("{}\zernike\zernike_positive\zernike_total.dat".format(respath), delimiter=" ", skiprows=1)
    with open("{}\COSINE_step_{}_Rs_{}_FAST.txt".format(respath, step, Rs_select)) as f:
        cosine = [float(x) for x in f.read().split()]
    min_cosine = np.min(cosine)
    max_cosine = np.max(cosine)
    #limiti (incluso il primo e l'ultimo -> per 10 divisioni sono 11 -> il loro indice va da 0 a 10)
    cosine_limits = np.arange(min_cosine,max_cosine+(max_cosine-min_cosine)/n_bins,(max_cosine-min_cosine)/n_bins)
    cosine_range = []
    len_range = []
    mean_bin_values = []
    for bin in np.arange(0,n_bins):
        elements_range = np.where((cosine >= cosine_limits[bin])&(cosine<cosine_limits[bin+1]))[0]
        mean_bin_values.append('%.2f' % float((cosine_limits[bin+1]+cosine_limits[bin])/2))
        len_range.append(len(elements_range))
        cosine_range.append(elements_range)

    n_test_points = min(len_range)
    bin_sampling = []
    bin_random = []
    for bin in np.arange(0,n_bins):
        test_index = sorted(random.sample(list(cosine_range[bin]), n_test_points))
        zernike_total = zernike[:, [i for i in test_index]]

        surf_sampling = SF.Surface(surface_sampling[:, :], patch_num=0, r0=R_zernike, theta_max=45)
        surf_random = SF.Surface(surface_random[:, :], patch_num=0, r0=R_zernike, theta_max=45)

        P = 0 #confronto tra patch random/samplin/totale per diverse roughness
        if P == 1:
            point_s = my_functions.find_nearest_vector(surface_sampling[:,:3],surf[test_index[10],:3])
            point_r = my_functions.find_nearest_vector(surface_random[:,:3],surf[test_index[10],:3])
            patch_s, _ = surf_sampling.BuildPatch(point_pos=point_s, Dmin=.5)
            patch_r, _ = surf_random.BuildPatch(point_pos=point_r, Dmin=.5)
            patch_t, _ = surf_obj.BuildPatch(point_pos=test_index[10], Dmin=.5)

            #res_s, c_s = SF.ConcatenateFigPlots([patch_t[:,:3],patch_s[:, :3]])
            res_s, c_s = SF.ConcatenateFigPlots([patch_s[:, :3]])
            SF.Plot3DPoints(res_s[:, 0], res_s[:, 1], res_s[:, 2], c_s, 0.3)

            res_r, c_r = SF.ConcatenateFigPlots([patch_r[:, :3]])
            SF.Plot3DPoints(res_r[:, 0], res_r[:, 1], res_r[:, 2], c_r, 0.3)

            res_t, c_t = SF.ConcatenateFigPlots([patch_t[:, :3]])
            SF.Plot3DPoints(res_t[:, 0], res_t[:, 1], res_t[:, 2], c_t, 0.3)

        zernike_sampling = np.zeros((121, n_test_points))
        zernike_random = np.zeros((121, n_test_points))
        j = 0
        for i in test_index:
            sys.stderr.write("\r Processing %i out of %i point for bin =%i (alpha=%i, beta=%i, gamma=%i, delta=%i)" % (
            j, n_test_points, bin, a, b, g, d))
            sys.stderr.flush()
            # finding antigen patch, plane and zernike descriptors..
            patch, mask = surf_sampling.BuildMixedPatch(point_pos=surf[i, :], Dmin=.5)
            surf_sampling.real_br = mask
            rot_patch, rot_ag_patch_nv = surf_sampling.PatchReorientNew(patch, verso)
            z = surf_sampling.FindOrigin(rot_patch)
            plane, weights, dist_plane, thetas = surf_sampling.CreatePlane(patch=rot_patch, z_c=z, Np=Npixel)
            new_plane = plane.copy()
            new_plane___ = plane.copy()
            if (np.shape(rot_patch)[1] == 4):
                new_plane_re = surf_sampling.FillTheGap_everywhere(plane_=np.real(plane))
                new_plane_im = surf_sampling.FillTheGap_everywhere(plane_=np.imag(plane))
                new_plane_re_ = surf_sampling.EnlargePixels(new_plane_re)
                new_plane_im_ = surf_sampling.EnlargePixels(new_plane_im)
                new_plane_ = new_plane_re_ + 1j * new_plane_im_ / np.max(np.abs(new_plane_im_))
            else:
                new_plane = surf_sampling.FillTheGap_everywhere(plane_=plane)
                ## enlarging plane..
                new_plane_ = surf_sampling.EnlargePixels(new_plane)
            try:
                zernike_env.img = new_plane_
            except:
                zernike_env = ZF.Zernike2d(new_plane_)
            # br_recon, br_coeff = zernike_env.ZernikeReconstruction(order=ZOrder, PLOT=0)
            br_coeff = zernike_env.ZernikeDecomposition(order=ZOrder)
            zernike_sampling[:, j] = np.absolute(br_coeff)

                                #############   RANDOM #############
            patch_r, mask_r = surf_random.BuildMixedPatch(point_pos=surf[i, :], Dmin=.5)
            surf_random.real_br = mask_r
            rot_patch_r, rot_ag_patch_nv_r = surf_random.PatchReorientNew(patch_r, verso)
            z_r = surf_random.FindOrigin(rot_patch_r)
            plane_r, weights_r, dist_plane_r, thetas_r = surf_random.CreatePlane(patch=rot_patch_r, z_c=z_r,
                                                                                 Np=Npixel)
            new_plane_r = plane_r.copy()
            new_plane___r = plane_r.copy()
            if (np.shape(rot_patch_r)[1] == 4):
                new_plane_re_r = surf_random.FillTheGap_everywhere(plane_=np.real(plane_r))
                new_plane_im_r = surf_random.FillTheGap_everywhere(plane_=np.imag(plane_r))
                new_plane_re_r = surf_random.EnlargePixels(new_plane_re_r)
                new_plane_im_r = surf_random.EnlargePixels(new_plane_im_r)
                new_plane_r = new_plane_re_r + 1j * new_plane_im_r / np.max(np.abs(new_plane_im_r))
            else:
                new_plane_r = surf_random.FillTheGap_everywhere(plane_=plane_r)
                ## enlarging plane..
                new_plane_r = surf_random.EnlargePixels(new_plane_r)
            try:
                zernike_env_r.img = new_plane_r
            except:
                zernike_env_r = ZF.Zernike2d(new_plane_r)
            # br_recon, br_coeff = zernike_env.ZernikeReconstruction(order=ZOrder, PLOT=0)
            br_coeff_r = zernike_env_r.ZernikeDecomposition(order=ZOrder)
            zernike_random[:, j] = np.absolute(br_coeff_r)
            j += 1

        res_inv = np.row_stack([test_index, zernike_sampling])
        res_inv_r = np.row_stack([test_index, zernike_random])

        bin_sampling.append(np.linalg.norm(np.array(zernike_total) - np.array(zernike_sampling), axis=0))
        bin_random.append(np.linalg.norm(np.array(zernike_total) - np.array(zernike_random), axis=0))

    np.savetxt("{}\EDOTTIA\\BOX_R_a{}_b{}_g{}_d{}.txt".format(respath, a, b, g, d), bin_random)
    np.savetxt("{}\EDOTTIA\\BOX_S_a{}_b{}_g{}_d{}.txt".format(respath, a, b, g, d), bin_sampling)

    ticks = mean_bin_values
    def set_box_color(bp, color):
        plt.setp(bp['boxes'], color=color)
        plt.setp(bp['whiskers'], color=color)
        plt.setp(bp['caps'], color=color)
        plt.setp(bp['medians'], color=color)

    bps = plt.boxplot(bin_sampling, positions=np.array(range(len(bin_sampling))) * 2.0 - 0.4, sym='', widths=0.6)
    bpr = plt.boxplot(bin_random, positions=np.array(range(len(bin_random))) * 2.0 + 0.4, sym='', widths=0.6)
    set_box_color(bps, '#D7191C')  # colors are from http://colorbrewer2.org/
    set_box_color(bpr, '#2C7BB6')
    # draw temporary red and blue lines and use them to create a legend
    plt.plot([], c='#D7191C', label='From total ({} points) and sampled surface ({} points)'.format(lag,number_points))
    plt.plot([], c='#2C7BB6', label='From total ({} points) and random surface ({} points)'.format(lag,number_points))
    plt.title(
        "Distance between Zernike vectors for $\\alpha$={}, $\\beta$={}, $\\gamma$={}, $\\delta$={} for different ranges of roughness ({} for each interval)".format(a, b, g, d, len(bin_sampling[1]), lag),fontsize=22)
    plt.xlabel("Mean cosine of the considered patches",fontsize=22)
    plt.ylabel("Norm of the distance between Zernike vectors",fontsize=22)

    plt.legend(fontsize=16)

    plt.xticks(range(0, len(ticks) * 2, 2), ticks,fontsize=22)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.show()

#############################   ROUGHNESS DELLE PATCH  ###########
print("Vuoi calcolare la MEDIA DEL COSENO di ciascuna possibile patch? y o n?")
ooo = input()
if ooo == "y":
    cosine = my_functions.RapidCos(ltmp, surf, surf_obj_scan, respath, Rs_select, step)


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
    ax1.set_ylabel("Number of sampled points over {} points".format(lag))
    plt.plot(alpha, alpha_points_plot)
    plt.show()



print('Vuoi vedere i coseni nella superifice?')
o = input()
if o == 'y':
    res, c = SF.ConcatenateFigPlots(list_=[surf[:, :3]])
    c1 = np.loadtxt("{}\COSINE_step_{}_Rs_{}_FAST.txt".format(respath, step, Rs_select))
    SF.Plot3DPoints(res[:, 0], res[:, 1], res[:, 2], color=c1)


            ################## ZERNIKE PER OGNI SINGOLO PUNTO   #####################
print("Vuoi fare ZERNIKE per OGNI {} PUNTI? y o n?".format(Npoint))
o = input()
if o == "y":
    index_all_possible_area = np.arange(ltmp)[::Npoint]
    ltmp = len(index_all_possible_area)
    zernike_sampling_inv_a = np.zeros((121, ltmp))

    all_patches = []
    for i in range(ltmp):
        sys.stderr.write("\r Processing %i out of %i point" % (i, ltmp))
        sys.stderr.flush()
        # finding antigen patch, plane and zernike descriptors..

        patch, mask = surf_obj.BuildPatch(point_pos=index_all_possible_area[i], Dmin=.5)
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
        # br_recon, br_coeff = zernike_env.ZernikeReconstruction(order=ZOrder, PLOT=0)
        br_coeff = zernike_env.ZernikeDecomposition(order=ZOrder)

        zernike_sampling_inv_a[:, i] = np.absolute(br_coeff)


    res_inv_ = np.row_stack([index_all_possible_area, zernike_sampling_inv_a])
    if verso == 1:
        np.savetxt("{}/zernike/zernike_positive/zernike_total.dat".format(respath), res_inv_, fmt="%.4e")

    else:
        np.savetxt("{}/zernike/zernike_negative/zernike_total.dat".format(respath), res_inv_, fmt="%.4e")



print("Vuoi vedere la media dei coseni per ciascuna patch? y o n?")
o = input()
if o == "y":
    with open("{}\COSINE_step_{}_Rs_{}_FAST.txt".format(respath, step, Rs_select)) as f:
        mean_cos = [float(x) for x in f.read().split()]

    x, bins, p = plt.hist(mean_cos, density=True,bins='auto', color='powderblue')
    #plt.hist(mean_cos, bins='auto', densitiy=True, color='powderblue')  # arguments are passed to np.histogram
    plt.title("Mean cosine value of the surface's {} patches".format(len(mean_cos)),fontsize=40)
    plt.xlabel('$R_i$', size=40)
    plt.ylabel('Counts density', size=40)
    plt.yticks(fontsize=40)
    plt.xticks(fontsize=40)
    plt.show()

    ############### SCREENING PARTENDO DA COSINE_FAST.txt    #################
print("Vuoi fare lo SCREENING partendo da cosine.txt per diversi valori di alpha? y o n?")
o = input()
if o == "y":
    with open("{}\COSINE_step_{}_Rs_{}_FAST.txt".format(respath, step, Rs_select)) as f:
        #mean_cos = [2 * (float(x)) - 1 for x in f.read().split()]
        mean_cos = [float(x) for x in f.read().split()]
    for i in alpha:
       index_possible_area = my_functions.ContinuousDistributionFast(ltmp,mean_cos, surf, surf_obj_scan, Rs_select, i, step, respath, screening)
       #index_possible_area = my_functions.ContinuousDistribution(mean_cos, surf, surf_obj_scan, Rs_select, i, step, respath, screening)
