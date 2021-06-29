#!/usr/bin/env python
# coding: utf-8

import os, sys
import numpy as np
import scipy as sp
import pandas as pd
import csv
import matplotlib.pyplot as plt

    ########    PARAMETERS  #######
fragment1 = 208
fragment2 = 220
cluster1 = 5
cluster2 = 2
sign1 = 1
sign2 = -1
R_zernike = 6
R_s = 4
Npoint1 = 1
Npoint2 = 1
#########   END PARAMETERS  #########

if sign1 == 1:
    verso1 = "positive"
else:
    verso1 = "negative"
if sign2 == 1:
    verso2 = "positive"
else:
    verso2= "negative"

pdb_file1 = "..\\{}\cluster{}.dms".format(fragment1, cluster1)
pdb_file2 = "..\\{}\cluster{}.dms".format(fragment2, cluster2)

path1  = "..\\{}\cluster{}\R_zernike_{}\R_s_{}\zernike".format(fragment1, cluster1, R_zernike, R_s)
path2  = "..\\{}\cluster{}\R_zernike_{}\R_s_{}\zernike".format(fragment2, cluster2, R_zernike, R_s)

#with open("{}\index_alpha.txt".format(path1)) as f:
#    index_possible_points1 = [int(float(x)) for x in f.read().split()]

surf_total1_ = pd.read_csv(pdb_file1)
lag1 = len(surf_total1_["x"])
surf_total1 = np.zeros((lag1, 6))
surf_total1[:, :] = surf_total1_[["x", "y", "z", "Nx", "Ny", "Nz"]]
print("Original dim of surf of cluster {} from fragment {}: {}".format(cluster1, fragment1, np.shape(surf_total1)[0]))


#surf_alpha_1 = surf_total1[index_possible_points1, :]
#rr1 = surf_total1_["Res"][index_possible_points1]
surf_alpha_1 = surf_total1[::Npoint1, :]
rr1 = surf_total1_["Res"][::Npoint1]


#with open("{}\index_alpha.txt".format(path2)) as f:
#    index_possible_points2 = [int(float(x)) for x in f.read().split()]


surf_total2_ = pd.read_csv(pdb_file2)
lag2 = len(surf_total2_["x"])
surf_total2 = np.zeros((lag2, 6))
surf_total2[:, :] = surf_total2_[["x", "y", "z", "Nx", "Ny", "Nz"]]
print("Original dim of surf of cluster {} from fragment {}: {}".format(cluster2, fragment2, np.shape(surf_total2)[0]))

#surf_alpha_2 = surf_total2[index_possible_points2, :]
#rr2 = surf_total2_["Res"][index_possible_points2]
surf_alpha_2 = surf_total2[::Npoint2, :]
rr2 = surf_total2_["Res"][::Npoint2]

try:
   # zern_1 = np.loadtxt("{}\zernike_{}\zernike_alpha.dat".format(path1, verso1))
   # zern_2 = np.loadtxt("{}\zernike_{}\zernike_alpha.dat".format(path2, verso2))

    zern_1 = np.loadtxt("{}\zernike_{}\zernike_total.dat".format(path1, verso1))
    zern_2 = np.loadtxt("{}\zernike_{}\zernike_total.dat".format(path2, verso2))
except:
    print("Invariants for complex are not present!")
    exit()



l1 = np.shape(zern_1)[1]
l2 = np.shape(zern_2)[1]


MAT = 1
SAVE_FIG = 0
PLOT = 0

res_dir = "."


def EuclideanDistance(a,b):    
    d = np.sqrt(np.sum( (a-b)**2 ) )
    return(d)


def Zscore(x):
    z = (x-np.mean(x))/np.std(x)
    return(z)

surf_1_Coord = surf_alpha_1[:,0:3]
surf_2_Coord = surf_alpha_2[:,0:3]

print('Vuoi calcolare la binding propensity?')
o =input()
if o =='y':

    if(l1*l2>1e8):
        MAT = 0



    if(MAT):
        d_all_theo = -2*np.dot(np.transpose(zern_1[1:,:]), zern_2[1:,:])
        d1_square = np.sum( zern_1[1:,:]*zern_1[1:,:], axis=0)
        d2_square = np.sum( zern_2[1:,:]*zern_2[1:,:], axis=0)


        for i in range(l1):
            d_all_theo[i,:] += d2_square
        for i in range(l2):
            d_all_theo[:,i] += d1_square


        d_all_theo = np.sqrt(d_all_theo)

        color_1 = np.min(d_all_theo, axis=1)
        color_2 = np.min(d_all_theo, axis=0)

        pos_min_1 = np.argmin(d_all_theo, axis=1)
        pos_min_2 = np.argmin(d_all_theo, axis=0)

    else:
        color_1 = np.zeros(l1)
        color_2 = np.zeros(l2)


        tmp = np.transpose(zern_2[1:,:])
        for i in range(l1):
            ddd =  np.sqrt(np.sum( (tmp-  zern_1[1:,i])**2 , axis=1))
            sys.stderr.write("\r %d of %d"%(i, l1))
            sys.stderr.flush()
            color_1[i] = np.min(ddd)

        tmp = np.transpose(zern_1[1:,:])
        for i in range(l2):
            sys.stderr.write("\r %d of %d"%(i, l2))
            sys.stderr.flush()
            ddd = np.sqrt(np.sum( (tmp-  zern_2[1:,i])**2 , axis=1))
            color_2[i] = np.min(ddd)

    len_1 = np.shape(surf_alpha_1)[0]
    len_2 = np.shape(surf_alpha_2)[0]




    df_1 = pd.DataFrame(surf_alpha_1[:,:3], columns= ["x", "y","z"]) #color_1_, rr1])
    df_1.insert(3,"c", color_1, True)
    df_1.insert(4,"bs", np.zeros(len_1), True)
    df_1.insert(5,"res", list(rr1), True)

    df_2 = pd.DataFrame(surf_alpha_2[:,:3], columns= ["x", "y","z"])
    df_2.insert(3,"c", color_2, True)
    df_2.insert(4,"bs", np.zeros(len_2), True)
    df_2.insert(5,"res", list(rr2), True)

    df_2.to_csv("{}\BindProp_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path2, fragment2, cluster2, verso2, fragment1, cluster1, verso1))



    df_1.to_csv("{}\BindProp_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path1, fragment1, cluster1, verso1, fragment2, cluster2, verso2))

print('Vuoi fare la media delle binding propensity?')
o=input()
if o =='y':
    file1 ="{}\BindProp_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path1, fragment1, cluster1, verso1, fragment2, cluster2, verso2)
    color_1 = pd.read_csv(file1, usecols=['c'], squeeze=True)

    smoothing_1 = []
    for i in range(l1):
        sys.stderr.write("\r Processing %i out of %i point" % (i, l1))
        sys.stderr.flush()
        d2 = (surf_alpha_1[:, 0] - surf_alpha_1[i, 0]) ** 2 + (
                            surf_alpha_1[:, 1] - surf_alpha_1[i, 1]) ** 2 + (
                                 surf_alpha_1[:, 2] - surf_alpha_1[i, 2]) ** 2

        mask = d2 <= R_zernike ** 2
        patch = color_1[mask]
        smoothing_1.append(np.mean(patch))

    smoothing1 = pd.DataFrame(smoothing_1, columns= ["smooth"])
    smoothing1.to_csv("{}\smooth_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path1, fragment1, cluster1, verso1, fragment2, cluster2, verso2))

    file2 ="{}\BindProp_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path2, fragment2, cluster2,
                                                                                         verso2, fragment1,
                                                                                       cluster1, verso1)
    color_2 = pd.read_csv(file2, usecols=['c'], squeeze=True)

    smoothing_2 = []
    for i in range(l2):
        sys.stderr.write("\r Processing %i out of %i point" % (i, l2))
        sys.stderr.flush()
        d2 = (surf_alpha_2[:, 0] - surf_alpha_2[i, 0]) ** 2 + (
                surf_alpha_2[:, 1] - surf_alpha_2[i, 1]) ** 2 + (
                     surf_alpha_2[:, 2] - surf_alpha_2[i, 2]) ** 2
        mask = d2 <= R_zernike ** 2
        patch = color_2[mask]
        smoothing_2.append(np.mean(patch))

    smoothing2 = pd.DataFrame(smoothing_2, columns=["smooth"])
    smoothing2.to_csv(
        "{}\smooth_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path2, fragment2,
                                                                                              cluster2, verso2,
                                                                                              fragment1, cluster1,
                                                                                              verso1))

print('Vuoi passare ai resiudi?')
o = input()
if o == 'y':
    smooth_point = pd.DataFrame()
    file_point ="{}\BindProp_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path1, fragment1, cluster1, verso1, fragment2, cluster2, verso2)
    file_smooth ="{}\smooth_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path1, fragment1, cluster1, verso1, fragment2, cluster2, verso2)

    smooth_point['res'] = pd.read_csv(file_point, usecols=['res'], squeeze=True)
    smooth_point['smooth'] = pd.read_csv(file_smooth, usecols=['smooth'], squeeze=True)

    #smooth_res = smooth_point.groupby('res',as_index=False)['smooth'].mean()
    smooth_point['average'] = smooth_point.groupby('res', as_index=False)['smooth'].transform('mean')
    smooth_res = smooth_point.drop_duplicates(subset='res',keep="first", inplace=False)

    smooth_res.to_csv(
        "{}\smooth_RES_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path1, fragment1,
                                                                                              cluster1, verso1,
                                                                                              fragment2, cluster2,
                                                                                              verso2))

    smooth_point = pd.DataFrame()
    file_point = "{}\BindProp_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path2,
                                                                                                         fragment2,
                                                                                                         cluster2,
                                                                                                         verso2,
                                                                                                         fragment1,
                                                                                                         cluster1,
                                                                                                         verso1)
    file_smooth = "{}\smooth_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path2,
                                                                                                        fragment2,
                                                                                                        cluster2,
                                                                                                        verso2,
                                                                                                        fragment1,
                                                                                                        cluster1,
                                                                                                        verso1)

    smooth_point['res'] = pd.read_csv(file_point, usecols=['res'], squeeze=True)
    smooth_point['smooth'] = pd.read_csv(file_smooth, usecols=['smooth'], squeeze=True)

    smooth_point['average'] = smooth_point.groupby('res', as_index=False)['smooth'].transform('mean')
    smooth_res = smooth_point.drop_duplicates(subset='res', keep="first", inplace=False)

    smooth_res.to_csv(
        "{}\smooth_RES_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path2, fragment2,
                                                                                                  cluster2, verso2,
                                                                                                  fragment1,
                                                                                                  cluster1,
                                                                                                  verso1))




print('Vuoi vedere il grafico della b p s?')
o = input()
if o == 'y':
    file1 ="{}\smooth_RES_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path1, fragment1, cluster1, verso1, fragment2, cluster2, verso2)
    smooth = pd.read_csv(file1, usecols=['average'], squeeze=True)
    res = pd.read_csv(file1, usecols=['res'], squeeze=True)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Index of the residues of the cluster {} of fragment {}'. format(cluster1-1, fragment1))
    ax1.set_ylabel("Smoothed binding propensity")
    plt.plot(smooth)
    plt.title('Binding propensity between cluster {} of fragment {} and cluster {} of fragment {}'.format(cluster1-1,fragment1,cluster2-1, fragment2), fontsize=10)
    plt.show()

    file2 ="{}\smooth_RES_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path2, fragment2, cluster2,
                                                                                         verso2, fragment1,
                                                                                 cluster1, verso1)
    smooth = pd.read_csv(file2, usecols=['average'], squeeze=True)
    res = pd.read_csv(file2, usecols=['res'], squeeze=True)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Index of the residues of the cluster {} of fragment {}'.format(cluster2 - 1, fragment2))
    ax1.set_ylabel("Smoothed binding propensity")
    plt.plot(smooth)
    plt.title(
        'Binding propensity between cluster {} of fragment {} and cluster {} of fragment {}'.format(cluster2 - 1,
                                                                                                    fragment2,
                                                                                                    cluster1 - 1,
                                                                                                    fragment1),
        fontsize=10)
    plt.show()
