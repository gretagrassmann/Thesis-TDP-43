#!/usr/bin/env python


import os, sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as mpl

import pandas as pd
from mayavi import mlab



import SurfaceFunc as SF
import matplotlib.pyplot as plt
import openpyxl

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

if fragment1 == 220:
    fragment_1 = 'B'
else:
    fragment_1 = 'A'
if fragment2 == 220:
    fragment_2 = 'B'
else:
    fragment_2 = 'A'

if sign1 == 1:
    verso1 = "positive"
else:
    verso1 = "negative"
if sign2 == 1:
    verso2 = "positive"
else:
    verso2= "negative"

path1 = "..\\{}\cluster{}\R_zernike_{}\R_s_{}\zernike".format(fragment1, cluster1, R_zernike, R_s)
path2 = "..\\{}\cluster{}\R_zernike_{}\R_s_{}\zernike".format(fragment2, cluster2, R_zernike, R_s)

def ColNorm(c):
    nc = (c-np.min(c))/(np.max(c)-np.min(c))
    return(nc)


def Smoothing(data, col, radius):
    '''
    This function performs a smoothing procedure of the 'col' vector according to spatial information contained in data.
    Each col value acquires a new value defined as the mean of the values of points that are closer than 'radius' in space.
    '''
    l = len(col)
    newcol = np.zeros(l)
    
    for i in range(l):
        
        d = np.sum((data - data[i,:])**2, axis = 1)
        
        mask = d <= radius**2
        
        newcol[i] = np.mean(col[mask])
        
    return(newcol)

print('Vuoi vedere il grafico per ciascun frammento?')
o = input()
if o =='y':
    for i in [220, 208]:
        if i == 220:
            j = [1, 2]
            frag = 'B'
        else:
            j = [1,2,3,4,5]
            frag = 'A'
        for ii in j:
            pdb_file = "..\\{}\cluster{}.dms".format(i, ii)
            res_original = pd.DataFrame()
            res_original['res'] = pd.read_csv(pdb_file, usecols=['Res'], squeeze=True)
            res_original = res_original.drop_duplicates(subset='res', keep="first", inplace=False).set_index('res')
            res_original['BP'] = 0

            file = '..\..\..\complementarity_regions\\{}.csv'.format(i)
            fragment_results = pd.DataFrame()
            fragment_results['res'] = pd.read_csv(file, usecols=['res'], squeeze=True)
            fragment_results['BP'] = pd.read_csv(file, usecols=['BP'], squeeze=True)
            fragment_results['cluster'] = pd.read_csv(file, usecols=['cluster'])
            fragment_results = fragment_results[fragment_results.cluster.isin([ii])]
            fragment_bp = pd.DataFrame()
            fragment_bp = fragment_results.groupby('res', as_index=False)['BP'].sum()
            interacting_res = list(fragment_bp['res'])
            fragment_bp = fragment_bp.set_index('res')


            for r in interacting_res:
               res_original.loc[r, 'BP'] = fragment_bp.loc[r, 'BP']

            res_original.to_csv('..\..\..\complementarity_regions\\{}_{}.csv'.format(i, ii))

            res_original.plot.bar()
            plt.tight_layout()
            plt.title('Results for cluster {} of fragment {}'.format(ii-1, frag))
        plt.show()

        cluster_mean = pd.DataFrame()
        cluster_mean['res'] = res_original.index
        cluster_mean['BP'] = 0
        cluster_mean = cluster_mean.set_index('res')
        for ii in j:
            file = '..\..\..\complementarity_regions\\{}_{}.csv'.format(i, ii)
            cluster_results = pd.DataFrame()
            cluster_results['res'] = pd.read_csv(file, usecols=['res'], squeeze=True)
            cluster_results['BP'] = pd.read_csv(file, usecols=['BP'], squeeze=True)
            interacting_res = list(cluster_results['res'])
            cluster_results = cluster_results.set_index('res')

            for r in interacting_res:
                cluster_mean.loc[r, 'BP'] = cluster_mean.loc[r, 'BP'] + cluster_results.loc[r, 'BP']

        cluster_mean.plot.bar()
        plt.tight_layout()
        plt.title('Results for fragment {}'.format(frag))
        plt.show()

#print('Vuoi vedere il grafico delle superifici migliori?')
#o = input()
o = 'n'
if o == 'y':

    comb = ['220_1 vs 220_1',
            '220_2 vs 220_2',
            '220_1 vs 220_2',
            '208_1 vs 208_1',
            '208_2 vs 208_2',
            '208_3 vs 208_3',
            '208_4 vs 208_4',
            '208_5 vs 208_5',
            '208_1 vs 208_2',
            '208_1 vs 208_3',
            '208_1 vs 208_4',
            '208_1 vs 208_5',
            '208_2 vs 208_3',
            '208_2 vs 208_4',
            '208_2 vs 208_5',
            '208_3 vs 208_4',
            '208_3 vs 208_5',
            '208_4 vs 208_5',
            '208_1 vs 220_1',
            '208_1 vs 220_2',
            '208_2 vs 220_1',
            '208_2 vs 220_2',
            '208_3 vs 220_1',
            '208_3 vs 220_2',
            '208_4 vs 220_1',
            '208_4 vs 220_2',
            '208_5 vs 220_1',
            '208_5 vs 220_2'
            ]
    pb_o = np.loadtxt('..\..\..\\results.txt')
    l = len(pb_o)

    bp = ColNorm(pb_o)
    bp = np.abs(bp -1.)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xticks(np.arange(l))
    ax1.set_xticklabels(comb, rotation='vertical', fontsize=8)
    ax1.set_xlabel('Combination')
    ax1.set_ylabel("Best binding propensity for each combination")
    plt.plot(bp)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xticks(np.arange(l))
    ax1.set_xticklabels(comb, rotation='vertical', fontsize=8)
    ax1.set_xlabel('Combination')
    ax1.set_ylabel("Best binding propensity for each combination")
    plt.plot(pb_o)

    plt.show()





PLOT = 0

Rsmooth = 6  #A, the radius of the sphere to be used in the smoothing procedure..
D = 1.       #A, the threshold distance to discriminate the 'islands', i.e. the groups of points belonging to the same patch..



file1 = "{}\BindProp_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path1, fragment1, cluster1, verso1, fragment2, cluster2, verso2)
file2 = "{}\BindProp_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path2, fragment2, cluster2, verso2, fragment1, cluster1, verso1)

Thres = .94 # the threshold for patch point selection based on rescaled binding propensity..


### Loading Binding propensity data...
prot1_ = pd.read_csv(file1)
prot2_ = pd.read_csv(file2)


prot1 = np.zeros((np.shape(prot1_)[0],4))
prot2 = np.zeros((np.shape(prot2_)[0],4))
prot1[:,:] = prot1_[["x","y","z","c"]]
prot2[:,:] = prot2_[["x","y","z","c"]]
    
res1 = np.array(prot1_["res"])
res2 = np.array(prot2_["res"])

print('Vuoi vedere l\' istogramma dell bp rough?')
o = input()
if o == 'y':
    fragment1 = 220
    cluster1 = 2
    total = []
    for fragment2 in [208, 220]:
        if fragment2 == 208:
            clusters = [1,2,3,4,5]
        else:
            clusters = [1,2]
        for cluster2 in clusters:
            verso1 = 'positive'
            verso2 = 'negative'
            path1 = "..\\{}\cluster{}\R_zernike_{}\R_s_{}\zernike".format(fragment1, cluster1, R_zernike, R_s)
            file1 = "{}\BindProp_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path1,
                                                                                                            fragment1,
                                                                                                            cluster1,
                                                                                                            verso1,
                                                                                                            fragment2,
                                                                                                            cluster2,
                                                                                                            verso2)

            if os.path.exists(file1):
                with open(file1, 'r') as f:
                    ### Loading Binding propensity data...
                    prot1_ = pd.read_csv(f)
                    prot1 = np.zeros((np.shape(prot1_)[0], 4))
                    prot1[:, :] = prot1_[["x", "y", "z", "c"]]
                    c1_ = prot1[:, -1]
                    total.extend(c1_)
            # Do Stuff with file
            else:
                verso1 = 'negative'
                verso2 = 'positive'
                path1 = "..\\{}\cluster{}\R_zernike_{}\R_s_{}\zernike".format(fragment1, cluster1, R_zernike, R_s)
                file1 = "{}\BindProp_fragment{}_cluster{}_verso_{}_VS_fragment{}_cluster{}_verso_{}.csv".format(path1,
                                                                                                                fragment1,
                                                                                                                cluster1,
                                                                                                                verso1,
                                                                                                                fragment2,
                                                                                                                cluster2,
                                                                                                                verso2)
                ### Loading Binding propensity data...
                prot1_ = pd.read_csv(file1)
                prot1 = np.zeros((np.shape(prot1_)[0], 4))
                prot1[:, :] = prot1_[["x", "y", "z", "c"]]
                c1_ = prot1[:, -1]
                total.extend(c1_)

    plt.hist(total)
    plt.xlabel("Rough binding propensity")
    plt.ylabel('Counts')
    plt.title("Rough binding propensity values for cluster {} of fragment {}".format(cluster1 - 1,fragment_1) )
    plt.show()

### rescaling of binding propensities: values are rescaled in such a way that all values are in the interval [0,1];
###                            
c1_ = prot1[:,-1]
c2_ = prot2[:,-1]

c1 = ColNorm(c1_)
c2 = ColNorm(c2_)

### flipping of binding propensities: values are flipped  in such a way that points with high complementarity have high values
###                                   of binding propensity and viceversa.
c1 = np.abs(c1 - 1.)
c2 = np.abs(c2 - 1.)



print('Plot before the smoothing')
if(PLOT):
    res, c = SF.ConcatenateFigPlots(list_=[prot1[:,:3]])
    SF.Plot3DPoints(res[:,0], res[:,1], res[:,2], color=c1)

    res, c = SF.ConcatenateFigPlots(list_=[prot2[:, :3]])
    SF.Plot3DPoints(res[:, 0], res[:, 1], res[:, 2], color=c2)

print('Plot histogram before the smoothing')
if(PLOT):
    plt.hist(c1)
    plt.xlabel("Binding propensity")
    plt.ylabel('Points in the surface')
    plt.title("Binding propensity between cluster {} from fragment {} and cluster {} from fragment {}".format(cluster1-1,fragment_1, cluster2-1, fragment_2))
    plt.show()

    plt.hist(c1)
    plt.xlabel("Binding propensity")
    plt.ylabel('Points in the surface')
    plt.title("Binding propensity between cluster {} from fragment {} and cluster {} from fragment {}".format(cluster2 - 1,fragment_2, cluster1-1, fragment_1))
    plt.show()

### Binding propensities are smoothed

newcol_1 = ColNorm(Smoothing(prot1, c1, Rsmooth))
newcol_2 = ColNorm(Smoothing(prot2, c2, Rsmooth))

original1 = ColNorm(Smoothing(prot1, c1_, Rsmooth))
original2 = ColNorm(Smoothing(prot2, c2_, Rsmooth))

print('Plot protein {} cluster {} (protein 1) after the smoothing'.format(fragment1,cluster1))
if(PLOT):
    res, c = SF.ConcatenateFigPlots(list_=[prot1[:,:3]])
    SF.Plot3DPoints(res[:,0], res[:,1], res[:,2], color = newcol_1)

print('Plot protein {} cluster {} (protein 2) after the smoothing'.format(fragment2,cluster2))
if(PLOT):
    res, c = SF.ConcatenateFigPlots(list_=[prot2[:,:3]])
    SF.Plot3DPoints(res[:,0], res[:,1], res[:,2], color = newcol_2)



### Interacting regions are found selecting all points whose binding prop is higher than Thres
mask1 = newcol_1 > Thres
res1 = res1[mask1]
original1 = original1[mask1]
## BP of the interacting points
bp_points1 = newcol_1[mask1]

mask2 = newcol_2 > Thres
res2 = res2[mask2]
original2 = original2[mask2]
## BP of the interacting points
bp_points2 = newcol_2[mask2]

### Isolating different regions... 
label_1 = SF.IsolateSurfaces(prot1[mask1,:3], minD= D)
label_2 = SF.IsolateSurfaces(prot2[mask2,:3], minD= D)

### Finding sizes of each region... 
labs1, counts1 = np.unique(label_1, return_counts=True)
labs2, counts2 = np.unique(label_2, return_counts=True)


### Finding relative sizes of each region (dividing by the total number of found points)...
freqs1 = counts1/np.sum(counts1)
freqs2 = counts2/np.sum(counts2)



print("prot1:",freqs1)
print("prot2:",freqs2)


### Discarding regions with small number of points..
labs_1 = labs1[freqs1 > 0.20]
print("Prot 1:")
if(len(labs_1) == 0):
    print("Warning! higly fragmented binding regions (no patch has more than 20% of the total idenitified points!")
    print("Biggest patch has been selected")
    labs_1 = [labs1[0]]
else:
    print("Number of found patches: %d"%len(labs_1))

    
labs_2 = labs2[freqs2 > 0.20]
print("Prot 2:")
if(len(labs_2) == 0):
    print("Warning! higly fragmented binding regions (no patch has more than 20% of the total idenitified points!")
    print("Biggest patch has been selected")
    labs_2 = [labs2[0]]
else:
    print("Number of found patches: %d"%len(labs_2))



### Finding residues of each region...
res1_set = []
bp1 = []
bp_residues1 = []
res1_set_repeated = []
for i in range(len(labs_1)):
    mask = label_1 == labs_1[i]
    res1_set_repeated.extend(res1[mask])
    res1_set.append(np.unique(res1[mask]))
    bp1.append(np.mean(original1[mask]))
    bp_residues1.extend(bp_points1[mask])

#bp_residues1 = [item for sublist in bp_residues1 for item in sublist]
#res1_set_repeated = [item for sublist in res1_set_repeated for item in sublist]
bp_mean_residues1 = pd.DataFrame(np.array(bp_residues1).T, columns=['BP'])
bp_mean_residues1['res'] = np.array(res1_set_repeated).T
bp_mean_residues1['BP'] = bp_mean_residues1.groupby('res', as_index=False)['BP'].transform('mean')
bp_mean_residues1 = bp_mean_residues1.drop_duplicates(subset='res', keep="first", inplace=False)
bp_mean_residues1['fragment'] = fragment1
bp_mean_residues1['cluster'] = cluster1
bp_mean_residues1['zernike'] = verso1
bp_mean_residues1['vs fragment'] = fragment2
bp_mean_residues1['vs cluster'] = cluster2
bp_mean_residues1['vs zernike'] = verso2


res2_set = []
bp2 = []
bp_residues2 = []
res2_set_repeated = []
for i in range(len(labs_2)):
    mask = label_2 == labs_2[i]
    res2_set_repeated.extend(res2[mask])
    res2_set.append(np.unique(res2[mask]))
    bp2.append(np.mean(original2[mask]))
    bp_residues2.extend(bp_points2[mask])
bp_mean_residues2 = pd.DataFrame(np.array(bp_residues2).T, columns=['BP'])
bp_mean_residues2['res'] = np.array(res2_set_repeated).T
bp_mean_residues2['BP'] = bp_mean_residues2.groupby('res', as_index=False)['BP'].transform('mean')
bp_mean_residues2 = bp_mean_residues2.drop_duplicates(subset='res', keep="first", inplace=False)
bp_mean_residues2['fragment'] = fragment2
bp_mean_residues2['cluster'] = cluster2
bp_mean_residues2['zernike'] = verso2
bp_mean_residues2['vs fragment'] = fragment1
bp_mean_residues2['vs cluster'] = cluster1
bp_mean_residues2['vs zernike'] = verso1



print("################### Prot 1  #################")
print(bp_mean_residues1)
for i in range(len(labs_1)):
    print("Residues patch %d:"%i)
    print(res1_set[i])
    print("")



print("################### Prot 2  ##################")
print(bp_mean_residues2)
for i in range(len(labs_2)):
    print("Residues patch %d:"%i)
    print(res2_set[i])
    print("")


bp_mean_residues1.to_csv('..\..\..\complementarity_regions\\{}.csv'.format(fragment1), mode='a', header=False)
bp_mean_residues2.to_csv('..\..\..\complementarity_regions\\{}.csv'.format(fragment2), mode='a', header=False)



#print('Choose the lowest mean value')
#print(bp1)
#print(bp2)