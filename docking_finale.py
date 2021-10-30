import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import cdist
from matplotlib import pyplot
import statistics
import seaborn as sb
from matplotlib.ticker import MaxNLocator

cutoff_int = 9.
cutoff_patch = 3.
n_model = 20 #20
pairings = 10  #10

b_35_res_B = [['PHE_229', 'IHS_256', 'ILE_257', 'SER_258'],
              ['null'],
              ['PHE_229', 'HIS_256', 'ILE_257'],
              ['PHE_229', 'PHE_231', 'HIS_256', 'ILE_257', 'SER_258'],
              ['PHE_229', 'PHE_231']
              ]
b_35_res_A = [['PHE_231'],
              ['null'],
              ['SER_254', 'VAL_255', 'IHS_256'],
              ['ALA_228', 'PHE_229', 'PHE_231', 'THR_233', 'SER_254', 'VAL_255', 'HIS_256', 'ILE_257',
               'SER_258'],
              ['THR_233', 'ILE_257']
              ]

b_35_res_B = [['PHE_229', 'IHS_256', 'ILE_257', 'SER_258'],
              ['null'],
              ['PHE_229', 'HIS_256', 'ILE_257'],
              ['PHE_229', 'PHE_231', 'HIS_256', 'ILE_257', 'SER_258'],
              ['PHE_229', 'PHE_231']
              ]
b_35_res_N_A = [[231],
                ['null'],
                [254, 255, 256],
                [228, 229, 231, 233, 254, 255, 256, 257,
                 258],
                [233, 257]
                ]
b_35_res_N_B = [[229, 256, 257, 258],
                ['null'],
                [229, 256, 257],
                [229, 231, 256, 257, 258],
                [229, 231]
                ]
N1 = [208, 208, 208, 208, 208, 208, 208, 208, 220, 208, 208, 208, 208, 208, 208]
C1 = [1, 1, 1, 1, 3, 1, 3, 4, 1, 3, 5, 4, 3, 4, 5]
N2 = [208, 208, 208, 208, 208, 220, 208, 208, 220, 220, 220, 220, 208, 208, 208]
C2 = [3, 1, 5, 4, 3, 1, 4, 4, 1, 1,1,1,5,5,5]

cut = [486, 486,486,486,486,486,486,486,393, 486, 486, 486, 486, 486, 486 ]

for i in range(0,pairings):
    num1 = N1[i]
    conf1 = C1[i]
    num2 = N2[i]
    conf2 = C2[i]

    if num1 == 208:
        namef1='A'
        b_35_conf1 = b_35_res_A[conf1-1]
        b_35_N_conf1 = b_35_res_N_A[conf1 - 1]
    else:
        namef1='B'
        b_35_conf1 = b_35_res_B[conf1 - 1]
        b_35_N_conf1 = b_35_res_N_B[conf1 - 1]
    if num2 == 208:
        namef2='A'
        b_35_conf2 = b_35_res_A[conf2 - 1]
        b_35_N_conf2 = b_35_res_N_A[conf2 - 1]
    else:
        namef2='B'
        b_35_conf2 = b_35_res_B[conf2 - 1]
        b_35_N_conf2 = b_35_res_N_B[conf2 - 1]
    name1 = '{}{}'.format(namef1,conf1)
    name2 = '{}{}'.format(namef2,conf2)

    path = '..\\..\\..\\HDOCK\\{}_{}'.format(name1,name2)

    BP1 = pd.DataFrame()
    BP1 = pd.read_csv('..\..\..\complementarity_regions\\WITHOUT_T_{}.csv'.format(num1))
    BP1 = BP1[(BP1['fragment'] == num1) & (BP1['cluster'] == conf1) & (BP1['vs fragment'] == num2) & (
            BP1['vs cluster'] == conf2)]
    BP1[['residue name', 'residue number']] = BP1['res'].str.split('_', 1, expand=True)
    BP1['BP'] = (BP1['BP']-BP1['BP'].mean())/BP1['BP'].std()
    #BP_1 = pd.DataFrame()
    #BP_1[['residue name', 'residue number']] = BP1['res'].str.split('_', 1, expand=True)
    #BP_1 = BP_1['residue number'].astype(int).to_list()

    BP2 = pd.DataFrame()
    BP2 = pd.read_csv('..\..\..\complementarity_regions\\WITHOUT_T_{}.csv'.format(num2))
    BP2 = BP2[(BP2['fragment'] == num2) & (BP2['cluster'] == conf2) & (BP2['vs fragment'] == num1) & (
            BP2['vs cluster'] == conf1)]
    BP2[['residue name', 'residue number']] = BP2['res'].str.split('_', 1, expand=True)
    BP2['BP'] = (BP2['BP']-BP2['BP'].mean())/BP2['BP'].std()
    #BP_2 = pd.DataFrame()
    #BP_2[['residue name', 'residue number']] = BP2['res'].str.split('_', 1, expand=True)
    #BP_2 = BP_2['residue number'].astype(int).to_list()

    b3_b5 = []
    for model in range(1, n_model + 1):
        ppdb = PandasPdb()
        ppdb.read_pdb('{}\model_{}.pdb'.format(path, model))
        ppdb.df['ATOM'].loc[:cut[i], 'chain_id'] = 1
        ppdb.df['ATOM'].loc[cut[i] + 1:, 'chain_id'] = 2
        CA = PandasPdb()
        CA = ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] == 'CA']

        dist = pd.DataFrame(
            cdist(CA[CA['chain_id'] == 1].iloc[:, 11:14], CA[CA['chain_id'] == 2].iloc[:, 11:14], metric='euclidean'))
        small = dist[dist < cutoff_int]

        res_1_index = small.loc[small.notna().any(axis=1)].index.to_list()  # INDICE DEI RESIDUI DI 1 (<-RIGHE) VICINI
        res_2_index = small.loc[:,
                      small.notna().any(axis=0)].columns.to_list()  # INDICE DEI RESIDUI DI 2 (<-COLONNE) VICINI
        CA_int_1 = CA[CA['chain_id'] == 1].iloc[res_1_index, :]  # RESIDUI DI 1 (<-RIGHE) VICINI
        CA_int_2 = CA[CA['chain_id'] == 2].iloc[res_2_index, :]  # RESIDUI DI 2 (<-COLONNE) VICINI

        dist_in_patch_1 = pd.DataFrame(
            cdist(ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == 1].iloc[:, 11:14], CA_int_1.iloc[:, 11:14],
                  metric='euclidean'))
        small_in_patch_1 = dist_in_patch_1[dist_in_patch_1 < cutoff_patch]
        res_patch_1_index = small_in_patch_1.loc[
            small_in_patch_1.notna().any(axis=1)].index.to_list()  # INDICE DEI RESIDUI NELLA PACTH 1
        PATCHES1 = ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == 1].iloc[res_patch_1_index, :]

        dist_in_patch_2 = pd.DataFrame(
            cdist(ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == 2].iloc[:, 11:14], CA_int_2.iloc[:, 11:14],
                  metric='euclidean'))
        small_in_patch_2 = dist_in_patch_2[dist_in_patch_2 < cutoff_patch]
        res_patch_2_index = small_in_patch_2.loc[
            small_in_patch_2.notna().any(axis=1)].index.to_list()  # INDICE DEI RESIDUI NELLA PACTH 2
        PATCHES2 = ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == 2].iloc[res_patch_2_index, :]

        #QUANTI RESIDUI B3 B5 COMPAIONO NEI LEGAMI CHE HO APPENA TROVATO
        count_1 = PATCHES1[PATCHES1['residue_number'].isin(b_35_N_conf1)]
        count_2 = PATCHES2[PATCHES2['residue_number'].isin(b_35_N_conf2)]
        b3_b5.append((count_1['residue_number'].nunique() + count_2['residue_number'].nunique())/(len(b_35_N_conf1)+len(b_35_N_conf2)))

        #BP DEI RESIDUI CHE COMPAIONO NEI LEGAMI CHE HO APPENA TROVATO E DI QUELLI NON INTERAGENTI
        BP1[BP1['residue number'].astype(int).isin(CA_int_1['residue_number'].astype(int).to_list())].to_csv('{}\{}_interactingBP_1.csv'.format(path,model))
        BP1[(~BP1['residue number'].astype(int).isin(CA_int_1['residue_number'].astype(int).to_list()))].to_csv('{}\\{}_non_interactingBP_1.csv'.format(path,model))
        BP2[BP2['residue number'].astype(int).isin(CA_int_2['residue_number'].astype(int).to_list())].to_csv('{}\{}_interactingBP_2.csv'.format(path,model))
        BP2[(~BP2['residue number'].astype(int).isin(CA_int_2['residue_number'].astype(int).to_list()))].to_csv('{}\\{}_non_interactingBP_2.csv'.format(path,model))

    j = 1
    selected_complexes = []
    colours = []
    for i in b3_b5:
        if i == max(b3_b5):
            selected_complexes.append(j)
            colours.append('tomato')
        else:
            colours.append('white')
        j += 1
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.set_axis_off()
    table = ax.table(
        cellText=[['{:.1%}'.format(i) for i in b3_b5]],
        fontsize=20,
        colLabels=["{}".format(i) for i in range(1,n_model+1)],
        colColours=["lightgray"] * n_model,
        cellColours=[colours],
        cellLoc='center',
        loc='upper left')
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    ax.set_title('{}-{}'.format(name1,name2), fontsize=20)
    plt.tight_layout()
    plt.show()

    for i in selected_complexes:
        print(b3_b5[i-1])
        complex1_interacting = pd.read_csv('{}\{}_interactingBP_1.csv'.format(path,i))
        complex1_NON_interacting = pd.read_csv('{}\{}_non_interactingBP_1.csv'.format(path,i))
        complex2_interacting = pd.read_csv('{}\{}_interactingBP_2.csv'.format(path, i))
        complex2_NON_interacting = pd.read_csv('{}\{}_non_interactingBP_2.csv'.format(path, i))

        interacting = complex1_interacting['BP'].to_list()+ complex2_interacting['BP'].to_list()
        non_interacting = complex1_NON_interacting['BP'].to_list()+ complex2_NON_interacting['BP'].to_list()
        sb.kdeplot(interacting, bw=.5, fill=True, color='powderblue')#,label='Interacting')
        sb.kdeplot(non_interacting, bw=0.5, fill=True, color='orange')#,label='Non interacting')
        #x, bins, p = plt.hist(interacting, density=True, bins='auto', color='powderblue',label='Interacting')  # bins = 'auto'
        #x, bins, p = plt.hist(non_interacting, density=True, bins='auto', color='orange',alpha=.4, label='Non interacting')  # bins = 'auto'
        plt.axvline(statistics.mean(interacting), color='b', linestyle='dashed', linewidth=2)#, label='Interacting mean Z-score')
        plt.axvline(statistics.mean(non_interacting), color='r', linestyle='dashed', linewidth=2)#, label='Non interacting mean Z-score')
        plt.title("Prediction {} for {}-{}".format(i, name1,name2),
                  fontsize=30)
        # plt.xlabel('BP', size=20)
        plt.ylabel('Density', size=30)
        plt.xlabel('Z-score', size=30)
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)
        plt.locator_params(axis='y', nbins=6)
        plt.locator_params(axis='x', nbins=6)
        #plt.legend(fontsize=20)
        plt.tight_layout()
        plt.show()