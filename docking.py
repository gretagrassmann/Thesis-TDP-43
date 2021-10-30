import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import cdist
from matplotlib import pyplot

cutoff_int = 9.
cutoff_patch = 3.
n_model = 20

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

for i in range(0,10): #0,10
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

    print('Vuoi fare i calcoli per {}_{}? y o n?'.format(name1,name2))
    o=input()
    if o == 'y':
        BP1 = pd.DataFrame()
        BP1 = pd.read_csv('..\..\..\complementarity_regions\\{}.csv'.format(num1))
        BP1 = BP1[(BP1['fragment']==num1) & (BP1['cluster']==conf1) & (BP1['vs fragment']==num2) & (BP1['vs cluster']==conf2)]
        BP_1 = pd.DataFrame()
        BP_1[['residue name', 'residue number']] = BP1['res'].str.split('_', 1, expand=True)
        BP_1= BP_1['residue number'].astype(int).to_list()
        BP2 = pd.DataFrame()
        BP2 = pd.read_csv('..\..\..\complementarity_regions\\{}.csv'.format(num2))
        BP2 = BP2[(BP2['fragment']==num2) & (BP2['cluster']==conf2) & (BP2['vs fragment']==num1) & (BP2['vs cluster']==conf1)]
        BP_2 = pd.DataFrame()
        BP_2[['residue name', 'residue number']] = BP2['res'].str.split('_', 1, expand=True)
        BP_2 = BP_2['residue number'].astype(int).to_list()

        found_residues_1 = []
        found_residues_2 = []
        found_count_BP1 = []
        found_count_BP2 = []
        score = []
        for model in range(1,n_model+1):
            ppdb = PandasPdb()
            ppdb.read_pdb('{}\model_{}.pdb'.format(path,model))
            ppdb.df['ATOM'].loc[:cut[i],'chain_id'] = 1
            ppdb.df['ATOM'].loc[cut[i]+1:,'chain_id'] = 2
            score.append(ppdb.df['OTHERS'].loc[3,'entry'][8:])
            CA = PandasPdb()
            CA= ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] == 'CA'] #.head()

            dist = pd.DataFrame(cdist(CA[CA['chain_id']==1].iloc[:,11:14], CA[CA['chain_id']==2].iloc[:,11:14], metric='euclidean'))
            small = dist[dist < cutoff_int]

            res_1_index = small.loc[small.notna().any(axis=1)].index.to_list()#INDICE DEI RESIDUI DI 1 (<-RIGHE) VICINI
            res_2_index = small.loc[:,small.notna().any(axis=0)].columns.to_list()#INDICE DEI RESIDUI DI 2 (<-COLONNE) VICINI
            CA_int_1 = CA[CA['chain_id']==1].iloc[res_1_index,:] #RESIDUI DI 1 (<-RIGHE) VICINI
            CA_int_2 = CA[CA['chain_id']==2].iloc[res_2_index,:] #RESIDUI DI 2 (<-COLONNE) VICINI

            dist_in_patch_1 = pd.DataFrame(cdist(ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == 1].iloc[:,11:14], CA_int_1.iloc[:,11:14], metric='euclidean'))
            small_in_patch_1 = dist_in_patch_1[dist_in_patch_1 < cutoff_patch]
            res_patch_1_index = small_in_patch_1.loc[small_in_patch_1.notna().any(axis=1)].index.to_list()#INDICE DEI RESIDUI NELLA PACTH 1
            PATCHES1 = ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == 1].iloc[res_patch_1_index,:]

            dist_in_patch_2 = pd.DataFrame(cdist(ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == 2].iloc[:,11:14], CA_int_2.iloc[:,11:14], metric='euclidean'))
            small_in_patch_2 = dist_in_patch_2[dist_in_patch_2 < cutoff_patch]
            res_patch_2_index = small_in_patch_2.loc[small_in_patch_2.notna().any(axis=1)].index.to_list()#INDICE DEI RESIDUI NELLA PACTH 2
            PATCHES2 = ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == 2].iloc[res_patch_2_index,:]

            count_1 = PATCHES1[PATCHES1['residue_number'].isin(b_35_N_conf1)]
            found_residues_1.append(count_1['residue_number'].nunique())
            count_2 = PATCHES2[PATCHES2['residue_number'].isin(b_35_N_conf2)]
            found_residues_2.append(count_2['residue_number'].nunique())

            count_BP_1 = PATCHES1[PATCHES1['residue_number'].isin(BP_1)]
            found_count_BP1.append('{:.2}'.format(count_BP_1['residue_number'].nunique()/PATCHES1['residue_number'].nunique()))

            count_BP_2 = PATCHES2[PATCHES2['residue_number'].isin(BP_2)]
            found_count_BP2.append('{:.2}'.format(count_BP_2['residue_number'].nunique()/PATCHES2['residue_number'].nunique()))

        final = pd.DataFrame({'BP_1': found_count_BP1,
                           'BP_2': found_count_BP2,
                              'b3b5_1': found_residues_1,
                              'b3b5_2': found_residues_2,
                              'SCORE': score}
                             )
        final.to_csv("{}\{}_{}.csv".format(path,name1,name2), mode = 'w')

    file = "{}\{}_{}.csv".format(path,name1,name2)
    results = pd.DataFrame()
    results = pd.read_csv(file)
    found_residues_1 = results['b3b5_1'].to_list()
    found_residues_2 = results['b3b5_2'].to_list()

    found_count_BP1 =results['BP_1'].to_list()
    found_count_BP2 = results['BP_2'].to_list()

    x = range(n_model)
    y1=range(1,len(b_35_N_conf1)+1)
    y2=range(len(b_35_N_conf2)+1)
    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True)
    #fig.suptitle('Number of interacting $\\beta 3$, $\\beta 5$ residues',fontsize=30)
    plt.subplot(2,1,1)
    pop = plt.bar(x,found_residues_1,color='powderblue')
    plt.ylabel('{}'.format(name1),fontsize=20)
    plt.yticks(y1,fontsize=20)
    plt.xticks(x,fontsize=20)
    plt.subplot(2,1,2)
    plt.yticks(y2,fontsize=20)
    gdp =plt.bar(x,found_residues_2, color='cornflowerblue')
    plt.ylabel('{}'.format(name2),fontsize=20)
    #plt.xlabel('Top {} complexes'.format(n_model),fontsize=20)
    for ax in axs:
        ax.label_outer()
    plt.xticks([])
    plt.tight_layout()
    #plt.show()


    #fig = plt.figure()
    #gs = fig.add_gridspec(2, hspace=0)
    #axs = gs.subplots(sharex=True)
    #fig.suptitle('Percentage of $\\beta 3$, $\\beta 5$ residues with high BP',fontsize=30)
    #plt.subplot(2,1,1)
    #pop = plt.bar(x,found_count_BP1,color='powderblue')
    #plt.ylabel('{}'.format(name1),fontsize=20)
    #plt.xticks(x,fontsize=20)
    #The below code will create the second plot.
    #plt.subplot(2,1,2)
    #gdp =plt.bar(x,found_count_BP2, color='cornflowerblue')
    #plt.ylabel('{}'.format(name2),fontsize=20)
    #plt.xlabel('Top {} complexes'.format(n_model),fontsize=20)
    # Hide x labels and tick labels for all but bottom plot.
    #for ax in axs:
    #    ax.label_outer()
    #plt.tight_layout()
    #plt.show()



    df=pd.DataFrame((found_count_BP1,found_count_BP2))
    vals = np.around(df.values,2)
    normal = plt.Normalize(vals.min(), vals.max()+0.1)
    colours = plt.cm.Blues(normal(vals))
    fig, ax = plt.subplots(figsize=(10,2))
    ax.set_axis_off()
    table = ax.table(
        cellText=[['{:.1%}'.format(i) for i in found_count_BP1],['{:.1%}'.format(i) for i in found_count_BP2]],fontsize=20,
        rowLabels=[name1,name2],
        colLabels=["{}".format(i) for i in range(n_model)],
        rowColours=["yellowgreen"] * 2,
        colColours=["lightgray"] * n_model,
        cellColours=colours,
        cellLoc='center',
        loc='upper left')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    ax.set_title('Interacting residues with high BP for the top {} predictions'.format(n_model),fontweight="bold")
    plt.tight_layout()
    #plt.show()

    df=pd.DataFrame((found_residues_1,found_residues_2))
    vals = np.around(df.values,2)
    normal = plt.Normalize(vals.min()-1, vals.max()+1)
    colours = plt.cm.Oranges(normal(vals))
    fig, ax = plt.subplots(figsize=(10,2))
    ax.set_axis_off()
    table = ax.table(
        cellText=[['{:.1%}'.format(i/len(b_35_N_conf1)) for i in found_residues_1],['{:.1%}'.format(i/len(b_35_N_conf2)) for i in found_residues_2]],fontsize=20,
        rowLabels=[name1,name2],
        #colLabels=["{}".format(i) for i in range(n_model)],
        rowColours=["yellowgreen"] * 4,
        #colColours=["lightgray"] * n_model,
        cellColours=colours,
        cellLoc='center',
        loc='upper left')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    #ax.set_title('Interacting $\\beta 3, \\beta 5$ residues for the top {} predictions'.format(n_model),fontweight="bold")
    plt.tight_layout()
    #plt.show()

    mean_BP = [(x+y)/2 for x, y in zip(found_count_BP1,found_count_BP2)]
    mean_b3b5 = [((x/len(b_35_N_conf1))+(y/len(b_35_N_conf2)))/2 for x, y in zip(found_residues_1,found_residues_2)]
    order_BP = sorted(range(len(mean_BP)), key=lambda k: mean_BP[k])
    order_b3b5 = sorted(range(len(mean_b3b5)), key=lambda k: mean_b3b5[k])

    fig = plt.figure()
    ind = np.arange(n_model)
    plt.plot(ind, order_BP, color='cornflowerblue', label='High BP residues')
    plt.plot(ind, order_b3b5,color='orange', label='$\\beta 3, \\beta 5$ residues')
    plt.xticks([])
    plt.yticks(np.arange(0,n_model,2), fontsize=20)
    plt.ylabel('Complex number',fontsize=20)
    plt.title("Complexes ordered according to the % of interacting residues ", fontsize=20)
    plt.legend(fontsize=15)
    plt.show()