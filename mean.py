import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)


first_fram = [208]
second_fram = [208]
prima_conf = [4]
num = 4-1  # = PRIMA CONF -1
second_conf = [4]

o='y'
if o == 'y':
    #HO DISGUSTOSAMENTE COPIATO SU FILE I NUMERI UNO ALLA VOLTA, QUESTA COSA VA AUTOMATIZZATA
    b3b5 = np.loadtxt("..\..\..\senzaUNFOLDING\\b3b52.txt")
    b = np.loadtxt("..\..\..\senzaUNFOLDING\\b2.txt")
    b3b5 = b3b5
    b = b
    with open("..\..\..\senzaUNFOLDING\\pairs2.txt") as fh:
        pairs = fh.read().split('\n')


    results = pd.DataFrame()
    results['b3b5'] = b3b5
    results['b'] = b
    results['pairs'] = pairs
    results.sort_values('b3b5', inplace=True)

    b3b5 = results['b3b5'].to_numpy()
    b = results['b'].to_numpy()
    pairs = results['pairs']

    #fig = plt.figure(figsize=(6,4))
    fig, ax = plt.subplots()

    ind = np.arange(len(pairs))
    width = 0.15

    plt.bar(ind, b3b5, width, color='brown', label='$m_{\\beta 3, \\beta 5}$')
    plt.bar(ind + width, b, width,color='burlywood', label='$m_\\beta$')
    plt.xticks(ind + width / 2, results['pairs'],rotation =90)
    plt.title("Mean Binding Propensity", fontsize=30)
    plt.legend(fontsize=20)


    #ax1 = fig.add_subplot(111)
    #ax1.set_title("BP", fontsize=30)
    #ax1.set_xticklabels(results['pairs'],rotation =90)
    #ax1.set_xticks(range(len(pairs)))
    #ax1.plot(b3b5, c='powderblue',linewidth=3, label='$\\beta$3 and $\\beta$5 residues')
    #ax1.plot(b, c='brown', linewidth=3, label='$\\beta$ residues')
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    #leg = ax1.legend(fontsize=20)
    plt.tight_layout()
    plt.show()



for i in first_fram:
    if i == 220:
        frag = 'B'
        b_res = [['PHE_229', 'IHS_256', 'ILE_257', 'SER_258'],
                 ['null'],
                 ['PHE_229', 'HIS_256', 'ILE_257'],
                 ['PHE_229', 'PHE_231', 'HIS_256', 'ILE_257', 'SER_258'],
                 ['PHE_229', 'PHE_231']
        ]
        b_35_res = [['PHE_229', 'IHS_256', 'ILE_257', 'SER_258'],
                 ['null'],
                 ['PHE_229', 'HIS_256', 'ILE_257'],
                 ['PHE_229', 'PHE_231', 'HIS_256', 'ILE_257', 'SER_258'],
                 ['PHE_229', 'PHE_231']
                 ]
    else:
        frag = 'A'
        b_res = [['PHE_221', 'PHE_231'],
                 ['GLN_213', 'TYR_214', 'GLY_215', 'ASP_216', 'VAL_217', 'MET_218', "ASP_219", 'VAL_220', 'ILE_222'],
                 ['SER_254', 'VAL_255', 'IHS_256'],
                 ['PHE_221', 'ALA_228', 'PHE_229', 'PHE_231', 'THR_233', 'SER_254', 'VAL_255', 'HIS_256', 'ILE_257', 'SER_258'],
                 ['VAL_217', 'MET_218', 'ASP_219', 'VAL_220', 'PHE_221', "ILE_222", 'THR_233', 'PHE_234', 'ILE_257']
        ]
        b_35_res = [['PHE_231'],
                 ['null'],
                 ['SER_254', 'VAL_255', 'IHS_256'],
                 ['ALA_228', 'PHE_229', 'PHE_231', 'THR_233', 'SER_254', 'VAL_255', 'HIS_256', 'ILE_257',
                  'SER_258'],
                 ['THR_233', 'ILE_257']
                 ]

    file = '..\..\..\complementarity_regions\\WITHOUT_T_{}.csv'.format(i)
    fragment_results = pd.DataFrame()
    fragment_results['res'] = pd.read_csv(file, usecols=['res'], squeeze=True)
    fragment_results['BP'] = pd.read_csv(file, usecols=['BP'], squeeze=True)
    #disattivare la prossima riga se non vuoi le BP riscalate
    fragment_results['BP'] = (fragment_results['BP']-fragment_results['BP'].mean())/fragment_results['BP'].std()
    fragment_results['cluster'] = pd.read_csv(file, usecols=['cluster'])
    fragment_results['fragmentvs'] = pd.read_csv(file, usecols=['vs fragment'])
    fragment_results['clustervs'] = pd.read_csv(file, usecols=['vs cluster'])

    for j in prima_conf:
        pdb_file = "..\\{}\cluster{}.dms".format(i, j)
        cluster_results = pd.DataFrame()
        cluster_results = fragment_results[fragment_results.cluster.isin([j])]

        n = []
        for ii in second_fram:
            if ii == 220:
                fragvs = 'B'
            else:
                fragvs = 'A'
            final_fragment_results = pd.DataFrame()
            final_fragment_results = cluster_results[cluster_results.fragmentvs.isin([ii])]
            for jj in second_conf:
                comp_region = pd.DataFrame()
                comp_region['Residues name and specifier'] = pd.read_csv(pdb_file, usecols=['Res'], squeeze=True)
                #comp_region['beta'] = comp_region['Residues name and specifier'].apply(lambda x: any([k in x for k in b_res[num]]))
                comp_region = comp_region.drop_duplicates(subset='Residues name and specifier', keep="first", inplace=False).set_index('Residues name and specifier')
                comp_region['BP'] = 0

                final_final_fragment_results = pd.DataFrame()
                final_final_fragment_results = final_fragment_results[final_fragment_results.clustervs.isin([jj])]
                        ### attiva se vuoi vedere solo beta
                total_final_final_fragment_results = final_final_fragment_results[final_final_fragment_results.res.isin([k for k in b_res[num]])]

                total_beta_mean = total_final_final_fragment_results['BP'].mean()
                print('{}_{} vs {}_{} mean TOTAL beta = %.3f'.format(frag,j, fragvs,jj) %total_beta_mean)


                beta_results = pd.DataFrame()
                beta_results = final_final_fragment_results[final_final_fragment_results.res.isin([k for k in b_35_res[num]])].drop_duplicates(subset='res', keep="first", inplace=False).set_index('res')
                beta_mean = beta_results['BP'].mean()
                print('{}_{} vs {}_{} mean beta 3 and 5 =%.3f'.format(frag,j, fragvs,jj) %beta_mean)

                all_res = final_final_fragment_results['res'].values.tolist()
                non_beta_results = pd.DataFrame()
                non_beta_results = final_final_fragment_results[final_final_fragment_results.res.isin([k for k in all_res if k not in b_res[num]])].drop_duplicates(subset='res', keep="first",
                                                                                                    inplace=False).set_index('res')
                non_beta_mean = non_beta_results['BP'].mean()
                print('{}_{} vs {}_{} mean NON beta = %.3f'.format(frag, j, fragvs, jj) % non_beta_mean)


                res = list(final_final_fragment_results['res'])
                final_final_fragment_results = final_final_fragment_results.drop_duplicates(subset='res', keep="first", inplace=False).set_index('res')

                for r in res:
                    comp_region.loc[r, 'BP'] = final_final_fragment_results.loc[r, 'BP']

                #comp_region['beta'] = np.where(comp_region.index.isin([k for k in b_res[num]]), 'blue','red' )
                #print(comp_region)
                x, bins, p = plt.hist(comp_region.BP, density=False, bins='auto', color='powderblue') #bins = 'auto'
                plt.title("BP values of all {} residues ({}{} vs {}{})".format(len(comp_region.BP),frag, j, fragvs, jj), fontsize=30)
                #plt.xlabel('BP', size=20)
                plt.ylabel('Counts', size=20)
                plt.yticks(fontsize=15)
                plt.xticks(fontsize=15)
                plt.show()


                comp_region.BP.rolling(5).mean()

                #SE USI BP ORGINALE, CON TUTTI I VALORI POSTIVI
                #ax = comp_region.plot.area(color='powderblue', legend=False)
                #ax.set_yticklabels([])
                #ax.set_yticks([])
                #ax.yaxis.set_tick_params(labelsize=30)
                #ax.xaxis.set_tick_params(labelsize=30)
                #ax.xaxis.label.set_visible(False)
                #ax.yaxis.label.set_visible(False)
                #ax.set_xlabel("Residues", fontsize=60)
                #ax.set_ylabel("BP", fontsize=60)
                #ax.text(1.5, 0.2, 'BP={}'.format((round(beta_mean, 3))), color='black',
                #        bbox=dict(facecolor='white', edgecolor='black'), fontsize=85)
                #plt.title('{}{} vs {}{}'.format(frag,j,fragvs,jj), fontsize=85)
                #n.append(len(final_final_fragment_results))

                #SE INVECE RISCALI E HAI VALORI POSITIVI E NEGATIVI
                fig, ax = plt.subplots()
                set_size(5, 2)
                ax.set_yticklabels([])
                ax.set_yticks([])
                ax.set_xticklabels([])
                ax.set_xticks([])
                ax.xaxis.label.set_visible(False)
                # split dataframe df into negative only and positive only values
                df_neg, df_pos = comp_region.clip(upper=0), comp_region.clip(lower=0)
                # stacked area plot of positive values
                df_pos.plot.area(ax=ax, stacked=True, linewidth=0., color='powderblue', legend=False)
                # reset the color cycle
                ax.set_prop_cycle(None)
                # stacked area plot of negative values, prepend column names with '_' such that they don't appear in the legend
                df_neg.rename(columns=lambda x: '_' + x).plot.area(ax=ax, stacked=True, linewidth=0., color='powderblue', legend=False)
                # rescale the y axis
                ax.set_ylim([df_neg.sum(axis=1).min(), df_pos.sum(axis=1).max()])
                #ax.text(.0, -1.8, '$m_{\\beta 3, \\beta 5}=0.887$', color='black', #.format((round(beta_mean, 3)))
                #bbox=dict(facecolor='white', edgecolor='cornflowerblue', linewidth=4), fontsize=30)
                plt.title('{}{} vs {}{}'.format(frag,j,fragvs,jj), fontsize=30)
        plt.show()





