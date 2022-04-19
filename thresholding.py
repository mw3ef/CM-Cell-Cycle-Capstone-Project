# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 16:30:40 2022

@author: Michelle Wu
"""
def dapi_and_ki67_analysis ():
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    
    ## pre-processing
    # load dapi and ki67 data file
    df = pd.read_excel (r'C:\Users\mwu\OneDrive\raw dapi data.xlsx')
    
    # split data according to the experimental conditions and convert them to numpy
    dmso = df.iloc[2:1533, [0, 5]].to_numpy()
    dmso_n = 1531 # total number of dmso cells
    
    nrg = df.iloc[2:1724, [1, 6]].to_numpy()
    nrg_n = 1722 # total number of nrg cells
    
    nrg_p38 = df.iloc[2:1689, [2, 7]].to_numpy()
    nrg_p38_n = 1687 # total number of nrg + p38 cells
    
    tt10 = df.iloc[2:1774, [3, 8]].to_numpy()
    tt10_n = 1772 # total number of tt10 cells
    
    # scatter plot 
    
    
    ## dmso condition
    # initialize the counts for each phase for dmso
    dmso_G0 = 0
    dmso_G1 = 0
    dmso_S = 0
    dmso_G2M = 0
    # dmso_aneu = 0
    # dmso_G04c = 0
    
    # loop through the data and threshold to get cell phase counts 
    for i in range(0, len(dmso)):
        dapi = dmso[i, 0]
        ki67 = dmso[i, 1]
        if(dapi > 2597558 and dapi < 3794223 and ki67 < 3.141):
            dmso_G0 += 1
        if(dapi > 2597558 and dapi < 3794223 and ki67 > 3.141):
            dmso_G1 += 1
        if(dapi > 3794223 and dapi < 6459892 and ki67 > 3.141):
            dmso_S += 1
        if(dapi > 6459892 and dapi < 6927843 and ki67 > 3.141):
            dmso_G2M += 1
        # if(dapi > 3794223 and dapi < 6459892 and ki67 < 3.141):
        #     dmso_aneu += 1
        # if(dapi > 6459892 and dapi < 6927843 and ki67 < 3.141):
        #     dmso_G04c += 1
    
    # made the final cell phase counts into an array
    dmso_height = np.array([dmso_G0, dmso_G1, dmso_S, dmso_G2M])
    dmso_height = (dmso_height/dmso_n)*100 # converting to percentages
    # print("dmso", dmso_height)
    
    ## nrg
    # initialize the counts for each phase for nrg
    nrg_G0 = 0
    nrg_G1 = 0
    nrg_S = 0
    nrg_G2M = 0
    # nrg_aneu = 0
    # nrg_G04c = 0
    
    # loop through the data and threshold to get cell phase counts
    for j in range(0, len(nrg)):
        dapi = nrg[j, 0]
        ki67 = nrg[j, 1]
        if(dapi > 2848206 and dapi < 3887269 and ki67 < 3.141):
            nrg_G0 += 1
        if(dapi > 2848206 and dapi < 3887269 and ki67 > 3.141):
            nrg_G1 += 1
        if(dapi > 3887269 and dapi < 6786273 and ki67 > 3.141):
            nrg_S += 1
        if(dapi > 6786273 and dapi < 7268209 and ki67 > 3.141):
            nrg_G2M += 1
        # if(dapi > 3887269 and dapi < 6786273 and ki67 < 3.141):
        #     nrg_aneu += 1
        # if(dapi > 6786273 and dapi < 7268209 and ki67 < 3.141):
        #     nrg_G04c += 1
            
    # made the final cell phase counts into an array        
    nrg_height = np.array([nrg_G0, nrg_G1, nrg_S, nrg_G2M])
    nrg_height = (nrg_height/nrg_n)*100 # converting to percentages
    # print("nrg", nrg_height)
    
    ## nrg + p38
    # initialize the counts for each phase for nrg + p38
    nrg_p38_G0 = 0
    nrg_p38_G1 = 0
    nrg_p38_S = 0
    nrg_p38_G2M = 0
    # nrg_p38_aneu = 0
    # nrg_p38_G04c = 0
    
    # loop through the data and threshold to get cell phase counts
    for k in range(0, len(nrg_p38)):
        dapi = nrg_p38[k, 0]
        ki67 = nrg_p38[k, 1]
        if(dapi > 2902792 and dapi < 3809748 and ki67 < 3.141):
            nrg_p38_G0 += 1
        if(dapi > 2902792 and dapi < 3809748 and ki67 > 3.141):
            nrg_p38_G1 += 1
        if(dapi > 3809748 and dapi < 6768766 and ki67 > 3.141):
            nrg_p38_S += 1
        if(dapi > 6768766 and dapi < 7220127 and ki67 > 3.141):
            nrg_p38_G2M += 1
        # if(dapi > 3809748 and dapi < 6768766 and ki67 < 3.141):
        #     nrg_p38_aneu += 1
        # if(dapi > 6768766 and dapi < 7220127 and ki67 < 3.141):
        #     nrg_p38_G04c += 1
    
    # made the final cell phase counts into an array
    nrg_p38_height = np.array([nrg_p38_G0, nrg_p38_G1, nrg_p38_S, nrg_p38_G2M])
    nrg_p38_height = (nrg_p38_height/nrg_p38_n)*100 # converting to percentages
    # print("nrg + p38", nrg_p38_height)
    
    ## tt10
    # initialize the counts for each phase for tt10
    tt10_G0 = 0
    tt10_G1 = 0
    tt10_S = 0
    tt10_G2M = 0
    # tt10_aneu = 0
    # tt10_G04c = 0
    
    # loop through the data and threshold to get cell phase counts
    for m in range(0, len(tt10)):
        dapi = tt10[m, 0]
        ki67 = tt10[m, 1]
        if(dapi > 2574537 and dapi < 4003839 and ki67 < 3.141):
            tt10_G0 += 1
        if(dapi > 2574537 and dapi < 4003839 and ki67 > 3.141):
            tt10_G1 += 1
        if(dapi > 4003839 and dapi < 6772490 and ki67 > 3.141):
            tt10_S += 1
        if(dapi > 6772490 and dapi < 6981477 and ki67 > 3.141):
            tt10_G2M += 1
        # if(dapi > 4003839 and dapi < 6772490 and ki67 < 3.141):
        #     tt10_aneu += 1
        # if(dapi > 6772490 and dapi < 6981477 and ki67 < 3.141):
        #     tt10_G04c += 1
            
    # made the final cell phase counts into an array
    tt10_height = np.array([tt10_G0, tt10_G1, tt10_S, tt10_G2M])
    tt10_height = (tt10_height/tt10_n)*100 # converting to percentages
    # print("tt10", tt10_height)
    
    # plotting bar plots
    stages = ["G0", "G1", "S", "G2/M"] # bar graph labels
    x_axis = np.arange(len(stages))
    
    w = 0.2 # bar widths
    plt.bar(x_axis -1.5*w, dmso_height, width=w, label = 'dmso', color = '#d7191c')
    plt.bar(x_axis -w/2, nrg_height, width=w, label = 'nrg', color = '#fdae61')
    plt.bar(x_axis +w/2, nrg_p38_height, width=w, label = 'nrg+p38', color = '#CF9FFF')
    plt.bar(x_axis +1.5*w, tt10_height, width=w, label = 'tt10', color = '#2c7bb6')
    
    plt.xticks(x_axis, stages)
    plt.legend()
    plt.title("Cell Counts in Each Stage Using DAPI and Ki67 Data")
    plt.xlabel("Cell Cycle Stages")
    plt.ylabel("Cell Counts Percentage")
    
    plt.show()
    
    return dmso_height, nrg_height, nrg_p38_height, tt10_height
    



