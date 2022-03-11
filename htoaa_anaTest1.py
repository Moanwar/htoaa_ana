#!/usr/bin/env python3


import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import pandas as pd
import xgboost as xgb
import pickle
import random
import os
from sklearn.model_selection import GridSearchCV
from pathlib import Path
import glob

from dataManager import processData, ggHPaths, BGenPaths, bEnrPaths, allVars, trainVars, disc, TTJetsPaths, WJetsPaths, ZJetsPaths

from htoaaRootFilesLoc1 import QCD_HT100to200_sample


from optparse import OptionParser
parser = OptionParser()
parser.add_option('--root' , dest='root', default=True)
(options, args) = parser.parse_args()



##############
## reload the data from ROOT, can just open pickle
## read determines if you want to load data from a root or pikl

root = options.root

## if you are doing dataVSMC using JetHT, set true
## maybe in the future, can read this from the where you read file names
#JetHT = True


if root==True:
    #print("Read data from root files", flush=True)
    print("Read data from root files")


    QCDData = pd.DataFrame()

    QCD_Samples = QCD_HT100to200_sample
    
    QCD_Paths = []
    for path in QCD_Samples:
        #print("path: {}".format(path))
        if "*" in path:
           QCD_Paths += glob.glob(path)
        else:
            QCD_Paths.append(path)
    print("QCD_Paths: {}".format(QCD_Paths))        
    
    for QCDDataPath in QCD_Paths:
    #for QCDDataPath in QCD_Paths[:5]:
        #print("\n\nhtoaa_anaTest1: QCDDataPath: ",QCDDataPath)
        QCDData = QCDData.append(processData(QCDDataPath, 'BGen', False), ignore_index=True, sort=False)
        #print("htoaa_anaTest1: QCDData: {}".format(QCDData))


QCDDf = QCDData
print("htoaa_anaTest1: QCDDf: {}".format(QCDDf))

dfdict = {
    'QCD': QCDDf,
}



cols = list(QCDDf.columns)
print("cols: ",cols)

print("htoaa_anaTest1: Siddh here 1")
## remove things i don't want plotted
toremove = ['final_weights', 'LHE_HT', 'QCD_corrections', 'lumi_weights',
            'PU_weights', 'LHE_weights']

print("htoaa_anaTest1: Siddh here 2")
for i in toremove:
    if i in cols: cols.remove(i)
print("htoaa_anaTest1: Siddh here 3")


for var in cols:
    if 'pt' in var:
        nbins = 40
    else:
        nbins = 20

    print("var",var)

    fig, (ax0, ax1) = plt.subplots(nrows=2, gridspec_kw={'height_ratios': [3, 1]},
                                   sharex=True)
    ax0.set_ylabel('events')
    ax1.set_ylabel('ratio')

    ## get the get range for histograms
    xmin = list()
    xmax = list()
    for dfkey, df in dfdict.items():
        xmintmp, xmaxtmp = np.percentile(df[var], [0,99.8])
        xmin.append(xmintmp)
        xmax.append(xmaxtmp)


    ## h41 and bdt score needs its own range tho
    if 'H4qvsQCD' in var:
        #range_local = (0, max(maxBGen, maxData))
        range_local = (0, max(xmax))
    elif 'BDTScore'==var:
        range_local = (0,0.8)
    else:
        #range_local = (min(minBGen, minData), max(maxBGen, maxData))
        range_local = (min(xmin), max(xmax))


    density = False
    hist_params = {'density': density, 'histtype': 'bar', 'range' : range_local, 'bins':nbins, 'stacked':True}

    ## prep backgroudn MC for plotting
    toplot = pd.DataFrame()
    toplotweights = pd.DataFrame()
    toplotlabel = list()
    for dfkey, df in dfdict.items():
        #toplot[dfkey] = df[var]
        #toplotweights[dfkey] = df['final_weights']
        toplot = pd.concat([toplot, df[var]], ignore_index=False, axis=1)
        toplotweights = pd.concat([toplotweights, df['final_weights']], ignore_index=False,
                  axis=1)

        toplotweights = toplotweights.rename({'final_weights': dfkey})
        toplotlabel.append(f'{dfkey} ({round(np.sum(df.final_weights))})')

    ## making color palette for the QCD stakcs
    pal = ['#603514', '#b940f2','#ec6f38', '#6acaf8', '#82f759']

    ## plotting background MC
    '''
    bgvals, bgbins, _ = ax0.hist(toplot.values.astype(float),
                                weights=toplotweights.values.astype(float),
                                label=toplotlabel,
                                color=pal, **hist_params)
    '''
    bgvals, bgbins, _ = ax0.hist(toplot.values.astype(float),
                                weights=toplotweights.values.astype(float),
                                label=toplotlabel,
                                color=pal[:1], **hist_params)


    ax0.set_title(var )# + ' JetHT')
    ax0.legend(loc='best', frameon=True)
    ax0.grid()

    filedest = 'dataVsMCDist_tmp/Parked/{}.png'
    plt.savefig(filedest.format(var))
    plt.show()
    plt.close()



print("htoaa_anaTest1 done")
