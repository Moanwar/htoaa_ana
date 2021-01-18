import pickle
import numpy as np
import pandas as pd
from analib import PhysObj, Event
from info import BGenFileNames, bEnrFileNames

import dataVsMC_DataManager as DM #import processData, jetVars, muonVars, PVVars, allVars, dataPath, ggHPath, bEnrPaths, BGenPaths, TTJetsPaths,
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import mplhep as hep
import uproot

## get that sweet CMS style plots
plt.style.use(hep.style.CMS)

## monte carlo ggH signal
#ggHDf = processData(ggHPath, 'ggH')

## monte carlo backgrounds
#BGenDf = processData(BGenPath, 'BGen')
## uncomment below if using multiple bg MC files

root = True

if root:
    BGenDf = pd.DataFrame()
    for fileName in DM.BGenPaths:
        tmpDf = DM.processData(fileName, 'BGen')
        BGenDf = BGenDf.append(tmpDf, ignore_index=True, sort=False)
    
    #bEnrDf = processData(bEnrPath, 'bEnr')
    ## uncommend below if using multiple bg MC files
    bEnrDf = pd.DataFrame()
    for fileName in DM.bEnrPaths:
        tmpDf = DM.processData(fileName, 'bEnr')
        bEnrDf = bEnrDf.append(tmpDf, ignore_index=True, sort=False)
    
    TTJetsDf = pd.DataFrame()
    for fileName in DM.TTJetsPaths:
        tmpdf = DM.processData(fileName, 'TTJets')
        TTJetsDf = TTJetsDf.append(tmpdf, ignore_index=True, sort=False)
    
    ZJetsDf = pd.DataFrame()
    for fileName in DM.ZJetsPaths:
        tmpdf = DM.processData(fileName, 'ZJets')
        ZJetsDf = ZJetsDf.append(tmpdf, ignore_index=True, sort=False)
    
    WJetsDf = pd.DataFrame()
    for fileName in DM.WJetsPaths:
        tmpdf = DM.processData(fileName, 'WJets')
        WJetsDf = WJetsDf.append(tmpdf, ignore_index=True, sort=False)
    
    ## datafile
    if not DM.JetHT:
        dataDf = DM.processData(DM.dataPath, 'data')
    
    if DM.JetHT:
        JetHTDf = pd.DataFrame()
        for fileName in DM.JetHTPaths:
            tmpdf = DM.processData(fileName, 'JetHT')
            JetHTDf = JetHTDf.append(tmpdf, ignore_index=True, sort=False)
    
    
    
    
    #pickle.dump(ggHDf, open('dataVsMC/ggHDf.pkl', 'wb'))
    pickle.dump(BGenDf, open('dataVsMC/BGenDf.pkl', 'wb'))
    pickle.dump(bEnrDf, open('dataVsMC/bEnrDf.pkl', 'wb'))
    #pickle.dump(dataDf, open('dataVsMC/dataDf.pkl', 'wb'))
    pickle.dump(TTJetsDf, open('dataVsMC/TTJetsDf.pkl', 'wb'))
    pickle.dump(ZJetsDf, open('dataVsMC/ZJetsDf.pkl', 'wb'))
    pickle.dump(WJetsDf, open('dataVsMC/WJetsDf.pkl', 'wb'))
    pickle.dump(JetHTDf, open('dataVsMC/JetHTDf.pkl', 'wb'))

else:
    #ggHDf = pickle.load(open('dataVsMC/ggHDf.pkl', 'rb'))
    BGenDf = pickle.load(open('dataVsMC/BGenDf.pkl', 'rb'))
    bEnrDf = pickle.load(open('dataVsMC/bEnrDf.pkl', 'rb'))
    #dataDf = pickle.load(open('dataVsMC/dataDf.pkl', 'rb'))
    TTJetsDf = pickle.load(open('dataVsMC/TTJetsDf.pkl', 'rb'))
    ZJetsDf = pickle.load(open('dataVsMC/ZJetsDf.pkl', 'rb'))
    WJetsDf = pickle.load(open('dataVsMC/WJetsDf.pkl', 'rb'))
    JetHTDf = pickle.load(open('dataVsMC/JetHTDf.pkl', 'rb'))


dfdict = {'BGenDf': BGenDf,
          'bEnrDf': bEnrDf,
          'TTJetsDf': TTJetsDf,
          'ZJetsDf': ZJetsDf,
          'WJetsDf': WJetsDf,
          'JetHTDf': JetHTDf}

print('BGenDf ', BGenDf.columns)
print('bEnrDf', bEnrDf.columns)
print('TTJetsDf', TTJetsDf.columns)
print('ZJetsDf', ZJetsDf.columns)
print('WJetsDf', WJetsDf.columns)
print('JetHTDf', JetHTDf.columns)

for var in JetHTDf.columns:
    if 'pt' in var:
        nbins = 80
    else:
        nbins = 40

    fig, ax = plt.subplots(figsize=(10,8))

    ## get the min/max value for hists
    xmin = list()
    xmax = list()

    ## old way. when only doing 1 background
    #minBGen, maxBGen = np.percentile(BGenDf[var], [0, 99.8])
    #minData, maxData = np.percentile(dataDf[var], [0, 99.8])

    for dfkey, df in dfdict.items():
        xmintmp, xmaxtmp = np.percentile(df[var], [0,99.8])
        xmin.append(xmintmp)
        xmax.append(xmaxtmp)

    if 'H4qvsQCD' in var:
        #range_local = (0, max(maxBGen, maxData))
        range_local = (0, max(xmax))
    else:
        #range_local = (min(minBGen, minData), max(maxBGen, maxData))
        range_local = (min(xmin), max(xmax))


    hist_params = {'density': True, 'histtype': 'barstacked', 'range' : range_local, 'bins':nbins}


    toplot = pd.DataFrame()
    toplotweights = pd.DataFrame()
    toplotlabel = list()
    JetHTDf['final_weights'] = 1
    for dfkey, df in dfdict.items():
        #toplot = np.append(toplot, df[var], 1)
        #toplotweights = np.append(toplotweights, df['final_weights'], 1)
        toplot[dfkey] = df[var]
        toplotweights[dfkey] = df['final_weights']
        toplotlabel.append(dfkey)

    ax.hist(toplot.values, weights=toplotweights.values,label=toplotlabel, **hist_params)

    #ax.hist(BGenDf[var], label='BGen', weights=BGenDf['final_weights'],
            #**hist_params)
    #ax.hist(dataDf[var], label='Data',
    ax.set_title(var + ' JetHT')
    ax.legend(loc='best', frameon=True)
    ax.grid()
    plt.savefig('dataVsMCDist/JetHT/{}.png'.format(var))
    plt.show()


    ## when used to do bEnr and BGen separately
    # fig, ax = plt.subplots(figsize=(10,8))
    # minbEnr, maxbEnr = np.percentile(bEnrDf[var], [0, 99.8])
    # if 'H4qvsQCD' in var:
    #     range_local = (0, max(maxbEnr, maxData))
    # else:
    #     range_local = (min(minbEnr, minData), max(maxbEnr, maxData))
    # ax.hist(bEnrDf[var], bins=nbins, label='bEnr', color='b', range=range_local,
    #         weights=bEnrDf['final_weights'], histtype=chart, density=d)
    # ax.hist(dataDf[var], bins=nbins, label='Data', color='r', range = range_local,
    #         histtype=chart, density=d)
    # ax.legend(loc='best', frameon=True)
    # ax.set_title(var + ' bEnr vs Data')
    # ax.grid()
    # plt.savefig('dataVsMCDist_fixed/bEnr/{}_bEnrVsData.png'.format(var))
    # plt.close()

# nbins=80
'''fig, ax = plt.subplots(figsize=(10,8))

del dfdict['JetHTDf']
del dfdict['TTJetsDf']

for dfkey, dfvalue in dfdict.items():
    range_local = np.percentile(dfvalue['LHE_HT'], [0,99.8])
    ax.hist(dfvalue['LHE_HT'], weights=dfvalue['LHE_weights'], label=dfkey,
            **hist_params)
    ax.set_title('LHE_HT ' + dfkey)
    ax.grid()
    plt.savefig(f'dataVsMCDist/JetHT/LHE_HT_{dfkey}.png')
    plt.show()'''



## old way when bgen benr plotted separate
# minbEnr, maxbEnr = np.percentile(bEnrDf['LHE_HT'], [0, 99.8])
# minBGen, maxBGen = np.percentile(BGenDf['LHE_HT'], [0, 99.8])
# range_local = (min(minbEnr, minBGen), max(maxbEnr, maxBGen))
# ax.hist(bEnrDf['LHE_HT'], weights=bEnrDf['LHE_weights'], label=bEnr, **hist_params)
# ax.hist(BGenDf['LHE_HT'], weights=BGenDf['LHE_weights'], label=bEnr, **hist_params)
# ax.legend(loc='best', frameon=True)
# ax.set_title( 'LHE_HT bEnr vs BGen')
# ax.grid()
# plt.savefig('dataVsMCDist/JetHT/LHE_HT.png')
# plt.close()





