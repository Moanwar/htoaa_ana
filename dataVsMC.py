import pickle
import numpy as np
import pandas as pd
from analib import PhysObj, Event
from info import BGenFileNames, bEnrFileNames

from data_managerDataVsMC import processData, BGenWeight, bEnrWeight, trainVars, extraVars, allVars, plotVars, dataPath, ggHPath, bEnrPath, BGenPath
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import mplhep as hep
import uproot

## get that sweet CMS style plots
plt.style.use(hep.style.CMS)


## uncomment when upload to kodiak

## monte carlo ggH signal
#ggHDf = processData(ggHPath, 'ggH')

## monte carlo backgrounds
#BGenDf = pd.DataFrame()
BGenDf = processData(BGenPath, 'BGen')
#for fileName in BGenFileNames:
#    tmpDf = processData(fileName, 'BGen')
#    BGenDf = BGenDf.append(tmpDf, ignore_index=True, sort=False)

#bEnrDf = pd.DataFrame()
bEnrDf = processData(bEnrPath, 'bEnr')
#for fileName in bEnrFileNames:
#    tmpDf = processData(fileName, 'bEnr')
#    bEnrDf = bEnrDf.append(tmpDf, ignore_index=True, sort=False)
#
## datafile
dataDf = processData(dataPath, 'data')


'''
pickle.dump(ggHDf, open('dataVsMC/ggHDf.pkl', 'wb'))
pickle.dump(BGenDf, open('dataVsMC/BGenDf.pkl', 'wb'))
pickle.dump(bEnrDf, open('dataVsMC/bEnrDf.pkl', 'wb'))
pickle.dump(dataDf, open('dataVsMC/dataDf.pkl', 'wb'))

ggHDf = pickle.load(open('dataVsMC/ggHDf.pkl', 'rb'))
BGenDf = pickle.load(open('dataVsMC/BGenDf.pkl', 'rb'))
bEnrDf = pickle.load(open('dataVsMC/bEnrDf.pkl', 'rb'))
dataDf = pickle.load(open('dataVsMC/dataDf.pkl', 'rb'))
'''

d = True
chart='step'

for var in trainVars + plotVars:
    print(var)
    if 'pt' in var or 'LHE' in var:
        nbins = 80
    else:
        nbins = 40

    fig, ax = plt.subplots(figsize=(10,8))
    minBGen, maxBGen = np.percentile(BGenDf[var], [0, 99.8])
    minData, maxData = np.percentile(dataDf[var], [0, 99.8])
    if 'H4qvsQCD' in var:
        range_local = (0, max(maxBGen, maxData))
    else:
        range_local = (min(minBGen, minData), max(maxBGen, maxData))

    ax.hist(BGenDf[var], bins=nbins, label='BGen', color='b', range = range_local,
             weights=BGenDf['final_weights'], histtype=chart, density=d)
    ax.hist(dataDf[var], bins=nbins, label='Data', color='r', range = range_local,
            histtype=chart, density=d)
    ax.set_title(var + ' BGen vs Data')
    ax.legend(loc='best', frameon=True)
    ax.grid()
    plt.savefig('dataVsMCDist/BGen/{}_BGenVsData.png'.format(var))
    plt.close()

    fig, ax = plt.subplots(figsize=(10,8))
    minbEnr, maxbEnr = np.percentile(bEnrDf[var], [0, 99.8])
    if 'H4qvsQCD' in var:
        range_local = (0, max(maxbEnr, maxData))
    else:
        range_local = (min(minbEnr, minData), max(maxbEnr, maxData))
    ax.hist(bEnrDf[var], bins=nbins, label='bEnr', color='b', range=range_local,
            weights=bEnrDf['final_weights'], histtype=chart, density=d)
    ax.hist(dataDf[var], bins=nbins, label='Data', color='r', range = range_local,
            histtype=chart, density=d)
    ax.legend(loc='best', frameon=True)
    ax.set_title(var + ' bEnr vs Data')
    ax.grid()
    plt.savefig('dataVsMCDist/bEnr/{}_bEnrVsData.png'.format(var))
    plt.close()

nbins=80
fig, ax = plt.subplots(figsize=(10,8))
minbEnr, maxbEnr = np.percentile(bEnrDf['LHE_HT'], [0, 99.8])
minBGen, maxBGen = np.percentile(BGenDf['LHE_HT'], [0, 99.8])
range_local = (min(minbEnr, minBGen), max(maxbEnr, maxBGen))
ax.hist(bEnrDf['LHE_HT'], bins=nbins, label='bEnr', color='b', #range=range_local,
        weights=bEnrDf['final_weights'], histtype=chart, density=d)
ax.hist(BGenDf['LHE_HT'], bins=nbins, label='BGen', color='r', #range = range_local,
        weights=BGenDf['final_weights'], histtype=chart, density=d)
ax.legend(loc='best', frameon=True)
ax.set_title( 'LHE_HT bEnr vs BGen')
ax.grid()
plt.savefig('dataVsMCDist/LHE_HT_bEnrVsBGen.png')
plt.close()




# maxJetIdx = np.argmax(jetPt, axis=1)
# maxMuonIdx = np.argmax(muonPt, axis=1)
# maxJetIdx = jetPt.idxmax(axis=1)
# maxMuonIdx = muonPt.idxmax(axis=1)
# maxJetIdx = maxJetIdx.fillna(value=4)
# maxMuonIdx = maxMuonIdx.fillna(value=4)

# new = maxJetIdx.combine(maxMuonIdx,  how = 'inner' ,indicator=False)

