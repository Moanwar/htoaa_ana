import pickle
import numpy as np
import pandas as pd
from analib import PhysObj, Event
from info import BGenFileNames, bEnrFileNames

import dataVsMC_DataManager as DM #import processData, jetVars, muonVars, PVVars, allVars, dataPath, ggHPath, bEnrPaths, BGenPaths, TTJetsPaths,
from dataManager import trainVars # this is so we can run our thing through the BDT XML correctly
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import mplhep as hep
import uproot

#np.seterr(divide='ignore', invalid='ignore')
# testesttset

## get that sweet CMS style plots
plt.style.use(hep.style.CMS)


## if reading from rootfiles, set true. if already have pickle of dataframe, set false
root = True

## function to add a BDTScore column to each of the background/signal/data DF
loadedModel = pickle.load(open('Htoaa_BDThigh disc.pkl', 'rb'))
def analyze(dataDf):
    prediction = loadedModel.predict_proba(dataDf[trainVars])
    dataDf = dataDf.assign(BDTScore=prediction[:,1])
    return dataDf

if root:
    ## monte carlo ggH signal
    ggHDf = DM.processData(DM.ggHPaths, 'ggH')

    ############## BGen ##############
    BGenDf = pd.DataFrame()
    for fileName in DM.BGenPaths:
        tmpDf = DM.processData(fileName, 'BGen')
        BGenDf = BGenDf.append(tmpDf, ignore_index=True, sort=False)
    ## calculate Xsec weights for each event
    #LHE_QCD_corrections = BGenDf.LHE_weights * BGenDf.QCD_correction
    #BGenDf['final_weights'] = LHE_QCD_corrections * 21.56/np.sum(LHE_QCD_corrections) * 1043278
    ###################################

    ############## bEnr ##############
    #bEnrDf = processData(bEnrPath, 'bEnr')
    ## uncommend below if using multiple bg MC files
    bEnrDf = pd.DataFrame()
    for fileName in DM.bEnrPaths:
        tmpDf = DM.processData(fileName, 'bEnr')
        bEnrDf = bEnrDf.append(tmpDf, ignore_index=True, sort=False)
    #LHE_QCD_corrections = bEnrDf.LHE_weights * bEnrDf.QCD_correction
    #bEnrDf['final_weights'] = LHE_QCD_corrections * 8.20/np.sum(LHE_QCD_corrections)
    ###################################

    ############## TTJets ##############
    TTJetsDf = pd.DataFrame()
    for fileName in DM.TTJetsPaths:
        tmpdf = DM.processData(fileName, 'TTJets')
        TTJetsDf = TTJetsDf.append(tmpdf, ignore_index=True, sort=False)

    #length = len(TTJetsDf)
    #TTJetsDf = TTJetsDf.assign(final_weights=831760.0/10244307/length)
    ###################################

    ############## ZJets ##############
    ZJetsDf = pd.DataFrame()
    for fileName in DM.ZJetsPaths:
        tmpdf = DM.processData(fileName, 'ZJets')
        ZJetsDf = ZJetsDf.append(tmpdf, ignore_index=True, sort=False)

    # length = len(ZJetsDf.loc[(ZJetsDf['LHE_HT']>=400) & (ZJetsDf['LHE_HT']<600)])
    # if length != 0:
    #     ZJetsDf.loc[(ZJetsDf['LHE_HT']>=400) & (ZJetsDf['LHE_HT']<600), 'LHE_weights'] = 145400/16704355/length

    # length = len(ZJetsDf.loc[(ZJetsDf['LHE_HT']>=600) & (ZJetsDf['LHE_HT']<800)])
    # if length != 0:
    #     ZJetsDf.loc[(ZJetsDf['LHE_HT']>=600) & (ZJetsDf['LHE_HT']<800), 'LHE_weights'] = 34000/14642701/length

    # length = len(ZJetsDf.loc[(ZJetsDf['LHE_HT']>=800)])
    # if length != 0:
    #     ZJetsDf.loc[(ZJetsDf['LHE_HT']>=800), 'LHE_weights'] = 18670/10561192/length

    # ZJetsDf = ZJetsDf.assign(final_weights = ZJetsDf['LHE_weights'])
    ###################################

    ############## WJets ##############
    WJetsDf = pd.DataFrame()
    for fileName in DM.WJetsPaths:
        tmpdf = DM.processData(fileName, 'WJets')
        WJetsDf = WJetsDf.append(tmpdf, ignore_index=True, sort=False)

    # length = len(WJetsDf.loc[(WJetsDf['LHE_HT']>=400) & (WJetsDf['LHE_HT']<600)])
    # if length != 0:
    #     WJetsDf.loc[(WJetsDf['LHE_HT']>=400) & (WJetsDf['LHE_HT']<600), 'LHE_weights'] = 315600/10071273/length

    # length =len(WJetsDf.loc[(WJetsDf['LHE_HT']>=600) & (WJetsDf['LHE_HT']<800)])
    # if length != 0:
    #     WJetsDf.loc[(WJetsDf['LHE_HT']>=600) & (WJetsDf['LHE_HT']<800), 'LHE_weights'] = 68570/15298056/length

    # length = len(WJetsDf.loc[(WJetsDf['LHE_HT']>=800)])
    # if length != 0:
    #     WJetsDf.loc[(WJetsDf['LHE_HT']>=800), 'LHE_weights'] = 34900/14627242/length

    # WJetsDf = WJetsDf.assign(final_weights = WJetsDf['LHE_weights'])
    ######################################



    # datafile
    if not DM.JetHT:
        dataDf = DM.processData(DM.ParkedDataPaths[0], 'data')
        pickle.dump(dataDf, open('dataVsMC/dataDf.pkl', 'wb'))


    if DM.JetHT:
        JetHTDf = pd.DataFrame()
        for fileName in DM.JetHTPaths:
            tmpdf = DM.processData(fileName, 'JetHT')
            JetHTDf = JetHTDf.append(tmpdf, ignore_index=True, sort=False)
        JetHTDf['final_weights'] = 1
        pickle.dump(JetHTDf, open('dataVsMC/JetHTDf.pkl', 'wb'))


    pickle.dump(ggHDf, open('dataVsMC/ggHDf.pkl', 'wb'))
    pickle.dump(BGenDf, open('dataVsMC/BGenDf.pkl', 'wb'))
    pickle.dump(bEnrDf, open('dataVsMC/bEnrDf.pkl', 'wb'))
    pickle.dump(TTJetsDf, open('dataVsMC/TTJetsDf.pkl', 'wb'))
    pickle.dump(ZJetsDf, open('dataVsMC/ZJetsDf.pkl', 'wb'))
    pickle.dump(WJetsDf, open('dataVsMC/WJetsDf.pkl', 'wb'))


else:
    ggHDf = pickle.load(open('dataVsMC/ggHDf.pkl', 'rb'))
    BGenDf = pickle.load(open('dataVsMC/BGenDf.pkl', 'rb'))
    bEnrDf = pickle.load(open('dataVsMC/bEnrDf.pkl', 'rb'))
    TTJetsDf = pickle.load(open('dataVsMC/TTJetsDf.pkl', 'rb'))
    ZJetsDf = pickle.load(open('dataVsMC/ZJetsDf.pkl', 'rb'))
    WJetsDf = pickle.load(open('dataVsMC/WJetsDf.pkl', 'rb'))
    if DM.JetHT:
        JetHTDf = pickle.load(open('dataVsMC/JetHTDf.pkl', 'rb'))
    else:
        dataDf = pickle.load(open('dataVsMC/dataDf.pkl', 'rb'))


ggHDf = analyze(ggHDf)
BGenDf = analyze(BGenDf)
bEnrDf = analyze(bEnrDf)
TTJetsDf = analyze(TTJetsDf)
ZJetsDf = analyze(ZJetsDf)
WJetsDf = analyze(WJetsDf)
if DM.JetHT:
    JetHTDf = analyze(JetHTDf)
else:
    dataDf = analyze(dataDf)



dfdict = { 'WJets': WJetsDf,
           'ZJets': ZJetsDf,
           'TTJets': TTJetsDf,
           'bEnr': bEnrDf,
           'BGen': BGenDf,
           }

if DM.JetHT:
    cols = list(JetHTDf.columns)
else:
    cols = list(dataDf.columns)

cols.remove('final_weights')
for var in cols:
    if 'pt' in var:
        nbins = 40
    else:
        nbins = 20

    print(var)

    fig, ax = plt.subplots(figsize=(10,8))

    ## get the min/max value for hists
    xmin = list()
    xmax = list()

    for dfkey, df in dfdict.items():
        xmintmp, xmaxtmp = np.percentile(df[var], [0,99.8])
        xmin.append(xmintmp)
        xmax.append(xmaxtmp)


    if 'H4qvsQCD' in var:
        #range_local = (0, max(maxBGen, maxData))
        range_local = (0, max(xmax))
    elif 'FatJet_pt'==var:
        range_local = (240, 700)
    elif 'BDTScore'==var:
        range_local = (0,0.7)
    else:
        #range_local = (min(minBGen, minData), max(maxBGen, maxData))
        range_local = (min(xmin), max(xmax))


    density = False
    hist_params = {'density': density, 'histtype': 'bar', 'range' : range_local, 'bins':nbins, 'stacked':True}


    toplot = pd.DataFrame()
    toplotweights = pd.DataFrame()
    toplotlabel = list()
    for dfkey, df in dfdict.items():
        #toplot = np.append(toplot, df[var], 1)
        #toplotweights = np.append(toplotweights, df['final_weights'], 1)
        toplot[dfkey] = df[var]
        toplotweights[dfkey] = df['final_weights']
        toplotlabel.append(f'{dfkey} ({round(np.sum(df.final_weights))})')

    if DM.JetHT:
        ax.hist(JetHTDf[var].values, label=f'JetHT ({round(np.sum(JetHTDf.final_weights))})', histtype='step',
                density=density, bins=nbins, linewidth=3, color='k',
                range=range_local, weights=JetHTDf.final_weights)
    else:
        ax.hist(dataDf[var].values, label=f'parkedData ({round(np.sum(dataDf.final_weights))})',
                histtype='step',
                density=density, bins=nbins, linewidth=3, color='k',
                range=range_local, weights=dataDf.final_weights)

    ax.hist(ggHDf[var].values, label=f'GGH ({round(np.sum(ggHDf.final_weights))})', histtype='step',
            density=density, bins=nbins, linewidth=3, color='r',
            range=range_local, weights=ggHDf.final_weights)

    ## making color palette for the QCD stakcs
    pal = ['#603514', '#b940f2','#ec6f38', '#6acaf8', '#82f759']
    ax.hist(toplot.values, weights=toplotweights.values,label=toplotlabel, color=pal, **hist_params)
    #ax.hist(BGenDf[var].values, weights=BGenDf['final_weights'].values, label='BGen', **hist_params)
    #ax.hist(bEnrDf[var].values, weights=bEnrDf['final_weights'].values, label='bEnr', **hist_params)
    #ax.hist(TTJetsDf[var].values, weights=TTJetsDf['final_weights'].values, label='TTJets', **hist_params)



    #ax.hist(dataDf[var], label='Data',
    ax.set_title(var )# + ' JetHT')
    ax.legend(loc='best', frameon=True)
    ax.grid()
    if DM.JetHT:
        filedest = 'dataVsMCDist/JetHT/{}.png'
    else:
        filedest = 'dataVsMCDist/Parked/{}.png'
    plt.savefig(filedest.format(var))
    plt.show()
    plt.close()


########################### Sensitivity Plots ##########################
sortedggH = ggHDf.sort_values(by='BDTScore', axis=0, kind='mergesort')
histHeight = sortedggH.final_weights.sum()/10
edges = [0]
cumuSum = 0
edgesloc = []
for i in range(sortedggH.BDTScore.size):
    cumuSum = cumuSum + sortedggH.final_weights.iloc[i]
    if cumuSum > histHeight:
        edges.append(sortedggH.BDTScore.iloc[i-1])
        cumuSum = sortedggH.final_weights.iloc[i]
        edgesloc.append(i-1)
edges.pop()
edges.append(sortedggH.BDTScore.iloc[-1])

fig, ax = plt.subplots(figsize=(8,8))
sighist = ax.hist(sortedggH.BDTScore, bins=edges,
                  weights=sortedggH.final_weights, histtype='step',
                  label='ggH', log=True, color='r', linewidth=3)
datahist = ax.hist(dataDf.BDTScore, bins=edges, histtype='step',
                   label='data', log=True, color='k', linewidth=3)

toplot = pd.DataFrame()
toplotweights = pd.DataFrame()
toplotlabel = list()
for dfkey, df in dfdict.items():
    toplot[dfkey] = df['BDTScore']
    toplotweights[dfkey] = df['final_weights']
    toplotlabel.append(dfkey)


bghist = ax.hist(toplot.values, weights=toplotweights.values, bins=edges,
        label=toplotlabel, color=pal, log=True, histtype='bar', stacked=True)


# !!! TODO
## figure out how to compute the sensitivity for the stacked plot bitch

ax.legend(loc='best')
ax.set_title('sensitivty')
plt.savefig(f'dataVsMCDist/Parked/Sensitivity.png')
plt.show()
plt.close()


########################################################################


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

nbins=40

#del dfdict['JetHTDf']



for dfkey, dfvalue in dfdict.items():
    fig, ax = plt.subplots(figsize=(10,8))
    range_local = np.percentile(dfvalue['LHE_HT'], [0,99.8])
    hist_params = {'density': False, 'histtype': 'bar', 'range' : range_local, 'bins':nbins, 'stacked':True}
    ax.hist(dfvalue['LHE_HT'], weights=dfvalue['final_weights'], label=dfkey,
            **hist_params)
    ax.set_title('LHE_HT ' + dfkey)
    ax.grid()
    plt.savefig(f'dataVsMCDist/Parked/LHE_HT_{dfkey}.png')
    plt.show()
    plt.close()
    #plt.clf()


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





