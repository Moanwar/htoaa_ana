#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 17:46:30 2021

@author: si_sutantawibul1
"""
import pandas as pd
import uproot
from htoaaRootFilesLoc import ParkedDataPaths
import dataVsMC_DataManager as DM
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.special import ndtri

hist_params = {'bins':30, 'range':(0,1500)}
append_params = {'ignore_index':True, 'sort':False}
pickledir = 'JetHTTrigEff/pickles/'



def normalError(total, passed, totalWeights):
#------------- definitions
## "passed"    = Weighted number of events passed in bin
## "total"     = Weighted number of total events in bin
## "totalErr"  = Uncertainty on "total" = sqrt(sum of squares of weights of "total" events in bin)


    ## in case the dataframe doesn't have a weight
    if not isinstance(totalWeights, pd.DataFrame):
        totalWeights = np.sqrt(len(total))

    level = .68
    alpha = (1-level)/2
    totalErr = np.sqrt(np.sum(totalWeights*totalWeights))
    avgWeight = pow(totalErr, 2)/total
    jitter = avgWeight/total
    average = passed/total
    sigma = np.sqrt(avgWeight * (average+jitter)*(1+jitter-average)/total)
    delta = -sigma*ndtri(1-alpha)

    ## upper is upper uncertainty, lower is lower uncertainty
    upper = (average + delta) if ((average + delta) < 1) else 1
    lower = (average - delta) if ((average - delta) > 0) else 0

    return upper, lower


root = False


#------------------ L1_SingleJet180 and HLT_AK8PFJet500  ----------------------
#-------------------------- ParkingBPH data -----------------------------------
#%%
if root:
    parkedDf = DM.processData(ParkedDataPaths[0], 'data', 'Parked', False, False)
    parkedDftrig = DM.processData(ParkedDataPaths[0], 'data', 'Parked', False, True)
    
    pickle.dump(parkedDf, open(pickledir + 'parkedDf.pkl', 'wb'))
    pickle.dump(parkedDftrig, open(pickledir + 'parkedDftrig.pkl', 'wb'))
else:
    parkedDf = pickle.load(open(pickledir + 'parkedDf.pkl', 'rb'))
    parkedDftrig = pickle.load(open(pickledir + 'parkedDftrig.pkl', 'rb'))

demParked, demParkedEdge = np.histogram(parkedDf.FatJet_pt,
                                        weights=parkedDf.final_weights,
                                        **hist_params)
numParked, numParkedEdge = np.histogram(parkedDftrig.FatJet_pt,
                                        weights=parkedDftrig.final_weights,
                                        **hist_params)
#ParkedErrup, ParkedErrdown = normalError(total = demParked,
#                                         passed = numParked,
#                                         totalWeights = parkedDf.final_weights)
quotParked = np.divide(numParked, demParked, where=demParked!=0)

#%%

#---------------------------- MuonEG data ------------------------------------
if root:
    MuonEGDf = pd.DataFrame()
    for fileName in DM.MuonEGPaths:
        tmpDf = DM.processData(fileName, 'MuonEG', 'MuonEG', MC=False, trigger=False)
        MuonEGDf = MuonEGDf.append(tmpDf, **append_params)
    MuonEGDftrig = pd.DataFrame()
    for fileName in DM.MuonEGPaths:
        tmpDf = DM.processData(fileName, 'MuonEG', 'MuonEG', MC=False, trigger=True)
        MuonEGDftrig = MuonEGDftrig.append(tmpDf, **append_params)
    
    pickle.dump(MuonEGDf, open(pickledir + 'MuonEGDf.pkl', 'wb'))
    pickle.dump(MuonEGDftrig, open(pickledir + 'MuonEGDftrig.pkl', 'wb'))
else:
    MuonEGDf = pickle.load(open(pickledir + 'MuonEGDf.pkl', 'rb'))
    MuonEGDftrig = pickle.load(open(pickledir + 'MuonEGDftrig.pkl', 'rb'))

demMuonEG, demMuonEGEdge = np.histogram(MuonEGDf.FatJet_pt,
                                        #weights=MuonEGDf.final_weights,
                                        **hist_params)
numMuonEG, numMuonEGEDge = np.histogram(MuonEGDftrig.FatJet_pt,
                                        #weights=MuonEGDftrig.final_weights,
                                        **hist_params)
quotMuonEG = np.divide(numMuonEG, demMuonEG, where=demMuonEG!=0)


#--------------------- QCD MC with ParkingBPH selection and weights -----------
if root:
    QCDParkingDf = pd.DataFrame()
    QCDParkingDftrig = pd.DataFrame()
    for fileName in DM.BGenPaths:
        tmpDf = DM.processData(fileName, 'BGen', 'Parked', MC=True, trigger=False)
        QCDParkingDf = QCDParkingDf.append(tmpDf, **append_params)
    
        tmpDf = DM.processData(fileName, 'BGen', 'Parked', MC=True, trigger=True)
        QCDParkingDftrig = QCDParkingDftrig.append(tmpDf, **append_params)
    
    for fileName in DM.bEnrPaths:
        tmpDf = DM.processData(fileName, 'bEnr', 'Parked', MC=True, trigger=False)
        QCDParkingDf = QCDParkingDf.append(tmpDf, **append_params)
    
        tmpDf = DM.processData(fileName, 'bEnr', 'Parked', MC=True, trigger=True)
        QCDParkingDftrig = QCDParkingDftrig.append(tmpDf, **append_params)
    
    pickle.dump(QCDParkingDf, open(pickledir + 'QCDParkingDf.pkl', 'wb'))
    pickle.dump(QCDParkingDftrig, open(pickledir + 'QCDParkingDftrig.pkl', 'wb'))
else:
    QCDParkingDf = pickle.load(open(pickledir + 'QCDParkingDf.pkl', 'rb'))
    QCDParkingDftrig = pickle.load(open(pickledir + 'QCDParkingDftrig.pkl', 'rb'))

demQCDParking, demQCDParkingEdge = np.histogram(QCDParkingDf.FatJet_pt,
                                                weights=QCDParkingDf.final_weights,
                                                **hist_params)
numQCDParking, numQCDParkingEdge = np.histogram(QCDParkingDftrig.FatJet_pt,
                                                weights=QCDParkingDftrig.final_weights,
                                                **hist_params)
quotQCDParking = np.divide(numQCDParking, demQCDParking, where=demQCDParking!=0)

#--------------- QCD MC with no offline muon or pt/IP/PU weights -------------

if root:
    QCDBaseDf = pd.DataFrame()
    QCDBaseDftrig = pd.DataFrame()
    for fileName in DM.BGenPaths:
        tmpDf = DM.processData(fileName, 'BGen', 'Base', MC=True, trigger=False)
        QCDBaseDf = QCDBaseDf.append(tmpDf, **append_params)

        tmpDf = DM.processData(fileName, 'BGen', 'Base', MC=True, trigger=True)
        QCDBaseDftrig = QCDBaseDftrig.append(tmpDf, **append_params)
    for fileName in DM.bEnrPaths:
        tmpDf = DM.processData(fileName, 'bEnr', 'Base', MC=True, trigger=False)
        QCDBaseDf = QCDBaseDf.append(tmpDf, **append_params)
    
        tmpDf = DM.processData(fileName, 'bEnr', 'Base', MC=True, trigger=True)
        QCDBaseDftrig = QCDBaseDftrig.append(tmpDf, **append_params)
    pickle.dump(QCDBaseDf, open(pickledir + 'QCDBaseDf.pkl', 'wb'))
    pickle.dump(QCDBaseDftrig, open(pickledir + 'QCDBaseDftrig.pkl', 'wb'))
else:
    QCDBaseDf = pickle.load(open(pickledir + 'QCDBaseDf.pkl', 'rb'))
    QCDBaseDftrig = pickle.load(open(pickledir + 'QCDBaseDftrig.pkl', 'rb'))

demQCDBase, demQCDBaseEdge = np.histogram(QCDBaseDf.FatJet_pt,
                                          weights=QCDBaseDf.final_weights,
                                          **hist_params)
numQCDBase, numQCDBaseEdge = np.histogram(QCDBaseDftrig.FatJet_pt,
                                          weights=QCDBaseDftrig.final_weights,
                                          **hist_params)
quotQCDBase = np.divide(numQCDBase, demQCDBase, where=demQCDBase!=0)



#---------------- ggH with parking selections and weights ---------------------

if root==True:
    ggHParkingDf = DM.processData(DM.ggHPaths, 'ggH', 'Parked', MC=True, trigger=False)
    ggHParkingDftrig = DM.processData(DM.ggHPaths, 'ggH', 'Parked', MC=True, trigger=True)
    pickle.dump(ggHParkingDf, open(pickledir + 'ggHParkingDf.pkl', 'wb'))
    pickle.dump(ggHParkingDftrig, open(pickledir + 'ggHParkingDftrig.pkl', 'wb'))
else:
    ggHParkingDf = pickle.load(open(pickledir + 'ggHParkingDf.pkl', 'rb'))
    ggHParkingDftrig = pickle.load(open(pickledir + 'ggHParkingDftrig.pkl', 'rb'))

demggHParked, demggHParkedEdge = np.histogram(ggHParkingDf.FatJet_pt,
                                              weights=ggHParkingDf.final_weights,
                                              **hist_params)
numggHParked, numggHParkedEdge = np.histogram(ggHParkingDftrig.FatJet_pt,
                                              weights=ggHParkingDftrig.final_weights,
                                              **hist_params)
quotggHParked = np.divide(numggHParked, demggHParked, where=demggHParked!=0)
#-----------------------------------------------------------------------------

#---------------- ggH without muon selection or pT/IP/PU weights --------------
if root==True:
    ggHBaseDf = DM.processData(DM.ggHPaths, 'ggH', 'Base', MC=True, trigger=False)
    ggHBaseDftrig = DM.processData(DM.ggHPaths, 'ggH', 'Base', MC=True, trigger=True)
    pickle.dump(ggHBaseDf, open(pickledir + 'ggHBaseDf.pkl', 'wb'))
    pickle.dump(ggHBaseDftrig, open(pickledir + 'ggHBaseDftrig.pkl', 'wb'))
else:
    ggHBaseDf = pickle.load(open(pickledir + 'ggHBaseDf.pkl', 'rb'))
    ggHBaseDftrig = pickle.load(open(pickledir + 'ggHBaseDftrig.pkl', 'rb'))

demggHBase, demggHBaseEdge = np.histogram(ggHBaseDf.FatJet_pt,
                                          weights=ggHBaseDf.final_weights,
                                          **hist_params)
numggHBase, numggHBaseEdge = np.histogram(ggHBaseDftrig.FatJet_pt,
                                          weights=ggHBaseDftrig.final_weights,
                                          **hist_params)
quotggHBase = np.divide(numggHBase, demggHBase, where=demggHBase!=0)
#------------------------------------------------------------------------------

#---------------------- ttbar MC with muon selection --------------------------
if root==True:
    TTBarMuonEG = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'MuonEG', MC=True, trigger=False)
    TTBarMuonEGtrig = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'MuonEG', MC=True, trigger=True)
    pickle.dump(TTBarMuonEG, open(pickledir + 'TTBarMuonEG.pkl', 'wb'))
    pickle.dump(TTBarMuonEGtrig, open(pickledir + 'TTBarMuonEGtrig.pkl', 'wb'))
else:
    TTBarMuonEG = pickle.load(open(pickledir + 'TTBarMuonEG.pkl', 'rb'))
    TTBarMuonEGtrig = pickle.load(open(pickledir + 'TTBarMuonEGtrig.pkl', 'rb'))

demTTBarMuonEG, demTTBarMuonEGEdge = np.histogram(TTBarMuonEG.FatJet_pt,
                                                  weights=TTBarMuonEG.final_weights,
                                                  **hist_params)
numTTBarMuonEG, numTTBarMuonEGEdge = np.histogram(TTBarMuonEGtrig.FatJet_pt,
                                                  weights=TTBarMuonEGtrig.final_weights,
                                                  **hist_params)
quotTTBarMuonEG = np.divide(numTTBarMuonEG, demTTBarMuonEG, where=demTTBarMuonEG!=0)
#------------------------------------------------------------------------------

#----------------------- ttbar mc without muon selection ---------------------
if root==True:
    TTBarBaseDf = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'Base', MC=True, trigger=False)
    TTBarBaseDftrig = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'Base', MC=True, trigger=True)
    pickle.dump(TTBarBaseDf, open(pickledir + 'TTBarBaseDf.pkl', 'wb'))
    pickle.dump(TTBarBaseDftrig, open(pickledir + 'TTBarBaseDftrig.pkl', 'wb'))
else:
    TTBarBaseDf = pickle.load(open(pickledir + 'TTBarBaseDf.pkl', 'rb'))
    TTBarBaseDftrig = pickle.load(open(pickledir + 'TTBarBaseDftrig.pkl', 'rb'))

demTTBarBase, demTTBarBaseEdge = np.histogram(TTBarBaseDf.FatJet_pt,
                                              weights=TTBarBaseDf.final_weights,
                                              **hist_params)
numTTBarBase, numTTBarBaseEdge = np.histogram(TTBarBaseDftrig.FatJet_pt,
                                              weights=TTBarBaseDftrig.final_weights,
                                              **hist_params)
quotTTBarBase = np.divide(numTTBarBase, demTTBarBase, where=demTTBarBase!=0)
#-----------------------------------------------------------------------------




## plotting all df to see the rang.....there's so many

## parkedDf, parkedDftrig, MuonEGDf, MuonEGDftrig, QCDParkingDf, QCDParkingDftrig
## QCDBaseDf, QCDBaseDftrig, ggHParkingDf, ggHParkingDftrig, ggHBaseDf, ggHBaseDftrig
## TTBarMuonEG, TTBarMuonEGtrig, TTBarBaseDf, TTBarBaseDftrig



toplot = [(parkedDf, parkedDftrig),
          (MuonEGDf, MuonEGDftrig),
          (QCDParkingDf, QCDParkingDftrig),
          (QCDBaseDf, QCDBaseDftrig),
          (ggHParkingDf, ggHParkingDftrig),
          (ggHBaseDf, ggHBaseDftrig),
          (TTBarMuonEG, TTBarMuonEGtrig),
          (TTBarBaseDf, TTBarBaseDftrig)]

denoms = [demParked, ]


numers = [numParked]
## convert toplot into the only col we want
for i in toplot:
    for j in i:
        j = j['FatJet_pt']



scatterlist = [quotParked, quotMuonEG, quotQCDParking, quotQCDBase,
               quotggHParked, quotggHBase, quotTTBarMuonEG, quotTTBarBase]
scatterlabel = ['parking data', 'MuonEG', 'QCD (parking)', 'QCD',
                'ggH (parking)', 'ggH', 'ttbar (muonEG)', 'ttbar']

plotdir = 'JetHTTrigEff/plots/'
edge = demTTBarBaseEdge[1:]
for i, val in enumerate(scatterlist):
    plt.scatter(edge, val, label=scatterlabel[i])

    print('edge: ', edge)
    print('val: ', val)
    #for j, edval in enumerate(edge):





    plt.title(scatterlabel[i])
    plt.savefig(f'{plotdir}{scatterlabel[i]}.png')
    plt.show()
#plt.legend(loc='best')
#plt.title('Efficiency plots')

#plt.show()


# hist_params = {'bins':40}
# #fig, ax = plt.subplots()


# fig, ax = plt.subplots()
# for i in toplot:
#     val, edge, _ = ax.hist(i['FatJet_pt'], **hist_params)
#     print(edge[-1])
# plt.show()













