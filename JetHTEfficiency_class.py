#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 19:40:46 2021

@author: si_sutantawibul1
"""
import numpy as np
import pandas as pd
import sys
import pickle
from scipy.special import ndtri
import matplotlib.pyplot as plt
import dataVsMC_DataManager as DM

class EfficiencyInfo(object) :
    def __init__(self, demDf, numDf, name):
        if not (isinstance(demDf, pd.DataFrame) or isinstance(numDf, pd.DataFrame)):
            print('init error: input needs to be dataframes')
            sys.exit()

        self.hasWeights = ('final_weights' in demDf)
        self.name = name
        self.demDf = demDf
        self.dem, self.demEdge = self.hist(demDf)
        self.num, self.numEdge = self.hist(numDf)

        #print(len(self.demEdge))


        #self.demEdge = self.getBinCenter(self.demEdge)


        #print(len(self.demEdge))
        #self.numEdge = self.getBinCenter(self.numEdge)
        self.quot = np.divide(self.num, self.dem, where=self.dem!=0)
        self.upErr, self.lowErr = self.computeError()


        self.plot()


    ## makes the pt into histograms
    def hist(self, df):
        if self.hasWeights:
            return np.histogram(df['FatJet_pt'], weights=df['final_weights'],
                                bins = 30, range = (0,1500))
        else:
            return np.histogram(df['FatJet_pt'], bins=30, range=(0,1500))

    def getBinCenter(self, arr):
        arrCen = list()
        for i in range(len(arr)-1):
            arrCen.append((arr[i+1]+arr[i])/2)
        return arrCen

    ## compute normal error for 1 bin
    ## returns np array of error. Code found at:
    ## https://root.cern.ch/doc/master/TEfficiency_8cxx_source.html#l02744
    def normalError(self, total, passed, weights):
        if 0 == total:
            return 0,0
        if self.hasWeights:
            totalErr = np.sqrt((weights*weights).sum())
        else:
            totalErr =  np.sqrt(len(total))


        level = 0.68
        alpha = (1-level)/2
        avgWgt = np.power(totalErr,2)/total
        jitter = avgWgt/total
        average = passed/total
        sigma = np.sqrt(avgWgt*(average+jitter)*(1+jitter-average)/total)
        delta = -sigma*ndtri(1-alpha)

        upper = min(delta,1-average)#(average + delta) if ((average + delta) < 1) else 1
        lower = min(delta,average)#(average - delta) if ((average - delta) > 0) else 0

        return upper, lower

    ## compute error for whole histogram
    def computeError(self, ):
        edges = self.demEdge
        df = self.demDf
        upperError = list()
        lowerError = list()
        for i in range(len(edges)-1):
            if self.hasWeights:
                wgts = df.final_weights.loc[(df.FatJet_pt > edges[i])
                                            & (df.FatJet_pt < edges[i+1])]
            else:
                wgts = 1

            tmpUp, tmpLow = self.normalError(self.dem[i], self.num[i], wgts)
            upperError.append(tmpUp)
            lowerError.append(tmpLow)

        return np.array(upperError), np.array(lowerError)

    def plot(self, ):
        edge = self.demEdge[1:]
        #edge = self.demEdge
        fig, ax = plt.subplots()
        ax.grid()
        #ax.scatter(edge, self.quot, linestyle='None')
        ax.set_ylim([0,1.05])
        ax.set_title(self.name)
        ax.errorbar(edge, self.quot, yerr=(self.upErr, self.lowErr),
                    linestyle='None',fmt='ok', capsize=5)

        plotdir = 'JetHTTrigEff/plots/'
        plt.savefig(f'{plotdir}{self.name}.png')





pickledir = 'JetHTTrigEff/pickles/'
root = True
hist_params = {'bins':30, 'range':(0,1500)}
append_params = {'ignore_index':True, 'sort':False}

#-------------------------- ParkingBPH data -----------------------------------

if root:
    parkedDf = DM.processData(DM.ParkedDataPaths[0], 'data', 'Parked', False, False)
    parkedDftrig = DM.processData(DM.ParkedDataPaths[0], 'data', 'Parked', False, True)
    
    pickle.dump(parkedDf, open(pickledir + 'parkedDf.pkl', 'wb'))
    pickle.dump(parkedDftrig, open(pickledir + 'parkedDftrig.pkl', 'wb'))
else:
    parkedDf = pickle.load(open(pickledir + 'parkedDf.pkl', 'rb'))
    parkedDftrig = pickle.load(open(pickledir + 'parkedDftrig.pkl', 'rb'))
parked = EfficiencyInfo(parkedDf, parkedDftrig, 'Parking data')
#------------------------------------------------------------------------------
#------------------------------ MuonEG data -----------------------------------
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
MuonEGDf['final_weights'] = 1
MuonEGDftrig['final_weights'] = 1
muonEG = EfficiencyInfo(demDf=MuonEGDf, numDf=MuonEGDftrig, name='MuonEG')
#-------------------------------------------------------------------------
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
QCD_parking = EfficiencyInfo(demDf=QCDParkingDf, numDf=QCDParkingDftrig, name='QCD (parking selection)')
#-------------------------------------------------------------------------
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
QCD_base = EfficiencyInfo(demDf=QCDBaseDf, numDf=QCDBaseDftrig, name='QCD (base selection)')
#------------------------------------------------------------------------------------
#---------------- ggH with parking selections and weights ---------------------

if root==True:
    ggHParkingDf = DM.processData(DM.ggHPaths, 'ggH', 'Parked', MC=True, trigger=False)
    ggHParkingDftrig = DM.processData(DM.ggHPaths, 'ggH', 'Parked', MC=True, trigger=True)
    pickle.dump(ggHParkingDf, open(pickledir + 'ggHParkingDf.pkl', 'wb'))
    pickle.dump(ggHParkingDftrig, open(pickledir + 'ggHParkingDftrig.pkl', 'wb'))
else:
    ggHParkingDf = pickle.load(open(pickledir + 'ggHParkingDf.pkl', 'rb'))
    ggHParkingDftrig = pickle.load(open(pickledir + 'ggHParkingDftrig.pkl', 'rb'))
ggH_parking = EfficiencyInfo(demDf=ggHParkingDf, numDf=ggHParkingDftrig, name='GGH (parking selection)')
#--------------------------------------------------------------------------------
#---------------- ggH without muon selection or pT/IP/PU weights --------------


if root==True:
    ggHBaseDf = DM.processData(DM.ggHPaths, 'ggH', 'Base', MC=True, trigger=False)
    ggHBaseDftrig = DM.processData(DM.ggHPaths, 'ggH', 'Base', MC=True, trigger=True)
    pickle.dump(ggHBaseDf, open(pickledir + 'ggHBaseDf.pkl', 'wb'))
    pickle.dump(ggHBaseDftrig, open(pickledir + 'ggHBaseDftrig.pkl', 'wb'))

else:
    ggHBaseDf = pickle.load(open(pickledir + 'ggHBaseDf.pkl', 'rb'))
    ggHBaseDftrig = pickle.load(open(pickledir + 'ggHBaseDftrig.pkl', 'rb'))
ggH_base = EfficiencyInfo(demDf=ggHBaseDf, numDf=ggHBaseDftrig, name='GGH (base selection)')

#--------------------------------------------------------------------------------
#---------------------- ttbar MC with muon selection --------------------------
if root==True:
    TTBarMuonEG = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'MuonEG', MC=True, trigger=False)
    TTBarMuonEGtrig = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'MuonEG', MC=True, trigger=True)
    pickle.dump(TTBarMuonEG, open(pickledir + 'TTBarMuonEG.pkl', 'wb'))
    pickle.dump(TTBarMuonEGtrig, open(pickledir + 'TTBarMuonEGtrig.pkl', 'wb'))
else:
    TTBarMuonEG = pickle.load(open(pickledir + 'TTBarMuonEG.pkl', 'rb'))
    TTBarMuonEGtrig = pickle.load(open(pickledir + 'TTBarMuonEGtrig.pkl', 'rb'))
TTBar_muonEG = EfficiencyInfo(demDf=TTBarMuonEG, numDf=TTBarMuonEGtrig, name='TTBar (MuonEG selection)')
#----------------------------------------------------------------------------------------
#----------------------- ttbar mc without muon selection ---------------------
if root==True:
    TTBarBaseDf = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'Base', MC=True, trigger=False)
    TTBarBaseDftrig = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'Base', MC=True, trigger=True)
    pickle.dump(TTBarBaseDf, open(pickledir + 'TTBarBaseDf.pkl', 'wb'))
    pickle.dump(TTBarBaseDftrig, open(pickledir + 'TTBarBaseDftrig.pkl', 'wb'))
else:
    TTBarBaseDf = pickle.load(open(pickledir + 'TTBarBaseDf.pkl', 'rb'))
    TTBarBaseDftrig = pickle.load(open(pickledir + 'TTBarBaseDftrig.pkl', 'rb'))
TTBar_base = EfficiencyInfo(demDf=TTBarBaseDf, numDf=TTBarBaseDftrig, name='TTBar (base selection)')
#--------------------------------------------------------------------------------

pickle.dump(ggH_base, open('JetHTTrigEff/pickles/histpickls/ggH_base.pkl', 'wb'))
