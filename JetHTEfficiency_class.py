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
from scipy.stats import norm
import matplotlib.pyplot as plt
import dataVsMC_DataManager as DM

class EfficiencyInfo(object) :
    def __init__(self, demDf, name):
        if not (isinstance(demDf, pd.DataFrame)):
            print('init error: input needs to be dataframes')
            sys.exit()

        self.hasWeights = ('final_weights' in demDf)
        self.name = name
        self.demDf = demDf
        self.numDf = self.getNumDf(demDf)
        self.dem, self.demEdge = self.hist(demDf)
        self.num, self.numEdge = self.hist(self.numDf)
        self.quot = np.divide(self.num, self.dem, where=self.dem!=0)
        self.upErr, self.lowErr = self.computeError()


        #self.plot()


    ## makes the pt into histograms
    def hist(self, df):
        if self.hasWeights:
            return np.histogram(df['FatJet_pt'], weights=df['final_weights'],
                                bins = 27, range = (150,1500))
        else:
            return np.histogram(df['FatJet_pt'], bins=27, range=(150,1500))

    def getNumDf(self, df):
        return df.loc[(df.L1_SingleJet180==True) & (df.HLT_AK8PFJet500==True)]


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
        #else:
        #    totalErr =  np.sqrt(len(total))

        level = 0.68
        alpha = (1-level)/2
        avgWgt = np.power(totalErr,2)/total
        jitter = avgWgt/total
        average = passed/total
        sigma = np.sqrt((avgWgt*(average+jitter)*(1+jitter-average))/total)
        delta = norm.ppf(1-alpha,0,sigma)
        #delta = -sigma*ndtri(1-alpha)

        upper = min(delta,1-average)#(average + delta) if ((average + delta) < 1) else 1
        lower = min(delta,average)#(average - delta) if ((average - delta) > 0) else 0

        return upper, lower

    ## compute error for whole histogram
    def computeError(self, ):
        edges = self.demEdge
        demdf = self.demDf
        upperError = list()
        lowerError = list()

        for i in range(len(edges)-1):
            if self.hasWeights:
                wgts = demdf.final_weights.loc[(demdf.FatJet_pt >= edges[i])
                                            & (demdf.FatJet_pt < edges[i+1])]
            else:
                wgts = 0

            tmpUp, tmpLow = self.normalError(self.dem[i], self.num[i], wgts)
            upperError.append(tmpUp)
            lowerError.append(tmpLow)


        return np.array(upperError), np.array(lowerError)

    def plot(self, ):
        edge = self.getBinCenter(self.demEdge)

        fig, ax = plt.subplots(figsize=(10,6))
        ax.grid()
        ax.set_ylim([-0.05,1.05])
        ax.set_title(self.name)

        ax.errorbar(edge, self.quot, yerr=(self.lowErr, self.upErr),
                    linestyle='None',fmt='ok', capsize=5)
        xerr = np.ones((2, len(self.quot)))*26
        ax.errorbar(edge, self.quot, xerr=xerr, linestyle='None', fmt='k')

        plotdir = 'JetHTTrigEff/plots/'
        plt.savefig(f'{plotdir}{self.name}.png')





pickledir = 'JetHTTrigEff/pickles/'
root = True
hist_params = {'bins':27, 'range':(150,1500)}
append_params = {'ignore_index':True, 'sort':False}



#-------------------------- ParkingBPH data -----------------------------------

if root:
    parkedDf = DM.processData(DM.ParkedDataPaths[0], 'data', 'Parked', MC=False, trigger=True)
    
    #pickle.dump(parkedDf, open(pickledir + 'parkedDf.pkl', 'wb'))

else:
    parkedDf = pickle.load(open(pickledir + 'parkedDf.pkl', 'rb'))

parked = EfficiencyInfo(parkedDf, 'Parking data')
parked.plot()


#------------------------------------------------------------------------------
#------------------------------ MuonEG data -----------------------------------

if root:
    MuonEGDf = pd.DataFrame()
    for fileName in DM.MuonEGPaths:
        tmpDf = DM.processData(fileName, 'MuonEG', 'MuonEG', MC=False, trigger=True)
        MuonEGDf = MuonEGDf.append(tmpDf, **append_params)

    MuonEGDf['final_weights'] = 1
    pickle.dump(MuonEGDf, open(pickledir + 'MuonEGDf.pkl', 'wb'))

else:
    MuonEGDf = pickle.load(open(pickledir + 'MuonEGDf.pkl', 'rb'))


muonEG = EfficiencyInfo(demDf=MuonEGDf, name='MuonEG')
muonEG.plot()

#-------------------------------------------------------------------------
#--------------------- QCD MC with ParkingBPH selection and weights -----------
#root = True
if root:
    QCDParkingDf = pd.DataFrame()
    QCDParkingDftrig = pd.DataFrame()
    for fileName in DM.BGenPaths:
        tmpDf = DM.processData(fileName, 'BGen', 'Parked', MC=True, trigger=True)
        QCDParkingDf = QCDParkingDf.append(tmpDf, **append_params)

    for fileName in DM.bEnrPaths:
        tmpDf = DM.processData(fileName, 'bEnr', 'Parked', MC=True, trigger=True)
        QCDParkingDf = QCDParkingDf.append(tmpDf, **append_params)

    pickle.dump(QCDParkingDf, open(pickledir + 'QCDParkingDf.pkl', 'wb'))

else:
    QCDParkingDf = pickle.load(open(pickledir + 'QCDParkingDf.pkl', 'rb'))

QCD_parking = EfficiencyInfo(demDf=QCDParkingDf, name='QCD (parking selection)')
QCD_parking.plot()

#-------------------------------------------------------------------------
#--------------- QCD MC with no offline muon or pt/IP/PU weights -------------

if root:
    QCDBaseDf = pd.DataFrame()
    QCDBaseDftrig = pd.DataFrame()
    for fileName in DM.BGenPaths:
        tmpDf = DM.processData(fileName, 'BGen', 'Base', MC=True, trigger=True)
        QCDBaseDf = QCDBaseDf.append(tmpDf, **append_params)

    for fileName in DM.bEnrPaths:
        tmpDf = DM.processData(fileName, 'bEnr', 'Base', MC=True, trigger=True)
        QCDBaseDf = QCDBaseDf.append(tmpDf, **append_params)

    pickle.dump(QCDBaseDf, open(pickledir + 'QCDBaseDf.pkl', 'wb'))

else:
    QCDBaseDf = pickle.load(open(pickledir + 'QCDBaseDf.pkl', 'rb'))

QCD_base = EfficiencyInfo(demDf=QCDBaseDf, name='QCD (base selection)')
QCD_base.plot()

#------------------------------------------------------------------------------------
#---------------- ggH with parking selections and weights ---------------------

if root==True:
    ggHParkingDf = DM.processData(DM.ggHPaths, 'ggH', 'Parked', MC=True, trigger=True)
    pickle.dump(ggHParkingDf, open(pickledir + 'ggHParkingDf.pkl', 'wb'))

else:
    ggHParkingDf = pickle.load(open(pickledir + 'ggHParkingDf.pkl', 'rb'))

ggH_parking = EfficiencyInfo(demDf=ggHParkingDf, name='GGH (parking selection)')
ggH_parking.plot()

#--------------------------------------------------------------------------------
#---------------- ggH without muon selection or pT/IP/PU weights --------------


if root==True:
    ggHBaseDf = DM.processData(DM.ggHPaths, 'ggH', 'Base', MC=True, trigger=True)
    pickle.dump(ggHBaseDf, open(pickledir + 'ggHBaseDf.pkl', 'wb'))
else:
    ggHBaseDf = pickle.load(open(pickledir + 'ggHBaseDf.pkl', 'rb'))
ggH_base = EfficiencyInfo(demDf=ggHBaseDf, name='GGH (base selection)')
ggH_base.plot()

#--------------------------------------------------------------------------------
#---------------------- ttbar MC with muon selection --------------------------


if root==True:
    TTBarMuonEG = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'MuonEG', MC=True, trigger=True)
    pickle.dump(TTBarMuonEG, open(pickledir + 'TTBarMuonEG.pkl', 'wb'))
else:
    TTBarMuonEG = pickle.load(open(pickledir + 'TTBarMuonEG.pkl', 'rb'))

TTBar_muonEG = EfficiencyInfo(demDf=TTBarMuonEG, name='TTBar (MuonEG selection)')
TTBar_muonEG.plot()

#----------------------------------------------------------------------------------------
#----------------------- ttbar mc without muon selection ---------------------
if root==True:
    TTBarBaseDf = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'Base', MC=True, trigger=True)
    pickle.dump(TTBarBaseDf, open(pickledir + 'TTBarBaseDf.pkl', 'wb'))

else:
    TTBarBaseDf = pickle.load(open(pickledir + 'TTBarBaseDf.pkl', 'rb'))

TTBar_base = EfficiencyInfo(demDf=TTBarBaseDf, name='TTBar (base selection)')
TTBar_base.plot()

#--------------------------------------------------------------------------------



#------------------ overlays - parking + QCD --------------------------------------
figdir = 'JetHTTrigEff/plots/overlays'
plot_params = {'fmt':'o', 'capsize':5, 'linestyle':'None',}# 'xerr':'xerr'}
fig, ax = plt.subplots(figsize=(12,7.2))
xerr = np.ones((2, 27))*26
edge = parked.getBinCenter(parked.demEdge)
ax.grid()
ax.set_ylim([-0.05, 1.05])
ax.errorbar(edge, parked.quot, yerr=(parked.lowErr, parked.upErr),
            xerr=xerr,label='Parking BPH',**plot_params)
ax.errorbar(edge, QCD_base.quot, yerr=(QCD_base.lowErr, QCD_base.upErr),
            xerr=xerr, label='QCD (base)', **plot_params)
ax.errorbar(edge, QCD_parking.quot, yerr=(QCD_parking.lowErr, QCD_parking.upErr),
            xerr=xerr, label='QCD (parking)', **plot_params)
ax.set_title('Efficiency Parking BPH + QCD')
ax.set_xlabel('L1_SingleJet180 and HLT_AK8PFJet500')
ax.set_ylabel('ratio')
ax.legend(loc='upper left')
plt.savefig(f'{figdir}/Parking_QCD_Overlay.png')
plt.show()
plt.close()

#------------------- overlays - parking + ggH --------------------------------
fig, ax = plt.subplots(figsize=(12,7.2))
ax.errorbar(edge, parked.quot, yerr=(parked.lowErr, parked.upErr),
            xerr=xerr,label='Parking BPH',**plot_params)
ax.errorbar(edge, ggH_base.quot, yerr=(ggH_base.lowErr, ggH_base.upErr),
            xerr=xerr, label='GGH', **plot_params)
ax.errorbar(edge, ggH_parking.quot, yerr=(ggH_parking.lowErr, ggH_parking.upErr),
            xerr=xerr, label='GGH (parking)', **plot_params)
ax.set_title('Efficiency Parking BPH + GGH')
ax.set_xlabel('L1_SingleJet180 and HLT_AK8PFJet500')
ax.set_ylabel('ratio')
ax.legend(loc='upper left')
plt.savefig(f'{figdir}/Parking_ggH_Overlay.png')
plt.show()
plt.close()


#------------------- overlays = muonEG + ttbar ------------------------------
fig, ax = plt.subplots(figsize=(12,7.2))
ax.errorbar(edge, muonEG.quot, yerr=(muonEG.lowErr, muonEG.upErr),
            xerr=xerr, label='MuonEG', **plot_params)
ax.errorbar(edge, TTBar_base.quot, yerr=(TTBar_base.lowErr, TTBar_base.upErr),
            xerr=xerr, label='TTBar', **plot_params)
ax.errorbar(edge, TTBar_muonEG.quot, yerr=(TTBar_muonEG.lowErr, TTBar_muonEG.upErr),
            xerr=xerr, label='TTBar (muonEG)', **plot_params)
ax.set_title('Efficiency MuonEG + TTBar')
ax.set_xlabel('L1_SingleJet180 and HLT_AK8PFJet500')
ax.set_ylabel('ratio')
ax.legend(loc='upper left')
plt.savefig(f'{figdir}/muonEG_TTBar_Overlay.png')
plt.show()
plt.close()

#---------------------
fig, ax = plt.subplots(figsize=(12,7.2))
ax.errorbar(edge, QCD_base.quot, yerr=(QCD_base.lowErr, QCD_base.upErr),
            xerr=xerr, label='QCD', **plot_params)
ax.errorbar(edge, TTBar_base.quot, yerr=(TTBar_base.lowErr, TTBar_base.upErr),
            xerr=xerr, label='TTBar', **plot_params)
ax.errorbar(edge, ggH_base.quot, yerr=(ggH_base.lowErr, ggH_base.upErr),
            xerr=xerr, label='ggH', **plot_params)
ax.set_title('Efficiency QCD + TTBar + ggH')
ax.set_xlabel('L1_SingleJet180 and HLT_AK8PFJet500')
ax.set_ylabel('ratio')
ax.legend(loc='upper left')
plt.savefig(f'{figdir}/QCD_TTBar_ggH_Overlay.png')
plt.show()
plt.close()

#----------------------------------------------------------------------------

#pickle.dump(muonEG, open('JetHTTrigEff/pickles/histpickls/muonEG.pkl', 'wb'))
