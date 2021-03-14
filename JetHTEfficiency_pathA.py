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
from EfficiencyInfo import EfficiencyInfo






pickledir = 'JetHTTrigEff/pickles/'
root = True
hist_params = {'bins':27, 'range':(150,1500)}
append_params = {'ignore_index':True, 'sort':False}
trigger='A'

#%%
#-------------------------- ParkingBPH data -----------------------------------

if root:
    parkedDf = DM.processData(DM.ParkedDataPaths[0], 'data', 'Parked', MC=False, trigger=trigger)
    
    pickle.dump(parkedDf, open(pickledir + 'parkedDfB.pkl', 'wb'))

else:
    parkedDf = pickle.load(open(pickledir + 'parkedDfB.pkl', 'rb'))

parked = EfficiencyInfo(parkedDf, 'Parking data')
parked.plotpt()
#%%

#------------------------------------------------------------------------------
#------------------------------ MuonEG data -----------------------------------

if root:
    MuonEGDf = pd.DataFrame()
    for fileName in DM.MuonEGPaths:
        tmpDf = DM.processData(fileName, 'MuonEG', 'MuonEG', MC=False, trigger=trigger)
        MuonEGDf = MuonEGDf.append(tmpDf, **append_params)

    MuonEGDf['final_weights'] = 1
    pickle.dump(MuonEGDf, open(pickledir + 'MuonEGDfB.pkl', 'wb'))

else:
    MuonEGDf = pickle.load(open(pickledir + 'MuonEGDfB.pkl', 'rb'))


muonEG = EfficiencyInfo(demDf=MuonEGDf, name='MuonEG')
muonEG.plotpt()

#-------------------------------------------------------------------------
#--------------------- QCD MC with ParkingBPH selection and weights -----------
#%%
if root:
    QCDParkingDf = pd.DataFrame()
    QCDParkingDftrig = pd.DataFrame()
    for fileName in DM.BGenPaths:
        tmpDf = DM.processData(fileName, 'BGen', 'Parked', MC=True, trigger=trigger)
        QCDParkingDf = QCDParkingDf.append(tmpDf, **append_params)

    for fileName in DM.bEnrPaths:
        tmpDf = DM.processData(fileName, 'bEnr', 'Parked', MC=True, trigger=trigger)
        QCDParkingDf = QCDParkingDf.append(tmpDf, **append_params)

    pickle.dump(QCDParkingDf, open(pickledir + 'QCDParkingDfB.pkl', 'wb'))

else:
    QCDParkingDf = pickle.load(open(pickledir + 'QCDParkingDfB.pkl', 'rb'))

QCD_parking = EfficiencyInfo(demDf=QCDParkingDf, name='QCD (parking selection)')
QCD_parking.plotpt()
#%%
#-------------------------------------------------------------------------
#--------------- QCD MC with no offline muon or pt/IP/PU weights -------------

if root:
    QCDBaseDf = pd.DataFrame()
    QCDBaseDftrig = pd.DataFrame()
    for fileName in DM.BGenPaths:
        tmpDf = DM.processData(fileName, 'BGen', 'Base', MC=True, trigger=trigger)
        QCDBaseDf = QCDBaseDf.append(tmpDf, **append_params)

    for fileName in DM.bEnrPaths:
        tmpDf = DM.processData(fileName, 'bEnr', 'Base', MC=True, trigger=trigger)
        QCDBaseDf = QCDBaseDf.append(tmpDf, **append_params)

    pickle.dump(QCDBaseDf, open(pickledir + 'QCDBaseDfB.pkl', 'wb'))

else:
    QCDBaseDf = pickle.load(open(pickledir + 'QCDBaseDfB.pkl', 'rb'))

QCD_base = EfficiencyInfo(demDf=QCDBaseDf, name='QCD (base selection)')
QCD_base.plotpt()

#------------------------------------------------------------------------------------
#---------------- ggH with parking selections and weights ---------------------

if root==True:
    ggHParkingDf = DM.processData(DM.ggHPaths, 'ggH', 'Parked', MC=True, trigger=trigger)
    pickle.dump(ggHParkingDf, open(pickledir + 'ggHParkingDfB.pkl', 'wb'))

else:
    ggHParkingDf = pickle.load(open(pickledir + 'ggHParkingDfB.pkl', 'rb'))

ggH_parking = EfficiencyInfo(demDf=ggHParkingDf, name='GGH (parking selection)')
ggH_parking.plotpt()

#--------------------------------------------------------------------------------
#---------------- ggH without muon selection or pT/IP/PU weights --------------


if root==True:
    ggHBaseDf = DM.processData(DM.ggHPaths, 'ggH', 'Base', MC=True, trigger=trigger)
    pickle.dump(ggHBaseDf, open(pickledir + 'ggHBaseDfB.pkl', 'wb'))
else:
    ggHBaseDf = pickle.load(open(pickledir + 'ggHBaseDfB.pkl', 'rb'))
ggH_base = EfficiencyInfo(demDf=ggHBaseDf, name='GGH (base selection)')
ggH_base.plotpt()

#--------------------------------------------------------------------------------
#---------------------- ttbar MC with muon selection --------------------------


if root==True:
    TTBarMuonEG = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'MuonEG', MC=True, trigger=trigger)
    pickle.dump(TTBarMuonEG, open(pickledir + 'TTBarMuonEGB.pkl', 'wb'))
else:
    TTBarMuonEG = pickle.load(open(pickledir + 'TTBarMuonEGB.pkl', 'rb'))

TTBar_muonEG = EfficiencyInfo(demDf=TTBarMuonEG, name='TTBar (MuonEG selection)')
TTBar_muonEG.plotpt()

#----------------------------------------------------------------------------------------
#----------------------- ttbar mc without muon selection ---------------------
if root==True:
    TTBarBaseDf = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'Base', MC=True, trigger=trigger)
    pickle.dump(TTBarBaseDf, open(pickledir + 'TTBarBaseDfB.pkl', 'wb'))

else:
    TTBarBaseDf = pickle.load(open(pickledir + 'TTBarBaseDfB.pkl', 'rb'))

TTBar_base = EfficiencyInfo(demDf=TTBarBaseDf, name='TTBar (base selection)')
TTBar_base.plotpt()

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

pickle.dump(muonEG, open('JetHTTrigEff/pickles/histpickls/muonEGB.pkl', 'wb'))




#%%
# ---------------------------- scale factor ---------------------------------
from ScaleFactor import ScaleFactor


def truncateDf(df, cut = 550):
    return df.loc[df.FatJet_pt > cut]


aparktruncDf = truncateDf(parkedDf, cut=550)
aqcdtruncDf= truncateDf(QCDParkingDf, cut=550)

aparkEff = EfficiencyInfo(demDf=aparktruncDf, name='Parked',
                           var='FatJet_pt', nbins=50, histrange=(550,1000))
aqcdEff = EfficiencyInfo(demDf=aqcdtruncDf, name='QCD(Parking selection)',
                          var='FatJet_pt', nbins=50, histrange=(550,1000))

Bscale = ScaleFactor(aparkEff, aqcdEff)
Bscale.plot(title='Scale Factor (AK8PFJet500)', xlabel='FatJet_pt', figdir='JetHTTrigEff/plots/A/')
