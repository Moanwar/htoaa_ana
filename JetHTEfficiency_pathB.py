#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 20:05:41 2021

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
#from EfficiencyInfo import EfficiencyInfo
import mplhep as hep
#hep.set_style(hep.style.ROOT)
plt.style.use(hep.style.ROOT)
from EfficiencyInfoB import EfficiencyInfoB

def truncateDf(df, cut = 400):
    return df.loc[df.FatJet_pt > cut]


pickledir = 'JetHTTrigEff/pickles/'
root = True
hist_params = {'bins':50, 'range':(1000)}
append_params = {'ignore_index':True, 'sort':False}
trigger='B'

ptparams = {'var':'FatJet_pt', 'nbins':50, 'histrange':(0,1000)}
msoftparams = {'var':'FatJet_msoftdrop', 'nbins':11, 'histrange':(90,200)}
ddbvlparams = {'var':'FatJet_btagDDBvL', 'nbins':20, 'histrange':(0.8,1)}


# ------------------------- parking -------------------

if root:
    parkedDf = DM.processData(DM.ParkedDataPaths[0], 'data', 'Parked', MC=False, trigger=trigger)
    
    pickle.dump(parkedDf, open(pickledir + 'parkedDfB.pkl', 'wb'))

else:
    parkedDf = pickle.load(open(pickledir + 'parkedDfB.pkl', 'rb'))

parkedpt = EfficiencyInfoB(parkedDf, name='Parked data', **ptparams)
parkedpt.plot()

truncDf = truncateDf(parkedDf)
parked_msoft = EfficiencyInfoB(demDf=truncDf, name='Parked Data', **msoftparams)
parked_msoft.plot(ylim=[0.8,1.0])

parked_ddbvl = EfficiencyInfoB(demDf=truncDf, name='Parked Data', **ddbvlparams)
parked_ddbvl.plot(ylim=[0.8,1.0])

#%%

# -------------- MuonEG --------------------
if root:
    MuonEGDf = pd.DataFrame()
    for fileName in DM.MuonEGPaths:
        tmpDf = DM.processData(fileName, 'MuonEG', 'MuonEG', MC=False, trigger=trigger)
        MuonEGDf = MuonEGDf.append(tmpDf, **append_params)

    MuonEGDf['final_weights'] = 1
    pickle.dump(MuonEGDf, open(pickledir + 'MuonEGDfB.pkl', 'wb'))

else:
    MuonEGDf = pickle.load(open(pickledir + 'MuonEGDfB.pkl', 'rb'))

muonEG_pt = EfficiencyInfoB(demDf=MuonEGDf, name='MuonEG', **ptparams)
muonEG_pt.plot()

truncDf = truncateDf(MuonEGDf)
muonEG_msoft = EfficiencyInfoB(demDf=truncDf, name='MuonEG', **msoftparams)
muonEG_msoft.plot()

muonEG_ddbvl = EfficiencyInfoB(demDf=truncDf, name='MuonEG', **ddbvlparams)
muonEG_ddbvl.plot()

#--------------------- QCD MC with ParkingBPH selection and weights -----------
#%%
root=True
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

QCD_parking = EfficiencyInfoB(demDf=QCDParkingDf, name='QCD (parking selection)',
                              **ptparams)
QCD_parking.plot()

truncDf = truncateDf(QCDParkingDf)
QCD_parking_msoft = EfficiencyInfoB(demDf=truncDf, name='QCD (parking selection)',
                                    **msoftparams)
QCD_parking_msoft.plot(ylim=[0.8,1.0])

QCD_parking_ddbvl = EfficiencyInfoB(demDf=truncDf, name='QCD (parking selection)',
                                    **ddbvlparams)
QCD_parking_ddbvl.plot(ylim=[0.8,1.0])

#%%

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

QCD_base_pt = EfficiencyInfoB(demDf=QCDBaseDf, name='QCD (base selection)',
                          **ptparams)
QCD_base_pt.plot()

truncDf = truncateDf(QCDBaseDf)
QCD_base_msoft = EfficiencyInfoB(demDf=truncDf, name='QCD (base selection)',
                          **msoftparams)
QCD_base_msoft.plot(ylim=[0.8,1.0])

QCD_base_ddbvl = EfficiencyInfoB(demDf=truncDf, name='QCD (base selection)',
                          **ddbvlparams)
QCD_base_ddbvl.plot(ylim=[0.8,1.0])


#---------------- ggH with parking selections and weights ---------------------

if root==True:
    ggHParkingDf = DM.processData(DM.ggHPaths, 'ggH', 'Parked', MC=True, trigger=trigger)
    pickle.dump(ggHParkingDf, open(pickledir + 'ggHParkingDfB.pkl', 'wb'))

else:
    ggHParkingDf = pickle.load(open(pickledir + 'ggHParkingDfB.pkl', 'rb'))

ggH_parking_pt = EfficiencyInfoB(demDf=ggHParkingDf, name='GGH (parking selection)',
                              **ptparams)
ggH_parking_pt.plot()

truncDf = truncateDf(ggHParkingDf)

ggH_parking_msoft = EfficiencyInfoB(demDf=truncDf, name='GGH (parking selection)',
                              **msoftparams)
ggH_parking_msoft.plot(ylim=[0.8,1.0])

ggH_parking_ddbvl = EfficiencyInfoB(demDf=truncDf, name='GGH (parking selection)',
                              **ddbvlparams)
ggH_parking_ddbvl.plot(ylim=[0.8,1.0])


#---------------- ggH without muon selection or pT/IP/PU weights --------------


if root==True:
    ggHBaseDf = DM.processData(DM.ggHPaths, 'ggH', 'Base', MC=True, trigger=trigger)
    pickle.dump(ggHBaseDf, open(pickledir + 'ggHBaseDfB.pkl', 'wb'))
else:
    ggHBaseDf = pickle.load(open(pickledir + 'ggHBaseDfB.pkl', 'rb'))

ggH_base_pt = EfficiencyInfoB(demDf=ggHBaseDf, name='GGH (base selection)',
                              **ptparams)
ggH_base_pt.plot()

truncDf = truncateDf(ggHBaseDf)
ggH_base_msoft = EfficiencyInfoB(demDf=truncDf, name='GGH (base selection)',
                              **msoftparams)
ggH_base_msoft.plot(ylim=[0.8,1.0])

ggH_base_ddbvl = EfficiencyInfoB(demDf=truncDf, name='GGH (base selection)',
                              **ddbvlparams)
ggH_base_ddbvl.plot(ylim=[0.8,1.0])



#---------------------- ttbar MC with muon selection --------------------------

if root==True:
    TTBarMuonEG = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'MuonEG', MC=True, trigger=trigger)
    pickle.dump(TTBarMuonEG, open(pickledir + 'TTBarMuonEGB.pkl', 'wb'))
else:
    TTBarMuonEG = pickle.load(open(pickledir + 'TTBarMuonEGB.pkl', 'rb'))

TTBar_muonEG_pt = EfficiencyInfoB(demDf=TTBarMuonEG, name='TTBar (MuonEG selection)',
                               **ptparams)
TTBar_muonEG_pt.plot()

truncDf = truncateDf(TTBarMuonEG)
TTBar_muonEG_msoft = EfficiencyInfoB(demDf=truncDf, name='TTBar (MuonEG selection)',
                               **msoftparams)
TTBar_muonEG_msoft.plot()

TTBar_muonEG_ddbvl = EfficiencyInfoB(demDf=truncDf, name='TTBar (MuonEG selection)',
                               **ddbvlparams)
TTBar_muonEG_ddbvl.plot()



#----------------------- ttbar mc without muon selection ---------------------
if root==True:
    TTBarBaseDf = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'Base', MC=True, trigger=trigger)
    pickle.dump(TTBarBaseDf, open(pickledir + 'TTBarBaseDfB.pkl', 'wb'))

else:
    TTBarBaseDf = pickle.load(open(pickledir + 'TTBarBaseDfB.pkl', 'rb'))

TTBar_base_pt = EfficiencyInfoB(demDf=TTBarBaseDf, name='TTBar (base selection)',
                               **ptparams)
TTBar_base_pt.plot()

truncDf = truncateDf(TTBarBaseDf)
TTBar_base_msoft = EfficiencyInfoB(demDf=truncDf, name='TTBar (base selection)',
                               **msoftparams)
TTBar_base_msoft.plot(ylim=[0.8,1.0])

TTBar_base_ddbvl = EfficiencyInfoB(demDf=truncDf, name='TTBar (base selection)',
                               **ddbvlparams)
TTBar_base_ddbvl.plot(ylim=[0.8,1.0])





#------------------------------- overlays -----------------------------------
#------------------ overlays - parking + QCD --------------------------------------
#%%
if 'var' in ptparams:
    ptparams.pop('var')
    msoftparams.pop('var')
    ddbvlparams.pop('var')
def makeOverlay(effInfo, labels, title, xlabel, figname, ylim=[-0.05,1.05]):
    edge = effInfo[0].getBinCenter(effInfo[0].demEdge)
    halfbinlen = (effInfo[0].demEdge[1] - effInfo[0].demEdge[0])/2
    xerr = np.ones((2, len(effInfo[0].quot)))*halfbinlen
    fig, ax = plt.subplots(figsize=(12, 7.2))
    ax.grid()
    ax.set_ylim(ylim)
    plot_params = {'fmt':'o', 'capsize':5, 'linestyle':'None',}

    for i, eff in enumerate(effInfo):
        ax.errorbar(edge, eff.quot, yerr=(eff.lowErr, eff.upErr),
                    xerr = xerr, label=labels[i], **plot_params)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.legend(loc='lower right')
    plt.savefig(f'JetHTTrigEff/plots/B/{figname}')
    plt.show()
    plt.close()




## ---------------- ParkingBPH data + QCD parking + QCD ----------------------
parkingQCDlabels =['Parking', 'QCD (parking selection)', 'QCD (base selection)']
parkingQCDtitle = 'Parking + QCD (trigger B)'
pt = 'FatJet_pt'
msoft = 'FatJet_msoftdrop'
ddbvl = 'FatJet_btagDDBvL'
makeOverlay(effInfo=[parkedpt, QCD_parking, QCD_base_pt],
            labels= parkingQCDlabels,
            title=parkingQCDtitle,
            xlabel=pt,
            figname = 'Parking_QCD_ptOverlay')

makeOverlay(effInfo=[parked_msoft, QCD_parking_msoft, QCD_base_msoft],
            labels=parkingQCDlabels,
            title=parkingQCDtitle,
            xlabel=msoft,
            figname = 'Parking_QCD_msoftOverlay',
            ylim=[0.8,1])

makeOverlay(effInfo=[parked_ddbvl, QCD_parking_ddbvl, QCD_base_ddbvl],
            labels = parkingQCDlabels,
            title=parkingQCDtitle,
            xlabel=ddbvl,
            figname = 'Parking_QCD_ddbvlOverlay')

## ---------------- parking + ggH (parking) + ggH (base) -------------------
parkingggHlabels = ['Parking', 'ggH (parking selection)' ,'ggH (base selection)']
parkingggHtitle = 'Parking + ggH (trigger B)'
makeOverlay(effInfo=[parkedpt, ggH_parking_pt, ggH_base_pt],
            labels=parkingggHlabels,
            title=parkingggHtitle,
            xlabel=pt,
            figname = 'Parking_ggH_ptOverlay')
makeOverlay(effInfo=[parked_msoft, ggH_parking_msoft, ggH_base_msoft],
            labels=parkingggHlabels,
            title=parkingggHtitle,
            xlabel=msoft,
            figname = 'Parking_ggH_msoftOverlay',
            ylim=[0.8,1])
makeOverlay(effInfo=[parked_ddbvl, ggH_parking_ddbvl, ggH_base_ddbvl],
            labels=parkingggHlabels,
            title=parkingggHtitle,
            xlabel=ddbvl,
            figname = 'Parking_ggH_ddbvlOverlay')

## -------------- MuonEG + ttBar (MuonEG) + ttBar (base) -----------------
muonegttbarlabels = ['MuonEG', 'TTBar (MuonEG selection)', 'TTBar (base selection)']
muonegttbartitle = 'MuonEG + TTBar (trigger B)'
makeOverlay(effInfo=[muonEG_pt, TTBar_muonEG_pt, TTBar_base_pt],
            labels=muonegttbarlabels,
            title=muonegttbartitle,
            xlabel=pt,
            figname='MuonEG_TTBar_ptOverlay')
makeOverlay(effInfo=[muonEG_msoft, TTBar_muonEG_msoft, TTBar_base_msoft],
            labels=muonegttbarlabels,
            title=muonegttbartitle,
            xlabel=msoft,
            figname='MuonEG_TTbar_msoftOverlay')
makeOverlay(effInfo=[muonEG_ddbvl, TTBar_muonEG_ddbvl, TTBar_base_ddbvl],
            labels=muonegttbarlabels,
            title=muonegttbartitle,
            xlabel=ddbvl,
            figname='MuonEG_TTBar_ddbvlOverlay')

## ----------------- qcd + ttbar + ggh -----------------------------
qcdttbargghlabels = ['QCD', 'TTBar', 'GGH']
qcdttbargghtitle = 'QCD + TTBar + GGH (trigger B)'
makeOverlay(effInfo=[QCD_base_pt, TTBar_base_pt, ggH_base_pt],
            labels=qcdttbargghlabels,
            title=qcdttbargghtitle,
            xlabel=pt,
            figname='QCD_TTBar_GGH_ptOverlay')
makeOverlay(effInfo=[QCD_base_msoft, TTBar_base_msoft, ggH_base_msoft],
            labels=qcdttbargghlabels,
            title=qcdttbargghtitle,
            xlabel=msoft,
            figname='QCD_TTBar_GGH_msoftOverlay',
            ylim=[0.8,1])
makeOverlay(effInfo=[QCD_base_ddbvl, TTBar_base_ddbvl, ggH_base_ddbvl],
            labels=qcdttbargghlabels,
            title=qcdttbargghtitle,
            xlabel=ddbvl,
            figname='QCD_TTBar_GGH_ddbvlOverlay',
            ylim=[0.8,1])



#%%
## ----------------------- scale factor plots --------------------------
from ScaleFactor import ScaleFactor

bparktruncDf = truncateDf(parkedDf, cut=400)
bqcdtruncDf= truncateDf(QCDParkingDf, cut=400)

bparkEff = EfficiencyInfoB(demDf=bparktruncDf, name='Parked',
                           var='FatJet_pt', nbins=50, histrange=(400,1000))
bqcdEff = EfficiencyInfoB(demDf=bqcdtruncDf, name='QCD(Parking selection)',
                          var='FatJet_pt', nbins=50, histrange=(400,1000))

Bscale = ScaleFactor(bparkEff, bqcdEff)
Bscale.plot(title='Scale Factor(BoostedDoubleB)', xlabel='FatJet_pt', figdir='JetHTTrigEff/plots/B/')



