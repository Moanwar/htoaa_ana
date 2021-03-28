#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 14:12:04 2021

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

def truncateDf(df):
    return df.loc[df.FatJet_pt > 400]

class EfficiencyInfoC(EfficiencyInfoB):
    def __init__(self, demDf, name, var, nbins, histrange):
        #EfficiencyInfoB.__init__(demDf, name, var, nbins, histrange)
        super().__init__(demDf, name, var, nbins, histrange)
        self.plotdir = 'JetHTTrigEff/plots/C/'

    def getNumDf(self, df):
        trig = np.logical_and(np.logical_or(df.L1_DoubleJet112er2p3_dEta_Max1p6,
                                            df.L1_DoubleJet150er2p5),
                              df.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71)
        # trig = ((df.L1_DoubleJet112er2p3_dEta_Max1p6 |
        #                 df.L1_DoubleJet150er2p5) &
        #                df.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71)
        ret = df[trig]
        return ret
        #return df.loc[(df.L1_DoubleJet112er2p3_dEta_Max1p6 ==True |
        #                df.L1_DoubleJet150er2p5==True) &
         #              df.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71==True]

# L1_DoubleJet112er2p3_dEta_Max1p6 or L1_DoubleJet150er2p5) and HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71

def btagCut(df):
    return df[(df.Jet_btagDeepB1 > 0.4184) & (df.Jet_btagDeepB2 > 0.4184)]

def ptCut(df):
    return df[(df.Jet_pt1 > 150) & (df.Jet_pt2 > 150)]

def bothCut(df):
    return df[(df.Jet_btagDeepB1 > 0.4184) &
              (df.Jet_btagDeepB2 > 0.4184) &
              (df.Jet_pt1 > 140) &
              (df.Jet_pt2 > 140)]

def reassignBtag(df):
    btags = pd.DataFrame(df.Jet_btagDeepB1, columns=['Jet_btagDeepB1'])
    btags = btags.assign(Jet_btagDeepB2=df.Jet_btagDeepB2)
    df = df.assign(Jet_btagDeepB1=btags.max(axis=1))
    df = df.assign(Jet_btagDeepB2=btags.min(axis=1))
    return df

#--------------- parkeddata -----------------------

pickledir = 'JetHTTrigEff/pickles/C/'
root = False
hist_params = {'bins':50, 'range':(600)}
append_params = {'ignore_index':True, 'sort':False}
trigger='C'

ptparams = {'var':'Jet_pt2', 'nbins':61, 'histrange':(30,1050)}
#msoftparams = {'var':'FatJet_msoftdrop', 'nbins':11, 'histrange':(90,200)}
btagparams = {'var':'Jet_btagDeepB2', 'nbins':20, 'histrange':(0,1)}
fatjetparams = {'var':'FatJet_pt', 'nbins':27, 'histrange':(150,1500)}

#%%
if root:
    parkedDf = DM.processData(DM.ParkedDataPaths[0], 'data', 'Parked', MC=False, trigger=trigger)
    
    pickle.dump(parkedDf, open(pickledir + 'parkedDf.pkl', 'wb'))

else:
    parkedDf = pickle.load(open(pickledir + 'parkedDf.pkl', 'rb'))

parkedDf_btagCut = btagCut(parkedDf)
parkedpt = EfficiencyInfoC(parkedDf_btagCut, name='Parked data', **ptparams)
parkedpt.plot(title='trigger C')


## reassign btag so that instead of btag associated with ptjets, it's higher, and lower btag
## si make it a function
parkedDf_ptCut = ptCut(parkedDf)
parkedDf_ptCut = reassignBtag(parkedDf_ptCut)
parkedbtag = EfficiencyInfoC(parkedDf_ptCut, name='Parked data', **btagparams)
parkedbtag.plot(title='trigger C')

parkedDf_ptbtagCut = bothCut(parkedDf)
parkedfatpt_cut = EfficiencyInfoC(parkedDf_ptbtagCut, name = 'Parked data ( pt, btag cut)',
                                 **fatjetparams)
parkedfatpt_cut.plot(title='trigger C')

 #%%

# -------------- MuonEG --------------------
root = False
if root:
    MuonEGDf = pd.DataFrame()
    for fileName in DM.MuonEGPaths:
        tmpDf = DM.processData(fileName, 'MuonEG', 'MuonEG', MC=False, trigger=trigger)
        MuonEGDf = MuonEGDf.append(tmpDf, **append_params)

    MuonEGDf['final_weights'] = 1
    pickle.dump(MuonEGDf, open(pickledir + 'MuonEGDfB.pkl', 'wb'))

else:
    MuonEGDf = pickle.load(open(pickledir + 'MuonEGDfB.pkl', 'rb'))

muonEG_btagCut = btagCut(MuonEGDf)
muonEGpt = EfficiencyInfoC(demDf=muonEG_btagCut, name='MuonEG', **ptparams)
muonEGpt.plot(title='trigger C')

muonEG_ptCut = ptCut(MuonEGDf)
muonEG_ptCut = reassignBtag(muonEG_ptCut)
muonEGbtag = EfficiencyInfoC(demDf=muonEG_ptCut, name='MuonEG', **btagparams)
muonEGbtag.plot(title='trigger C')


muonEG_ptbtagCut = bothCut(MuonEGDf)
muonEGfatpt_cut = EfficiencyInfoC(muonEG_ptbtagCut, name = 'MuonEG ( pt, btag cut)',
                                 **fatjetparams)
muonEGfatpt_cut.plot(title='trigger C')
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

QCDParkingDf_btagCut = btagCut(QCDParkingDf)
QCD_parking = EfficiencyInfoC(demDf=QCDParkingDf_btagCut, name='QCD (parking selection)',
                              **ptparams)
QCD_parking.plot(title='trigger C')

QCDParkingDf_ptCut = ptCut(QCDParkingDf)
QCDParkingDf_ptCut = reassignBtag(QCDParkingDf_ptCut)
QCD_parking_btag = EfficiencyInfoC(demDf=QCDParkingDf_ptCut, name='QCD (parking selection)',
                                    **btagparams)
QCD_parking_btag.plot(title='trigger C')

QCDParkingDf_ptbtagCut = bothCut(QCDParkingDf)
QCD_parkingfatpt_cut = EfficiencyInfoC(QCDParkingDf_ptbtagCut, name = 'QCD (parking selection) ( pt, btag cut)',
                                 **fatjetparams)
QCD_parkingfatpt_cut.plot(title='trigger C')
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

QCDBaseDf_btagCut = btagCut(QCDBaseDf)
QCD_base_pt = EfficiencyInfoC(demDf=QCDBaseDf_btagCut, name='QCD (base selection)',
                          **ptparams)
QCD_base_pt.plot(title='trigger C')

QCDBaseDf_ptCut = ptCut(QCDBaseDf)
QCDBaseDf_ptCut = reassignBtag(QCDBaseDf_ptCut)
QCD_base_btag = EfficiencyInfoC(demDf=QCDBaseDf_ptCut, name='QCD (base selection)',
                          **btagparams)
QCD_base_btag.plot(title='trigger C')

QCDBaseDf_ptbtagCut = bothCut(QCDBaseDf)
QCD_parkingfatpt_cut = EfficiencyInfoC(QCDBaseDf_ptbtagCut, name = 'QCD (base selection) ( pt, btag cut)',
                                 **fatjetparams)
QCD_parkingfatpt_cut.plot(title='trigger C')
#---------------- ggH with parking selections and weights ---------------------
#%%
root=False
if root==True:
    ggHParkingDf = DM.processData(DM.ggHPaths, tag='ggH', dataSet='Parked', MC=True, trigger=trigger)
    pickle.dump(ggHParkingDf, open(pickledir + 'ggHParkingDfB.pkl', 'wb'))

else:
    ggHParkingDf = pickle.load(open(pickledir + 'ggHParkingDfB.pkl', 'rb'))

ggHParkingDf_btagCut = btagCut(ggHParkingDf)
ggH_parking_pt = EfficiencyInfoC(demDf=ggHParkingDf_btagCut, name='GGH (parking selection)',
                              **ptparams)
ggH_parking_pt.plot(title='trigger C')

ggHParkingDf_ptCut = ptCut(ggHParkingDf)
#ggHParkingDf_ptCut = reassignBtag(ggHParkingDf_ptCut)
btags = pd.DataFrame(ggHParkingDf_ptCut.Jet_btagDeepB1)
btags = btags.assign(two= ggHParkingDf_ptCut.Jet_btagDeepB2)
ggHParkingDf_ptCut = ggHParkingDf_ptCut.assign(Jet_btagDeepB2 = btags.min(axis=1))
ggH_parking_btag = EfficiencyInfoC(demDf=ggHParkingDf_ptCut, name='GGH (parking selection)',
                              **btagparams)
ggH_parking_btag.plot(title='trigger C')

ggHParkingDf_ptbtagCut = bothCut(ggHParkingDf)
ggH_parkingfatpt_cut = EfficiencyInfoC(ggHParkingDf_ptbtagCut, name = 'GGH (parking selection) ( pt, btag cut)',
                                 **fatjetparams)
ggH_parkingfatpt_cut.plot(title='trigger C')
#---------------- ggH without muon selection or pT/IP/PU weights --------------
#%%
if root==True:
    ggHBaseDf = DM.processData(DM.ggHPaths, 'ggH', 'Base', MC=True, trigger=trigger)
    pickle.dump(ggHBaseDf, open(pickledir + 'ggHBaseDf.pkl', 'wb'))
else:
    ggHBaseDf = pickle.load(open(pickledir + 'ggHBaseDf.pkl', 'rb'))

ggHBaseDf_btagCut = btagCut(ggHBaseDf)
ggH_base_pt = EfficiencyInfoC(demDf=ggHBaseDf_btagCut, name='GGH (base selection)',
                              **ptparams)
ggH_base_pt.plot(title='trigger C')

ggHBaseDf_ptCut = ptCut(ggHBaseDf)
ggHBaseDf_ptCut = reassignBtag(ggHBaseDf_ptCut)
# btags = pd.DataFrame(ggHBaseDf_ptCut.Jet_btagDeepB1, columns=['Jet_btagDeepB1'])
# btags.loc[:,'Jet_btagDeepB2'] = ggHBaseDf_ptCut.Jet_btagDeepB2
# ggHBaseDf_ptCut = ggHBaseDf_ptCut.assign(Jet_btagDeepB1=btags.max(axis=1))
# ggHBaseDf_ptCut = ggHBaseDf_ptCut.assign(Jet_btagDeepB2=btags.min(axis=1))

ggH_base_btag = EfficiencyInfoC(demDf=ggHBaseDf_ptCut, name='GGH (base selection)',
                              **btagparams)
ggH_base_btag.plot(title='trigger C')

ggHBaseDf_ptbtagCut = bothCut(ggHBaseDf)
ggH_basefatpt_cut = EfficiencyInfoC(ggHBaseDf_ptbtagCut, name = 'GGH (base selection) (pt, btag cut)',
                                    var='FatJet_pt',nbins=27, histrange=(150,1500))#**fatjetparams)
ggH_basefatpt_cut.plot(title='trigger C')
#%%
#---------------------- ttbar MC with muon selection --------------------------

if root==True:
    TTBarMuonEG = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'MuonEG', MC=True, trigger=trigger)
    pickle.dump(TTBarMuonEG, open(pickledir + 'TTBarMuonEGB.pkl', 'wb'))
else:
    TTBarMuonEG = pickle.load(open(pickledir + 'TTBarMuonEGB.pkl', 'rb'))

TTBarMuonEG_btagCut = btagCut(TTBarMuonEG)
TTBar_muonEG_pt = EfficiencyInfoC(demDf=TTBarMuonEG_btagCut, name='TTBar (MuonEG selection)',
                               **ptparams)
TTBar_muonEG_pt.plot(title='trigger C')

TTBarMuonEG_ptCut = ptCut(TTBarMuonEG)
TTBarMuonEG_ptCut = reassignBtag(TTBarMuonEG_ptCut)
TTBar_muonEG_btag = EfficiencyInfoC(demDf=TTBarMuonEG_ptCut, name='TTBar (MuonEG selection)',
                               **btagparams)
TTBar_muonEG_btag.plot(title='trigger C')

TTBarMuonEG_ptbtagCut = bothCut(TTBarMuonEG)
TTBar_muonEGfatpt_cut = EfficiencyInfoC(TTBarMuonEG_ptbtagCut, name = 'TTBar (MuonEG selection) ( pt, btag cut)',
                                 **fatjetparams)
TTBar_muonEGfatpt_cut.plot(title='trigger C')

#----------------------- ttbar mc without muon selection ---------------------
if root==True:
    TTBarBaseDf = DM.processData(DM.TTJetsPaths[0], 'TTJets', 'Base', MC=True, trigger=trigger)
    pickle.dump(TTBarBaseDf, open(pickledir + 'TTBarBaseDfB.pkl', 'wb'))

else:
    TTBarBaseDf = pickle.load(open(pickledir + 'TTBarBaseDfB.pkl', 'rb'))

TTBarBaseDf_btagCut = btagCut(TTBarBaseDf)
TTBar_base_pt = EfficiencyInfoC(demDf=TTBarBaseDf_btagCut, name='TTBar (base selection)',
                               **ptparams)
TTBar_base_pt.plot(title = 'trigger C')

TTBarBaseDf_ptCut = ptCut(TTBarBaseDf)
TTBarBaseDf_ptCut = reassignBtag(TTBarBaseDf_ptCut)
TTBar_base_btag = EfficiencyInfoC(demDf=TTBarBaseDf_ptCut, name='TTBar (base selection)',
                               **btagparams)
TTBar_base_btag.plot(title='trigger C')


TTBarBaseDf_ptbtagCut = bothCut(TTBarBaseDf)
TTBarBasefatpt_cut = EfficiencyInfoC(TTBarBaseDf_ptbtagCut, name = 'TTBar (base selection) ( pt, btag cut)',
                                 **fatjetparams)
TTBarBasefatpt_cut.plot(title='trigger C')

#------------------------ overlays ----------------------------
#------------------ overlays - parking + QCD ----------------------------------
#%%
if 'var' in ptparams:
    ptparams.pop('var')
    btagparams.pop('var')

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
    plt.savefig(f'JetHTTrigEff/plots/C/{figname}.png')
    plt.show()
    plt.close()

parkingQCDlabels =['Parking', 'QCD (parking selection)', 'QCD (base selection)']
parkingQCDtitle = 'Parking + QCD (trigger C)'
pt = 'Jet_pt'
btag = 'Jet_btagDeepB'

makeOverlay(effInfo=[parkedpt, QCD_parking, QCD_base_pt],
            labels= parkingQCDlabels,
            title=parkingQCDtitle,
            xlabel=pt,
            figname = 'Parking_QCD_ptOverlay')
makeOverlay(effInfo=[parkedbtag, QCD_parking_btag, QCD_base_btag],
            labels= parkingQCDlabels,
            title=parkingQCDtitle,
            xlabel=btag,
            figname = 'Parking_QCD_btagOverlay')


## ---------------- parking + ggH (parking) + ggH (base) -------------------
parkingggHlabels = ['Parking', 'ggH (parking selection)' ,'ggH (base selection)']
parkingggHtitle = 'Parking + ggH (trigger C)'
makeOverlay(effInfo=[parkedpt, ggH_parking_pt, ggH_base_pt],
            labels=parkingggHlabels,
            title=parkingggHtitle,
            xlabel=pt,
            figname = 'Parking_ggH_ptOverlay')
makeOverlay(effInfo=[parkedbtag, ggH_parking_btag, ggH_base_btag],
            labels=parkingggHlabels,
            title=parkingggHtitle,
            xlabel=btag,
            figname = 'Parking_ggH_btagOverlay',
            )

## -------------- MuonEG + ttBar (MuonEG) + ttBar (base) -----------------
muonegttbarlabels = ['MuonEG', 'TTBar (MuonEG selection)', 'TTBar (base selection)']
muonegttbartitle = 'MuonEG + TTBar (trigger C)'
makeOverlay(effInfo=[muonEGpt, TTBar_muonEG_pt, TTBar_base_pt],
            labels=muonegttbarlabels,
            title=muonegttbartitle,
            xlabel=pt,
            figname='MuonEG_TTBar_ptOverlay')
makeOverlay(effInfo=[muonEGbtag, TTBar_muonEG_btag, TTBar_base_btag],
            labels=muonegttbarlabels,
            title=muonegttbartitle,
            xlabel=btag,
            figname='MuonEG_TTbar_msoftOverlay')


## ----------------- qcd + ttbar + ggh -----------------------------
qcdttbargghlabels = ['QCD', 'TTBar', 'GGH']
qcdttbargghtitle = 'QCD + TTBar + GGH (trigger C)'
makeOverlay(effInfo=[QCD_base_pt, TTBar_base_pt, ggH_base_pt],
            labels=qcdttbargghlabels,
            title=qcdttbargghtitle,
            xlabel=pt,
            figname='QCD_TTBar_GGH_ptOverlay')
makeOverlay(effInfo=[QCD_base_btag, TTBar_base_btag, ggH_base_btag],
            labels=qcdttbargghlabels,
            title=qcdttbargghtitle,
            xlabel=btag,
            figname='QCD_TTBar_GGH_msoftOverlay')