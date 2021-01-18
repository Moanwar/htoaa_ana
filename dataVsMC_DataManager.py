#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 19:33:18 2020

@author: si_sutantawibul1
"""

import uproot
import pandas as pd
from analib import PhysObj, Event
import sys
import os
import pickle
import numpy as np
from dataManager import getSubJetData, getnSVCounts


JetHT = True

dataPath = 'data/2018D_Parked_promptD-v1_200218_214714_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root'
ggHPath = 'MC/nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId_999k.root'
#ggHPath = 'GGH_HPT.root'
BGenPath = 'MC/QCD_BGen_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root'
bEnrPath = 'MC/QCD_bEnriched_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root'

BGenPaths = ['QCD_BGenFilter/QCD_HT200to300.root',
             'QCD_BGenFilter/QCD_HT300to500.root']#,
             # 'QCD_BGenFilter/QCD_HT500to700.root',
             # 'QCD_BGenFilter/QCD_HT700to1000.root',
             # 'QCD_BGenFilter/QCD_HT1000to1500.root',
             # 'QCD_BGenFilter/QCD_HT1500to2000.root',
             # 'QCD_BGenFilter/QCD_HT2000toInf.root']

bEnrPaths = ['QCD_bEnriched/QCD_HT200to300.root',
             'QCD_bEnriched/QCD_HT300to500.root']#,
             # 'QCD_bEnriched/QCD_HT500to700.root',
             # 'QCD_bEnriched/QCD_HT700to1000.root',
             # 'QCD_bEnriched/QCD_HT1000to1500.root',
             # 'QCD_bEnriched/QCD_HT1500to2000.root',
             # 'QCD_bEnriched/QCD_HT2000toInf.root']

TTJetsPaths = ['MC/TTJets/TTJets_Skim_nFat1_doubB_0p8_deepB_Med_massH_70_msoft_70.root',
              'MC/TTJets/TTJets_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240.root']#,
              #'MC/TTJets/TTJets_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root']

WJetsPaths = ['MC/WJets/WJets_HT-800toInf_Skim_nFat1_doubB_0p8_deepB_Med_massH_70_msoft_70.root',
              'MC/WJets/WJets_HT-800toInf_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240.root']#,
              # 'MC/WJets/WJets_HT400to600_Skim_nFat1_doubB_0p8_deepB_Med_massH_70_msoft_70.root',
              # 'MC/WJets/WJets_HT400to600_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240.root',
              # 'MC/WJets/WJets_HT600to800_Skim_nFat1_doubB_0p8_deepB_Med_massH_70_msoft_70.root',
              # 'MC/WJets/WJets_HT600to800_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240.root',
              # 'MC/WJets/WJets_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root',]

ZJetsPaths = ['MC/ZJets/ZJets_HT-800toInf_Skim_nFat1_doubB_0p8_deepB_Med_massH_70_msoft_70.root',
              'MC/ZJets/ZJets_HT-800toInf_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240.root']#,
              # 'MC/ZJets/ZJets_HT400to600_Skim_nFat1_doubB_0p8_deepB_Med_massH_70_msoft_70.root',
              # 'MC/ZJets/ZJets_HT400to600_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240.root',
              # 'MC/ZJets/ZJets_HT600to800_Skim_nFat1_doubB_0p8_deepB_Med_massH_70_msoft_70.root',
              # 'MC/ZJets/ZJets_HT600to800_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240.root',
              # 'MC/ZJets/ZJets_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root',
              # ]

JetHTPaths = ['data/JetHT/2018A_JetHT_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240.root',
              'data/JetHT/2018B_JetHT_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240.root']#,
              #'data/JetHT/2018C_JetHT_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240.root',
              #'data/JetHT/2018D_JetHT_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240.root',
              #]


BGenWeight = [1, 0.259, 0.0515, 0.01666, 0.00905, 0.003594, 0.001401]
bEnrWeight =[ 1, 0.33, 0.034, 0.034, 0.0024, 0.00024, 0.00044]

jetVars = ['FatJet_pt',
           'FatJet_eta',
           'FatJet_mass',
           'FatJet_msoftdrop',
           'FatJet_btagCSVV2',
           'FatJet_btagDeepB',
           'FatJet_msoftdrop',
           'FatJet_btagDDBvL',
           'FatJet_deepTagMD_H4qvsQCD']

muonVars = ['Muon_pt',
            'Muon_eta',
            'Muon_ip3d',
            'Muon_softId']

PVVars = ['PV_npvs', 'PV_npvsGood']

allVars = list(jetVars + muonVars + PVVars + ['LHE_HT'])
allVars.sort()

muonR = pickle.load(open('muontensor/MuonRtensor.p', 'rb'))
muonL = pickle.load(open('muontensor/MuonLtensor.p', 'rb'))

ptkeys = list(muonL.keys()) + [999999]
ptkeys.remove('meta')
ipkeys = list(muonL[ptkeys[0]].keys()) + [999999]
npvsGkeys = muonR[6][2]['H']

tagslist = ['bEnr', 'BGen', 'data', 'JetHT', 'WJets', 'TTJets', 'ZJets']


def getMaxPt(physobj, col, varlist):
    colidx = physobj[col].idxmax(axis=1).to_numpy()
    rowidx = list(range(len(colidx)))
    maxPtData = pd.DataFrame()

    for var in varlist:
        npArr = physobj[var].to_numpy()

        maxPtData[var] = npArr[rowidx, colidx]
    #maxPtData = maxPtData.reindex(physobj[col].index)
    return maxPtData

## I don't have the getSubjetData here becuase i imported the dataManager one because yeah

def processData (filePath, tag): #JetHT=False):
    ## open file, get events
    fileName, fileExtension = os.path.splitext(filePath)
    
    print(filePath)

    if fileExtension != '.root':
        print('this program only supports .root  files')
        sys.exit()

    if tag not in tagslist:
        print('check yo tags')
        sys.exit()

    f = uproot.open(fileName + '.root')
    events = f.get('Events')

    jets = PhysObj('jets' + fileName)
    muons = PhysObj('muons' + fileName)
    other = PhysObj('other' + fileName)

    ## data doens't have LHE_HT
    if tag == 'data':
        allVars.remove('LHE_HT')

    ## fill the PhysObjs with data from the root file
    ## fatjets vars
    jets['FatJet_pt'] = pd.DataFrame(events.array('FatJet_pt'))
    jets['FatJet_eta'] = pd.DataFrame(np.abs(events.array('FatJet_eta')))
    jets['FatJet_mass'] = pd.DataFrame(events.array('FatJet_mass'))
    jets['FatJet_btagCSVV2'] = pd.DataFrame(events.array('FatJet_btagCSVV2'))
    jets['FatJet_btagDeepB'] = pd.DataFrame(events.array('FatJet_btagDeepB'))
    jets['FatJet_msoftdrop'] = pd.DataFrame(events.array('FatJet_msoftdrop'))
    jets['FatJet_btagDDBvL'] = pd.DataFrame(events.array('FatJet_btagDDBvL'))
    jets['FatJet_deepTagMD_H4qvsQCD'] = pd.DataFrame(events.array('FatJet_deepTagMD_H4qvsQCD'))
    jets['FatJet_n2b1'] = pd.DataFrame(events.array('FatJet_n2b1'))
    jets['SubJet_mass(1)'] = getSubJetData(1,'SubJet_mass', events)
    jets['SubJet_mass(2)'] = getSubJetData(2, 'SubJet_mass', events)
    jets['SubJet_tau1(1)'] = getSubJetData(1, 'SubJet_tau1', events)
    jets['FatJet_n3b1'] = pd.DataFrame(events.array('FatJet_n3b1'))
    jets['FatJet_tau2'] = pd.DataFrame(events.array('FatJet_tau2'))
    jets['SubJet_n2b1(1)'] = getSubJetData(1, 'SubJet_n2b1', events)
    jets['SubJet_pt(1)/FatJet_pt'] = getSubJetData(1, 'SubJet_pt', events)/jets.FatJet_pt
    jets['SubJet_pt(2)/FatJet_pt'] = getSubJetData(2, 'SubJet_pt', events)/jets.FatJet_pt
    jets['SubJet_btagDeepB(2)'] = getSubJetData(2, 'SubJet_btagDeepB', events)
    jets['SubJet_tau1(2)'] = getSubJetData(2, 'SubJet_tau1', events)

    ## JetHT triggers
    ## add them all up, so passing at least one is just see if event has trig > 0
    if JetHT:
        trig1 = events.array('HLT_AK8PFJet500')
        trig2 = events.array('HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4')
        trig3 = events.array('HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71')
        trig4 = events.array('HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5')
        other['HLT_trigger'] = pd.DataFrame(trig1+trig2+trig3+trig4)

    ## muons
    if not JetHT:
        muons['Muon_pt'] = pd.DataFrame(events.array('Muon_pt'))
        muons['Muon_eta'] = pd.DataFrame(np.abs(events.array('Muon_eta')))
        muons['Muon_ip3d'] = pd.DataFrame(events.array('Muon_ip3d'))
        muons['Muon_softId'] = pd.DataFrame(events.array('Muon_softId'))
        muons['Muon_IP'] = pd.DataFrame(events.array('Muon_dxy')/events.array('Muon_dxyErr')).abs()




    ## other vars
    if tag != 'data' and tag != 'JetHT':
        other['LHE_HT'] = pd.DataFrame(events.array('LHE_HT'))
    other['PV_npvs'] = pd.DataFrame(events.array('PV_npvs'))
    other['PV_npvsGood'] = pd.DataFrame(events.array('PV_npvsGood'))
    '''for var in allVars:
        df = pd.DataFrame(events.array(var))
        if 'eta' in var:
            df = df.abs()
        if 'FatJet' in var:
            jets[var] = df
        elif 'Muon' in var:
            muons[var] = df
        else:
            other[var] = df'''

    #muons['Muon_IP'] = (muons['Muon_dxy']/muons['Muon_dxyErr']).abs()


    ## make Event object
    ev = Event(jets, muons, other)

    ## cutting events
    ## jet cuts
    jets.cut(jets['FatJet_pt'] > 170)
    jets.cut(jets['FatJet_eta'] < 2.4)
    jets.cut(jets['FatJet_btagDDBvL'] > 0.8)
    jets.cut(jets['FatJet_btagDeepB'] > 0.4184)
    jets.cut(jets['FatJet_msoftdrop'] > 90)
    jets.cut(jets['FatJet_msoftdrop'] <= 200)
    jets.cut(jets['FatJet_mass'] > 90)
    jets.cut(jets['FatJet_mass'] <= 200)
    ## muon cuts
    if not JetHT:
        muons.cut((muons['Muon_softId'] > 0))
        muons.cut(muons['Muon_eta'] < 2.4)
        muons.cut(muons['Muon_pt'] > 7)
        muons.cut(muons['Muon_IP'] > 2)
        muons.cut(muons['Muon_ip3d'] < 0.5)
    if JetHT:
        other.cut(other['HLT_trigger'] > 0)


    ## sync so all events cut to same events after apply individual cuts
    ev.sync()


    ## rename the columns of LHE_HT, PV_npvs, PV_npvsGood to match the ones that get
    ## passed into getMaxPt
    if tag != 'data' and tag != 'JetHT':
        other.LHE_HT = other.LHE_HT.rename({0:'LHE_HT'}, axis='columns')
    other.PV_npvs = other.PV_npvs.rename({0:'PV_npvs'}, axis='columns')
    other.PV_npvsGood =other.PV_npvsGood.rename({0:'PV_npvsGood'}, axis='columns')


    ## if nothing's left after cut, return empty dataframe
    if (jets.FatJet_pt.empty):
        maxPtData = pd.DataFrame()

    else:
        maxPtJets = getMaxPt(jets, 'FatJet_pt', jetVars)#.reindex(jets.FatJet_pt.index)
        if JetHT:
            maxPtData = maxPtJets
        if not JetHT:

            maxPtMuons = getMaxPt(muons, 'Muon_pt', muonVars + ['Muon_IP'])#.reindex(muons.Muon_pt.index)
            maxPtData = pd.concat([maxPtJets,
                               maxPtMuons], axis=1)

        maxPtData = maxPtData.assign(PV_npvs=other.PV_npvs.to_numpy())
        maxPtData = maxPtData.assign(PV_npvsGood=other.PV_npvsGood.to_numpy())

        ## secondary vertex
        maxPtData['FatJet_nSV'] = getnSVCounts(jets, events)


        if tag != 'data' and tag != 'JetHT':
            maxPtData = maxPtData.assign(LHE_HT=other.LHE_HT.to_numpy())#pd.concat([maxPtData, other.LHE_HT], axis=1)
            #maxPtData = maxPtData.assign(LHE_HT=other.LHE_HT)

            ## LHE_weights
            if tag == 'ggH':
                 maxPtData['LHE_weights'] = 1
            elif tag == 'BGen':
                maxPtData.loc[(maxPtData['LHE_HT']>200) & (maxPtData['LHE_HT']<=300),
                              'LHE_weights'] = BGenWeight[0]
                maxPtData.loc[(maxPtData['LHE_HT']>300) & (maxPtData['LHE_HT']<=500),
                              'LHE_weights'] = BGenWeight[1]
                maxPtData.loc[(maxPtData['LHE_HT']>500) & (maxPtData['LHE_HT']<=700),
                              'LHE_weights'] = BGenWeight[2]
                maxPtData.loc[(maxPtData['LHE_HT']>700) & (maxPtData['LHE_HT']<=1000),
                              'LHE_weights'] = BGenWeight[3]
                maxPtData.loc[(maxPtData['LHE_HT']>1000) & (maxPtData['LHE_HT']<=1500),
                              'LHE_weights'] = BGenWeight[4]
                maxPtData.loc[(maxPtData['LHE_HT']>1500) & (maxPtData['LHE_HT']<=2000),
                              'LHE_weights'] = BGenWeight[5]
                maxPtData.loc[maxPtData['LHE_HT']>2000,
                              'LHE_weights'] = BGenWeight[6]
                wgt = 4.346 - 0.356*np.log(maxPtData.LHE_HT)/np.log(2)
                wgt[wgt<0.1] = 0.1
                maxPtData['QCD_correction'] = wgt

                maxPtData = maxPtData.assign(final_weights=
                                             maxPtData['LHE_weights']*
                                             maxPtData['QCD_correction']
                                             *21.56) # this is Xsec weight
            elif tag == 'bEnr':
                maxPtData.loc[(maxPtData['LHE_HT']>200) & (maxPtData['LHE_HT']<=300),
                              'LHE_weights'] = bEnrWeight[0]
                maxPtData.loc[(maxPtData['LHE_HT']>300) & (maxPtData['LHE_HT']<=500),
                              'LHE_weights'] = bEnrWeight[1]
                maxPtData.loc[(maxPtData['LHE_HT']>500) & (maxPtData['LHE_HT']<=700),
                              'LHE_weights'] = bEnrWeight[2]
                maxPtData.loc[(maxPtData['LHE_HT']>700) & (maxPtData['LHE_HT']<=1000),
                              'LHE_weights'] = bEnrWeight[3]
                maxPtData.loc[(maxPtData['LHE_HT']>1000) & (maxPtData['LHE_HT']<=1500),
                              'LHE_weights'] = bEnrWeight[4]
                maxPtData.loc[(maxPtData['LHE_HT']>1500) & (maxPtData['LHE_HT']<=2000),
                              'LHE_weights'] = bEnrWeight[5]
                maxPtData.loc[maxPtData['LHE_HT']>2000,
                              'LHE_weights'] = bEnrWeight[6]

                wgt = 4.346 - 0.356*np.log(maxPtData.LHE_HT)/np.log(2)
                wgt[wgt<0.1] = 0.1
                maxPtData['QCD_correction'] = wgt

                maxPtData = maxPtData.assign(final_weights=
                                             maxPtData['LHE_weights']*
                                             maxPtData['QCD_correction']*
                                             8.20)

            elif tag == 'WJets':
                maxPtData.loc[(maxPtData['LHE_HT']>=400) & (maxPtData['LHE_HT']<600),
                              'LHE_weights'] = 315600/10071273
                maxPtData.loc[(maxPtData['LHE_HT']>=600) & (maxPtData['LHE_HT']<800),
                              'LHE_weights'] = 68570/15298056
                maxPtData.loc[(maxPtData['LHE_HT']>=800),
                              'LHE_weights'] = 34900/14627242

                maxPtData = maxPtData.assign(final_weights = maxPtData['LHE_weights'])
            elif tag == 'ZJets':
                maxPtData.loc[(maxPtData['LHE_HT']>=400) & (maxPtData['LHE_HT']<600),
                              'LHE_weights'] = 145400/16704355
                maxPtData.loc[(maxPtData['LHE_HT']>=600) & (maxPtData['LHE_HT']<800),
                              'LHE_weights'] = 34000/14642701
                maxPtData.loc[(maxPtData['LHE_HT']>=800),
                              'LHE_weights'] = 18670/10561192

                maxPtData = maxPtData.assign(final_weights = maxPtData['LHE_weights'])

            elif tag == 'TTJets':
                maxPtData = maxPtData.assign(final_weights=831760.0/10244307)

            if not JetHT:
                ## npvs Ratio weights
                for i in range(len(ptkeys)-1):
                    for j in range(len(ipkeys)-1):
                        for k in range(len(npvsGkeys)):
                            maxPtData.loc[(maxPtData.Muon_pt >= ptkeys[i]) &
                                          (maxPtData.Muon_pt < ptkeys[i+1]) &
                                          (maxPtData.Muon_IP >= ipkeys[j]) &
                                          (maxPtData.Muon_IP < ipkeys[j+1]) &
                                          (maxPtData.Muon_eta < 1.5) &
                                          (maxPtData.PV_npvsGood == k+1),
                                          'PU_weights'] = muonR[ptkeys[i]][ipkeys[j]]['L'][k]
                            maxPtData.loc[(maxPtData.Muon_pt >= ptkeys[i]) &
                                          (maxPtData.Muon_pt < ptkeys[i+1]) &
                                          (maxPtData.Muon_IP >= ipkeys[j]) &
                                          (maxPtData.Muon_IP < ipkeys[j +1]) &
                                          (maxPtData.Muon_eta >= 1.5) &
                                          (maxPtData.PV_npvsGood == k+1),
                                          'PU_weights'] = muonR[ptkeys[i]][ipkeys[j]]['H'][k]
                # lumi weights
                for i in range(len(ptkeys)-1):
                    for j in range(len(ipkeys)-1):
                        maxPtData.loc[(maxPtData.Muon_pt >= ptkeys[i]) &
                                      (maxPtData.Muon_pt < ptkeys[i+1]) &
                                      (maxPtData.Muon_IP >= ipkeys[j]) &
                                      (maxPtData.Muon_IP < ipkeys[j+1]) &
                                      (maxPtData.Muon_eta < 1.5),
                                      'lumi_weights'] = muonL[ptkeys[i]][ipkeys[j]]['L']
                        maxPtData.loc[(maxPtData.Muon_pt >= ptkeys[i]) &
                                      (maxPtData.Muon_pt < ptkeys[i+1]) &
                                      (maxPtData.Muon_IP >= ipkeys[j]) &
                                      (maxPtData.Muon_IP  < ipkeys[j +1]) &
                                      (maxPtData.Muon_eta >= 1.5),
                                      'lumi_weights'] = muonL[ptkeys[i]][ipkeys[j]]['H']
        
                maxPtData = maxPtData.assign(final_weights =
                                             maxPtData['lumi_weights']*
                                             maxPtData['PU_weights']*
                                             maxPtData['final_weights'])



    maxPtData = maxPtData.dropna(how='all')
    maxPtData = maxPtData.fillna(0)

    return maxPtData

def main():
    processData(BGenPath, 'BGen')

if __name__ == "__main__":
    main()

