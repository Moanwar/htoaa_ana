#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 16:01:35 2020

@author: si_sutantawibul1
"""

import uproot
import pandas as pd
from analib import PhysObj, Event
import sys
import os
import pickle


dataPath = 'data/2018D_Parked_promptD-v1_200218_214714_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root'
# ggHPath = 'MC/nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId_999k.root'
# BGenPath = 'MC/QCD_BGen_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root'
# bEnrPath = 'MC/QCD_bEnriched_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root'

ggHPath = 'GGH_HPT.root'
setnames = ['QCD_HT200to300', 'QCD_HT300to500', 'QCD_HT500to700',
            'QCD_HT700to1000','QCD_HT1000to1500', 'QCD_HT1500to2000',
            'QCD_HT2000toInf']
BGenPaths = ['QCD_BGenFilter/' + i + '.root' for i in setnames]
bEnrPaths = ['QCD_bEnriched/' + i + '.root' for i in setnames]

BGenWeight = [1, 0.259, 0.0515, 0.01666, 0.00905, 0.003594, 0.001401]
bEnrWeight =[ 1, 0.33, 0.034, 0.034, 0.024, 0.0024, 0.00044]

trainVars = [
    'FatJet_pt', 
    'FatJet_eta', 
    'FatJet_mass', 
    'FatJet_btagCSVV2', 
    'FatJet_btagDeepB', 
    'FatJet_msoftdrop', 
    'FatJet_btagDDBvL',
    'FatJet_deepTagMD_H4qvsQCD'
]
trainVars.sort()

otherVars = [
    'Muon_softId',
    'Muon_eta',
    'Muon_pt',
    'Muon_dxy',
    'Muon_dxyErr',
    'Muon_ip3d',
    'FatJet_pt',
    'FatJet_eta',
    'FatJet_btagDDBvL',
    'FatJet_btagDeepB',
    'FatJet_msoftdrop',
    'FatJet_mass',
    'PV_npvsGood',
    'PV_npvs',
    'LHE_HT'
]
otherVars.sort()

allVars = list(set(trainVars + otherVars))
allVars.sort()

muonR = pickle.load(open('muontensor/MuonRtensor.p', 'rb'))
muonL = pickle.load(open('muontensor/MuonLtensor.p', 'rb'))

## also added on a 999999 treated as infinity
ptkeys = list(muonL.keys()) + [999999]
ptkeys.remove('meta')
ipkeys = list(muonL[ptkeys[0]].keys()) + [999999]
npvsGkeys = muonR[6][2]['H']

tags = ['ggH', 'BGen', 'bEnr', 'data']

def processData(filePath, tag):

    if tag not in tags:
        print("please enter valid tag: 'ggH', 'BGen', 'bEnr', or 'data'")
        sys.exit()

    ## open file, get events
    fileName, fileExtension = os.path.splitext(filePath)

    print(fileName)

    ## check to make sure it is root file
    if fileExtension != '.root':
        print('this program only supports .root  files')
        sys.exit()

    f = uproot.open(fileName + '.root')
    events = f.get('Events')
    ## make PhysObj of the event
    data = PhysObj('data_' + fileName)

    if tag == 'data':
        allVars.remove('LHE_HT')
    for var in allVars: 
        data[var] = pd.DataFrame(events.array(var))
        ## makes eta positive only
        if 'eta' in var: 
            data[var] = data[var].abs()
    data['Muon_IP'] = (data['Muon_dxy']/data['Muon_dxyErr']).abs()

    ## make event object
    ev = Event(data)

    ## apply cuts
    # data.cut(data[cutVar] > cutDict[cutVar])
    data.cut((data['Muon_softId'] == True))
    data.cut(data['Muon_eta'] < 2.4)
    data.cut(data['Muon_pt'] > 7)
    data.cut(data['Muon_IP'] > 2)
    data.cut(data['Muon_ip3d'] < 0.5)
    data.cut(data['FatJet_pt'] > 240)
    data.cut(data['FatJet_eta'] < 2.4)
    data.cut(data['FatJet_btagDDBvL'] > 0.8)
    data.cut(data['FatJet_btagDeepB'] > 0.4184)
    data.cut(data['FatJet_msoftdrop'] > 90)
    data.cut(data['FatJet_msoftdrop'] <= 200)
    data.cut(data['FatJet_mass'] > 90)
    data.cut(data['FatJet_mass'] <= 200)
    data.cut(data['PV_npvsGood'] > 0)

    ev.sync()


    ## return none if dataframe is empty
    if data['FatJet_pt'].empty:
        return


    ## keep only max Pt jet of event
    colidx = data['FatJet_pt'].idxmax(axis = 1).to_numpy()
    rowidx = list(range(len(colidx)))


    maxPtData = pd.DataFrame()

    for var in allVars + ['Muon_IP']:
        npArr = data[var].to_numpy()
        maxPtData[var] = npArr[rowidx, colidx]


    ## define target to distinguish between signal and background
    ## during training
    if tag == 'data':
        maxPtData['target'] = None
    elif tag == 'ggH':
        maxPtData['target'] = 1
    else:
        maxPtData['target'] = 0




    ## LHE weights
    #if tag == 'data':
    #    maxPtData['LHE_weights'] = 1
    if tag == 'ggH':
        maxPtData['LHE_weights'] = 43920/999000
    elif tag == 'BGen':
        maxPtData.loc[(maxPtData['LHE_HT']>200) & (maxPtData['LHE_HT']<300),
                      'LHE_weights'] = BGenWeight[0]
        maxPtData.loc[(maxPtData['LHE_HT']>300) & (maxPtData['LHE_HT']<500),
                      'LHE_weights'] = BGenWeight[1]
        maxPtData.loc[(maxPtData['LHE_HT']>500) & (maxPtData['LHE_HT']<700),
                      'LHE_weights'] = BGenWeight[2]
        maxPtData.loc[(maxPtData['LHE_HT']>700) & (maxPtData['LHE_HT']<1000),
                      'LHE_weights'] = BGenWeight[3]
        maxPtData.loc[(maxPtData['LHE_HT']>1000) & (maxPtData['LHE_HT']<1500),
                      'LHE_weights'] = BGenWeight[4]
        maxPtData.loc[(maxPtData['LHE_HT']>1500) & (maxPtData['LHE_HT']<2000),
                      'LHE_weights'] = BGenWeight[5]
        maxPtData.loc[maxPtData['LHE_HT']>2000,
                      'LHE_weights'] = BGenWeight[6]

    elif tag == 'bEnr':
        maxPtData.loc[(maxPtData['LHE_HT']>200) & (maxPtData['LHE_HT']<300),
                      'LHE_weights'] = bEnrWeight[0]
        maxPtData.loc[(maxPtData['LHE_HT']>300) & (maxPtData['LHE_HT']<500),
                      'LHE_weights'] = bEnrWeight[1]
        maxPtData.loc[(maxPtData['LHE_HT']>500) & (maxPtData['LHE_HT']<700),
                      'LHE_weights'] = bEnrWeight[2]
        maxPtData.loc[(maxPtData['LHE_HT']>700) & (maxPtData['LHE_HT']<1000),
                      'LHE_weights'] = bEnrWeight[3]
        maxPtData.loc[(maxPtData['LHE_HT']>1000) & (maxPtData['LHE_HT']<1500),
                      'LHE_weights'] = bEnrWeight[4]
        maxPtData.loc[(maxPtData['LHE_HT']>1500) & (maxPtData['LHE_HT']<2000),
                      'LHE_weights'] = bEnrWeight[5]
        maxPtData.loc[maxPtData['LHE_HT']>2000,
                      'LHE_weights'] = bEnrWeight[6]

    ## npvs Ratio weights
    if tag != 'data':
        for i in range(len(ptkeys)-1):
            for j in range(len(ipkeys)-1):
                for k in range(len(npvsGkeys)):
                    #print('{} {} {}'.format(ptkeys[i],ipkeys[j],npvsGkeys[k]))
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


        ## lumi weights
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
                             (maxPtData.Muon_IP < ipkeys[j +1]) &
                             (maxPtData.Muon_eta >= 1.5),
                             'lumi_weights'] = muonL[ptkeys[i]][ipkeys[j]]['H']

        ## ultimate weight
        maxPtData = maxPtData.assign(final_weights =
                                     maxPtData['lumi_weights']*
                                     maxPtData['PU_weights']*
                                     maxPtData['LHE_weights'])


    maxPtData = maxPtData.dropna(axis = 1, how = 'all')
    maxPtData = maxPtData.dropna(how = 'all')
    maxPtData = maxPtData.fillna(0)



    return maxPtData







## cuts
## On top of the regular selection, you should require the existence of at
## least one muon passing softId, with |eta| < 2.4, ip3d < 0.5, pT > 6, and
## (dxy/dxyErr) > 2, in both data and MC.

            # muons.cut(muons.softId > 0.9)
            # muons.cut(abs(muons.eta) < 2.4)
            # muons.cut(muons.pt > 7)
            # muons.cut(muons.ip > 2)
            #     jets.cut(jets.pt > 170)#240)#170)
            #     jets.cut(abs(jets.eta)<2.4)
            #     jets.cut(jets.DDBvL > 0.8)#0.8)#0.6)
            #     jets.cut(jets.DeepB > 0.4184)
            #     jets.cut(jets.msoft > 90)#90)#0.25)
            #     #
            #     jets.cut(jets.mass > 90)
            #     jets.cut(jets.msoft < 200)
            #     jets.cut(jets.npvsG >= 1)




