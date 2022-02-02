#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import uproot
import pandas as pd
from analib import PhysObj, Event
import sys
import os
import pickle
import numpy as np
from htoaaRootFilesLoc import TTJetsPaths, WJetsPaths, bEnrPaths, BGenPaths, ZJetsPaths, ParkedDataPaths, JetHTPaths, ggHPaths

BGenWeight = [1, 0.259, 0.0515, 0.01666, 0.00905, 0.003594, 0.001401]
bEnrWeight =[ 1, 0.33, 0.034, 0.034, 0.024, 0.0024, 0.00044]
ZJetsWeight = [ 145400. / 16704355, 34000. / 14642701, 18670. / 10561192]
WJetsWeight = [315600. / 10071273, 68570./ 15298056, 34900. / 14627242]
TJetsWeight = [831760.0 / 10244307]
ParkedDataWeight = [7.1055]

BGenDict = dict(zip(BGenPaths, BGenWeight))
bEnrDict = dict(zip(bEnrPaths, bEnrWeight))
ZJetsDict = dict(zip(ZJetsPaths, ZJetsWeight))
WJetsDict = dict(zip(WJetsPaths, WJetsWeight))
TTJetsDict = dict(zip(TTJetsPaths, TJetsWeight))
ParkedDataDict = dict(zip(ParkedDataPaths, ParkedDataWeight))

## this is to choose which variable sets to use
## 'b' = baselien
## 'h' = baseline + high discriminatory
## 'm' = baseline + high + medium
disc = 'h'

## if you are doing dataVSMC using JetHT, set true
## maybe in the future, can read this from the where you read file names
JetHT = False

## vars to grab maxPt from jet physobj
jetVars = ['FatJet_pt',
           'FatJet_eta',
           'FatJet_mass',
           'FatJet_msoftdrop',
           'FatJet_btagCSVV2',
           'FatJet_btagDeepB',
           'FatJet_msoftdrop',
           'FatJet_btagDDBvL',
           'FatJet_deepTagMD_H4qvsQCD',
           'FatJet_n2b1',
           'SubJet_mass(1)',
           'SubJet_mass(2)',
           'SubJet_tau1(1)',
           'FatJet_n3b1',
           'FatJet_tau2',
           'SubJet_n2b1(1)',
           'SubJet_pt(1)|FatJet_pt',
           'SubJet_pt(2)|FatJet_pt',
           'SubJet_btagDeepB(2)',
           'SubJet_tau1(2)',
           'FatJet_nSV']

if disc == 'h' or disc == 'b':
    for i in ['FatJet_n3b1','FatJet_tau2','SubJet_n2b1(1)','SubJet_pt(1)|FatJet_pt',
              'SubJet_pt(2)|FatJet_pt','SubJet_btagDeepB(2)','SubJet_tau1(2)']:
        jetVars.remove(i)
    if disc == 'b':
        for j in ['FatJet_n2b1','SubJet_mass(1)', 'SubJet_mass(2)', 'SubJet_tau1(1)',
                  'FatJet_nSV']:
            jetVars.remove(j)


## vars to grab maxpt from muon physobj
muonVars = ['Muon_pt',
            'Muon_eta',
            'Muon_ip3d',
            'Muon_softId']

## others
otherVars = [ 'PV_npvs', 'PV_npvsGood', 'LHE_HT']

allVars = list(jetVars + muonVars + otherVars)
allVars.sort()

## define training variables
standardVars = [
    'FatJet_pt',
    'FatJet_eta',
    'FatJet_mass',
    'FatJet_btagCSVV2',
    'FatJet_btagDeepB',
    'FatJet_msoftdrop',
    'FatJet_btagDDBvL',
    'FatJet_deepTagMD_H4qvsQCD'
    ]
highDiscVars = [
    'FatJet_n2b1',
    'SubJet_mass(1)',
    'SubJet_mass(2)',
    'SubJet_tau1(1)',
    'FatJet_nSV'
    ]
mediumDiscVars = [
    'FatJet_n3b1',
    'FatJet_tau2',
    'SubJet_n2b1(1)',
    'SubJet_pt(1)|FatJet_pt',
    'SubJet_pt(2)|FatJet_pt',
    'SubJet_btagDeepB(2)',
    'SubJet_tau1(2)'
    ]

if disc == 'b':
    trainVars = standardVars
elif disc == 'h':
    trainVars = standardVars + highDiscVars
elif disc == 'm':
    trainVars = standardVars + highDiscVars + mediumDiscVars

## lumi and npvs ratio tensors
muonR = pickle.load(open('weights/MuonRtensor.p', 'rb'))
muonL = pickle.load(open('weights/MuonLtensor.p', 'rb'))

ptkeys = list(muonL.keys()) + [999999]
ptkeys.remove('meta')
ipkeys = list(muonL[ptkeys[0]].keys()) + [999999]
npvsGkeys = muonR[6][2]['H']

## to make sure I sent in the right kind of stuff
tagslist = ['bEnr', 'BGen', 'data', 'JetHT', 'WJets', 'TTJets', 'ZJets', 'ggH']

## returns dataframe containing maxpt jet/muon and corresponding variable
## takes 1) PhysObj to be extracted from, 2) column of maxPt, 3) list of vars
## in the PhysObj (will also become column names)
def getMaxPt(physobj, col, varlist):
    varlistcopy = varlist
    if 'FatJet_nSV' in varlist:
        varlistcopy.remove('FatJet_nSV')
    colidx = physobj[col].idxmax(axis=1).to_numpy()
    rowidx = list(range(len(colidx)))
    maxPtData = pd.DataFrame()

    for var in varlistcopy:
        npArr = physobj[var].to_numpy()
        maxPtData[var] = npArr[rowidx, colidx]
    return maxPtData

## returns dataframe of subjet variable
## takes 1) the number associated with subjet, 2) subjet variable name as string
## 3) f.get('Events') from the rootfile
def getSubJetData(subjetnum, subjetvarname, events):
    subjetidx = events.array(f'FatJet_subJetIdx{subjetnum}')
    subjetvar = events.array(subjetvarname)
    padto = subjetvar.counts.max() + 1
    subjetvar = subjetvar.pad(padto).fillna(0)

    return pd.DataFrame(subjetvar[subjetidx])

## returns dataframe of secondary vertex <0.8 counts
## takes 1) physobj of jets, 2) f.get('Events') from root
## takes the idx of the max pt fatjet and index idx. get eta and phi of the corresponding
## fatjet and SV
def getnSVCounts(data, events):
    eventidx = data.FatJet_pt.index
    maxptjetidx = data['FatJet_pt'].idxmax(axis = 1).to_numpy() #same as colidx

    jeteta = events.array('FatJet_eta')[eventidx, maxptjetidx]
    jetphi = events.array('FatJet_phi')[eventidx, maxptjetidx]
    ## sveta and svphi are arrays of values corresponding to the event
    sveta = events.array('SV_eta')[eventidx]
    svphi = events.array('SV_phi')[eventidx]

    dr = np.sqrt(np.power(jeteta - sveta, 2) + np.power(dphi(jetphi, svphi), 2))
    dr = pd.DataFrame(dr)
    nSVcounts = (dr < 0.8).sum(axis=1)
    return nSVcounts

def dphi(jetphi, svphi):
    jetphi[jetphi < 0 ] = jetphi[jetphi<0] + 2*np.pi
    svphi[svphi < 0] = svphi[svphi<0] + 2*np.pi
    return jetphi - svphi

def processData (filePath, tag, BDT):
    if tag == 'JetHT' or tag == 'data':
        MC = False
    else:
        MC = True
    if BDT: # this is to cover my ass
        JethT = False

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

    ## fille the physObjs
    ## still break into disc and train or dataVSMC

    ## baseline
    jets['FatJet_pt'] = pd.DataFrame(events.array('FatJet_pt'))
    jets['FatJet_eta'] = pd.DataFrame(np.abs(events.array('FatJet_eta')))
    jets['FatJet_mass'] = pd.DataFrame(events.array('FatJet_mass'))
    jets['FatJet_btagCSVV2'] = pd.DataFrame(events.array('FatJet_btagCSVV2'))
    jets['FatJet_btagDeepB'] = pd.DataFrame(events.array('FatJet_btagDeepB'))
    jets['FatJet_msoftdrop'] = pd.DataFrame(events.array('FatJet_msoftdrop'))
    jets['FatJet_btagDDBvL'] = pd.DataFrame(events.array('FatJet_btagDDBvL'))
    jets['FatJet_deepTagMD_H4qvsQCD'] = pd.DataFrame(events.array('FatJet_deepTagMD_H4qvsQCD'))

    if 'h' == disc or 'm' == disc:
        jets['FatJet_n2b1'] = pd.DataFrame(events.array('FatJet_n2b1'))
        jets['SubJet_mass(1)'] = getSubJetData(1,'SubJet_mass', events)
        jets['SubJet_mass(2)'] = getSubJetData(2, 'SubJet_mass', events)
        jets['SubJet_tau1(1)'] = getSubJetData(1, 'SubJet_tau1', events)
    if 'm' == disc:
        jets['FatJet_n3b1'] = pd.DataFrame(events.array('FatJet_n3b1'))
        jets['FatJet_tau2'] = pd.DataFrame(events.array('FatJet_tau2'))
        jets['SubJet_n2b1(1)'] = getSubJetData(1, 'SubJet_n2b1', events)
        jets['SubJet_pt(1)|FatJet_pt'] = getSubJetData(1, 'SubJet_pt', events)/jets.FatJet_pt
        jets['SubJet_pt(2)|FatJet_pt'] = getSubJetData(2, 'SubJet_pt', events)/jets.FatJet_pt
        jets['SubJet_btagDeepB(2)'] = getSubJetData(2, 'SubJet_btagDeepB', events)
        jets['SubJet_tau1(2)'] = getSubJetData(2, 'SubJet_tau1', events)



    ## if the thing is not being used to train BDT, assume it's for dataVsMC
    if not BDT:
        ## LHE triggers for JetHT
        if 'JetHT' == tag:
            trig1 = events.array('HLT_AK8PFJet500')
            trig2 = events.array('HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4')
            trig3 = events.array('HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71')
            trig4 = events.array('HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5')
            other['HLT_trigger'] = pd.DataFrame(trig1+trig2+trig3+trig4)
        ## if it's not JetHT, it's probably parked which needs the muons
        else:
            muons['Muon_pt'] = pd.DataFrame(events.array('Muon_pt'))
            muons['Muon_eta'] = pd.DataFrame(np.abs(events.array('Muon_eta')))
            muons['Muon_ip3d'] = pd.DataFrame(events.array('Muon_ip3d'))
            muons['Muon_softId'] = pd.DataFrame(events.array('Muon_softId'))
            muons['Muon_IP'] = pd.DataFrame(events.array('Muon_dxy')/events.array('Muon_dxyErr')).abs()

    if MC:
        other['LHE_HT'] = pd.DataFrame(events.array('LHE_HT'))
    other['PV_npvs'] = pd.DataFrame(events.array('PV_npvs'))
    other['PV_npvsGood'] = pd.DataFrame(events.array('PV_npvsGood'))


    ## make Event object
    ev = Event(jets, muons, other)

    ## cuts
    ## jet cuts
    jets.cut(jets['FatJet_pt'] > 170)
    jets.cut(jets['FatJet_eta'] < 2.4)
    jets.cut(jets['FatJet_btagDDBvL'] > 0.8)
    jets.cut(jets['FatJet_btagDeepB'] > 0.4184)
    jets.cut(jets['FatJet_msoftdrop'] > 90)
    jets.cut(jets['FatJet_msoftdrop'] <= 200)
    jets.cut(jets['FatJet_mass'] > 90)
    jets.cut(jets['FatJet_mass'] <= 200)
    other.cut(other['PV_npvsGood'] >= 1)

    if not BDT:
        if JetHT:
            other.cut(other['HLT_trigger'] > 0)
        else:
            muons.cut((muons['Muon_softId'] > 0))
            muons.cut(muons['Muon_eta'] < 2.4)
            muons.cut(muons['Muon_pt'] > 7)
            muons.cut(muons['Muon_IP'] > 2)
            muons.cut(muons['Muon_ip3d'] < 0.5)


    ## sync so all events to propagate the cuts
    ev.sync()

    ## if nothing's left after cut, return empty dataframe
    if (jets.FatJet_pt.empty):
        return pd.DataFrame()


    ## rename the columns of otherVars PhysObj so they can be passed to
    ## maxptData later
    if MC:
        other.LHE_HT = other.LHE_HT.rename({0:'LHE_HT'}, axis=1)
    other.PV_npvs = other.PV_npvs.rename({0:'PV_npvs'}, axis='columns')
    other.PV_npvsGood =other.PV_npvsGood.rename({0:'PV_npvsGood'}, axis='columns')

    maxPtJets = getMaxPt(jets, 'FatJet_pt', jetVars)
    if (not JetHT) and (not BDT):
        maxPtMuons = getMaxPt(muons, 'Muon_pt', muonVars + ['Muon_IP'])
        maxPtData = pd.concat([maxPtJets, maxPtMuons], axis=1)
    else:
        maxPtData = maxPtJets

    maxPtData = maxPtData.assign(PV_npvs=other.PV_npvs.to_numpy())
    maxPtData = maxPtData.assign(PV_npvsGood=other.PV_npvsGood.to_numpy())
    if MC:
        maxPtData = maxPtData.assign(LHE_HT=other.LHE_HT.to_numpy())

    ## secondary vertex
    if disc!='b':
        maxPtData['FatJet_nSV'] = getnSVCounts(jets, events)

    ## define target to distinguish between signal and background
    ## during training
    if BDT:
        if tag == 'ggH':
            maxPtData['target'] = 1
        else:
            maxPtData['target'] = 0

    if tag == 'data':
        maxPtData['final_weights'] = 7.1055


    if tag == 'ggH':
         maxPtData['LHE_weights'] = 1
         wgt = 3.9 - 0.4*np.log(maxPtData.FatJet_pt)/np.log(2)
         wgt[wgt<0.1] = 0.1
         maxPtData['ggH_weights'] = wgt
         maxPtData['final_weights'] = (maxPtData['LHE_weights'] *
                                       maxPtData['ggH_weights'])
    elif tag == 'BGen':
        maxPtData['LHE_weights'] = BGenDict[filePath]
        wgt = 4.346 - 0.356*np.log(maxPtData.LHE_HT)/np.log(2)
        wgt[wgt<0.1] = 0.1
        maxPtData['QCD_correction'] = wgt
        Xsec_wgt = 21.56

        maxPtData = maxPtData.assign(final_weights =
                                     maxPtData['LHE_weights']*
                                     maxPtData['QCD_correction']*
                                     Xsec_wgt)

    elif tag == 'bEnr':
        maxPtData['LHE_weights'] = bEnrDict[filePath]
        wgt = 4.346 - 0.356*np.log(maxPtData.LHE_HT)/np.log(2)
        wgt[wgt<0.1] = 0.1
        maxPtData['QCD_correction'] = wgt
        Xsec_wgt = 8.2

        maxPtData = maxPtData.assign(final_weights=
                                     maxPtData['LHE_weights']*
                                     maxPtData['QCD_correction']*
                                     Xsec_wgt)

    elif tag == 'WJets':
        maxPtData['LHE_weights'] = WJetsDict[filePath]
        maxPtData = maxPtData.assign(final_weights = maxPtData['LHE_weights'])

    elif tag == 'ZJets':
        maxPtData['LHE_weights'] = ZJetsDict[filePath]
        maxPtData = maxPtData.assign(final_weights = maxPtData['LHE_weights'])

    elif tag == 'TTJets':
        maxPtData['LHE_weights'] = TTJetsDict[filePath]
        maxPtData = maxPtData.assign(final_weights=maxPtData['LHE_weights'])

    if (not BDT) and (not JetHT) and (tag!='ggH') and (tag!='data'):
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


