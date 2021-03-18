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
import numpy as np
import math
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

JetVars = [
    'FatJet_pt',
    'FatJet_eta',
    'FatJet_mass',
    'FatJet_btagCSVV2',
    'FatJet_btagDeepB',
    'FatJet_msoftdrop',
    'FatJet_btagDDBvL',
    'FatJet_deepTagMD_H4qvsQCD',
    ## high discriminatroy
    'FatJet_n2b1',
    'SubJet_mass1',
    'SubJet_mass2',
    'SubJet_tau1',
    ## medium
    'FatJet_n3b1',
    'FatJet_tau2',
    'SubJet_n2b1',
    'SubJet_pt',
    'SubJet_btagDeepB',
    'SubJet_tau1',
]
JetVars.sort()

otherVars = [
    # 'Muon_softId',
    # 'Muon_eta',
    # 'Muon_pt',
    # 'Muon_dxy',
    # 'Muon_dxyErr',
    # 'Muon_ip3d',
    'FatJet_pt',
    'FatJet_eta',
    'FatJet_btagDDBvL',
    'FatJet_btagDeepB',
    'FatJet_msoftdrop',
    'FatJet_mass',
    # 'PV_npvsGood',
    # 'PV_npvs',
    'LHE_HT'
]
otherVars.sort()

allVars = list(set(JetVars + otherVars))
allVars.sort()

## for when plot
## disc is the discriminatory variables
##    None = standard vars (default for now)
##    h = standard + high
##    m = standard + high + medium
disc = 'h'
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


if disc == None:
    trainVars = standardVars
elif disc == 'h':
    trainVars = standardVars + highDiscVars
elif disc == 'm':
    trainVars = standardVars + highDiscVars + mediumDiscVars

#muonR = pickle.load(open('muontensor/MuonRtensor.p', 'rb'))
#muonL = pickle.load(open('muontensor/MuonLtensor.p', 'rb'))

## also added on a 999999 treated as infinity
#ptkeys = list(muonL.keys()) + [999999]
#ptkeys.remove('meta')
#ipkeys = list(muonL[ptkeys[0]].keys()) + [999999]
#npvsGkeys = muonR[6][2]['H']

tags = ['ggH', 'BGen', 'bEnr', 'data', 'TTJets', 'WJets', 'ZJets']

def getSubJetData(subjetnum, subjetvarname, events):
    '''idxa1 = events.array('FatJet_subJetIdx1')
    idxa2 = events.array('FatJet_subJetIdx2')
    idxa1f = pd.DataFrame(idxa1)
    idxa2f = pd.DataFrame(idxa2)
    submass = events.array('SubJet_mass')
    subtau = events.array('SubJet_tau1')
    data['SubJet_mass(1)'] = pd.DataFrame(submass[idxa1[idxa1!=-1]]).add(idxa1f[idxa1f==-1]*0,fill_value=0)
    data['SubJet_mass2'] = pd.DataFrame(submass[idxa2[idxa2!=-1]]).add(idxa2f[idxa2f==-1]*0,fill_value=0)
    data['SubJet_tau1'] = pd.DataFrame(subtau[ idxa1[idxa1!=-1]]).add(idxa1f[idxa1f==-1]*0,fill_value=0)
    '''
    subjetidx = events.array(f'FatJet_subJetIdx{subjetnum}')
    subjetidxDF = pd.DataFrame(subjetidx)
    subjetvar = events.array(subjetvarname)

    return pd.DataFrame(subjetvar[subjetidx[subjetidx!=-1]]).add(subjetidxDF[subjetidxDF==-1]*0,fill_value=0)


## get secondary vertex info
## takes the idx of the max pt fatjet and index idx. get eta and phi of the corresponding
## fatjet and SV. Takes the original PhysObj or jet PysObj and the opened root file
def getnSVCounts(data, events):
    eventidx = data.FatJet_pt.index
    maxptjetidx = data['FatJet_pt'].idxmax(axis = 1).to_numpy() #same as colidx

    jeteta = events.array('FatJet_eta')[eventidx, maxptjetidx]
    jetphi = events.array('FatJet_phi')[eventidx, maxptjetidx]
    ## sveta and svphi are arrays of values corresponding to the event
    sveta = events.array('SV_eta')[eventidx]
    svphi = events.array('SV_phi')[eventidx]

    dphi = np.arccos(np.cos(jetphi+3.1)) + np.arccos(np.cos(svphi+3.1))
    dr = np.sqrt(np.power(jeteta - sveta, 2) + np.power(dphi, 2))
    dr = pd.DataFrame(dr)
    nSVcounts = (dr < 0.8).sum(axis=1)
    return nSVcounts



## makes the DF for putting into the BDT
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
    data = PhysObj('Jet' + fileName)

    ############## filling ##############
    ## filling the dataframe with events
    '''if tag == 'data':
        allVars.remove('LHE_HT')
    for var in allVars:
        ## in try because some in trainVars don't exist in the root file
        ## but i need to have it when training
        try:
            data[var] = pd.DataFrame(events.array(var))
            ## makes eta positive only
            if 'eta' in var:
                data[var] = data[var].abs()
        except:
            continue'''

    ## standard
    data['FatJet_pt'] = pd.DataFrame(events.array('FatJet_pt'))
    data['FatJet_eta'] = pd.DataFrame(np.abs(events.array('FatJet_eta')))
    data['FatJet_mass'] = pd.DataFrame(events.array('FatJet_mass'))
    data['FatJet_btagCSVV2'] = pd.DataFrame(events.array('FatJet_btagCSVV2'))
    data['FatJet_btagDeepB'] = pd.DataFrame(events.array('FatJet_btagDeepB'))
    data['FatJet_msoftdrop'] = pd.DataFrame(events.array('FatJet_msoftdrop'))
    data['FatJet_btagDDBvL'] = pd.DataFrame(events.array('FatJet_btagDDBvL'))
    data['FatJet_deepTagMD_H4qvsQCD'] = pd.DataFrame(events.array('FatJet_deepTagMD_H4qvsQCD'))

    if tag != 'data':
        data['LHE_HT'] = pd.DataFrame(events.array('LHE_HT'))

    if disc == 'h' or disc == 'm':
        ## high-disc
        data['FatJet_n2b1'] = pd.DataFrame(events.array('FatJet_n2b1'))
        ## high-disc subjets
        data['SubJet_mass(1)'] = getSubJetData(1,'SubJet_mass', events)
        data['SubJet_mass(2)'] = getSubJetData(2, 'SubJet_mass', events)
        data['SubJet_tau1(1)'] = getSubJetData(1, 'SubJet_tau1', events)

    if disc == 'm':
        ## medium-disc
        data['FatJet_n3b1'] = pd.DataFrame(events.array('FatJet_n3b1'))
        data['FatJet_tau2'] = pd.DataFrame(events.array('FatJet_tau2'))
        ## midium subjets
        data['SubJet_n2b1(1)'] = getSubJetData(1, 'SubJet_n2b1', events)
        data['SubJet_pt(1)|FatJet_pt'] = getSubJetData(1, 'SubJet_pt', events)/data.FatJet_pt
        data['SubJet_pt(2)|FatJet_pt'] = getSubJetData(2, 'SubJet_pt', events)/data.FatJet_pt
        data['SubJet_btagDeepB(2)'] = getSubJetData(2, 'SubJet_btagDeepB', events)
        data['SubJet_tau1(2)'] = getSubJetData(2, 'SubJet_tau1', events)


    ## make event object
    ev = Event(data)

    ## apply cuts
    # data.cut(data[cutVar] > cutDict[cutVar])
    # data.cut((data['Muon_softId'] == True))
    # data.cut(data['Muon_eta'] < 2.4)
    # data.cut(data['Muon_pt'] > 7)
    # data.cut(data['Muon_IP'] > 2)
    # data.cut(data['Muon_ip3d'] < 0.5)
    data.cut(data['FatJet_pt'] > 170) # 240)
    data.cut(data['FatJet_eta'] < 2.4)
    data.cut(data['FatJet_btagDDBvL'] > 0.8)
    data.cut(data['FatJet_btagDeepB'] > 0.4184)
    data.cut(data['FatJet_msoftdrop'] > 90)
    data.cut(data['FatJet_msoftdrop'] <= 200)
    data.cut(data['FatJet_mass'] > 90)
    data.cut(data['FatJet_mass'] <= 200)
    # data.cut(data['PV_npvsGood'] > 0)

    ev.sync()


    ## return none if dataframe is empty
    if data['FatJet_pt'].empty:
        return


    ## keep only max Pt jet of event
    colidx = data['FatJet_pt'].idxmax(axis = 1).to_numpy()
    rowidx = list(range(len(colidx)))


    maxPtData = pd.DataFrame()
    toiter = (trainVars + ['LHE_HT'])


    if disc != None:
        toiter.remove('FatJet_nSV')

    for var in toiter:
        #print(var)
        npArr = data[var].to_numpy()
        maxPtData[var] = npArr[rowidx, colidx]



    ############# Secondary Vertex stuff ################
    ## get phi and eta for the situation
    if disc != None:
        # eventidx = data.FatJet_pt.index
        # maxptjetidx = colidx

        # jeteta = events.array('FatJet_eta')[eventidx, maxptjetidx]
        # jetphi = events.array('FatJet_phi')[eventidx, maxptjetidx]
        # sveta = events.array('SV_eta')[eventidx]
        # svphi = events.array('SV_phi')[eventidx]

        # dr = np.sqrt(np.power(jeteta - sveta, 2) + np.power(jetphi - svphi, 2))
        # dr = pd.DataFrame(dr)
        # nSVcounts = (dr < 0.8).sum(axis=1)
        maxPtData['FatJet_nSV'] = getnSVCounts(data, events)

    #######################################################

    ## define target to distinguish between signal and background
    ## during training
    if tag == 'data':
        maxPtData['target'] = None
    elif tag == 'ggH':
        maxPtData['target'] = 1
    else:
        maxPtData['target'] = 0




    ## assign (where applicable) LHE_weights, QCD_correction, Xsec weights
    #if tag == 'data':
    #    maxPtData['LHE_weights'] = 1
    if tag == 'ggH':
        maxPtData['LHE_weights'] = 1
        wgt = 3.9 - 0.4*np.log(maxPtData.FatJet_pt)/np.log(2)
        wgt[wgt<0.1] = 0.1
        maxPtData['ggH_weights'] = wgt
        maxPtData['final_weights'] = maxPtData['LHE_weights'] * maxPtData['ggH_weights']
    elif tag == 'BGen':
        maxPtData['LHE_weights'] = BGenDict[filePath]
        '''maxPtData.loc[(maxPtData['LHE_HT']>200) & (maxPtData['LHE_HT']<300),
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
                      'LHE_weights'] = BGenWeight[6]'''

        wgt = 4.346 - 0.356*np.log(maxPtData.LHE_HT)/np.log(2)
        wgt[wgt<0.1] = 0.1
        maxPtData['QCD_correction'] = wgt



        maxPtData = maxPtData.assign(final_weights=
                                     maxPtData['LHE_weights']*
                                     maxPtData['QCD_correction'])
        maxPtData['final_weights'] = maxPtData['final_weights']*(21.56/maxPtData['final_weights'].sum())


    elif tag == 'bEnr':
        maxPtData['LHE_weights'] = bEnrDict[filePath]
        '''maxPtData.loc[(maxPtData['LHE_HT']>200) & (maxPtData['LHE_HT']<300),
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
                      'LHE_weights'] = bEnrWeight[6]'''

        wgt = 4.346 - 0.356*np.log(maxPtData.LHE_HT)/np.log(2)
        wgt[wgt<0.1] = 0.1
        maxPtData['QCD_correction'] = wgt

        maxPtData = maxPtData.assign(final_weights=
                                     maxPtData['LHE_weights']*
                                     maxPtData['QCD_correction'])
        maxPtData['final_weights'] = maxPtData['final_weights']*(8.20/maxPtData['final_weights'].sum())

    elif tag == 'TTJets':
        maxPtData['LHE_weights'] = TTJetsDict[filePath]

    elif tag == 'ZJets':
        maxPtData['LHE_weights'] = ZJetsDict[filePath]
        # maxPtData.loc[(maxPtData['LHE_HT']>=400) & (maxPtData['LHE_HT']<600),
        #               'LHE_weights'] = 145400/16704355
        # maxPtData.loc[(maxPtData['LHE_HT']>=600) & (maxPtData['LHE_HT']<800),
        #               'LHE_weights'] = 34000/14642701
        # maxPtData.loc[(maxPtData['LHE_HT']>=800),
        #               'LHE_weights'] = 18670/10561192

        maxPtData = maxPtData.assign(final_weights = maxPtData['LHE_weights'])

    elif tag == 'WJets':
        maxPtData['LHE_weights'] = WJetsDict[filePath]
        # maxPtData.loc[(maxPtData['LHE_HT']>=400) & (maxPtData['LHE_HT']<600),
        #               'LHE_weights'] = 315600/10071273
        # maxPtData.loc[(maxPtData['LHE_HT']>=600) & (maxPtData['LHE_HT']<800),
        #               'LHE_weights'] = 68570/15298056
        # maxPtData.loc[(maxPtData['LHE_HT']>=800),
        #               'LHE_weights'] = 34900/14627242

        maxPtData = maxPtData.assign(final_weights = maxPtData['LHE_weights'])


    # tag the thing
    maxPtData['tag'] = tag

    maxPtData = maxPtData.dropna(axis = 1, how = 'all')
    maxPtData = maxPtData.dropna(how = 'all')
    maxPtData = maxPtData.fillna(0)


    return maxPtData





