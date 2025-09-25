#! /usr/bin/env python
## Counts events passing leptonic categories

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

MAX_EVT = -1     ## Maximum number of events to process per MC sample
PRT_EVT = 10000  ## Print every Nth event while processing
VERBOSE = False
DEBUG   = [] ## [luminosityBlock, event] to debug
YEAR = '2018'  ## 2016APV, 2016, 2017, 2018

IN_DIR = '/eos/cms/store/user/abrinke1/NanoPostv2/%s/' % YEAR
#SAMPS = ['GluGluH_M-55','VBFH_M-15','WH_M-15','ZH_M-30','TTH_M-55']
SAMPS = ['VBFH_M-15']

# SAMPS = ['SingleMuon']
# IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/%s/data/PNet_v2_2024_11_22/SingleMuon/r1_Run%sC/' % (YEAR,YEAR)

CUTS = ['all','nPV','noise',
        'fatH_presel','mu_veto','ele_veto',
        'bjet_veto_dR','bjet_veto','MET_veto',
        'fatTop_sel','fatWZ_sel',
        'fatTop_veto','fatWZ_veto',
        'VBFjj_sel','VBFjj_veto','VBFjjLo','VBFjjHi',
        'VBFjjLoPtLo','VBFjjLoPtHi','VBFjjHiPtLo','VBFjjHiPtHi',
        'gg0lHi','gg0lLo','VjjHi','VjjLo','tt0l0b','tt0l1b','ZvvHi','ZvvLo',
        'X4bSB','X4bSR','4GenB']
CATS = {}
CATS['gg0l'] = [['all'],['nPV'],['noise'],['fatH_presel'],
                ['trigJetHT','trigBTag'],['mu_veto'],['ele_veto'],
                ['VBFjj_veto'],['fatWZ_veto'],['MET_veto'],['fatTop_veto'],
                ['bjet_veto_dR'],['gg0lHi','gg0lLo'],
                ['X4bSB'],['X4bSR'],['4GenB']]
CATS['VBFjj'] = [['all'],['nPV'],['noise'],['fatH_presel'],
                 ['trigJetHT','trigBTag'],['mu_veto'],['ele_veto'],
                 ['fatTop_veto'],['fatWZ_veto'],['MET_veto'],['bjet_veto'],
                 ['VBFjj_sel'],['VBFjjLo','VBFjjHi'],
                 ['VBFjjLoPtLo','VBFjjLoPtHi','VBFjjHiPtLo','VBFjjHiPtHi'],
                 ['X4bSB'],['X4bSR'],['4GenB']]
CATS['Vjj'] = [['all'],['nPV'],['noise'],['fatH_presel'],['fatWZ_sel'],
               ['trigJetHT','trigBTag'],['mu_veto'],['ele_veto'],
               ['MET_veto'],['fatTop_veto'],
               ['bjet_veto_dR'],['VjjHi','VjjLo'],
               ['X4bSB'],['X4bSR'],['4GenB']]
CATS['tt0l'] = [['all'],['nPV'],['noise'],['fatH_presel'],['fatTop_sel'],
               ['trigJetHT','trigBTag'],['mu_veto'],['ele_veto'],
               ['MET_veto'],['tt0l0b','tt0l1b'],
               ['X4bSB'],['X4bSR'],['4GenB']]
CATS['Zvv'] = [['all'],['nPV'],['noise'],['fatH_presel'],
               ['trigMET'],['MET_pt'],['MET_dPhi'],['mu_veto'],['ele_veto'],
               ['bjet_veto_dR'],['ZvvHi','ZvvLo'],
               ['X4bSB'],['X4bSR'],['4GenB']]


## Function to iterate counters based on per-event flags
def fill_counts(counter, selected):
    for cat in CATS.keys():
        for iCuts in range(len(CATS[cat])):
            if len(CATS[cat][iCuts]) == 1:
                cutstr = CATS[cat][iCuts][0]
                key = '%d_%s' % (iCuts, cutstr)
                if selected[cutstr]:
                    counter[cat][key] = counter[cat][key] + 1
                else: ## Subsequent cuts depend on earlier cuts
                    break
            else:
                atleastone = False
                for jCut in range(len(CATS[cat][iCuts])):
                    cutstr = CATS[cat][iCuts][jCut]
                    key = '%d_%s' % (iCuts, cutstr)
                    if selected[cutstr]:
                        counter[cat][key] = counter[cat][key] + 1
                        atleastone = True
                if atleastone:
                    counter[cat]['%d_any' % iCuts] = counter[cat]['%d_any' % iCuts] + 1
                else: ## Subsequent cuts depend on earlier cuts
                    break
        ## End loop: for iCuts in range(len(CATS[cat]))
    ## End loop: for cat in CATS.keys()
    return counter
## End function: fill_counts(counter, selected)    


########################################################
##                  MAIN FUNCTION
## Count passing events for each sample in each category
########################################################
count = {}
for samp in SAMPS:

    print('\n\n*** Beginning to look at sample %s ***\n' % samp)

    count[samp] = {}
    for cat in CATS.keys():
        count[samp][cat] = {}
        for iCuts in range(len(CATS[cat])):
            if len(CATS[cat][iCuts]) == 1:
                count[samp][cat]['%d_%s' % (iCuts, CATS[cat][iCuts][0])] = 0
            else:
                count[samp][cat]['%d_any' % iCuts] = 0
                for jCut in range(len(CATS[cat][iCuts])):
                    count[samp][cat]['%d_%s' % (iCuts, CATS[cat][iCuts][jCut])] = 0
        ## End loop: for iCuts in range(len(CATS[cat]))
    ## End loop: for cat in CATS.keys()

    in_file_names = []
    # for fn in os.listdir(IN_DIR):
    #     if samp == 'SingleMuon':
    #         if fn.endswith('7.root'):
    #             in_file_names.append(IN_DIR+fn)
    in_file_names = [IN_DIR+samp+'_Skim.root']

    chains = {}
    chains['Events'] = 0

    ## Combine input files into a single "chain"
    for i in range(len(in_file_names)):
        print('Adding file %s' % in_file_names[i])
        for key in chains.keys():
            if i == 0: chains[key] = R.TChain(key)
            chains[key].Add( in_file_names[i] )

    ## Loop through events, select, and count
    ch = chains['Events']  ## Shortcut expression
    nEntries = ch.GetEntries()
    print('\nEntering loop over %d events\n' % (nEntries))
    
    ##################
    ## MAIN EVENT LOOP
    for iEvt in range(nEntries):

        if iEvt >= MAX_EVT and MAX_EVT > 0: break
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))

        ch.GetEntry(iEvt)

        if len(DEBUG) > 0:
            if ch.luminosityBlock != DEBUG[0] or ch.event != DEBUG[1]:
                continue
            print('Starting LS = %d, event = %d' % (ch.luminosityBlock, ch.event))

        ## Initialize per-event bools for each cut
        sel = {}
        for cut in CUTS:
           sel[cut] = False
        sel['all'] = True

        ## nPV (TODO: not currently documented? not implemented in Lep?)
        sel['nPV'] = (ch.PV_npvsGood > 0)
        if not sel['nPV']:
            count[samp] = fill_counts(count[samp], sel)
            continue
        ## Table 5 in Ch. 4 of AN2023_047_v4
        sel['noise'] = (ch.Flag_goodVertices and ch.Flag_globalSuperTightHalo2016Filter and
                        ch.Flag_HBHENoiseFilter and ch.Flag_HBHENoiseIsoFilter and
                        ch.Flag_eeBadScFilter and ch.Flag_BadPFMuonFilter and
                        ch.Flag_BadPFMuonDzFilter and ch.Flag_EcalDeadCellTriggerPrimitiveFilter and
                        ('2016' in YEAR or ch.Flag_ecalBadCalibFilter))
        if not sel['noise']:
            count[samp] = fill_counts(count[samp], sel)
            continue


        ## Higgs candidate AK8 jet selection from Table 1 in Ch. 4 of AN2023_047_v4
        xFatH,xFatH_X4b = -99,-99
        xFatTops,xFatWZs = [],[]
        xFatH_X4b,xFatTop_tag,xFatWZ_tag = -99,-99,-99
        for iFat in range(ch.nFatJet):
            if   ch.FatJet_pt_nom[iFat]  <= 250: continue
            if  abs(ch.FatJet_eta[iFat]) >= 2.4: continue
            if     ch.FatJet_jetId[iFat] != 6:   continue
            if ch.FatJet_msoftdrop_nom[iFat] <= 20: continue
            if (YEAR == '2018' and ch.FatJet_particleNet_TvsQCD[iFat] > 0.970) or \
               (YEAR == '2017' and ch.FatJet_particleNet_TvsQCD[iFat] > 0.970) or \
               ('2016' in YEAR and ch.FatJet_particleNet_TvsQCD[iFat] > 0.957 + 0.001*('ost' in YEAR)):
                if ch.FatJet_pt_nom[iFat] > 300:
                    if ch.FatJet_particleNet_TvsQCD[iFat] > xFatTop_tag:
                        xFatTop_tag = ch.FatJet_particleNet_TvsQCD[iFat]
                        xFatTops.insert(0, iFat)
                    else:
                        xFatTops.append(iFat)
            if (YEAR == '2018' and ch.FatJet_particleNet_WZvsQCD[iFat] > 0.9873) or \
               (YEAR == '2017' and ch.FatJet_particleNet_WZvsQCD[iFat] > 0.9858) or \
               ('2016' in YEAR and ch.FatJet_particleNet_WZvsQCD[iFat] > 0.9843):
                if ch.FatJet_particleNet_WZvsQCD[iFat] > xFatWZ_tag:
                    xFatWZ_tag = ch.FatJet_particleNet_WZvsQCD[iFat]
                    xFatWZs.insert(0, iFat)
                else:
                    xFatWZs.append(iFat)
            if ch.FatJet_particleNetMD_XbbvsQCD[iFat] <= 0.75: continue
            iFatH_X4b = 0.5*(ch.FatJet_PNet_X4b_v2a_Haa4b_score[iFat] + \
                             ch.FatJet_PNet_X4b_v2b_Haa4b_score[iFat])
            if iFatH_X4b > xFatH_X4b:
                xFatH     = iFat
                xFatH_X4b = iFatH_X4b
        ## End loop: for iFat in range(ch.nFatJet)
        if xFatH < 0:
            count[samp] = fill_counts(count[samp], sel)
            continue
        sel['fatH_presel'] = True

        ## Set 4-vectors for AK8 jets
        vFatH,vFatWZ,vFatTop = R.TLorentzVector(),None,None
        vFatH.SetPtEtaPhiM(ch.FatJet_pt_nom[xFatH], ch.FatJet_eta[xFatH],
                           ch.FatJet_phi[xFatH], ch.FatJet_mass_nom[xFatH])
        xFatTop,xFatWZ = -99,-99
        if len(xFatTops) > 0:
            if xFatTops[0] != xFatH:
                xFatTop = xFatTops[0]
            elif len(xFatTops) > 1:
                xFatTop = xFatTops[1]
        if len(xFatWZs) > 0:
            if xFatWZs[0] != xFatH:
                xFatWZ = xFatWZs[0]
            elif len(xFatWZs) > 1:
                xFatWZ = xFatWZs[1]
        if xFatTop >= 0:
            vFatTop = R.TLorentzVector()
            vFatTop.SetPtEtaPhiM(ch.FatJet_pt_nom[xFatTop], ch.FatJet_eta[xFatTop],
                                 ch.FatJet_phi[xFatTop], ch.FatJet_mass_nom[xFatTop])
        if xFatWZ >= 0:
            vFatWZ = R.TLorentzVector()
            vFatWZ.SetPtEtaPhiM(ch.FatJet_pt_nom[xFatWZ], ch.FatJet_eta[xFatWZ],
                                 ch.FatJet_phi[xFatWZ], ch.FatJet_mass_nom[xFatWZ])


        ## Trigger: Table 18 in Appendix C of AN2023_047_v4
        if YEAR == '2018':
            sel['trigJetHT'] = (ch.HLT_AK8PFJet500 or ch.HLT_PFHT1050 or ch.HLT_PFJet500 or \
                                ch.HLT_AK8PFHT800_TrimMass50 or ch.HLT_AK8PFJet400_TrimMass30)
            sel['trigBTag']  = (ch.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 or \
                                ch.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 or \
                                ch.HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2 or \
                                ch.HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 or \
                                ch.HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5)
            sel['trigMET'] = (ch.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or \
                              ch.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 or \
                              ch.HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1 or \
                              ch.HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned or \
                              ch.HLT_PFMETTypeOne140_PFMHT140_IDTight)
        elif YEAR == '2017':
            sel['trigJetHT'] = (ch.HLT_AK8PFJet500 or ch.HLT_PFHT1050 or ch.HLT_PFJet500 or \
                                ch.HLT_AK8PFHT800_TrimMass50 or ch.HLT_AK8PFJet400_TrimMass30 or \
                                ch.HLT_AK8PFJet360_TrimMass30)
            sel['trigBTag'] = (ch.HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2 or \
                               ch.HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2 or \
                               ch.HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5 or \
                               ch.HLT_AK8PFHT750_TrimMass50 or ch.HLT_AK8PFJet380_TrimMass30 or \
                               ch.HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33 or \
                               ch.HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0)
            sel['trigMET'] = (ch.HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1 or \
                              ch.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or \
                              ch.HLT_PFMET120_PFMHT120_IDTight or \
                              ch.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 or \
                              ch.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or \
                              ch.HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60 or \
                              ch.HLT_PFMETTypeOne120_PFMHT120_IDTight or \
                              ch.HLT_PFMET140_PFMHT140_IDTight or \
                              ch.HLT_PFMETTypeOne140_PFMHT140_IDTight or \
                              ch.HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned or \
                              ch.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight)
        elif '2016' in YEAR:
            sel['trigJetHT'] = (ch.HLT_PFJet450 or ch.HLT_DiCentralPFJet430 or \
                                ch.HLT_PFHT650_WideJetMJJ900DEtaJJ1p5 or ch.HLT_PFHT750_4JetPt50 or \
                                ch.HLT_PFHT800 or ch.HLT_PFHT900 or ch.HLT_AK8PFJet360_TrimMass30 or \
                                ch.HLT_AK8PFJet450 or ch.HLT_AK8PFHT650_TrimR0p1PT0p03Mass50 or \
                                ch.HLT_AK8PFHT700_TrimR0p1PT0p03Mass50)
            sel['trigBTag']  = (ch.HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20 or \
                                ch.HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20 or \
                                ch.HLT_PFHT400_SixJet30_DoubleBTagCSV_p056 or \
                                ch.HLT_PFHT450_SixJet40_BTagCSV_p056 or \
                                ch.HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20 or \
                                ch.HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6 or \
                                ch.HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160 or \
                                ch.HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6 or \
                                ch.HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172 or \
                                ch.HLT_DoubleJet90_Double30_TripleBTagCSV_p08 or \
                                ch.HLT_QuadJet45_TripleBTagCSV_p087 or \
                                ch.HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460 or \
                                ch.HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500 or \
                                ch.HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200 or \
                                ch.HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240)
            ## TODO: several HLT paths missing from AN, see:
            ## https://mattermost.web.cern.ch/cms-exp/pl/t1z5gutjr3bmdy9hk5ifoxcwpr
            sel['trigMET'] = (ch.HLT_MET200 or ch.HLT_PFMET110_PFMHT110_IDTight or \
                              ch.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or \
                              ch.HLT_PFMET120_PFMHT120_IDTight or \
                              ch.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or \
                              ch.HLT_PFMET170_HBHECleaned or \
                              ch.HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight)

        ## Muon selection from Table 2 in Ch. 4 of AN2023_047_v4
        hasTrgMu = False
        for iMu in range(ch.nMuon):
            if ch.Muon_pt[iMu] <= 26.0 + 3.0*(YEAR == '2017'): continue
            if abs(ch.Muon_eta[iMu]) >= 2.4:  continue
            if ch.Muon_mediumId[iMu] == 0 and \
               (ch.Muon_pt[iMu] <= 200 or
                ch.Muon_highPtId[iMu] == 0):  continue
            if ch.Muon_miniPFRelIso_all[iMu] >= 0.10: continue
            if abs(ch.Muon_dxy[iMu])         >= 0.02: continue
            if abs(ch.Muon_dz[iMu])          >= 0.10: continue
            hasTrgMu = True
            break
        ## End loop: for iMu in range(ch.nMuon)
        sel['mu_veto'] = (not hasTrgMu)

        ## Electron selection from Table 3 in Ch. 4 of AN2023_047_v4
        hasTrgEle = False
        for iEle in range(ch.nElectron):
            if ch.Electron_pt[iEle] <= 30.0 + 5.0*(not '2016' in YEAR): continue
            if         abs(ch.Electron_eta[iEle]) >= 2.5 or                                                             (1.44 <= abs(ch.Electron_eta[iEle]) <= 1.57): continue
            if ch.Electron_mvaFall17V2Iso_WP90[iEle] == 0: continue
            if ch.Electron_mvaFall17V2Iso_WP80[iEle] == 0 and \
               (ch.Electron_pt[iEle] <= 35 or
                ch.Electron_cutBased_HEEP[iEle] == 0): continue
            if abs(ch.Electron_dxy[iEle])     >= 0.02: continue
            if abs(ch.Electron_dz[iEle])      >= 0.10: continue
            hasTrgEle = True
            break
        ## End loop: for iEle in range(ch.nElectron)
        sel['ele_veto'] = (not hasTrgEle)

        ## AK4 jet selection from Table 4 in Ch. 4 of AN2023_047_v4
        nBJet = 0
        btagWPM = 0.2783 if YEAR == '2018' else (0.3040 if YEAR == '2017' else (0.2489 if 'ost' in YEAR else 0.2598))
        sel['bjet_veto_dR'] = True
        sel['bjet_veto'] = True
        vJets = []
        vJetsB = []
        vJetsLF = []
        vJetsBNonTop = []
        for iJet in range(ch.nJet):
            if   ch.Jet_pt_nom[iJet] <= 30: continue
            if    ch.Jet_jetId[iJet] != 6:  continue
            if   ch.Jet_puId[iJet]  < 4 and \
               ch.Jet_pt_nom[iJet] <= 50: continue
            vJet = R.TLorentzVector()
            vJet.SetPtEtaPhiM(ch.Jet_pt_nom[iJet], ch.Jet_eta[iJet],
                              ch.Jet_phi[iJet], ch.Jet_mass_nom[iJet])
            if vJet.DeltaR(vFatH) <= 0.8: continue
            vJets.append(vJet)
            ## Check for b-tagged AK4 jets
            if       abs(ch.Jet_eta[iJet]) < 2.4 and \
                ch.Jet_btagDeepFlavB[iJet] > btagWPM:
                sel['bjet_veto'] = False
                if vJet.DeltaR(vFatH) <= 1.2:
                    sel['bjet_veto_dR'] = False
                ## Store b-jets for e.g. VBFjj selection
                vJetsB.append(vJet)
                ## Store b-jets not overlapping AK8(top) for tt0l selection
                if xFatTop >= 0 and vJet.DeltaR(vFatTop) > 0.8:
                    vJetsBNonTop.append(vJet)
            else:
                ## Store non-b-jets for e.g. VBFjj selection
                vJetsLF.append(vJet)
        ## End loop: for iJet in range(ch.nJet)

        vMET = R.TLorentzVector()
        vMET.SetPtEtaPhiM(ch.MET_T1_pt, 0, ch.MET_T1_phi, 0)
        ## TODO: use MET_T1? or MET_T1_smear?

        sel['fatTop_veto'] = (xFatTop  < 0)
        sel['fatTop_sel']  = (xFatTop >= 0)
        sel['fatWZ_veto']  = (xFatWZ  < 0)
        sel['fatWZ_sel']   = (xFatWZ >= 0)
        sel['MET_veto'] = (vMET.Pt() <= 200)
        sel['MET_pt']   = (vMET.Pt()  > 200)
        sel['MET_dPhi'] = (abs(vMET.DeltaPhi(vFatH)) > (np.pi/2))
        
        if len(vJets) >= 2:
            j1,j2 = vJets[0],vJets[1]
        if len(vJetsLF) >= 2:
            q1,q2 = vJetsLF[0],vJetsLF[1]
        sel['VBFjj_veto'] = (len(vJets)    < 2  or abs(j1.Eta() - j2.Eta()) <= 2.2  or (j1+j2).M() <= 450)
        sel['VBFjj_sel']  = (len(vJetsLF) >= 2 and abs(q1.Eta() - q2.Eta())  > 2.2 and (q1+q2).M()  > 450)

        ## Category selection: prior cuts get applied in chain
        sel['gg0lHi'] = (vFatH.Pt()  < 400)
        sel['gg0lLo'] = (vFatH.Pt() >= 400)
        
        sel['VBFjjHi'] = (sel['VBFjj_sel'] and abs(q1.Eta() - q2.Eta()) > 3.0 and (q1+q2).M() > 900)
        sel['VBFjjLo'] = (sel['VBFjj_sel'] and not sel['VBFjjHi'])
        for xSel in ['VBFjjLo','VBFjjHi']:
            sel[xSel+'PtLo'] = (sel[xSel] and vFatH.Pt()  < 400)
            sel[xSel+'PtHi'] = (sel[xSel] and vFatH.Pt() >= 400)

        if sel['fatWZ_sel']:
            sel['VjjHi'] = (vFatWZ.Pt()  > 400)
            sel['VjjLo'] = (vFatWZ.Pt() <= 400)

        if sel['fatTop_sel']:
            sel['tt0l0b'] = (len(vJetsBNonTop) == 0)
            sel['tt0l1b'] = (len(vJetsBNonTop) >= 1)

        if sel['MET_pt'] and sel['MET_dPhi']:
            sel['ZvvHi'] = (vMET.Pt()  > 300)
            sel['ZvvLo'] = (vMET.Pt() <= 300)

        sel['X4bSB'] = (xFatH_X4b > 0.84)
        sel['X4bSR'] = (xFatH_X4b > 0.96)
        sel['4GenB'] = ch.FatJet_nBHadrons[xFatH] >= 4
        count[samp] = fill_counts(count[samp], sel)

    ## End loop: for iEvt in range(nEntries)

    print('\n*** Finished looking at sample %s ***\n' % samp)

## End loop: for samp in SAMPS

print('\n*** Finished looking at all samples! ***\n\n')

for samp in SAMPS:
    for cat in CATS.keys():
        print('\nCut-flow for %s in %s:' % (cat, samp))
        for iCuts in range(len(CATS[cat])):
            if len(CATS[cat][iCuts]) == 1:
                cutstr = CATS[cat][iCuts][0]
                key = '%d_%s' % (iCuts, cutstr)
                print('%s : %d' % (key, count[samp][cat][key]))
            else:
                for jCut in range(len(CATS[cat][iCuts])):
                    cutstr = CATS[cat][iCuts][jCut]
                    key = '%d_%s' % (iCuts, cutstr)
                    print('  * %s : %d' % (key, count[samp][cat][key]))
                print('%d_any : %d' % (iCuts, count[samp][cat]['%d_any' % iCuts]))
        ## End loop: for iCuts in range(len(CATS[cat]))
    ## End loop: for cat in CATS.keys()
## End loop: for samp in SAMPS

print('\n*** All done!!! ***\n\n')
