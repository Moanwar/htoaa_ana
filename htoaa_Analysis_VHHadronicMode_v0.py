#htoaa analysis main code

import os
import sys
from datetime import datetime
#import time
print(f"htoaa_Analysis_VHHadronicMode:: here1 {datetime.now() = }"); sys.stdout.flush()
import subprocess
import json
from urllib.request import urlopen
import glob
from collections import OrderedDict as OD
import time
import tracemalloc
import math
print(f"htoaa_Analysis_VHHadronicMode:: here2 {datetime.now() = }"); sys.stdout.flush()
import numpy as np
from copy import copy, deepcopy
print(f"htoaa_Analysis_VHHadronicMode:: here3 {datetime.now() = }"); sys.stdout.flush()
#import uproot
#import uproot3 as uproot
import uproot as uproot
print(f"htoaa_Analysis_VHHadronicMode:: here4 {datetime.now() = }"); sys.stdout.flush()
#import parse
from parse import *
print(f"htoaa_Analysis_VHHadronicMode:: here4.1 {datetime.now() = }"); sys.stdout.flush()
import logging
print(f"htoaa_Analysis_VHHadronicMode:: here5 {datetime.now() = }"); sys.stdout.flush()

# comment test3
'''
GGF -> H->aa->4b boosted analysis macro

References:
  * Coffea framework used for TTGamma analysis: https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py
* Coffea installation: /home/siddhesh/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea
'''
print(f"htoaa_Analysis_VHHadronicMode:: here6 {datetime.now() = }"); sys.stdout.flush()
#import coffea.processor as processor
from coffea import processor, util
from coffea.nanoevents import schemas
from coffea.nanoevents.methods import nanoaod, vector
from coffea.analysis_tools import PackedSelection, Weights
#from coffea.lookup_tools import extractor
from coffea.lookup_tools.dense_lookup import dense_lookup
#from coffea.lumi_tools import LumiMask
#import hist
from coffea import hist # /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea/hist/hist_tools.py
import awkward as ak
#import uproot
#from dask.distributed import Client
print(f"htoaa_Analysis_VHHadronicMode:: here7 {datetime.now() = }"); sys.stdout.flush()
from particle import Particle # For PDG particle listing https://github.com/scikit-hep/particle
print(f"htoaa_Analysis_VHHadronicMode:: here8 {datetime.now() = }"); sys.stdout.flush()


from htoaa_Settings import *
print(f"htoaa_Analysis_VHHadronicMode:: here9 {datetime.now() = }"); sys.stdout.flush()
from htoaa_CommonTools import (
    GetDictFromJsonFile, akArray_isin,
    selectRunLuminosityBlock,
    calculate_lumiScale, getLumiScaleForPhSpOverlapRewgtMode, getSampleHTRange, # update_crosssection, 
    getNanoAODFile, setXRootDRedirector,  xrdcpFile,
    selectMETFilters, selectAK4Jets, selectMuons, selectElectrons,
    selGenPartsWithStatusFlag,
    getHiggsPtRewgtForGGToHToAATo4B, getTopPtRewgt, getPURewgts, getHTReweight,
    calculateAverageOfArrays, calculateMaxOfTwoArrays, calculateMaxOfArrays,  array_PutLowerBound,
    printVariable, insertInListBeforeThisElement,
)
print(f"htoaa_Analysis_VHHadronicMode:: here10 {datetime.now() = }"); sys.stdout.flush()
from htoaa_Samples import (
    kData, kQCD_bEnrich, kQCD_bGen, kQCDIncl
)
print(f"htoaa_Analysis_VHHadronicMode:: here11 {datetime.now() = }"); sys.stdout.flush()

from inspect import currentframe, getframeinfo
print(f"htoaa_Analysis_VHHadronicMode:: here12 {datetime.now() = }"); sys.stdout.flush()
frameinfo = getframeinfo(currentframe())
print(f"htoaa_Analysis_VHHadronicMode:: here13 {datetime.now() = }"); sys.stdout.flush()


# use GOldenJSON

 
printLevel = 0
nEventToReadInBatch = 2*10**4 # 0.5*10**5 # 0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze = -1 # 1000 # 100000 # -1
flushStdout = False
#pd.set_option('display.max_columns', None)  

#print("".format())

sWeighted = "Wtd: "




    
# -----------------------------------------------------------------------------------
def get_GenPartDaughters(awkArray, index_GenPart):
    if printLevel >= 9:
        print(f"\n get_GenPartDaughters:: awkArray: {awkArray},   index_GenPart: {index_GenPart}"); sys.stdout.flush()
    
    return False


# -----------------------------------------------------------------------------------
class ObjectSelection:
    def __init__(self, era):
        self.era = era
        
        self.tagger_btagDeepB = 'DeepCSV'
        self.wp_btagDeepB = 'M'
        self.wp_ParticleNetMD_XbbvsQCD = 'L'
        self.wp_ParticleNetMD_Hto4b_Htoaa4bOverQCD = 'WP-40' # 'WP-40' 'WP-80' #'WP-60'

        self.FatJetPtThsh  = 250 # 170 # 400, 170
        self.FatJetEtaThsh = 2.4
        self.FatJetJetID   = int(JetIDs.tightIDPassingLeptonVeto)

        self.FatJetMSoftDropThshLow  = 50 # 20 # 90
        self.FatJetMSoftDropThshHigh = 9999 #200

        self.FatJetParticleNetMD_Xbb_Thsh       = 0.8
        self.FatJetParticleNetMD_XbbvsQCD_Thsh  = bTagWPs[self.era]['ParticleNetMD_XbbvsQCD'][self.wp_ParticleNetMD_XbbvsQCD]
        self.FatJetDeepTagMD_bbvsLight_Thsh     = 0.98
        self.FatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_Thsh   = bTagWPs[self.era]['ParticleNetMD_Hto4b_Htoaa4bOverQCD'][self.wp_ParticleNetMD_Hto4b_Htoaa4bOverQCD]
        self.FatJetZHbb_plus_Xbb_Thsh = 0.4
        self.FatJetZHbb_Xbb_avg_Thsh  = 0.4
        self.FatJetZHbb_Thsh          = 0.7

        self.nSV_matched_leadingFatJet_Thsh = 3

    
        self.NonHto4bFatJetPNet_WZvsQCD_Thsh = 0.98 # 0.94

        self.MuonPtThsh         = 10
        self.MuonMVAId          =  3 # (1=MvaLoose, 2=MvaMedium, 3=MvaTight, 4=MvaVTight, 5=MvaVVTight)
        self.MuonMiniIsoId      =  3 # (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)
        self.MuonMVATTHThsh     = 0.5
        self.ElectronPtThsh     = 10
        self.ElectronMVAId      = 'mvaFall17V2Iso_WP80' # 'mvaFall17V2Iso_WP80', 'mvaFall17V2Iso_WP90' 'mvaFall17V2Iso_WPL'
        self.ElectronMVATTHThsh = 0.3

        self.NLeptonsTight_MaxThsh    = 0
                
        self.Ak4JetDeepJetB_Thsh = bTagWPs[self.era]['AK4DeepJet']['M']

        self.JetPtThshForHT = 30.0
        self.JetEtaThshForHT = 2.4
        
        self.nFatJetMin = 1
        self.GenHTThsh  = 100.0
        self.LHEHTThsh  = 100.0
        

    def selectFatJets(self, events):
        '''
        Select FatJets
        '''
        '''
        jets.cut(jets['FatJet_pt'] > 170)
        jets.cut(jets['FatJet_eta'].abs() < 2.4)
        jets.cut(jets['FatJet_btagDDBvL'] > 0.8)
        jets.cut(jets['FatJet_btagDeepB'] > 0.4184)
        jets.cut(jets['FatJet_msoftdrop'] > 90)
        jets.cut(jets['FatJet_msoftdrop'] < 200)#<= 200)
        jets.cut(jets['FatJet_mass'] > 90)
        #jets.cut(jets['FatJet_mass'] <= 200)
        other.cut(other['PV_npvsGood'] >= 1)
        '''

        maskSelFatJets = (
            (events.FatJet.pt > self.FatJetPtThsh) &
            (abs(events.FatJet.eta) < self.FatJetEtaThsh) &
            (events.FatJet.btagDeepB > bTagWPs[self.era][self.tagger_btagDeepB][self.wp_btagDeepB])
        )
        if printLevel >= 15:
            #print(f"era: {self.era}, bTagWPs[self.era]: {bTagWPs[self.era]}")
            print(f"selectFatJets()::maskSelFatJets {len(maskSelFatJets)}: {maskSelFatJets.to_list()}")
        return events.FatJet[maskSelFatJets]




    def selectGenHiggs(self, events):
        maskGenHiggs = (
            (events.GenPart.pdgId  == 25) & # pdgId:: 25: H0
            (events.GenPart.status == 62)   # statu 62: outgoing subprocess particle with primordial kT included https://pythia.org/latest-manual/ParticleProperties.html
        )
        if printLevel >= 13:
            print(f"\n maskGenHiggs:  {maskGenHiggs.to_list()} ")
            print(f"\n events.GenPart[maskGenHiggs]:  {events.GenPart[maskGenHiggs].to_list()} ")
            print(f"\n events.GenPart[maskGenHiggs].pt:  {events.GenPart[maskGenHiggs].pt.to_list()} ")
            print(f"\n events.GenPart[maskGenHiggs].mass:  {events.GenPart[maskGenHiggs].mass.to_list()} ")
        return events.GenPart[maskGenHiggs]

    def selectGenABoson(self, events):
        maskGenA = (
            (events.GenPart.pdgId == 36)
        )
        if printLevel >= 15:
            print(f"\n maskGenA:  {maskGenA.to_list()} ")
            print(f"\n events.GenPart[maskGenA]:  {events.GenPart[maskGenA].to_list()} ")
            print(f"\n events.GenPart[maskGenA].mass:  {events.GenPart[maskGenA].mass.to_list()} ")
        return events.GenPart[maskGenA]


    def GenHT(self, events):
        maskForGenHT = (
            (events.GenJet.pt > self.JetPtThshForHT) &
            (abs(events.GenJet.eta) < self.JetEtaThshForHT)
        )
        selGenJetPt = events.GenJet[maskForGenHT].pt
        GenHT = ak.sum(selGenJetPt, axis=-1)
        
        if printLevel >= 15:
            print(f"\nevents.GenJet.fields: {events.GenJet.fields}")
            print(f"\nevents.GenJet: {events.GenJet.to_list()}")
            print(f"\nevents.GenJet.pt ({len(events.GenJet.pt)}): {events.GenJet.pt.to_list()} ")
            print(f"\nmaskForGenHT: {maskForGenHT.to_list()} ")
            print(f"\nselGenJetPt: {selGenJetPt.to_list()} ")
            print(f"\nGenHT ({len(GenHT)}): {GenHT.to_list()} ")
            
        return GenHT



'''
 OBJ: TBranch	nFatJet	slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis : 0 at: 0x2733fa0
 OBJ: TBranch	FatJet_CHS_mass	CHS mass : 0 at: 0x2734550
 OBJ: TBranch	FatJet_CHS_mpruned	CHS pruned mass : 0 at: 0x2734ae0
 OBJ: TBranch	FatJet_CHS_msoftdrop	CHS soft drop mass : 0 at: 0x27350d0
 OBJ: TBranch	FatJet_CHS_pt	CHS pt : 0 at: 0x2735710
 OBJ: TBranch	FatJet_CHS_tau1	CHS Nsubjettiness (1 axis) : 0 at: 0x2735ce0
 OBJ: TBranch	FatJet_CHS_tau2	CHS Nsubjettiness (2 axis) : 0 at: 0x2736320
 OBJ: TBranch	FatJet_CHS_tau3	CHS Nsubjettiness (3 axis) : 0 at: 0x2736960
 OBJ: TBranch	FatJet_CHS_tau4	CHS Nsubjettiness (4 axis) : 0 at: 0x2736fa0
 OBJ: TBranch	FatJet_area	jet catchment area, for JECs : 0 at: 0x27375e0
 OBJ: TBranch	FatJet_btagCSVV2	 pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2) : 0 at: 0x2737be0
 OBJ: TBranch	FatJet_btagDDBvLV2	DeepDoubleX V2(mass-decorrelated) discriminator for H(Z)->bb vs QCD : 0 at: 0x2738270
 OBJ: TBranch	FatJet_btagDDCvBV2	DeepDoubleX V2 (mass-decorrelated) discriminator for H(Z)->cc vs H(Z)->bb : 0 at: 0x2738900
 OBJ: TBranch	FatJet_btagDDCvLV2	DeepDoubleX V2 (mass-decorrelated) discriminator for H(Z)->cc vs QCD : 0 at: 0x2738f90
 OBJ: TBranch	FatJet_btagDeepB	DeepCSV b+bb tag discriminator : 0 at: 0x2739620
 OBJ: TBranch	FatJet_btagHbb	Higgs to BB tagger discriminator : 0 at: 0x2739c80
 OBJ: TBranch	FatJet_chEmEF	charged Electromagnetic Energy Fraction : 0 at: 0x273a290
 OBJ: TBranch	FatJet_chHEF	charged Hadron Energy Fraction : 0 at: 0x273a8a0
 OBJ: TBranch	FatJet_deepTagMD_H4qvsQCD	Mass-decorrelated DeepBoostedJet tagger H->4q vs QCD discriminator : 0 at: 0x273aea0
 OBJ: TBranch	FatJet_deepTagMD_HbbvsQCD	Mass-decorrelated DeepBoostedJet tagger H->bb vs QCD discriminator : 0 at: 0x273b540
 OBJ: TBranch	FatJet_deepTagMD_TvsQCD	Mass-decorrelated DeepBoostedJet tagger top vs QCD discriminator : 0 at: 0x273bbe0
 OBJ: TBranch	FatJet_deepTagMD_WvsQCD	Mass-decorrelated DeepBoostedJet tagger W vs QCD discriminator : 0 at: 0x273c280
 OBJ: TBranch	FatJet_deepTagMD_ZHbbvsQCD	Mass-decorrelated DeepBoostedJet tagger Z/H->bb vs QCD discriminator : 0 at: 0x273c910
 OBJ: TBranch	FatJet_deepTagMD_ZHccvsQCD	Mass-decorrelated DeepBoostedJet tagger Z/H->cc vs QCD discriminator : 0 at: 0x273cfb0
 OBJ: TBranch	FatJet_deepTagMD_ZbbvsQCD	Mass-decorrelated DeepBoostedJet tagger Z->bb vs QCD discriminator : 0 at: 0x273d650
 OBJ: TBranch	FatJet_deepTagMD_ZvsQCD	Mass-decorrelated DeepBoostedJet tagger Z vs QCD discriminator : 0 at: 0x273dcf0
 OBJ: TBranch	FatJet_deepTagMD_bbvsLight	Mass-decorrelated DeepBoostedJet tagger Z/H/gluon->bb vs light flavour discriminator : 0 at: 0x273e380
 OBJ: TBranch	FatJet_deepTagMD_ccvsLight	Mass-decorrelated DeepBoostedJet tagger Z/H/gluon->cc vs light flavour discriminator : 0 at: 0x273ea30
 OBJ: TBranch	FatJet_deepTag_H	DeepBoostedJet tagger H(bb,cc,4q) sum : 0 at: 0x273f0e0
 OBJ: TBranch	FatJet_deepTag_QCD	DeepBoostedJet tagger QCD(bb,cc,b,c,others) sum : 0 at: 0x273f750
 OBJ: TBranch	FatJet_deepTag_QCDothers	DeepBoostedJet tagger QCDothers value : 0 at: 0x273fdc0
 OBJ: TBranch	FatJet_deepTag_TvsQCD	DeepBoostedJet tagger top vs QCD discriminator : 0 at: 0x2740440
 OBJ: TBranch	FatJet_deepTag_WvsQCD	DeepBoostedJet tagger W vs QCD discriminator : 0 at: 0x2740ab0
 OBJ: TBranch	FatJet_deepTag_ZvsQCD	DeepBoostedJet tagger Z vs QCD discriminator : 0 at: 0x2741120
 OBJ: TBranch	FatJet_eta	eta : 0 at: 0x2741790
 OBJ: TBranch	FatJet_mass	mass : 0 at: 0x2741d60
 OBJ: TBranch	FatJet_msoftdrop	Corrected soft drop mass with PUPPI : 0 at: 0x2742330
 OBJ: TBranch	FatJet_msoftdrop_uncorr	Uncorrected soft drop mass with PUPPI : 0 at: 0x27429a0
 OBJ: TBranch	FatJet_muEF	muon Energy Fraction : 0 at: 0x2743020
 OBJ: TBranch	FatJet_n2b1	N2 with beta=1 : 0 at: 0x2743620
 OBJ: TBranch	FatJet_n2b2	N2 with beta=2 : 0 at: 0x2743bf0
 OBJ: TBranch	FatJet_n3b1	N3 with beta=1 : 0 at: 0x27441c0
 OBJ: TBranch	FatJet_n3b2	N3 with beta=2 : 0 at: 0x2744790
 OBJ: TBranch	FatJet_neEmEF	neutral Electromagnetic Energy Fraction : 0 at: 0x2744d60
 OBJ: TBranch	FatJet_neHEF	neutral Hadron Energy Fraction : 0 at: 0x2745370
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_Haa01b	Mass-decorrelated multiclass ParticleNet tagger for H->aa->bbbb : 0 at: 0x2745970
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_Haa2b	Mass-decorrelated multiclass ParticleNet tagger for H->aa->bbbb : 0 at: 0x2746020
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_Haa3b	Mass-decorrelated multiclass ParticleNet tagger for H->aa->bbbb : 0 at: 0x27466d0
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_Haa4b	Mass-decorrelated multiclass ParticleNet tagger for H->aa->bbbb : 0 at: 0x2746d80
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_QCD0b	Mass-decorrelated multiclass ParticleNet tagger for H->aa->bbbb : 0 at: 0x2747430
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_QCD1b	Mass-decorrelated multiclass ParticleNet tagger for H->aa->bbbb : 0 at: 0x2747ae0
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_QCD2b	Mass-decorrelated multiclass ParticleNet tagger for H->aa->bbbb : 0 at: 0x2748190
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_QCD3b	Mass-decorrelated multiclass ParticleNet tagger for H->aa->bbbb : 0 at: 0x2748840
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_QCD4b	Mass-decorrelated multiclass ParticleNet tagger for H->aa->bbbb : 0 at: 0x2748ef0
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_binaryLF_Haa4b	Mass-decorrelated binary ParticleNet tagger for H->aa->bbbb vs. QCD (<=2 b) : 0 at: 0x27495a0
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_binaryLF_QCDlf	Mass-decorrelated binary ParticleNet tagger for H->aa->bbbb vs. QCD (<= 2b) : 0 at: 0x2749c70
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_binary_Haa4b	Mass-decorrelated binary ParticleNet tagger for H->aa->bbbb vs. QCD : 0 at: 0x274a340
 OBJ: TBranch	FatJet_particleNetMD_Hto4b_binary_QCD	Mass-decorrelated binary ParticleNet tagger for H->aa->bbbb vs. QCD : 0 at: 0x274aa10
 OBJ: TBranch	FatJet_particleNetMD_QCD	Mass-decorrelated ParticleNet tagger raw QCD score : 0 at: 0x274b0d0
 OBJ: TBranch	FatJet_particleNetMD_Xbb	Mass-decorrelated ParticleNet tagger raw X->bb score. For X->bb vs QCD tagging, use Xbb/(Xbb+QCD) : 0 at: 0x274b760
 OBJ: TBranch	FatJet_particleNetMD_Xcc	Mass-decorrelated ParticleNet tagger raw X->cc score. For X->cc vs QCD tagging, use Xcc/(Xcc+QCD) : 0 at: 0x274be20
 OBJ: TBranch	FatJet_particleNetMD_Xqq	Mass-decorrelated ParticleNet tagger raw X->qq (uds) score. For X->qq vs QCD tagging, use Xqq/(Xqq+QCD). For W vs QCD tagging, use (Xcc+Xqq)/(Xcc+Xqq+QCD) : 0 at: 0x274c4e0
 OBJ: TBranch	FatJet_particleNet_H4qvsQCD	ParticleNet tagger H(->VV->qqqq) vs QCD discriminator : 0 at: 0x274cbd0
 OBJ: TBranch	FatJet_particleNet_HbbvsQCD	ParticleNet tagger H(->bb) vs QCD discriminator : 0 at: 0x274d260
 OBJ: TBranch	FatJet_particleNet_HccvsQCD	ParticleNet tagger H(->cc) vs QCD discriminator : 0 at: 0x274d8e0
 OBJ: TBranch	FatJet_particleNet_QCD	ParticleNet tagger QCD(bb,cc,b,c,others) sum : 0 at: 0x274df60
 OBJ: TBranch	FatJet_particleNet_TvsQCD	ParticleNet tagger top vs QCD discriminator : 0 at: 0x274e5d0
 OBJ: TBranch	FatJet_particleNet_WvsQCD	ParticleNet tagger W vs QCD discriminator : 0 at: 0x274ec50
 OBJ: TBranch	FatJet_particleNet_ZvsQCD	ParticleNet tagger Z vs QCD discriminator : 0 at: 0x274f2d0
 OBJ: TBranch	FatJet_particleNet_mass	ParticleNet mass regression : 0 at: 0x274f950
 OBJ: TBranch	FatJet_particleNet_massA_Hto4b_v0	ParticleNet a boson mass regression for H->aa->bbbb (mass regr, mode 3 [original]) : 0 at: 0x274ffc0
 OBJ: TBranch	FatJet_particleNet_massA_Hto4b_v1	ParticleNet a boson mass regression for H->aa->bbbb (mass regr, mode 3 [final]) : 0 at: 0x2750690
 OBJ: TBranch	FatJet_particleNet_massA_Hto4b_v2	ParticleNet a boson mass regression for H->aa->bbbb (mass regr, mode 0) : 0 at: 0x2750d50
 OBJ: TBranch	FatJet_particleNet_massA_Hto4b_v3	ParticleNet a boson mass regression for H->aa->bbbb (logMass regr, mode 0) : 0 at: 0x2751410
 OBJ: TBranch	FatJet_particleNet_massA_Hto4b_v4	ParticleNet a boson mass regression for H->aa->bbbb (mass ratio regr, mode 3) : 0 at: 0x2751ad0
 OBJ: TBranch	FatJet_particleNet_massH_Hto4b_v0	ParticleNet Higgs mass regression for H->aa->bbbb (mass regr, mode 3 [2nd try]) : 0 at: 0x2752190
 OBJ: TBranch	FatJet_particleNet_massH_Hto4b_v00	ParticleNet Higgs mass regression for H->aa->bbbb (mass regr, mode 3 [original]) : 0 at: 0x2752850
 OBJ: TBranch	FatJet_particleNet_massH_Hto4b_v1	ParticleNet Higgs mass regression for H->aa->bbbb (mass regr, mode 3 [final]) : 0 at: 0x2752f20
 OBJ: TBranch	FatJet_particleNet_massH_Hto4b_v2	ParticleNet Higgs mass regression for H->aa->bbbb (mass regr, mode 0) : 0 at: 0x27535e0
 OBJ: TBranch	FatJet_particleNet_massH_Hto4b_v3	ParticleNet Higgs mass regression for H->aa->bbbb (logMass regr, mode 0) : 0 at: 0x2753ca0
 OBJ: TBranch	FatJet_particleNet_massH_Hto4b_v4	ParticleNet Higgs mass regression for H->aa->bbbb (mass ratio regr, mode 3) : 0 at: 0x2754360
 OBJ: TBranch	FatJet_phi	phi : 0 at: 0x215e6a0
 OBJ: TBranch	FatJet_pt	pt : 0 at: 0x215ec70
 OBJ: TBranch	FatJet_rawFactor	1 - Factor to get back to raw pT : 0 at: 0x215f240
 OBJ: TBranch	FatJet_tau1	Nsubjettiness (1 axis) : 0 at: 0x215f8b0
 OBJ: TBranch	FatJet_tau2	Nsubjettiness (2 axis) : 0 at: 0x215feb0
 OBJ: TBranch	FatJet_tau3	Nsubjettiness (3 axis) : 0 at: 0x21604b0
 OBJ: TBranch	FatJet_tau4	Nsubjettiness (4 axis) : 0 at: 0x2160ab0
 OBJ: TBranch	FatJet_lsf3	Lepton Subjet Fraction (3 subjets) : 0 at: 0x21610b0
 OBJ: TBranch	FatJet_chMultiplicity	(Puppi-weighted) Number of charged particles in the jet : 0 at: 0x21616c0
 OBJ: TBranch	FatJet_charge	jet charge (integer) : 0 at: 0x2161d40
 OBJ: TBranch	FatJet_jetId	Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto : 0 at: 0x2162340
 OBJ: TBranch	FatJet_neMultiplicity	(Puppi-weighted) Number of neutral particles in the jet : 0 at: 0x2162990
 OBJ: TBranch	FatJet_subJetIdx1	index of first subjet : 0 at: 0x2163010
 OBJ: TBranch	FatJet_subJetIdx2	index of second subjet : 0 at: 0x2163670
 OBJ: TBranch	FatJet_electronIdx3SJ	index of electron matched to jet : 0 at: 0x2163cd0
 OBJ: TBranch	FatJet_muonIdx3SJ	index of muon matched to jet : 0 at: 0x27605b0
 OBJ: TBranch	FatJet_nConstituents	Number of particles in the jet : 0 at: 0x2760c10

'''

    
    
class HToAATo4bProcessor(processor.ProcessorABC):
    def __init__(self, datasetInfo={}):
        print(f"HToAATo4bProcessor::__init__():: {datasetInfo = }")
         
        global runMode_SignalGenChecks;       runMode_SignalGenChecks  = False; # True
        global runMode_QCDGenValidation;      runMode_QCDGenValidation = False; # True
        global runMode_GenLHEPlots;           runMode_GenLHEPlots      = False
        global runMode_SignificancsScan2D;    runMode_SignificancsScan2D = False
        global runMode_OptimizePNetTaggerCut; runMode_OptimizePNetTaggerCut = False
        global runMode_2018HEM1516IssueValidation; runMode_2018HEM1516IssueValidation = False
        global runMode_SignalGenCuts;         runMode_SignalGenCuts = True; # set False for final round. True for optimization studies.
        
        

        ak.behavior.update(nanoaod.behavior)

        self.datasetInfo = datasetInfo
        self.objectSelector = ObjectSelection(era=self.datasetInfo["era"])
        datasetName_part1 = self.datasetInfo['datasetNameFull'].split('/')[1]
        print(f"{datasetName_part1 = }")

        # Identify and lable samples --------------------------------------------------
        self.datasetInfo['isSignal'       ]  = False
        self.datasetInfo['isSignalGGH'    ]  = False
        self.datasetInfo['isSignalVBFH'   ]  = False
        self.datasetInfo['isSignalWH'     ]  = False
        self.datasetInfo['isSignalVH'     ]  = False
        self.datasetInfo['isSignalTTH'    ]  = False
        self.datasetInfo['isQCD'          ]  = False
        self.datasetInfo['isQCDIncl'      ]  = False
        self.datasetInfo['isQCD_bEnrich'  ]  = False
        self.datasetInfo['isQCD_bGen'     ]  = False
        self.datasetInfo['isTTbar'        ]  = False
        self.datasetInfo['isHToBB'        ]  = False
        self.datasetInfo['isPythiaTuneCP5']  = False        
        if self.datasetInfo['isMC']:
            self.datasetInfo['isSignalGGH']      = True if "SUSY_GluGluH_01J_HToAATo4B" in datasetName_part1 else False
            self.datasetInfo['isSignalVBFH']     = True if "SUSY_VBFH_HToAATo4B"        in datasetName_part1 else False
            self.datasetInfo['isSignalWH']       = True if "SUSY_WH_WToAll_HToAATo4B"   in datasetName_part1 else False
            self.datasetInfo['isSignalVH']       = True if "SUSY_ZH_ZToAll_HToAATo4B"   in datasetName_part1 else False
            self.datasetInfo['isSignalTTH']      = True if "SUSY_TTH_TTToAll_HToAATo4B" in datasetName_part1 else False
            self.datasetInfo['isSignal']         = (self.datasetInfo['isSignalGGH']   or \
                                                     self.datasetInfo['isSignalVBFH'] or \
                                                     self.datasetInfo['isSignalWH']   or \
                                                     self.datasetInfo['isSignalVH']   or \
                                                     self.datasetInfo['isSignalTTH'] )
            self.datasetInfo['isQCDIncl']        = True if kQCDIncl      in self.datasetInfo['sample_category'] else False
            self.datasetInfo['isQCD_bEnrich']    = True if kQCD_bEnrich  in self.datasetInfo['sample_category'] else False
            self.datasetInfo['isQCD_bGen']       = True if kQCD_bGen     in self.datasetInfo['sample_category'] else False
            self.datasetInfo['isQCD']            = (self.datasetInfo['isQCDIncl']     or \
                                                     self.datasetInfo['isQCD_bEnrich'] or \
                                                     self.datasetInfo['isQCD_bGen'])
            sample_HT_Min, sample_HT_Max = getSampleHTRange( self.datasetInfo["datasetNameFull"] )
            self.datasetInfo['sample_HT_Min']    = sample_HT_Min
            self.datasetInfo['sample_HT_Max']    = sample_HT_Max
            self.datasetInfo['isTTbar']          = True if datasetName_part1.startswith('TTTo') else False
            self.datasetInfo['isHToBB']          = True if "HToBB"   in datasetName_part1 else False
            self.datasetInfo['isPythiaTuneCP5']  = True if 'TuneCP5' in datasetName_part1 else False

            if self.datasetInfo['isQCD_bGen']:
                # 'Corrections' variable defined in htoaa_Settings.py
                fitFunctionFormat_  = Corrections["HTRewgt"]["QCD_bGen"][self.datasetInfo["era"]]["FitFunctionFormat"] 
                fitFunctionHTRange_ = ""
                for sHTBin in Corrections["HTRewgt"]["QCD_bGen"][self.datasetInfo["era"]]:
                    if "HT%dto" % (self.datasetInfo['sample_HT_Min']) in sHTBin:
                        fitFunctionHTRange_ = sHTBin
                        fitFunction_  = Corrections["HTRewgt"]["QCD_bGen"][self.datasetInfo["era"]][sHTBin]

                        
                self.datasetInfo['HTRewgt'] = {
                    "fitFunctionFormat":  fitFunctionFormat_,
                    "fitFunction":        fitFunction_,
                    "fitFunctionHTRange": fitFunctionHTRange_,                    
                }

            print(f"{self.datasetInfo['isSignal'       ] = }")
            print(f"{self.datasetInfo['isSignalGGH'    ] = }")
            print(f"{self.datasetInfo['isSignalVBFH'   ] = }")
            print(f"{self.datasetInfo['isSignalWH'     ] = }")
            print(f"{self.datasetInfo['isSignalVH'     ] = }")
            print(f"{self.datasetInfo['isSignalTTH'    ] = }")
            print(f"{self.datasetInfo['isQCD'          ] = }")
            print(f"{self.datasetInfo['isQCDIncl'      ] = }")
            print(f"{self.datasetInfo['isQCD_bEnrich'  ] = }")
            print(f"{self.datasetInfo['isQCD_bGen'     ] = }")
            print(f"{self.datasetInfo['isTTbar'        ] = }")
            print(f"{self.datasetInfo['isHToBB'        ] = }")
            print(f"{self.datasetInfo['isPythiaTuneCP5'] = }")
            

        ## List of all analysis selection condition ---------------------------------------------
        global HLT_AK8PFJet330_name
        HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4" 
        global sTrgSelection
        sTrgSelection = 'Trg_Combo_AK4AK8Jet_HT' # 'Trg_Combo_AK4AK8Jet_HT'
        
        # sel_names_all = dict of {"selection name" : [list of different cuts]}; for cut-flow table 
        self.sel_names_all = OD([
            ("Presel",                    [
                "nPV",
                "METFilters",
                "leadingFatJetPt",
                "leadingFatJetEta",
                "JetID",
                "leadingNonHto4bFatJetPt", # leadingNonHto4bFatJet leadingNonHto4bFatJet_asSingletons
                "leadingNonHto4bFatJetEta",
                "leadingNonHto4bFatJetJetID",    
                #"L1_SingleJet180",
                #HLT_AK8PFJet330_name,
                sTrgSelection,
                #"leadingFatJetBtagDeepB",
                "leadingFatJetMSoftDrop",
                "leadingNonHto4bFatMSoftDrop",             
                #"leadingFatJetZHbb_Xbb_avg",
                "leadingFatJetZHbb",
                #"leadingFatJetDeepTagMD_bbvsLight", #"leadingFatJetParticleNetMD_Xbb",
                #"leadingFatJetParticleNetMD_XbbvsQCD",
                #"leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD",
                #"leadingFatJet_nSV"
                #
                "PNetWZvsQCDLoose", 
                "nLeptonsTight",
            ]),
        ])

        
        if runMode_OptimizePNetTaggerCut:
            # Optimize ParticleNet tagger tagger cut by monitoring S/B
            self.sel_names_all = OD([
                ("Presel",                    [
                    "nPV",
                    "METFilters",
                    "leadingFatJetPt",
                    "leadingFatJetEta",
                    "JetID",
                    #"L1_SingleJet180",
                    #HLT_AK8PFJet330_name,
                    #"leadingFatJetBtagDeepB",
                    "leadingFatJetMSoftDrop",
                    "leadingFatJetZHbb_plus_Xbb",
                    #"leadingFatJetDeepTagMD_bbvsLight", #"leadingFatJetParticleNetMD_Xbb",
                    #"leadingFatJetParticleNetMD_XbbvsQCD",
                    #"leadingFatJet_nSV"
                    "leadingNonHto4bFatJetPt", # leadingNonHto4bFatJet leadingNonHto4bFatJet_asSingletons
                    "leadingNonHto4bFatJetEta",
                    "leadingNonHto4bFatJetJetID",                    
                ]),
            ])
            self.objectSelector.FatJetPtThsh = 170
            self.objectSelector.FatJetMSoftDropThshLow  =  60 # 90
            self.objectSelector.FatJetMSoftDropThshHigh = 160 # 140
            self.objectSelector.FatJetZHbb_plus_Xbb_Thsh = 0.4





        if not self.datasetInfo['isMC']: 
            self.sel_names_all["Presel"].insert(0, "run:ls")

        else:
            if self.datasetInfo['isQCD']: #self.sel_names_all["Presel"].append("QCDStitch")
                self.sel_names_all["Presel"] = insertInListBeforeThisElement(
                    list1                  = self.sel_names_all["Presel"], 
                    sConditionToAdd        = "QCDStitch", 
                    addBeforeThisCondition = "METFilters")                

        if self.datasetInfo["era"] == Era_2018:
            # 2018HEM1516Issue ----------------
            #self.sel_names_all["Presel"].append("2018HEM1516Issue")
            # Update self.sel_names_all["Presel"] by adding "2018HEM1516Issue" before "leadingFatJetMSoftDrop" in the list              
            self.sel_names_all["Presel"] = insertInListBeforeThisElement(
                list1                  = self.sel_names_all["Presel"], 
                sConditionToAdd        = "2018HEM1516Issue", 
                addBeforeThisCondition = "leadingFatJetMSoftDrop"
            )



        self.sel_names_all["SRWP40"] = self.sel_names_all["Presel"] + [
            "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP40"
        ]
        self.sel_names_all["SRWP60"] = self.sel_names_all["Presel"] + [
            "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP60"
        ]
        self.sel_names_all["SRWP80"] = self.sel_names_all["Presel"] + [
            "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP80"
        ]
        self.sel_names_all["SBWP80to40"] = self.sel_names_all["Presel"] + [
            "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP80to40"
        ]
        self.sel_names_all["SBWP95to60"] = self.sel_names_all["Presel"] + [
            "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP95to60"
        ]
        self.sel_names_all["SBWP99to80"] = self.sel_names_all["Presel"] + [
            "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP99to80"
        ]

        #for sSelName_ in ["SRWP40", "SRWP60", "SRWP80"]:
        for sSelName_ in []:
            # massA windows
            for sMassAWindowName in massPseudoscalarA_windows_dict:
                sSelWindowName_ = "%s_%s" % (sSelName_, sMassAWindowName)
                
                self.sel_names_all[sSelWindowName_] = self.sel_names_all[sSelName_] + [
                    "leadingFatJetParticleNet_massA_Hto4b_avg_%s" % (sMassAWindowName)
                ]

            # massH windows
            for sMassHiggsWindowName in massHiggs_windows_dict:
                sSelWindowName_ = "%s_%s" % (sSelName_, sMassHiggsWindowName)

                if 'Msoftdrop' in sMassHiggsWindowName:  # == 'MsoftdropHiggsWindow':
                    self.sel_names_all[sSelWindowName_] = self.sel_names_all[sSelName_] + [
                        "leadingFatJetMsoftdrop_%s" % (sMassHiggsWindowName)
                    ]
                elif 'PNet_massH_Hto4b' in sMassHiggsWindowName: # == 'PNet_massH_Hto4b_HiggsWindow':
                    self.sel_names_all[sSelWindowName_] = self.sel_names_all[sSelName_] + [
                        "leadingFatJetPNet_massH_Hto4b_%s" % (sMassHiggsWindowName)
                    ]



        

        self.sel_conditions_all_list = set()
        for sel_conditions_ in self.sel_names_all.values():
            self.sel_conditions_all_list.update( sel_conditions_ )
        print(f"{self.sel_conditions_all_list = }")
        
        # selection region addition each SR conditions successively
        #for iCondition in range(self.sel_names_all["Presel"].index(HLT_AK8PFJet330_name), len(self.sel_names_all["Presel"]) - 1):
        for iCondition in range(self.sel_names_all["Presel"].index("leadingFatJetMSoftDrop"), len(self.sel_names_all["Presel"]) - 1):
            conditionName = self.sel_names_all["Presel"][iCondition]
            self.sel_names_all["sel_%s" % conditionName] = self.sel_names_all["Presel"][0 : (iCondition+1)]
        
        print(f"self.sel_names_all: {json.dumps(self.sel_names_all, indent=4)}")
        




        self.histosExtensions = ['']
        dataLSSelGoldenJSON = None
        self.SFs_ParticleNetMD_XbbvsQCD = None
        
        if not self.datasetInfo['isMC']: ## Data
            # data LS selection Golden JSON
            dataLSSelGoldenJSON = None  
            print(f'{kData} {self.datasetInfo["era"]}: Reading {sFilesGoldenJSON[self.datasetInfo["era"]]} ')
            if 'https:' in sFilesGoldenJSON[self.datasetInfo["era"]]:
                with urlopen(sFilesGoldenJSON[self.datasetInfo["era"]]) as fDataGoldenJSON:
                    dataLSSelGoldenJSON = json.load(fDataGoldenJSON)
            else: 
                with open(sFilesGoldenJSON[self.datasetInfo["era"]]) as fDataGoldenJSON:
                    dataLSSelGoldenJSON = json.load(fDataGoldenJSON)
            if dataLSSelGoldenJSON == None:
                logging.critical(f'htoaa_Analysis_VHHadronicMode.py::main():: {sFilesGoldenJSON[self.datasetInfo["era"]] = } could not read.')
                exit(0) 

            # convert runNumber in str to int
            dataLSSelGoldenJSON = {int(k): v for k, v in dataLSSelGoldenJSON.items()} 
            self.datasetInfo['dataLSSelGoldenJSON'] = dataLSSelGoldenJSON
            #print(f"{dataLSSelGoldenJSON = }")

        else: ## MC

            # lumiScale --------------------------------------------------------------------------------------------------
            if sTrgSelection not in Luminosities_forGGFMode[self.datasetInfo["era"]]:
                logging.critical(f'htoaa_Analysis_VHHadronicMode.py::main():: {sTrgSelection = } not in {Luminosities_forGGFMode[self.datasetInfo["era"]] = }.')
                exit(0) 

            self.datasetInfo["lumiScale"] = calculate_lumiScale(
                luminosity   = Luminosities_forGGFMode[self.datasetInfo["era"]][sTrgSelection][0], 
                crossSection = self.datasetInfo["sample_crossSection"], 
                sumEvents    = self.datasetInfo["sample_sumEvents"])
            print(f'luminosity: {Luminosities_forGGFMode[self.datasetInfo["era"]][sTrgSelection][0] = }, \
                    crossSection: {self.datasetInfo["sample_crossSection"]}, \
                    sumEvents: {self.datasetInfo["sample_sumEvents"]}, \
                    lumiScale: {self.datasetInfo["lumiScale"] }')

            # MC PURewgt --------------------------------------------------------------------------------------------------
            print(f'MC {self.datasetInfo["era"]} PU reweighting:: ip file: {Corrections["PURewgt"][self.datasetInfo["era"]]["inputFile"]}, histogram: {Corrections["PURewgt"][self.datasetInfo["era"]]["histogramName"]} ')
            with uproot.open(Corrections["PURewgt"][self.datasetInfo["era"]]["inputFile"]) as f_:
                #print(f"{f_.keys() = }"); sys.stdout.flush() 
                self.hPURewgt = f_['%s' % Corrections["PURewgt"][self.datasetInfo["era"]]["histogramName"]].to_hist()
                
        
            # set self.pdgId_BHadrons for 'QCD_bGenFilter' sample requirement ---------------------------------------------
            self.pdgId_BHadrons = []
            bHadrons_ = Particle.findall(lambda p: p.pdgid.has_bottom) # Find all bottom hadrons
            print(f"List of B-hadrons for QCD B-GEN-filter ({len(bHadrons_)}):")
            print("%s %-20s %15s %15s" %(" "*4, "pdg name", "pdgId", "Mass in MeV"))
            for bHadron in bHadrons_:
                print("%s %-20s %15d %15s" % (" "*4, str(bHadron), bHadron.pdgid.abspid, str(bHadron.mass)))
                if bHadron.pdgid.abspid not in self.pdgId_BHadrons:
                    self.pdgId_BHadrons.append(bHadron.pdgid.abspid)
            print(f"self.pdgId_BHadrons ({len(self.pdgId_BHadrons)}): {self.pdgId_BHadrons}")
            #self.pdgId_BHadrons = list(set(self.pdgId_BHadrons))
            #print(f" after duplicate removal --> \nself.pdgId_BHadrons ({len(self.pdgId_BHadrons)}): {self.pdgId_BHadrons}")
            
            ## MC QCD
            if self.datasetInfo['isQCD'] and \
                self.datasetInfo["MCSamplesStitchOption"] == MCSamplesStitchOptions.PhSpOverlapRewgt and \
                SplitQCDInGENCats:
                # for QCD, make histograms in category of number of GEN b quarks matching to leading fat jet (AK8) 
                self.histosExtensions = HistogramNameExtensions_QCD 

            ## MC ParticleNetMD_XbbvsQCD SFs
            print(f" {Corrections['ParticleNetMD_XbbvsQCD'][self.datasetInfo['era']][self.objectSelector.wp_ParticleNetMD_XbbvsQCD]['SFs'] = } "); sys.stdout.flush()
            print(f" {Corrections['ParticleNetMD_XbbvsQCD'][self.datasetInfo['era']][self.objectSelector.wp_ParticleNetMD_XbbvsQCD]['pT_binEdges'] = } "); sys.stdout.flush()

            self.SFs_ParticleNetMD_XbbvsQCD = dense_lookup(
                np.array( Corrections['ParticleNetMD_XbbvsQCD'][self.datasetInfo["era"]][self.objectSelector.wp_ParticleNetMD_XbbvsQCD]['SFs'] ), # list of SFs
                [ np.array(Corrections['ParticleNetMD_XbbvsQCD'][self.datasetInfo["era"]][self.objectSelector.wp_ParticleNetMD_XbbvsQCD]['pT_binEdges']) ] # list of bin edges for all axes of SF histogram
                )

        
        
        #dataset_axis = hist.axis.StrCategory(name="dataset", label="", categories=[], growth=True)
        #muon_axis = hist.axis.Regular(name="massT", label="Transverse Mass [GeV]", bins=50, start=15, stop=250)
        dataset_axis    = hist.Cat("dataset", "Dataset")
        systematic_axis = hist.Cat("systematic", "Systematic Uncertatinty")

        

        cutFlow_axis          = hist.Bin("CutFlow",                r"Cuts",                       21,    -0.5,    20.5)
        cutFlow50_axis        = hist.Bin("CutFlow50",              r"Cuts",                       51,    -0.5,    50.5)
        nObject_axis          = hist.Bin("nObject",                r"No. of object",              21,    -0.5,    20.5)
        nObject10_axis        = hist.Bin("nObject10",              r"No. of object",              11,    -0.5,    10.5)
        nObject10_axis1       = hist.Bin("nObject10_1",            r"No. of object",              11,    -0.5,    10.5)
        nObject50_axis        = hist.Bin("nObject50",              r"No. of object",              51,    -0.5,    50.5)
        nObject200_axis       = hist.Bin("nObject200",             r"No. of object",             201,    -0.5,   200.5)
        pt4TeV_axis           = hist.Bin("Pt4TeV",                 r"$p_{T}$ [GeV]",             200,       0,    4000)
        pt_axis               = hist.Bin("Pt",                     r"$p_{T}$ [GeV]",             200,       0,    1000)
        ptLow_axis            = hist.Bin("PtLow",                  r"$p_{T}$ [GeV]",             400,       0,     200)
        ptUltraLow_axis       = hist.Bin("PtUltraLow",             r"$p_{T}$ [GeV]",             200,       0,     0.1)
        pt1to10_axis          = hist.Bin("Pt1to10",                r"$p_{T}$ [GeV]",             100,       0,      10)
        pt2TeV_axis           = hist.Bin("Pt2TeV",                 r"$p_{T}$ [GeV]",            2000,       0,    2000)
        log2Pt2TeV_axis       = hist.Bin("Log2Pt2TeV",             r"Log2($p_{T}$) [GeV]",      200,  math.log2(1),    math.log2(2000))
        eta_axis              = hist.Bin("Eta",                    r"$#eta$",                    100,      -6,       6)
        phi_axis              = hist.Bin("Phi",                    r"$\phi$",                    100,   -3.14,    3.13)
        #mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       200, 0, 600)
        #mass_axis             = hist.Bin("Mass",      r"$m$ [GeV]",       400, 0, 200)
        mass_axis             = hist.Bin("Mass",                   r"$m$ [GeV]",                 300,       0,     300)
        mass_axis1            = hist.Bin("Mass1",                  r"$m$ [GeV]",               20*70,       0,     70)
        mass_axis2            = hist.Bin("Mass2",                  r"$m$ [GeV]",                2*70,       0,     70)
        mass10_axis           = hist.Bin("Mass10",                 r"$m$ [GeV]",                 300,       0,      10)
        logMass3_axis         = hist.Bin("logMass3",               r"$m$ [GeV]",                 300,       0,       3)
        mlScore_axis          = hist.Bin("MLScore",                r"ML score",                  100,    -1.1,     1.1)
        mlScore_axis1         = hist.Bin("MLScore1",               r"ML score",                  100,    -1.1,     1.1)
        mlScore_axis1k        = hist.Bin("MLScore1k",              r"ML score",                 1100,     0.0,     1.1)
        mlScore_axis2k        = hist.Bin("MLScore2k",              r"ML score",                 2100,     0.0,     2.1)
        jetN2_axis            = hist.Bin("N2",                     r"N2b1",                      100,       0,       3)
        jetN3_axis            = hist.Bin("N3",                     r"N3b1",                      100,       0,       5)
        jetTau_axis           = hist.Bin("TauN",                   r"TauN",                      100,       0,       1)
        deltaR_axis           = hist.Bin("deltaR",                 r"$delta$ r ",                500,       0,       5)
        deltaPhi_axis         = hist.Bin("deltaPhi",               r"$delta$ phi ",             1000,       0,       3.14) # <<<<<<<<<
        #HT_axis               = hist.Bin("HT",                     r"HT",                       3000,       0,    3000)
        HT_axis               = hist.Bin("HT",                     r"HT",                       4000,       0,    4000)
        PytPartStatus_axis    = hist.Bin("PytPartStatus",          r"PytPartStatus",             421,  -210.5,   210.5)
        boolean_axis          = hist.Bin("Boolean",                r"Boolean",                     2,    -0.5,     1.5)
        pdgId_axis            = hist.Bin("PdgId",                  r"PdgId",                     101,    -0.5,   100.5)
        alphaS_axis           = hist.Bin("alphaS",                 r"alphaS",                    101,    0.01,     0.2)
        PU_axis               = hist.Bin("PU",                     r"PU",                         99,     0.0,    99.0)
        Ratio_axis            = hist.Bin("Ratio",                  r"Ratio",                     100,     0.0,    2.0)
        
        sXaxis      = 'xAxis'
        sXaxisLabel = 'xAxisLabel'
        sYaxis      = 'yAxis'
        sYaxisLabel = 'yAxisLabel'

        # General or GEN-level histograms ---------------------------------------------------------------------------------------------
        histos = OD([
            # ('histogram_name',  {sXaxis: hist.Bin() axis,  sXaxisLabel: "histogram axis label"})

            ('hPV_npvs_beforeSel',                        {sXaxis: PU_axis,                sXaxisLabel: r"No. of primary vertices - before selection"}),
            ('hPV_npvsGood_beforeSel',                    {sXaxis: PU_axis,                sXaxisLabel: r"No. of good primary vertices - before selection"}),            
        ])

        if self.datasetInfo['isMC'] and runMode_GenLHEPlots: 
            histos.update(OD([
            ('hGenLHE_HT_all',                            {sXaxis: HT_axis,                sXaxisLabel: r"LHE HT [GeV]"}),
            ('hGenLHE_HTIncoming_all',                    {sXaxis: HT_axis,                sXaxisLabel: r"LHE HTIncoming [GeV]"}),
            ('hGenLHE_Vpt_all',                           {sXaxis: HT_axis,                sXaxisLabel: r"LHE Vpt [GeV]"}),
            ('hGenLHE_AlphaS_all',                        {sXaxis: alphaS_axis,             sXaxisLabel: r"LHE AlphaS [GeV]"}),
            ('hGenLHE_Njets_all',                         {sXaxis: nObject50_axis,         sXaxisLabel: r"LHE Njets [GeV]"}),
            ('hGenLHE_Nb_all',                            {sXaxis: nObject50_axis,         sXaxisLabel: r"LHE Nb [GeV]"}),
            ('hGenLHE_Nc_all',                            {sXaxis: nObject50_axis,         sXaxisLabel: r"LHE Nc [GeV]"}),
            ('hGenLHE_Nuds_all',                          {sXaxis: nObject50_axis,         sXaxisLabel: r"LHE Nuds [GeV]"}),
            ('hGenLHE_Nglu_all',                          {sXaxis: nObject50_axis,         sXaxisLabel: r"LHE Nglu [GeV]"}),
            ('hGenLHE_NpNLO_all',                         {sXaxis: nObject200_axis,        sXaxisLabel: r"LHE NpNLO [GeV]"}),
            ('hGenLHE_NpLO_all',                          {sXaxis: nObject200_axis,        sXaxisLabel: r"LHE NpLO [GeV]"}),

            ('hPileup_nTrueInt',                          {sXaxis: PU_axis,                sXaxisLabel: r"Pile up"}),
            ('hPileup_nPU',                               {sXaxis: PU_axis,                sXaxisLabel: r"Pile up"}),

            ]))

        if (self.datasetInfo['isSignal'] or self.datasetInfo['isHToBB']) and runMode_SignalGenChecks:
            histos.update(OD([
                ('hGenHiggsPt_all',                           {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs (pdgId: 25, status=62))$ [GeV]"}),
                ('hGenHiggsLog2Pt_all',                       {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs (pdgId: 25, status=62))$ [GeV]"}),
                
            ]))

        if self.datasetInfo['isSignal'] and runMode_SignalGenChecks:
            histos.update(OD([
                ('hGenHiggsPt_GenHToAATo4B',                  {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(GEN Higgs (pdgId: 25, status=62))$ [GeV]"}),
                ('hGenHiggsPt_sel',                           {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(GEN Higgs (pdgId: 25, status=62))$ [GeV]"}),
                ('hGenHiggsPt_sel_wGenCuts',                  {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(GEN Higgs (pdgId: 25, status=62))$ [GeV]"}),

                ('hGenHiggsMass_all_0',                         {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN H) [GeV]"}),
                ('hMass_GenA_all_0',                            {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A) [GeV]"}),
                ('hGenHiggsMass_all',                         {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN H) [GeV]"}),
                ('hMass_GenA_all',                            {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A) [GeV]"}),
                ('hMass_GenAApair_all',                       {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN HToAA) [GeV]"}),
                ('hMass_GenAToBBbarpair_all',                 {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN AToBB) [GeV]"}),
                ('hMass_Gen4BFromHToAA_all',                  {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN HTOAATo4B) [GeV]"}),
                ('hMass_GenAToBBbarpair_all_1',               {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN AToBB) [GeV]"}),
                ('hMass_Gen4BFromHToAA_all_1',                {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN HTOAATo4B) [GeV]"}),
                ('hDeltaR_GenH_GenB_max',                     {sXaxis: deltaR_axis,     sXaxisLabel: r"$Delta$r (GEN H, GEN B)_{max}"}),



                
                # 2-D distribution
                ('hMass_GenA1_vs_GenA2_all',                       {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A1) [GeV]",
                                                                    sYaxis: mass_axis1,      sYaxisLabel: r"m (GEN A2) [GeV]"}),
                ('hMass_GenA1ToBBbar_vs_GenA2ToBBbar_all',         {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A1ToBBbar) [GeV]",
                                                                    sYaxis: mass_axis1,      sYaxisLabel: r"m (GEN A2ToBBbar) [GeV]"}),
                ('hMass_GenAHeavy_vs_GenALight_all',               {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A heavy) [GeV]",
                                                                    sYaxis: mass_axis1,      sYaxisLabel: r"m (GEN A light) [GeV]"}),
                ('hMass_GenH_vs_GenAHeavy_all',                    {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN H) [GeV]",
                                                                    sYaxis: mass_axis1,      sYaxisLabel: r"m (GEN A heavy) [GeV]"}),
                ('hMass_GenH_vs_GenALight_all',                    {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN H) [GeV]",
                                                                    sYaxis: mass_axis1,      sYaxisLabel: r"m (GEN A light) [GeV]"}),
                ('hMassGenH_vs_maxDRGenHGenB_all',                 {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN H) [GeV]",
                                                                    sYaxis: deltaR_axis,      sYaxisLabel: r"$Delta$r (GEN H, GEN B)_{max}"}),
                ('hMassGenAHeavy_vs_maxDRGenHGenB_all',            {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A heavy) [GeV]",
                                                                    sYaxis: deltaR_axis,      sYaxisLabel: r"$Delta$r (GEN H, GEN B)_{max}"}),
                ('hMassGenALight_vs_maxDRGenHGenB_all',            {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A light) [GeV]",
                                                                    sYaxis: deltaR_axis,      sYaxisLabel: r"$Delta$r (GEN H, GEN B)_{max}"}),

            ]))

        if self.datasetInfo['isQCD'] and runMode_QCDGenValidation:
            histos.update(OD([
                ('hCutFlow',                                  {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),
                ('hCutFlowWeighted',                          {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),

                ('hNEventsQCD',                               {sXaxis: cutFlow50_axis,  sXaxisLabel: 'Cuts'}),
                ('hNEventsQCDUnweighted',                     {sXaxis: cutFlow50_axis,  sXaxisLabel: 'Cuts'}),

                # QCD sample sticking            
                ('hGenLHE_HT_SelQCDbEnrich',                  {sXaxis: HT_axis,         sXaxisLabel: r"LHE HT [GeV]"}),
                ('hGenLHE_HT_SelQCDbGen',                     {sXaxis: HT_axis,         sXaxisLabel: r"LHE HT [GeV]"}),
                ('hGenLHE_HT_SelQCDbHadron',                  {sXaxis: HT_axis,         sXaxisLabel: r"LHE HT [GeV]"}),
                ('hGenLHE_HT_QCDStitchCutBQuarkPt',           {sXaxis: HT_axis,         sXaxisLabel: r"LHE HT [GeV]"}),
                ('hGenLHE_HT_QCDStitchCutBHadron',            {sXaxis: HT_axis,         sXaxisLabel: r"LHE HT [GeV]"}),
                ('hGenLHE_HT_QCDStitch',                      {sXaxis: HT_axis,         sXaxisLabel: r"LHE HT [GeV]"}),
                ('hGenLHE_HT_QCD_bEnrich_PhSp',               {sXaxis: HT_axis,         sXaxisLabel: r"LHE HT [GeV]"}),
                ('hGenLHE_HT_QCD_bGen_PhSp',                  {sXaxis: HT_axis,         sXaxisLabel: r"LHE HT [GeV]"}),
                ('hGenLHE_HT_QCD_Incl_Remnant_PhSp',          {sXaxis: HT_axis,         sXaxisLabel: r"LHE HT [GeV]"}),
                #('hGenBquark_Status_all',                     {sXaxis: PytPartStatus_axis, sXaxisLabel: r"GEN Bquark Pythia status"}),
                #('hGenBquark_first_Status_all',               {sXaxis: PytPartStatus_axis, sXaxisLabel: r"GEN first Bquark Pythia status"}),
                #('hGenBquark_first_PdgId_all',                {sXaxis: pdgId_axis,       sXaxisLabel: r"GEN first Bquark pdgId"}),

                #('hGenBquark_first_isPrompt_all',             {sXaxis: boolean_axis,    sXaxisLabel: r"GEN first Bquark isPrompt"}),

                
                ('hGenBquark_forthLeadingPt_UltraLow_all',                       {sXaxis: ptUltraLow_axis,      sXaxisLabel: r"pT (GEN B, forth leading pT) [GeV]"}),
                
                ('hGenBquark_leadingPt_all',                                     {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, leading pT) [GeV]"}),
                ('hGenBquark_subleadingPt_all',                                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, subleading pT) [GeV]"}),
                ('hGenBquark_thirdLeadingPt_all',                                {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, third leading pT) [GeV]"}),
                ('hGenBquark_forthLeadingPt_all',                                {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, forth leading pT) [GeV]"}),
                ('hGenBquark_leadingPt_QCDStitchCutBQuarkPt',                    {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, leading pT) [GeV]"}),
                ('hGenBquark_subleadingPt_QCDStitchCutBQuarkPt',                 {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, subleading pT) [GeV]"}),
                ('hGenBquark_thirdLeadingPt_QCDStitchCutBQuarkPt',               {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, third leading pT) [GeV]"}),
                ('hGenBquark_forthLeadingPt_QCDStitchCutBQuarkPt',               {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, forth leading pT) [GeV]"}),
                ('hGenBquark_leadingPt_QCDStitchCutBHadron',                     {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, leading pT) [GeV]"}),
                ('hGenBquark_subleadingPt_QCDStitchCutBHadron',                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, subleading pT) [GeV]"}),
                ('hGenBquark_thirdLeadingPt_QCDStitchCutBHadron',                {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, third leading pT) [GeV]"}),
                ('hGenBquark_forthLeadingPt_QCDStitchCutBHadron',                {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B, forth leading pT) [GeV]"}),

                # LeadingPt
                ('hLeadingPtGenBquark_pt_all',                                   {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, leading pT) [GeV]"}),
                ('hLeadingPtGenBquark_eta_all',                                  {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, leading pT)"}),
                ('hLeadingPtGenBquarkHardSctred_pt_all',                         {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, leading pT) [GeV]"}),
                ('hLeadingPtGenBquarkHardSctred_eta_all',                        {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, leading pT)"}),
                ('hLeadingPtGenBHadron_pt_all',                                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, leading pT) [GeV]"}),
                ('hLeadingPtGenBHadron_eta_all',                                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, leading pT)"}),
                ('hLeadingPtGenBHadronStatus2_pt_all',                           {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, leading pT) [GeV]"}),
                ('hLeadingPtGenBHadronStatus2_eta_all',                          {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, leading pT)"}),

                ('hLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt',                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, leading pT) [GeV]"}),
                ('hLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt',                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, leading pT)"}),
                ('hLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt',        {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, leading pT) [GeV]"}),
                ('hLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt',       {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, leading pT)"}),
                ('hLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt',                 {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, leading pT) [GeV]"}),
                ('hLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt',                {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, leading pT)"}),
                ('hLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt',          {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, leading pT) [GeV]"}),
                ('hLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt',         {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, leading pT)"}),

                ('hLeadingPtGenBquark_pt_QCDStitchCutBHadron',                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, leading pT) [GeV]"}),
                ('hLeadingPtGenBquark_eta_QCDStitchCutBHadron',                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, leading pT)"}),
                ('hLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron',        {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, leading pT) [GeV]"}),
                ('hLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron',       {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, leading pT)"}),
                ('hLeadingPtGenBHadron_pt_QCDStitchCutBHadron',                 {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, leading pT) [GeV]"}),
                ('hLeadingPtGenBHadron_eta_QCDStitchCutBHadron',                {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, leading pT)"}),
                ('hLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron',          {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, leading pT) [GeV]"}),
                ('hLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron',         {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, leading pT)"}),


                # SubleadingPt
                ('hSubleadingPtGenBquark_pt_all',                                   {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBquark_eta_all',                                  {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, sub-leading pT)"}),
                ('hSubleadingPtGenBquarkHardSctred_pt_all',                         {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBquarkHardSctred_eta_all',                        {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, sub-leading pT)"}),
                ('hSubleadingPtGenBHadron_pt_all',                                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBHadron_eta_all',                                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, sub-leading pT)"}),
                ('hSubleadingPtGenBHadronStatus2_pt_all',                           {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBHadronStatus2_eta_all',                          {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, sub-leading pT)"}),

                ('hSubleadingPtGenBquark_pt_QCDStitchCutBQuarkPt',                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBquark_eta_QCDStitchCutBQuarkPt',                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, sub-leading pT)"}),
                ('hSubleadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt',        {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt',       {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, sub-leading pT)"}),
                ('hSubleadingPtGenBHadron_pt_QCDStitchCutBQuarkPt',                 {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBHadron_eta_QCDStitchCutBQuarkPt',                {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, sub-leading pT)"}),
                ('hSubleadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt',          {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt',         {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, sub-leading pT)"}),

                ('hSubleadingPtGenBquark_pt_QCDStitchCutBHadron',                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBquark_eta_QCDStitchCutBHadron',                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, sub-leading pT)"}),
                ('hSubleadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron',        {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron',       {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, sub-leading pT)"}),
                ('hSubleadingPtGenBHadron_pt_QCDStitchCutBHadron',                 {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBHadron_eta_QCDStitchCutBHadron',                {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, sub-leading pT)"}),
                ('hSubleadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron',          {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, sub-leading pT) [GeV]"}),
                ('hSubleadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron',         {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, sub-leading pT)"}),


                # Third-LeadingPt
                ('hThirdLeadingPtGenBquark_pt_all',                                   {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBquark_eta_all',                                  {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, third-leading pT)"}),
                ('hThirdLeadingPtGenBquarkHardSctred_pt_all',                         {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBquarkHardSctred_eta_all',                        {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, third-leading pT)"}),
                ('hThirdLeadingPtGenBHadron_pt_all',                                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBHadron_eta_all',                                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, third-leading pT)"}),
                ('hThirdLeadingPtGenBHadronStatus2_pt_all',                           {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBHadronStatus2_eta_all',                          {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, third-leading pT)"}),

                ('hThirdLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt',                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt',                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, third-leading pT)"}),
                ('hThirdLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt',        {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt',       {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, third-leading pT)"}),
                ('hThirdLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt',                 {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt',                {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, third-leading pT)"}),
                ('hThirdLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt',          {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt',         {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, third-leading pT)"}),

                ('hThirdLeadingPtGenBquark_pt_QCDStitchCutBHadron',                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBquark_eta_QCDStitchCutBHadron',                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, third-leading pT)"}),
                ('hThirdLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron',        {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron',       {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, third-leading pT)"}),
                ('hThirdLeadingPtGenBHadron_pt_QCDStitchCutBHadron',                 {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBHadron_eta_QCDStitchCutBHadron',                {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, third-leading pT)"}),
                ('hThirdLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron',          {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, third-leading pT) [GeV]"}),
                ('hThirdLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron',         {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, third-leading pT)"}),

                
                # Fourth-LeadingPt
                ('hFourthLeadingPtGenBquark_pt_all',                                   {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBquark_eta_all',                                  {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, fourth-leading pT)"}),
                ('hFourthLeadingPtGenBquarkHardSctred_pt_all',                         {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBquarkHardSctred_eta_all',                        {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, fourth-leading pT)"}),
                ('hFourthLeadingPtGenBHadron_pt_all',                                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBHadron_eta_all',                                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, fourth-leading pT)"}),
                ('hFourthLeadingPtGenBHadronStatus2_pt_all',                           {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBHadronStatus2_eta_all',                          {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, fourth-leading pT)"}),

                ('hFourthLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt',                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt',                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, fourth-leading pT)"}),
                ('hFourthLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt',        {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt',       {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, fourth-leading pT)"}),
                ('hFourthLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt',                 {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt',                {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, fourth-leading pT)"}),
                ('hFourthLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt',          {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt',         {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, fourth-leading pT)"}),

                ('hFourthLeadingPtGenBquark_pt_QCDStitchCutBHadron',                  {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBquark_eta_QCDStitchCutBHadron',                 {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark, fourth-leading pT)"}),
                ('hFourthLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron',        {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B quark from hard subprocess, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron',       {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B quark from hard subprocess, fourth-leading pT)"}),
                ('hFourthLeadingPtGenBHadron_pt_QCDStitchCutBHadron',                 {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBHadron_eta_QCDStitchCutBHadron',                {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron, fourth-leading pT)"}),
                ('hFourthLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron',          {sXaxis: ptLow_axis,      sXaxisLabel: r"pT (GEN B hadron w/ status=2, fourth-leading pT) [GeV]"}),
                ('hFourthLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron',         {sXaxis: eta_axis,        sXaxisLabel: r"eta (GEN B hadron w/ status=2, fourth-leading pT)"}),

            ]))

        

        # RECO-level histograms --------------------------------------------------------------------------------------------------------------
        if self.datasetInfo['isSignal']:
            histos.update(OD([
                ('hIdxFatJetMatchedToGenBFromHToAATo4B',                   {sXaxis: nObject_axis,    sXaxisLabel: r"IdxFatJetMatchedToGenBFromHToAATo4B"}),
                ('hIdxFatJetMaxPNetMD_Hto4b_Haa4bOverQCD',                 {sXaxis: nObject_axis,    sXaxisLabel: r"IdxFatJetMaxPNetMD_Hto4b_Haa4bOverQCD"}),
                ('hIdxFatJetMaxPNetMD_Hto4b_Haa4bOverQCD_1',               {sXaxis: nObject_axis,    sXaxisLabel: r"IdxFatJetMaxPNetMD_Hto4b_Haa4bOverQCD"}),
                ('hIdxFatJetMaxZHbb_plus_Xbb',                             {sXaxis: nObject_axis,    sXaxisLabel: r"IdxFatJetMaxZHbb_plus_Xbb"}),
                ('hIdxFatJetMaxZHbb_plus_Xbb_1',                           {sXaxis: nObject_axis,    sXaxisLabel: r"IdxFatJetMaxZHbb_plus_Xbb_1"}),
                ('hLeadingBtagFatJetPtOverLeadingFatJetPt_Sig',            {sXaxis: Ratio_axis,      sXaxisLabel: r"LeadingBtagFatJetPtOverLeadingFatJetPt_Sig"}),
                    
            ]))
        
        for sel_name in self.sel_names_all.keys(): # loop of list of selections

            # for QCD, make histograms in category of number of GEN b quarks matching to leading fat jet (AK8) 
            for sHExt_0 in self.histosExtensions:    
                sHExt = "_%s" % (sel_name)
                if sHExt_0 != '':
                    sHExt += "_%s" % (sHExt_0)

                histos.update(OD([

                    ('hCutFlow'+sHExt,                                  {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),
                    ('hCutFlowWeighted'+sHExt,                          {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),

                    ('hPV_npvs'+sHExt,                               {sXaxis: PU_axis,                sXaxisLabel: r"No. of primary vertices - signal region"}),
                    ('hPV_npvsGood'+sHExt,                           {sXaxis: PU_axis,                sXaxisLabel: r"No. of good primary vertices - signal region"}),

                    ('nSelFatJet'+sHExt,                                {sXaxis: nObject_axis,    sXaxisLabel: 'No. of selected FatJets'}),
                    ('hLeadingFatJetPt'+sHExt,                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                    ('hLeadingFatJetEta'+sHExt,                         {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                    ('hLeadingFatJetPhi'+sHExt,                         {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
                    ('hLeadingFatJetMass'+sHExt,                        {sXaxis: mass_axis,       sXaxisLabel: r"m (leading FatJet) [GeV]"}),
                    ('hLeadingFatJetMSoftDrop'+sHExt,                   {sXaxis: mass_axis,       sXaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                    ('hLeadingFatJetId'+sHExt,                          {sXaxis: nObject_axis,    sXaxisLabel: r"jet Id (leading FatJet)"}),
                    ('hLeadingBtagFatJetPtOverLeadingFatJetPt'+sHExt,   {sXaxis: Ratio_axis,      sXaxisLabel: r"LeadingBtagFatJetPtOverLeadingFatJetPt"}),
                 
                    ('hLeadingFatJetBtagDeepB'+sHExt,                   {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagDeepB"}),
                    ('hLeadingFatJetBtagDDBvLV2'+sHExt,                 {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagDDBvLV2"}),

                    ('hLeadingFatJetBtagDDCvBV2'+sHExt,                 {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagDDCvBV2"}),
                    ('hLeadingFatJetBtagHbb'+sHExt,                     {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagHbb"}),
                    ('hLeadingFatJetDeepTagMD_H4qvsQCD'+sHExt,          {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD"}),
                    ('hLeadingFatJetDeepTagMD_HbbvsQCD'+sHExt,          {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD"}),
                    ('hLeadingFatJetDeepTagMD_ZHbbvsQCD'+sHExt,         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD"}),
                    ('hLeadingFatJetDeepTagMD_ZHccvsQCD'+sHExt,         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHccvsQCD"}),
                    
                    ('hLeadingFatJetDeepTagMD_ZbbvsQCD'+sHExt,          {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetdeepTagMD_ZbbvsQCD"}),
                    ('hLeadingFatJetDeepTagMD_ZvsQCD'+sHExt,            {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZvsQCD"}),
                    ('hLeadingFatJetDeepTagMD_bbvsLight'+sHExt,         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight"}),
                    ('hLeadingFatJetDeepTagMD_ccvsLight'+sHExt,         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ccvsLight"}),
                    ('hLeadingFatJetDeepTag_H'+sHExt,                   {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTag_H"}),
                    ('hLeadingFatJetDeepTag_QCD'+sHExt,                 {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTag_QCD"}),
                    ('hLeadingFatJetDeepTag_QCDothers'+sHExt,           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTag_QCDothers"}),
                    
                    ('hLeadingFatJetN2b1'+sHExt,                        {sXaxis: jetN2_axis,      sXaxisLabel: r"LeadingFatJetn2b1"}),
                    ('hLeadingFatJetN3b1'+sHExt,                        {sXaxis: jetN3_axis,      sXaxisLabel: r"LeadingFatJetn3b1"}),
                    ('hLeadingFatJetTau1'+sHExt,                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau1"}),
                    ('hLeadingFatJetTau2'+sHExt,                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau2"}),
                    ('hLeadingFatJetTau3'+sHExt,                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau3"}),
                    ('hLeadingFatJetTau4'+sHExt,                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau4"}),

                    ('hLeadingFatJetTau4by3'+sHExt,                     {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau4by3"}),
                    ('hLeadingFatJetTau3by2'+sHExt,                     {sXaxis: jetTau_axis,     sXaxisLabel: r"hLeadingFatJetTau3by2"}),
                    ('hLeadingFatJetTau2by1'+sHExt,                     {sXaxis: jetTau_axis,     sXaxisLabel: r"hLeadingFatJetTau2by1"}),
                    
                    ('hLeadingFatJetNBHadrons'+sHExt,                   {sXaxis: nObject_axis,    sXaxisLabel: r"LeadingFatJetNBHadrons"}),
                    ('hLeadingFatJetNCHadrons'+sHExt,                   {sXaxis: nObject_axis,    sXaxisLabel: r"LeadingFatJetNCHadrons"}),            
                    ('hLeadingFatJetNConstituents'+sHExt,               {sXaxis: nObject50_axis,  sXaxisLabel: r"LeadingFatJetNConstituents"}),
                    ('hLeadingFatJetNBHadronsFromHToAA'+sHExt,                   {sXaxis: nObject_axis,    sXaxisLabel: r"LeadingFatJetNBHadronsFromHToAA"}),
                    ('hLeadingFatJetNBHadrons_Sig'+sHExt,                   {sXaxis: nObject_axis,    sXaxisLabel: r"LeadingFatJetNBHadrons Sig"}),

                    ('hLeadingFatJetParticleNetMD_QCD'+sHExt,           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_QCD"}),
                    ('hLeadingFatJetParticleNetMD_Xbb'+sHExt,           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_Xbb"}),
                    ('hLeadingFatJetParticleNetMD_Xcc'+sHExt,           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_Xcc"}),
                    ('hLeadingFatJetParticleNetMD_Xqq'+sHExt,           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_Xqq"}),

                    ('hLeadingFatJetParticleNetMD_XbbOverQCD'+sHExt,    {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD Xbb/(Xbb + QCD)"}),
                    ('hLeadingFatJetParticleNetMD_XccOverQCD'+sHExt,    {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD Xbb/(Xcc + QCD)"}),
                    ('hLeadingFatJetParticleNetMD_XqqOverQCD'+sHExt,    {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD Xbb/(Xqq + QCD)"}),

                    ('hLeadingFatJetParticleNet_H4qvsQCD'+sHExt,        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_H4qvsQCD"}),
                    ('hLeadingFatJetParticleNet_HbbvsQCD'+sHExt,        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_HbbvsQCD"}),
                    ('hLeadingFatJetParticleNet_HccvsQCD'+sHExt,        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_HccvsQCD"}),
                    ('hLeadingFatJetParticleNet_QCD'+sHExt,             {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_QCD"}),
                    
                    ('hLeadingFatJetParticleNet_mass'+sHExt,            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_mass"}),

                    ('hLeadingFatJetZHbb_plus_Xbb'+sHExt,               {sXaxis: mlScore_axis2k,  sXaxisLabel: r"LeadingFatJetZHbb_plus_Xbb"}),
                    
                    # ParticleNetMD HToAATo4B
                    ('hLeadingFatJetParticleNetMD_Hto4b_Haa01b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Haa01b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_Haa2b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Haa2b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_Haa3b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Haa3b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_Haa4b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Haa4b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_QCD0b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_QCD0b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_QCD1b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_QCD1b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_QCD2b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_QCD2b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_QCD3b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_QCD3b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_QCD4b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_binaryLF_Haa4b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binaryLF_Haa4b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_binaryLF_QCDlf'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binaryLF_QCDlf"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_binary_Haa4b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binary_Haa4b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_binary_QCD'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_Haa34b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Haa34b"}),
                    
                    ('hLeadingFatJetParticleNetMD_Hto4b_binary_Haa4b_avg'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binary_Haa4b_avg"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_Haa4b_avg'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Haa4b_avg"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_QCD01234b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_QCD01234b"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_binary_QCD_avg'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binary_QCD_avg"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_QCD_avg'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_QCD_avg"}),
                    
                    ('hLeadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Htoaa4bOverQCD"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_Htoaa34bOverQCD'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Htoaa34bOverQCD"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_binaryLF_Htoaa4bOverQCD'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binaryLF_Htoaa4bOverQCD"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_binary_Htoaa4bOverQCD'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binary_Htoaa4bOverQCD"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_binary_Htoaa4bOverQCD_avg'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binary_Htoaa4bOverQCD_avg"}),
                    ('hLeadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_avg'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Htoaa4bOverQCD_avg"}),
                    
                    ('hLeadingFatJetParticleNet_massA_Hto4b_v0'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"LeadingFatJetParticleNet_massA_Hto4b_v0"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_v1'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"LeadingFatJetParticleNet_massA_Hto4b_v1"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_v2'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"LeadingFatJetParticleNet_massA_Hto4b_v2"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_v3'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"LeadingFatJetParticleNet_massA_Hto4b_v3"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_v4'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"LeadingFatJetParticleNet_massA_Hto4b_v4"}),
                    
                    ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v01'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v01"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v02'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v02"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v03'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v03"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v12'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v12"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v13'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v13"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v23'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v23"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v012'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v012"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v013'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v013"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v023'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v023"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v123'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v123"}),
                    ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v0123'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v0123"}),                                        

                    ('hLeadingFatJetParticleNet_massH_Hto4b_v0'+sHExt,            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v0"}),
                    ('hLeadingFatJetParticleNet_massH_Hto4b_v00'+sHExt,           {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v00"}),
                    ('hLeadingFatJetParticleNet_massH_Hto4b_v1'+sHExt,            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v1"}),
                    ('hLeadingFatJetParticleNet_massH_Hto4b_v2'+sHExt,            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v2"}),
                    ('hLeadingFatJetParticleNet_massH_Hto4b_v3'+sHExt,            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v3"}),
                    ('hLeadingFatJetParticleNet_massH_Hto4b_v4'+sHExt,            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v4"}),

                    ('hLeadingFatJetParticleNet_massH_Hto4b_avg_v0123'+sHExt,     {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_avg_v0123"}),
                    
                    #(''+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b"}),
                    
                    #('hLeadingFatJetParticleNetMD_Hto4b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b"}),                    
                    #('hLeadingFatJetParticleNetMD_Hto4b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b"}),
                    #('hLeadingFatJetParticleNetMD_Hto4b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b"}),
                    #('hLeadingFatJetParticleNetMD_Hto4b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b"}),
                    #('hLeadingFatJetParticleNetMD_Hto4b'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b"}),




                    ## SubJet corresponding to leadingFatJet
                    ('hLeadingFatJet_nSubJets'+sHExt,                   {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of subjets in leadingFatJet "}),
                    ('hLeadingFatJet_nSubJets_bTag_L'+sHExt,            {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of subjets (loose bTag WP) in leadingFatJet "}),
                    ('hLeadingFatJet_nSubJets_bTag_M'+sHExt,            {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of subjets (medium bTag WP) in leadingFatJet "}),


                    ## SV
                    ('hLeadingFatJet_nSV'+sHExt,                        {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of secondary vertices within leadingFatJet "}),
                    ('hLeadingFatJet_mass_SV_MaxdxySig'+sHExt,          {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]"}),
                    ('hLeadingFatJet_logMass_SV_MaxdxySig'+sHExt,       {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)"}),


                    ## MET
                    ('hMET_pT'+sHExt,                                   {sXaxis: pt_axis,         sXaxisLabel: r"MET pT [GeV]"}),
                    ('hMET_sumEt'+sHExt,                                {sXaxis: pt4TeV_axis,     sXaxisLabel: r"MET sumEt [GeV]"}),


                    ## nLeptons_matched_leadingFatJet
                    ('hLeadingFatJet_nLeptons'+sHExt,                   {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of iso-leptons within leadingFatJet "}),
                    ('hnLeptonsTight'+sHExt,                            {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of tight leptons "}),
                    ('hnLeptons_nonoverlap_leadingFatJet'+sHExt,        {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of tight leptons nonoverlaping leadingFatJet"}),
                    ('hnLeptons_nonoverlap_selFatJets'+sHExt,        {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of tight leptons nonoverlaping selected FatJets"}),
                    

                    ## AK4 jets
                    ('hnAK4Jets_NonoverlapLeadingFatJet'+sHExt,                     {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of AK4 jets non-overlap FatJet H->4b "}),
                    ('hPtLeadingAK4Jets_NonoverlapLeadingFatJet'+sHExt,             {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading AK4 jets non-overlap FatJet H->4b) [GeV]"}),
                    ('hnAK4Jets_bTag_NonoverlapLeadingFatJet'+sHExt,                {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of bTag AK4 jets non-overlap FatJet H->4b "}),
                    ('hPtLeadingAK4Jets_bTag_NonoverlapLeadingFatJet'+sHExt,        {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading bTag AK4 jets non-overlap FatJet H->4b) [GeV]"}),

                    ('hnAK4JetsCentral_NonoverlapLeadingFatJet'+sHExt,                     {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of central AK4 jets non-overlap FatJet H->4b "}),
                    ('hPtLeadingAK4JetsCentral_NonoverlapLeadingFatJet'+sHExt,             {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading central AK4 jets non-overlap FatJet H->4b) [GeV]"}),
                    ('hnAK4JetsCentral_bTag_NonoverlapLeadingFatJet'+sHExt,                {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of bTag central AK4 jets non-overlap FatJet H->4b "}),
                    ('hPtLeadingAK4JetsCentral_bTag_NonoverlapLeadingFatJet'+sHExt,        {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading bTag central AK4 jets non-overlap FatJet H->4b) [GeV]"}),

                    ('hnAK4Jets_NonoverlapSelFatJets'+sHExt,                     {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of AK4 jets non-overlap FatJet H->4b "}),
                    ('hPtLeadingAK4Jets_NonoverlapSelFatJets'+sHExt,             {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading AK4 jets non-overlap FatJet H->4b) [GeV]"}),
                    ('hnAK4Jets_bTag_NonoverlapSelFatJets'+sHExt,                {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of bTag AK4 jets non-overlap FatJet H->4b "}),
                    ('hPtLeadingAK4Jets_bTag_NonoverlapSelFatJets'+sHExt,        {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading bTag AK4 jets non-overlap FatJet H->4b) [GeV]"}),

                    ('hnAK4JetsCentral_NonoverlapSelFatJets'+sHExt,                     {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of central AK4 jets non-overlap FatJet H->4b "}),
                    ('hPtLeadingAK4JetsCentral_NonoverlapSelFatJets'+sHExt,             {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading central AK4 jets non-overlap FatJet H->4b) [GeV]"}),
                    ('hnAK4JetsCentral_bTag_NonoverlapSelFatJets'+sHExt,                {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of bTag central AK4 jets non-overlap FatJet H->4b "}),
                    ('hPtLeadingAK4JetsCentral_bTag_NonoverlapSelFatJets'+sHExt,        {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading bTag central AK4 jets non-overlap FatJet H->4b) [GeV]"}),


                    ## leadingNonHto4bFatJet 
                    ('hLeadingNonHto4bFatJetPt'+sHExt,                    {sXaxis: pt_axis,           sXaxisLabel: r"$p_{T}(leading NonHto4bFatJet)$ [GeV]"}),
                    ('hLeadingNonHto4bFatJetEta'+sHExt,                   {sXaxis: eta_axis,          sXaxisLabel: r"\eta (leading NonHto4bFatJet)"}),
                    ('hLeadingNonHto4bFatJetPhi'+sHExt,                   {sXaxis: phi_axis,          sXaxisLabel: r"\phi (leading NonHto4bFatJet)"}),
                    ('hLeadingNonHto4bFatJetMass'+sHExt,                  {sXaxis: mass_axis,         sXaxisLabel: r"m (leading NonHto4bFatJet) [GeV]"}),
                    ('hLeadingNonHto4bFatJetMSoftDrop'+sHExt,             {sXaxis: mass_axis,         sXaxisLabel: r"m_{soft drop} (leading NonHto4bFatJet) [GeV]"}),
                    ('hLeadingNonHto4bFatJetId'+sHExt,                    {sXaxis: nObject_axis,      sXaxisLabel: r"jet Id (leading NonHto4bFatJet)"}),                     
                    ('hdPhi_LeadingFJ_LeadingNonHto4bFJ'+sHExt,           {sXaxis: deltaPhi_axis,     sXaxisLabel: r"hdPhi_LeadingFJ_LeadingNonHto4bFJ"}),
                    ('hLeadingNonHto4bFatJetDeepTagMD_WvsQCD'+sHExt,      {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingFatJetDeepTagMD_WvsQCD"}),
                    ('hLeadingNonHto4bFatJetDeepTagMD_ZvsQCD'+sHExt,      {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZvsQCD"}),
                    ('hLeadingNonHto4bFatJetDeepTag_WvsQCD'+sHExt,        {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingFatJetDeepTag_WvsQCD"}),
                    ('hLeadingNonHto4bFatJetDeepTag_ZvsQCD'+sHExt,        {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingFatJetDeepTag_ZvsQCD"}),
                    ('hLeadingNonHto4bFatJetDeepTag_VvsQCD_max'+sHExt,    {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_DeepTag_VvsQCD_max"}),
                    ('hLeadingNonHto4bFatJetDeepTag_W'+sHExt,             {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_DeepTag_W"}),
                    ('hLeadingNonHto4bFatJetDeepTag_Z'+sHExt,             {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_DeepTag_Z"}),
                    ('hLeadingNonHto4bFatJetDeepTag_V_max'+sHExt,         {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_DeepTag_V_max"}),
                    ('hLeadingNonHto4bFatJetDeepTag_WZvsQCD'+sHExt,       {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_DeepTag_WZvsQCD"}),
                    ('hLeadingNonHto4bFatJetDeepTag_WZvsQCD2'+sHExt,      {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_DeepTag_WZvsQCD2"}),
                    ('hLeadingNonHto4bFatJetParticleNet_WvsQCD'+sHExt,    {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJetParticleNet_WvsQCD"}),
                    ('hLeadingNonHto4bFatJetParticleNet_ZvsQCD'+sHExt,    {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJetParticleNet_ZvsQCD"}),
                    ('hLeadingNonHto4bFatJetPNet_VvsQCD_max'+sHExt,       {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_PNet_VvsQCD_max"}),
                    ('hLeadingNonHto4bFatJetPNet_W'+sHExt,                {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_PNet_W"}),
                    ('hLeadingNonHto4bFatJetPNet_Z'+sHExt,                {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_PNet_Z"}),
                    ('hLeadingNonHto4bFatJetPNet_V_max'+sHExt,            {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_PNet_V_max"}),
                    ('hLeadingNonHto4bFatJetPNet_WZvsQCD'+sHExt,          {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_PNet_WZvsQCD"}),
                    ('hLeadingNonHto4bFatJetPNet_WZvsQCD2'+sHExt,         {sXaxis: mlScore_axis1k,    sXaxisLabel: r"LeadingNonHto4bFatJet_PNet_WZvsQCD2"}),
                    
                    
                ]))

                ### 2-D distribution --------------------------------------------------------------------------------------------------------
                histos.update(OD([
                    ('hLeadingFatJetEta_vs_Phi'+sHExt,             
                     {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)",
                      sYaxis: phi_axis,        sYaxisLabel: r"\phi (leading FatJet)"}),                    
    
                    ('hLeadingFatJetParticleNet_massH_Hto4b_avg_vs_massA_Hto4b_avg'+sHExt,     
                     {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_avg",
                      sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg"}),     

                    ('hLeadingFatJetMass_vs_massA_Hto4b_avg'+sHExt,     
                     {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMass",
                      sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg"}),

                    ('hLeadingFatJetMSoftDrop_vs_massA_Hto4b_avg'+sHExt,     
                     {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMSoftDrop",
                      sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg"}),                                                           
            ]))


                if runMode_2018HEM1516IssueValidation:
                    histos.update(OD([
                        ('hLeadingFatJetPt_DataPreHEM1516Issue'+sHExt,                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta_DataPreHEM1516Issue'+sHExt,                         {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi_DataPreHEM1516Issue'+sHExt,                         {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
                        ('hLeadingFatJetPt_DataWithHEM1516Issue'+sHExt,                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta_DataWithHEM1516Issue'+sHExt,                         {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi_DataWithHEM1516Issue'+sHExt,                         {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),

                        ('hLeadingFatJetPt_HEM1516IssueEtaPhiCut'+sHExt,    {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta_HEM1516IssuePhiCut'+sHExt,      {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi_HEM1516IssueEtaCut'+sHExt,      {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
                        ('hLeadingFatJetPt_HEM1516IssueEtaPhiCut_DataPreHEM1516Issue'+sHExt,    {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta_HEM1516IssuePhiCut_DataPreHEM1516Issue'+sHExt,      {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi_HEM1516IssueEtaCut_DataPreHEM1516Issue'+sHExt,      {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
                        ('hLeadingFatJetPt_HEM1516IssueEtaPhiCut_DataWithHEM1516Issue'+sHExt,    {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta_HEM1516IssuePhiCut_DataWithHEM1516Issue'+sHExt,      {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi_HEM1516IssueEtaCut_DataWithHEM1516Issue'+sHExt,      {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),

                        ('hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix'+sHExt,    {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix'+sHExt,      {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix'+sHExt,      {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
                        ('hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix_DataPreHEM1516Issue'+sHExt,    {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix_DataPreHEM1516Issue'+sHExt,      {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix_DataPreHEM1516Issue'+sHExt,      {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
                        ('hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix_DataWithHEM1516Issue'+sHExt,    {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix_DataWithHEM1516Issue'+sHExt,      {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix_DataWithHEM1516Issue'+sHExt,      {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),

                        # 2018 HEM15/16 issue validation
                        ('hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt'+sHExt,    {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt'+sHExt,      {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt'+sHExt,      {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
                        ('hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt_DataPreHEM1516Issue'+sHExt,    {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt_DataPreHEM1516Issue'+sHExt,      {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt_DataPreHEM1516Issue'+sHExt,      {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
                        ('hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt_DataWithHEM1516Issue'+sHExt,    {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt_DataWithHEM1516Issue'+sHExt,      {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt_DataWithHEM1516Issue'+sHExt,      {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
                    ]))


                if runMode_SignificancsScan2D:
                    histos.update(OD([

                        ## 2-D hLeadingFatJetDeepTagMD_H4qvsQCD
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetDeepTagMD_HbbvsQCD'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: mlScore_axis1,   sYaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetDeepTagMD_ZHbbvsQCD'+sHExt,           
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: mlScore_axis1,   sYaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetDeepTagMD_bbvsLight'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: mlScore_axis1,   sYaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_nSubJets'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (loose bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (medium bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_nSV'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of secondary vertices within leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: mass10_axis,     sYaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: logMass3_axis,  sYaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_MET_pT'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_H4qvsQCD_vs_MET_sumEt'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJetDeepTagMD_HbbvsQCD
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetDeepTagMD_ZHbbvsQCD'+sHExt,           
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: mlScore_axis1,   sYaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetDeepTagMD_bbvsLight'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: mlScore_axis1,   sYaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_nSubJets'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (loose bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (medium bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_nSV'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of secondary vertices within leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: mass10_axis,     sYaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: logMass3_axis,   sYaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_MET_pT'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_HbbvsQCD_vs_MET_sumEt'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJetDeepTagMD_ZHbbvsQCD
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetDeepTagMD_bbvsLight'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: mlScore_axis1,   sYaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_nSubJets'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (loose bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (medium bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_nSV'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of secondary vertices within leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: mass10_axis,     sYaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: logMass3_axis,  sYaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_MET_pT'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_MET_sumEt'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJetDeepTagMD_bbvsLight
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_nSubJets'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (loose bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (medium bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_nSV'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of secondary vertices within leadingFatJet "}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: mass10_axis,     sYaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: logMass3_axis,  sYaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_MET_pT'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJetDeepTagMD_bbvsLight_vs_MET_sumEt'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJetParticleNetMD_XbbOverQCD
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetParticleNetMD_XqqOverQCD'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: mlScore_axis1,   sYaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetParticleNet_mass'+sHExt,           
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: mass_axis,   sYaxisLabel: r"LLeadingFatJetParticleNet_mass"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_nSubJets'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets in leadingFatJet "}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (loose bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (medium bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_nSV'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of secondary vertices within leadingFatJet "}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: mass10_axis,     sYaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: logMass3_axis,  sYaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_MET_pT'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJetParticleNetMD_XbbOverQCD_vs_MET_sumEt'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XbbOverQCD",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJetParticleNetMD_XqqOverQCD
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetParticleNet_mass'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: mass_axis,   sYaxisLabel: r"LeadingFatJetParticleNet_mass"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_nSubJets'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets in leadingFatJet "}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (loose bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (medium bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_nSV'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of secondary vertices within leadingFatJet "}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: mass10_axis,     sYaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: logMass3_axis,  sYaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_MET_pT'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJetParticleNetMD_XqqOverQCD_vs_MET_sumEt'+sHExt,             
                        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_XqqOverQCD",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJetParticleNet_mass
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: mass_axis1,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: mass_axis1,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_nSubJets'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets in leadingFatJet "}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (loose bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of subjets (medium bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_nSV'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: nObject10_axis,  sYaxisLabel: r"No. of secondary vertices within leadingFatJet "}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: mass10_axis,     sYaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]"}),
                        ('hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: logMass3_axis,  sYaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)"}),
                        ('hLeadingFatJetParticleNet_mass_vs_MET_pT'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJetParticleNet_mass_vs_MET_sumEt'+sHExt,             
                        {sXaxis: mass_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_mass",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJet_nSubJets
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: nObject10_axis1,  sYaxisLabel: r"No. of subjets (loose bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: nObject10_axis1,  sYaxisLabel: r"No. of subjets (medium bTag WP) in leadingFatJet "}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJet_nSV'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: nObject10_axis1,  sYaxisLabel: r"No. of secondary vertices within leadingFatJet "}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: mass10_axis,     sYaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]"}),
                        ('hLeadingFatJet_nSubJets_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: logMass3_axis,   sYaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)"}),
                        ('hLeadingFatJet_nSubJets_vs_MET_pT'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJet_nSubJets_vs_MET_sumEt'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJet_nSubJets_bTag_L
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJet_nSV'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: nObject10_axis1,  sYaxisLabel: r"No. of secondary vertices within leadingFatJet "}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: mass10_axis,     sYaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: logMass3_axis,   sYaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_MET_pT'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJet_nSubJets_bTag_L_vs_MET_sumEt'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_L",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJet_nSubJets_bTag_M
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJet_nSV'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: nObject10_axis1,  sYaxisLabel: r"No. of secondary vertices within leadingFatJet "}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: mass10_axis,     sYaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: logMass3_axis,   sYaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_MET_pT'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJet_nSubJets_bTag_M_vs_MET_sumEt'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSubJets_bTag_M",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJet_nSV
                        ('hLeadingFatJet_nSV_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_nSV_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_nSV_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJet_nSV_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJet_nSV_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJet_nSV_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJet_nSV_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJet_nSV_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJet_nSV_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJet_nSV_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJet_nSV_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJet_nSV_vs_MET_pT'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJet_nSV_vs_MET_sumEt'+sHExt,             
                        {sXaxis: nObject10_axis,    sXaxisLabel: r"LeadingFatJet_nSV",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJet_mass_SV_MaxdxySig
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_MET_pT'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJet_mass_SV_MaxdxySig_vs_MET_sumEt'+sHExt,             
                        {sXaxis: mass10_axis,    sXaxisLabel: r"Mass of secondary vertices within leadingFatJet w/ max. dxySig [GeV]",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hLeadingFatJet_logMass_SV_MaxdxySig
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_MET_pT'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: pt_axis,         sYaxisLabel: r"MET pT [GeV]"}),
                        ('hLeadingFatJet_logMass_SV_MaxdxySig_vs_MET_sumEt'+sHExt,             
                        {sXaxis: logMass3_axis,  sXaxisLabel: r"log(Mass of secondary vertices within leadingFatJet w/ max. dxySig)",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hMET_pT
                        ('hMET_pT_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hMET_pT_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hMET_pT_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hMET_pT_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hMET_pT_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hMET_pT_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hMET_pT_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hMET_pT_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hMET_pT_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hMET_pT_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hMET_pT_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),
                        ('hMET_pT_vs_MET_sumEt'+sHExt,             
                        {sXaxis: pt_axis,    sXaxisLabel: r"MET_pT",
                        sYaxis: pt4TeV_axis,     sYaxisLabel: r"MET sumEt [GeV]"}),


                        ## 2-D hMET_sumEt
                        ('hMET_sumEt_vs_LeadingFatJetMass'+sHExt,             
                        {sXaxis: pt4TeV_axis,    sXaxisLabel: r"MET_sumEt",
                        sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hMET_sumEt_vs_LeadingFatJetMSoftDrop'+sHExt,             
                        {sXaxis: pt4TeV_axis,    sXaxisLabel: r"MET_sumEt",
                        sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hMET_sumEt_vs_LeadingFatJetN2b1'+sHExt,             
                        {sXaxis: pt4TeV_axis,    sXaxisLabel: r"MET_sumEt",
                        sYaxis: jetN2_axis,      sYaxisLabel: r"LeadingFatJetn2b1"}),
                        ('hMET_sumEt_vs_LeadingFatJetN3b1'+sHExt,             
                        {sXaxis: pt4TeV_axis,    sXaxisLabel: r"MET_sumEt",
                        sYaxis: jetN3_axis,      sYaxisLabel: r"LeadingFatJetn3b1"}),
                        ('hMET_sumEt_vs_LeadingFatJetTau1'+sHExt,             
                        {sXaxis: pt4TeV_axis,    sXaxisLabel: r"MET_sumEt",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau1"}),
                        ('hMET_sumEt_vs_LeadingFatJetTau2'+sHExt,             
                        {sXaxis: pt4TeV_axis,    sXaxisLabel: r"MET_sumEt",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau2"}),
                        ('hMET_sumEt_vs_LeadingFatJetTau3'+sHExt,             
                        {sXaxis: pt4TeV_axis,    sXaxisLabel: r"MET_sumEt",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau3"}),
                        ('hMET_sumEt_vs_LeadingFatJetTau4'+sHExt,             
                        {sXaxis: pt4TeV_axis,    sXaxisLabel: r"MET_sumEt",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4"}),
                        ('hMET_sumEt_vs_LeadingFatJetTau4by3'+sHExt,             
                        {sXaxis: pt4TeV_axis,    sXaxisLabel: r"MET_sumEt",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"LeadingFatJetTau4by3"}),
                        ('hMET_sumEt_vs_LeadingFatJetTau3by2'+sHExt,             
                        {sXaxis: pt4TeV_axis,    sXaxisLabel: r"MET_sumEt",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau3by2"}),
                        ('hMET_sumEt_vs_LeadingFatJetTau2by1'+sHExt,             
                        {sXaxis: pt4TeV_axis,    sXaxisLabel: r"MET_sumEt",
                        sYaxis: jetTau_axis,     sYaxisLabel: r"hLeadingFatJetTau2by1"}),

                        ## 2-D FatJetPt vs ML scores
                        ('hLeadingFatJetPt_vs_PNetMD_Hto4b_Haa3b'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis1k,  sYaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Haa3b"}),
                        ('hLeadingFatJetPt_vs_PNetMD_Hto4b_Haa4b'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis1k,  sYaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Haa4b"}),
                        ('hLeadingFatJetPt_vs_PNetMD_Hto4b_binaryLF_Haa4b'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis1k,  sYaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binaryLF_Haa4b"}),
                        ('hLeadingFatJetPt_vs_PNetMD_Hto4b_binary_Haa4b'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis1k,  sYaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binary_Haa4b"}),
                        ('hLeadingFatJetPt_vs_PNetMD_Hto4b_Haa34b'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis1k,  sYaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Haa34b"}),
                        ('hLeadingFatJetPt_vs_PNetMD_Hto4b_binary_Haa4b_avg'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis1k,  sYaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binary_Haa4b_avg"}),
                        ('hLeadingFatJetPt_vs_PNetMD_Hto4b_Haa4b_avg'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis1k,  sYaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Haa4b_avg"}),
                        ('hLeadingFatJetPt_vs_PNetMD_Hto4b_Htoaa4bOverQCD'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis1k,  sYaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Htoaa4bOverQCD"}),
                        ('hLeadingFatJetPt_vs_PNetMD_Hto4b_Htoaa34bOverQCD'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis1k,  sYaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Htoaa34bOverQCD"}),
                        ('hLeadingFatJetPt_vs_PNetMD_Hto4b_binary_Htoaa4bOverQCD_avg'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis1k,  sYaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binary_Htoaa4bOverQCD_avg"}),
                        ('hLeadingFatJetPt_vs_PNetMD_Hto4b_Htoaa4bOverQCD_avg'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis1k,  sYaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Htoaa4bOverQCD_avg"}),
                        ('hLeadingFatJetPt_vs_ZHbb_plus_Xbb'+sHExt,             
                        {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                        sYaxis: mlScore_axis2k,  sYaxisLabel: r"LeadingFatJetParticleNetMD ZHbb_plus_Xbb"}),

                    ]))
            
            
        for statusFlag_ in GENPART_STATUSFLAGS_LIST:
            histos['hGenBquark_first_%s_all' % (statusFlag_)] = {sXaxis: boolean_axis,    sXaxisLabel: r"GEN first Bquark %s"  % (statusFlag_)}
        
        self._accumulator = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int)
        })
        
        for histName, histAttributes in histos.items():
            #hXaxis = histAttributes[sXaxis].copy()
            hXaxis = deepcopy(histAttributes[sXaxis])
            hXaxis.label = histAttributes[sXaxisLabel]

            if sYaxis not in histAttributes.keys():
                # TH1
                self._accumulator.add({
                    histName: hist.Hist(
                        "Counts",
                        dataset_axis,
                        hXaxis, #nObject_axis,
                        systematic_axis,
                    )
                })
            else:
                # TH2
                hYaxis = deepcopy(histAttributes[sYaxis])
                hYaxis.label = histAttributes[sYaxisLabel]
                
                self._accumulator.add({
                    histName: hist.Hist(
                        "Counts",
                        dataset_axis,
                        hXaxis, #nObject_axis,
                        hYaxis,
                        systematic_axis,
                    )
                })

            
        print(f"HToAATo4bProcessor::__init__():: END", flush=flushStdout)

        

                

        '''
        self._accumulator = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int),
            #'hnFatJet_level0': hist.Hist.new.Reg(20, -0.5, 19.5, name="No. of FatJet").Weight(),
            #'hnFatJet_level1': hist.Hist.new.Reg(20, -0.5, 19.5, name="No. of FatJet").Weight(),
            'nSelFatJet': hist.Hist(
                "Counts",
                dataset_axis,
                nObject_axis,
                systematic_axis,
            ),
            'hLeadingFatJetPt': hist.Hist(
                "Counts",
                dataset_axis,
                pt_axis,
                systematic_axis,
            ),
            
        })

        self._accumulator.add({
             'hLeadingFatJetDeepTagMD_HbbvsQCD_1': hist.Hist(
                "Counts",
                dataset_axis,
                mlScore_axis,
                systematic_axis,
            ),
        })
        '''


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata["dataset"] # dataset label
        print(f"process():: {self.datasetInfo['sample_category'] = }, {dataset = }", flush=flushStdout)

        if printLevel >= 20:
            print(f"nEvents: {len(events)}")
        if printLevel >= 1:
            print(f"\n events.fields ({type(events.fields)}): {events.fields}"); sys.stdout.flush()
            #print(f"\n events.GenPart.fields: {events.GenPart.fields}")
            #print(f"\n events.L1.fields: {events.L1.fields}")
            #printVariable('\n events.L1.SingleJet180', events.L1.SingleJet180)
            #print(f"{len(events) = },  {ak.sum(events.L1.SingleJet180) = }")
            #print(f"\n events.HLT.fields: {events.HLT.fields}")
            #printVariable('\n events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4', events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4)
            #print(f"{len(events) = },  {ak.sum(events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4) = }")
            #print(f"\n events.FatJet.fields: {events.FatJet.fields}")
            #print(f"\n events.LHE.fields: {events.LHE.fields}")
            #print(f"\n events.LHE.HT: {events.LHE.HT.to_list()}")
            #print(f"{events.LHE.fields = } ")
            #print(f"{events.LHE = } ")
            #printVariable('events.LHE', events.LHE)
            #print(f"{events.run.fields = }")
            #printVariable('events.run', events.run)
            #print(f"{events.luminosityBlock.fields = }")
            #printVariable('events.luminosityBlock', events.luminosityBlock)
            #print(f"\n events.Pileup.fields: {events.Pileup.fields}")
            #print(f"\n events.Pileup.nTrueInt: {events.Pileup.nTrueInt}")
            #print(f"\n events.Pileup.nPU: {events.Pileup.nPU}")
            #print(f"\n events.PV.fields: {events.PV.fields}")
            #print(f"\n events.PV.npvs: {events.PV.npvs}")
            #print(f"\n events.PV.npvsGood: {events.PV.npvsGood}")
            #print(f"\n {events.Flag.fields = } ")
            #print(f"{events.Flag.goodVertices = }, \n{events.Flag.globalSuperTightHalo2016Filter = }, \n{events.Flag.HBHENoiseFilter = }, \n{events.Flag.HBHENoiseIsoFilter = }, \n{events.Flag.EcalDeadCellTriggerPrimitiveFilter = }, \n{events.Flag.BadPFMuonFilter = }, \n{events.Flag.BadPFMuonDzFilter = },")
            #print(f"{events.Flag.hfNoisyHitsFilter = }, \n{events.Flag.eeBadScFilter = }, \n{events.Flag.ecalBadCalibFilter = }, ")

            #print(f"\n events.FatJet.fields: {events.FatJet.fields}")
            #print(f"\n events.FatJet.pt: {events.FatJet.pt}")
            #print(f"\n events.FatJet.deepTagMD_bbvsLight: {events.FatJet.deepTagMD_bbvsLight}")
            #print(f"\n events.FatJet.particleNetMD_Xbb: {events.FatJet.particleNetMD_Xbb}")
            #print(f"\n events.FatJet.btagDeepB: {events.FatJet.btagDeepB}")

            #print(f"{events.SV.fields = }")
            #printVariable('\n ak.count(events.SV.x, axis=1)', ak.count(events.SV.x, axis=1))
            #printVariable('\n events.SV', events.SV)

            #printVariable('\n events.FatJet.pt', events.FatJet.pt)
            #printVariable('\n events.FatJet', events.FatJet)

            #printVariable('\n events.AssociatedSV', events.AssociatedSV)

            #printVariable('\n events.SV.p4', events.SV.p4)

            print(f"{events.FatJet.fields = } ")
            #print(f"{events.SV.fields = } ")
            #print(f"{events.FatJetSVs_sVIdx = } ")

            #printVariable('\n events.MET.pt[:10]', events.MET.pt[:10])
            #printVariable('\n events.MET.sumEt[:10]', events.MET.sumEt[:10])

            #print(f"{events.FatJet.fields = }")
            #printVariable('\n events.FatJet.subJetIdx1', events.FatJet.subJetIdx1)
            #printVariable('\n events.FatJet.subJetIdx2', events.FatJet.subJetIdx2)

            #printVariable('\n events.FatJet.subJetIdxG', events.FatJet.subJetIdxG)

            #printVariable('\n events.FatJet.subJetIdx1G', events.FatJet.subJetIdx1G)
            #printVariable('\n events.FatJet.subJetIdx2G', events.FatJet.subJetIdx2G)

             
        if nEventsToAnalyze != -1:
            #print(f"\n (run:ls:event): {ak.zip([events.run, events.luminosityBlock, events.event])}") 
            printVariable('\n (run:ls:event): ', ak.zip([events.run, events.luminosityBlock, events.event])); #sys.stdout.flush()           

        #if not self.datasetInfo['isMC']:
        if self.datasetInfo['isMC']:
            output = self.accumulator.identity()
            systematics_shift = [None]
            for _syst in systematics_shift:
                output += self.process_shift(events, _syst)
        else:
            print(f" {np.unique(events.run, return_counts=True) = } "); sys.stdout.flush()
            output = self.process_shift(events, None)


        return output


    
    def process_shift(self, events, shift_syst=None):
        
        output = self.accumulator.identity()
        dataset = events.metadata["dataset"] # dataset label
        #print(f"process_shift():: {shift_syst = } dataset: {dataset}", flush=flushStdout)

        

        ones_list  = np.ones(len(events))
        trues_list = np.ones(len(events), dtype=bool)
        falses_list = np.full(len(events), False)


          
        ##################
        # OBJECT SELECTION
        ##################


        # Gen-level selection ---------------------------------------------------------------------
        genHiggs = None
        genHT    = None
        idx_GenB_fromHToAA = None
        mask_SignalHToAATo4B_Boosted = None
        if self.datasetInfo['isSignal'] or self.datasetInfo['isHToBB']: 
            genHiggs  = self.objectSelector.selectGenHiggs(events)

        if self.datasetInfo['isSignal']: 
            genHT     = self.objectSelector.GenHT(events)

            # m(bbar from A) and m(4b from HToAA) ----------
            '''
            genHiggsCollection = events.GenPart[(
                (events.GenPart.pdgId  == 25) & # pdgId:: 25: H0
                (events.GenPart.status == 62)   # statu 62: outgoing subprocess particle with primordial kT included https://pythia.org/latest-manual/ParticleProperties.html
            )]
            '''

            genACollection = self.objectSelector.selectGenABoson(events)
            genA_First  = genACollection[:, 0]
            genA_Second = genACollection[:, 1]
            
            idxGenA_sortByMass = ak.argsort(genACollection.mass, axis=-1, ascending=False)
            # genACollection[idxGenA_sortByMass[0]] : Leading mass GenA
            # genACollection[idxGenA_sortByMass[1]] : Subleading mass GenA
            
            
            genBBar_pairs_all = ak.argcombinations(events.GenPart, 2, fields=['b', 'bbar'])
            genBBar_pairs = genBBar_pairs_all[(
                (abs(events.GenPart[genBBar_pairs_all['b'   ]].pdgId) == 5) &
                (abs(events.GenPart[genBBar_pairs_all['bbar']].pdgId) == 5) &
                ((events.GenPart[genBBar_pairs_all['b']].pdgId) == (-1*events.GenPart[genBBar_pairs_all['bbar']].pdgId)  ) &
                (events.GenPart[genBBar_pairs_all['b']].genPartIdxMother == events.GenPart[genBBar_pairs_all['bbar']].genPartIdxMother) &
                (events.GenPart[ events.GenPart[genBBar_pairs_all['b'   ]].genPartIdxMother ].pdgId == 36) &
                (events.GenPart[ events.GenPart[genBBar_pairs_all['bbar']].genPartIdxMother ].pdgId == 36) 
            )]
            idx_GenB_fromHToAA = ak.concatenate([genBBar_pairs['b'], genBBar_pairs['bbar']], axis=-1)

            # LorentVector of GenB quarks from HToAATo4b
            nEvents_11 = ak.num(events.GenPart[genBBar_pairs['b']][:, 0].pt, axis=0)
            #mass_bQuark = 4.18
            #print(f"\n np.full(nEvents_11, mass_bQuark): {np.full(nEvents_11, mass_bQuark)}")

            # https://coffeateam.github.io/coffea/modules/coffea.nanoevents.methods.vector.html
            # bQuark from 1st A
            LVGenB_0 = ak.zip(
                {
                    "pt"  : events.GenPart[genBBar_pairs['b']][:, 0].pt,
                    "eta" : events.GenPart[genBBar_pairs['b']][:, 0].eta,
                    "phi" : events.GenPart[genBBar_pairs['b']][:, 0].phi,
                    "mass": np.full(nEvents_11, MASS_BottomQuark),
                },
                with_name="PtEtaPhiMLorentzVector",
                behavior=vector.behavior,
            )

            # bbarQuark from 1st A
            LVGenBbar_0 = ak.zip(
                {
                    "pt"  : events.GenPart[genBBar_pairs['bbar']][:, 0].pt,
                    "eta" : events.GenPart[genBBar_pairs['bbar']][:, 0].eta,
                    "phi" : events.GenPart[genBBar_pairs['bbar']][:, 0].phi,
                    "mass": np.full(nEvents_11, MASS_BottomQuark),
                },
                with_name="PtEtaPhiMLorentzVector",
                behavior=vector.behavior,
            )

            # bQuark from 2nd A
            LVGenB_1 = ak.zip(
                {
                    "pt"  : events.GenPart[genBBar_pairs['b']][:, 1].pt,
                    "eta" : events.GenPart[genBBar_pairs['b']][:, 1].eta,
                    "phi" : events.GenPart[genBBar_pairs['b']][:, 1].phi,
                    "mass": np.full(nEvents_11, MASS_BottomQuark),
                },
                with_name="PtEtaPhiMLorentzVector",
                behavior=vector.behavior,
            )

            # bbarQuark from 2nd A
            LVGenBbar_1 = ak.zip(
                {
                    "pt"  : events.GenPart[genBBar_pairs['bbar']][:, 1].pt,
                    "eta" : events.GenPart[genBBar_pairs['bbar']][:, 1].eta,
                    "phi" : events.GenPart[genBBar_pairs['bbar']][:, 1].phi,
                    "mass": np.full(nEvents_11, MASS_BottomQuark),
                },
                with_name="PtEtaPhiMLorentzVector",
                behavior=vector.behavior,
            )            

            dr_GenH_GenB = ak.concatenate([genHiggs.delta_r(LVGenB_0), genHiggs.delta_r(LVGenBbar_0), genHiggs.delta_r(LVGenB_1), genHiggs.delta_r(LVGenBbar_1)], axis=-1)
            max_dr_GenH_GenB = ak.max(dr_GenH_GenB, axis=-1)    
            mask_SignalHToAATo4B_Boosted = (max_dr_GenH_GenB < 0.8)

            if printLevel >= 10:
                printVariable("\n genBBar_pairs['b']", genBBar_pairs['b'])  
                printVariable("\n [genBBar_pairs['bbar']", genBBar_pairs['bbar'])
                printVariable("\n events.GenPart[genBBar_pairs['b']]", events.GenPart[genBBar_pairs['b']])
                printVariable("\n events.GenPart[genBBar_pairs['bbar']]", events.GenPart[genBBar_pairs['bbar']])
                printVariable("\n ak.concatenate([genBBar_pairs['b'], genBBar_pairs['bbar']], axis=-1)", ak.concatenate([genBBar_pairs['b'], genBBar_pairs['bbar']], axis=-1))
                printVariable("\n events.GenPart[ ak.concatenate([genBBar_pairs['b'], genBBar_pairs['bbar']], axis=-1)]", events.GenPart[ ak.concatenate([genBBar_pairs['b'], genBBar_pairs['bbar']], axis=-1)] )
                printVariable("\n genB_fromHToAA", genB_fromHToAA)

            
                
        # QCD MC ----------------------------------------------
        mask_genBHadrons_status2_eventwise                            = None
        mask_genBQuarks_hardSctred_eventwise                          = None
        mask_genBHadrons_status2_and_noGenBQuarksHardSctred_eventwise = None
        mask_QCD_stitch_CutBHadron_eventwise                          = None
        mask_QCD_stitch_CutBQuarkPt_eventwise                         = None
        mask_QCD_stitch_eventwise                                     = None
        mask_QCD_bEnrich_PhSp                                         = None
        mask_QCD_bGen_PhSp                                            = None
        mask_QCD_Incl_Remnant_PhSp                                    = None
        mask_genBQuarksHardSctred_genBHadronsStatus2                  = None
        vGenBQuarksHardSctred_genBHadronsStatus2_sel                  = None
        if self.datasetInfo['isMC'] and self.datasetInfo['isQCD'] :
            mask_genLHEHTLt100 = (events.LHE.HT < 100)

            if printLevel >= 12:
                printVariable('\n events.LHE.HT', events.LHE.HT); sys.stdout.flush()
                printVariable('\n mask_genLHEHTLt100', mask_genLHEHTLt100); sys.stdout.flush()
            
            
            genBQuarks = events.GenPart[(
                (abs(events.GenPart.pdgId) == PDGID_BottomQuark )
            )]
            genBQuarks_pT = ak.sort(genBQuarks.pt, axis=-1, ascending=False)
            genBQuarks_first = ak.firsts(genBQuarks)
            #genBQuarks_first = ak.firsts(genBQuarks, axis=-1)
            mask_genBQuarks = (ak.count(genBQuarks.pdgId, axis=1) >= 1)

            mask_genBQuarks_pTAbvTrsh = ak.any((genBQuarks.pt > 15.0), axis=1)

            idx_genBQuarks_pTsort = ak.argsort(genBQuarks.pt, axis=-1, ascending=False)
            
            
            if printLevel >= 12:
                printVariable('\n genBQuarks', genBQuarks); sys.stdout.flush()
                printVariable(' genBQuarks.pt', genBQuarks.pt); sys.stdout.flush()
                printVariable(' (genBQuarks.pt > 15.0)', (genBQuarks.pt > 15.0)); sys.stdout.flush()
                printVariable(' mask_genBQuarks_pTAbvTrsh', mask_genBQuarks_pTAbvTrsh); sys.stdout.flush()

                printVariable('\n idx_genBQuarks_pTsort', idx_genBQuarks_pTsort); sys.stdout.flush()
                printVariable(' genBQuarks[idx_genBQuarks_pTsort]', genBQuarks[idx_genBQuarks_pTsort]); sys.stdout.flush()
                printVariable(' genBQuarks[idx_genBQuarks_pTsort].pt', genBQuarks[idx_genBQuarks_pTsort].pt); sys.stdout.flush()
                
                printVariable('\n events.GenPart[(events.GenPart.status == 2)].pdgId', events.GenPart[(events.GenPart.status == 2)].pdgId); sys.stdout.flush()


            # Check if event has B-hadron with pythia status==2, which are QCD-BGenFilter requirement -----------------
            mask_genBHadrons_status2 = None
            for ipdgId_tmp in range(len(self.pdgId_BHadrons)):
                pdgId_tmp = self.pdgId_BHadrons[ipdgId_tmp]
                mask_genBHadrons_status2_i = (
                    (events.GenPart.status == 2) &
                    (abs(events.GenPart.pdgId) == pdgId_tmp)
                )

                if ipdgId_tmp == 0:
                    mask_genBHadrons_status2 = mask_genBHadrons_status2_i
                else:
                    mask_genBHadrons_status2 = mask_genBHadrons_status2 | mask_genBHadrons_status2_i

                if printLevel >= 100:
                    printVariable('\n\n %d \t mask_genBHadrons_status2_i' % pdgId_tmp, mask_genBHadrons_status2_i); sys.stdout.flush()
                    printVariable('\n     \t mask_genBHadrons_status2  ', mask_genBHadrons_status2); sys.stdout.flush()
                    
            mask_genBHadrons_status2_eventwise = ak.any(mask_genBHadrons_status2, axis=1)

            genBHadrons_status2 = events.GenPart[mask_genBHadrons_status2]
            idx_genBHadrons_status2_pTsort = ak.argsort(genBHadrons_status2.pt, axis=-1, ascending=False)

            
            mask_genBHadrons = None
            for ipdgId_tmp in range(len(self.pdgId_BHadrons)):
                pdgId_tmp = self.pdgId_BHadrons[ipdgId_tmp]
                mask_genBHadrons_i = (
                    (abs(events.GenPart.pdgId) == pdgId_tmp)
                )

                if ipdgId_tmp == 0:
                    mask_genBHadrons = mask_genBHadrons_i
                else:
                    mask_genBHadrons = mask_genBHadrons | mask_genBHadrons_i                   

            genBHadrons = events.GenPart[mask_genBHadrons]
            idx_genBHadrons_pTsort = ak.argsort(genBHadrons.pt, axis=-1, ascending=False)
           
                    
            if printLevel >= 12:
                #genBHadrons_status2 = events.GenPart[mask_genBHadrons_status2]
                printVariable('\n genBHadrons_status2', genBHadrons_status2); sys.stdout.flush()
                printVariable(' genBHadrons_status2.pdgId', genBHadrons_status2.pdgId); sys.stdout.flush()
                printVariable(' ak.any(mask_genBHadrons_status2, axis=1)', ak.any(mask_genBHadrons_status2, axis=1)); sys.stdout.flush()


            # Check if events has b-quark outgoing from hard subprocess -----------------------------------------------
            mask_genBQuarks_hardSctred = (
                (abs(events.GenPart.pdgId) == PDGID_BottomQuark ) &
                (events.GenPart.status == 23)
            )
            mask_genBQuarks_hardSctred_eventwise = ak.any(mask_genBQuarks_hardSctred, axis=1)

            genBQuarks_hardSctred = events.GenPart[mask_genBQuarks_hardSctred]
            idx_genBQuarks_hardSctred_pTsort = ak.argsort(genBQuarks_hardSctred.pt, axis=-1, ascending=False)
            
            if printLevel >= 12:
                #genBQuarks_hardSctred = events.GenPart[mask_genBQuarks_hardSctred]
                printVariable('\n genBQuarks_hardSctred', genBQuarks_hardSctred); sys.stdout.flush()
                printVariable(' genBQuarks_hardSctred.hasFlags("isHardProcess")', genBQuarks_hardSctred.hasFlags("isHardProcess")); sys.stdout.flush()
                printVariable(' genBQuarks_hardSctred.hasFlags("fromHardProcess")', genBQuarks_hardSctred.hasFlags("fromHardProcess")); sys.stdout.flush()
                printVariable(' ak.any(mask_genBQuarks_hardSctred, axis=1)', ak.any(mask_genBQuarks_hardSctred, axis=1)); sys.stdout.flush()

                printVariable(' genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort].pt', genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort].pt); sys.stdout.flush()

            # QCD stitch cut-based: conditions -----------------------------------------------------------------------------------
            # option 1: GEN b-quark pT > 15 GeV for QCD BEnrich and QCD bGen samples.
            if self.datasetInfo['isQCD_bEnrich'] or self.datasetInfo['isQCD_bGen']:
                mask_QCD_stitch_CutBQuarkPt_eventwise = mask_genBQuarks_pTAbvTrsh
            elif self.datasetInfo['isQCDIncl']:
                mask_QCD_stitch_CutBQuarkPt_eventwise = (
                    (mask_genBQuarks_pTAbvTrsh == False) |
                    (mask_genLHEHTLt100 == True)
                )
                
            if printLevel >= 12:
                printVariable('\n mask_QCD_stitch_CutBQuarkPt_eventwise', mask_QCD_stitch_CutBQuarkPt_eventwise); sys.stdout.flush()


            # Option 2: B-Hadron for QCD bGEN. b-quark from hard subprocess for QCD bEnrich
            if self.datasetInfo['isQCD_bEnrich']:
                mask_QCD_stitch_CutBHadron_eventwise = trues_list
            elif self.datasetInfo['isQCD_bGen']:
                mask_QCD_stitch_CutBHadron_eventwise = (
                    (mask_genBQuarks_hardSctred_eventwise == False)
                )
            elif self.datasetInfo['isQCDIncl']:
                mask_QCD_stitch_CutBHadron_eventwise = (
                    (
                        (mask_genBQuarks_hardSctred_eventwise == False) &
                        (mask_genBHadrons_status2_eventwise == False)
                    ) |
                    (
                        (mask_genLHEHTLt100 == True)
                    )
                )
                
            if printLevel >= 12:
                printVariable('\n mask_QCD_stitch_CutBHadron_eventwise', mask_QCD_stitch_CutBHadron_eventwise); sys.stdout.flush()

            if printLevel >= 10:
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 4)
                #mask_tmp1 = (genBHadrons_status2[idx_genBHadrons_status2_pTsort].pt < 0.01)
                
                #printVariable('\n (genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp)].pt)', (genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp)].pt)); sys.stdout.flush()
                printVariable('\n (genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp)][:, 3].pt)', (genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp)][:, 2].pt)); sys.stdout.flush()

                genBHadrons_status2_fourth_pt = genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp)][:, 2].pt
                mask_tmp1_ = (genBHadrons_status2_fourth_pt < 0.01)
                printVariable('\n genBHadrons_status2_fourth_pt[mask_tmp1_]', genBHadrons_status2_fourth_pt[mask_tmp1_]); sys.stdout.flush()

            if printLevel >= 13:
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 4)
                genBQuarks_fourth_pT = genBQuarks_pT[mask_tmp][:, 3]
                mask_tmp1_ = (genBQuarks_fourth_pT < 0.01)

                printVariable('\n genBQuarks_pT[mask_tmp][:, 3]', genBQuarks_pT[mask_tmp][:, 3]); sys.stdout.flush()
                printVariable('\n genBQuarks_fourth_pT[mask_tmp1_]', genBQuarks_fourth_pT[mask_tmp1_]); sys.stdout.flush()

            if printLevel >= 10:
                mask_tmp1_ = (genBQuarks.pt < 1e-3)
                mask_tmp2_ = ak.any(mask_tmp1_, axis=1)

                # genBQuarks_pT[mask_tmp][:, 3]
                #printVariable('\n ', ); sys.stdout.flush()
                #printVariable('\n aw.count(genBQuarks.pt, axis=1)', ak.count(genBQuarks.pt, axis=1)); sys.stdout.flush()
                printVariable('\n aw.count(genBQuarks.pt[mask_tmp1_], axis=1)', ak.count(genBQuarks.pt[mask_tmp1_], axis=1)); sys.stdout.flush()
                printVariable('\n genBQuarks[mask_tmp2_].pt', genBQuarks[mask_tmp2_].pt); sys.stdout.flush()
                printVariable('\n genBQuarks_pT[mask_tmp2_]', genBQuarks_pT[mask_tmp2_]); sys.stdout.flush()

            if printLevel >= 10:
                mask_tmp1_ = (genBQuarks.pt < 1e-3)
                mask_tmp2_ = ak.any(mask_tmp1_, axis=1)

                # genBQuarks_pT[mask_tmp][:, 3]
                #printVariable('\n ', ); sys.stdout.flush()
                #printVariable('\n aw.count(genBQuarks.pt, axis=1)', ak.count(genBQuarks.pt, axis=1)); sys.stdout.flush()
                printVariable('\n aw.count(genBQuarks.pt[mask_tmp1_], axis=1)', ak.count(genBQuarks.pt[mask_tmp1_], axis=1)); sys.stdout.flush()
                printVariable('\n genBQuarks[mask_tmp2_].pt', genBQuarks[mask_tmp2_].pt); sys.stdout.flush()
                #printVariable('\n genBQuarks_pT[mask_tmp2_]', genBQuarks_pT[mask_tmp2_]); sys.stdout.flush()
                printVariable('\n genBQuarks[mask_tmp2_]', genBQuarks[mask_tmp2_]); sys.stdout.flush()


            # Phase space_ QCD_bEnrich, QCD_bGen, QCD_Incl_Remnant
            mask_QCD_bEnrich_PhSp = (mask_genBQuarks_hardSctred_eventwise == True)
            mask_QCD_bGen_PhSp = (
                (mask_genBQuarks_hardSctred_eventwise == False) &
                (mask_genBHadrons_status2_eventwise   == True)
            )
            mask_QCD_Incl_Remnant_PhSp = (
                (
                    (mask_genBQuarks_hardSctred_eventwise == False) &
                    (mask_genBHadrons_status2_eventwise   == False)
                ) |
                (
                    (mask_genLHEHTLt100 == True)
                )
            )
                
            if self.datasetInfo["MCSamplesStitchOption"] == MCSamplesStitchOptions.PhSpOverlapRewgt:
                mask_QCD_stitch_eventwise = trues_list # select all events
            else:
                mask_QCD_stitch_eventwise = mask_QCD_stitch_CutBHadron_eventwise

            mask_genBHadrons_status2_and_noGenBQuarksHardSctred_eventwise = (
                (mask_genBHadrons_status2_eventwise == True) &
                (mask_genBQuarks_hardSctred_eventwise == False)
            )

            if printLevel >= 100:
                mask_phsp_1 = (
                    mask_QCD_bEnrich_PhSp | mask_QCD_bGen_PhSp | mask_QCD_Incl_Remnant_PhSp
                )
                print(f" {ak.sum(mask_QCD_bEnrich_PhSp & mask_QCD_bGen_PhSp) = }")
                print(f" {ak.sum(mask_QCD_bGen_PhSp & mask_QCD_Incl_Remnant_PhSp) = }")
                print(f" {ak.sum(mask_QCD_bEnrich_PhSp & mask_QCD_Incl_Remnant_PhSp) = }")
                print(f" {ak.sum(   (mask_QCD_bEnrich_PhSp | mask_QCD_bGen_PhSp | mask_QCD_Incl_Remnant_PhSp) ) = } ") 
                print(f" {ak.sum( ~ (mask_QCD_bEnrich_PhSp | mask_QCD_bGen_PhSp | mask_QCD_Incl_Remnant_PhSp) ) = } ")

            '''
            #  genBQuarks collection to match to reco AK8 jet at the later stage ----------------------
            mask_genBQuarksHardSctred_genBHadronsStatus2 = mask_genBQuarks_hardSctred | mask_genBHadrons_status2
            
            vGenBQuarksHardSctred_genBHadronsStatus2 = ak.zip(
                {
                    "pt"  : events.GenPart[mask_genBQuarksHardSctred_genBHadronsStatus2].pt,
                    "eta" : events.GenPart[mask_genBQuarksHardSctred_genBHadronsStatus2].eta,
                    "phi" : events.GenPart[mask_genBQuarksHardSctred_genBHadronsStatus2].phi,
                    "mass": ak.where(
                        abs(events.GenPart[mask_genBQuarksHardSctred_genBHadronsStatus2].pdgId) == 5,
                        MASS_BottomQuark,
                        events.GenPart[mask_genBQuarksHardSctred_genBHadronsStatus2].mass
                    )
                },
                with_name="PtEtaPhiMLorentzVector",
                behavior=vector.behavior,
            )            

            # From genBQuarks or genBHadrons collection, remove children of genBQuarks or genBHadrons existed in the collection
            dR_parent_child_max = 0.1            
            idx_pairs_genBQuarksHardSctred_genBHadronsStatus2 = ak.argcombinations(vGenBQuarksHardSctred_genBHadronsStatus2, 2, fields=['b1', 'b2'])
            mask_nearbyPairs_genBQuarksHardSctred_genBHadronsStatus2 = vGenBQuarksHardSctred_genBHadronsStatus2[idx_pairs_genBQuarksHardSctred_genBHadronsStatus2['b1']].delta_r(
                vGenBQuarksHardSctred_genBHadronsStatus2[idx_pairs_genBQuarksHardSctred_genBHadronsStatus2['b2']]
            ) < dR_parent_child_max
            mask_distinct_genBQuarksHardSctred_genBHadronsStatus2 = ~ akArray_isin(
                testArray = ak.local_index(vGenBQuarksHardSctred_genBHadronsStatus2), 
                referenceArray = idx_pairs_genBQuarksHardSctred_genBHadronsStatus2[mask_nearbyPairs_genBQuarksHardSctred_genBHadronsStatus2]['b2'] # b2: children 
                )
            vGenBQuarksHardSctred_genBHadronsStatus2_sel = vGenBQuarksHardSctred_genBHadronsStatus2[mask_distinct_genBQuarksHardSctred_genBHadronsStatus2]
            '''
        # --------------------------------------------------------------------------------------------------
        

        # MC ttbar ----------------------------------------------
        #mask_1 = None
        if self.datasetInfo['isMC'] and self.datasetInfo['isTTbar'] :
            mask_genTopQuark = (
                (abs(events.GenPart.pdgId) == PDGID_TopQuark ) & 
                (events.GenPart.hasFlags("isLastCopy"))
            )   
            if printLevel >= 13:
                printVariable('\n mask_genTopQuark \n', mask_genTopQuark); sys.stdout.flush()
                printVariable('\n events.GenPart[mask_genTopQuark] \n', events.GenPart[mask_genTopQuark]); sys.stdout.flush()
                printVariable('\n events.GenPart[mask_genTopQuark] (pdgId, status) \n', 
                              ak.zip([
                                events.GenPart[mask_genTopQuark].pdgId,
                                events.GenPart[mask_genTopQuark].status
                              ])
                              ); sys.stdout.flush()
                printVariable('\n events.GenPart[mask_genTopQuark] (pt, eta, phi) \n', 
                              ak.zip([
                                events.GenPart[mask_genTopQuark].pt,
                                events.GenPart[mask_genTopQuark].eta,
                                events.GenPart[mask_genTopQuark].phi,
                              ])
                              ); sys.stdout.flush()
                print(f"{GENPART_STATUSFLAGS._member_map_ = }"); sys.stdout.flush()
                for GENPART_STATUSFLAG_name, GENPART_STATUSFLAG_number in GENPART_STATUSFLAGS._member_map_.items():
                    print(f"{GENPART_STATUSFLAG_name = }, {GENPART_STATUSFLAG_number = }"); sys.stdout.flush()
                    printVariable('\n events.GenPart[mask_genTopQuark].statusFlags %s \n'%GENPART_STATUSFLAG_name, 
                                  selGenPartsWithStatusFlag(events.GenPart[mask_genTopQuark].statusFlags, GENPART_STATUSFLAG_number)); sys.stdout.flush()
                    printVariable('events.GenPart[mask_genTopQuark].statusFlags %s << \n'%GENPART_STATUSFLAG_name, 
                                  events.GenPart[mask_genTopQuark].hasFlags([GENPART_STATUSFLAG_name]), ); sys.stdout.flush()
                    

                #printVariable('\n events.GenPart[mask_genTopQuark]', events.GenPart[mask_genTopQuark].); sys.stdout.flush()
                #printVariable('\n events.GenPart[mask_genTopQuark]', events.GenPart[mask_genTopQuark]); sys.stdout.flush()

            if printLevel >= 3:
                wgt_TopPt = getTopPtRewgt(
                    eventsGenPart = events.GenPart[mask_genTopQuark],
                    isPythiaTuneCP5 = self.datasetInfo['isPythiaTuneCP5']
                    )
                printVariable('\n wgt_TopPt', wgt_TopPt)

            














        # Reco-level -----------------------------------------------------------------------------------
        # FatJet selection
        #selFatJet = self.objectSelector.selectFatJets(events)

        '''
        mask_FatJetPt = (events.FatJet.pt > self.objectSelector.FatJetPtThsh)
        mask_FatJetEta = (abs(events.FatJet.eta) < self.objectSelector.FatJetEtaThsh)
        mask_FatJetBtagDeepB = (events.FatJet.btagDeepB > bTagWPs[self.objectSelector.era][self.objectSelector.tagger_btagDeepB][self.objectSelector.wp_btagDeepB])
        mask_FatJetMSoftDrop = (
            (events.FatJet.msoftdrop > self.objectSelector.FatJetMSoftDropThshLow) &
            (events.FatJet.msoftdrop < self.objectSelector.FatJetMSoftDropThshHigh)
        )
        
        selFatJet = events.FatJet[(
            mask_FatJetPt &
            mask_FatJetEta &
            mask_FatJetBtagDeepB #&
            #mask_FatJetMSoftDrop
        )]
        '''

        

        ################## 
        # EVENT VARIABLES
        ##################
        '''
        leadingFatJet = None
        if ( not self.datasetInfo['isMC']) and (self.datasetInfo["era"] == Era_2018):
            #  HEM15/16 issue in jets with -3.2<eta<-1.3 and -1.57<phi< -0.87 in 2018 data (runs>=319077, i.e. last certified run of 2018B, and all of 2018C+D) https://twiki.cern.ch/twiki/bin/view/CMS/JetMET#
            run_FatJetEta_FatJetPhi = ak.zip({'run': events.run, 'FatJetEta': events.FatJet.eta, 'FatJetPhi': events.FatJet.phi})   
            mask_jets_surviving_HEM15_16_issue = ~ (
                (run_FatJetEta_FatJetPhi.run >= 319077) &
                (run_FatJetEta_FatJetPhi.FatJetEta > -3.2 ) & (run_FatJetEta_FatJetPhi.FatJetEta < -1.3 ) &
                (run_FatJetEta_FatJetPhi.FatJetPhi > -1.57) & (run_FatJetEta_FatJetPhi.FatJetPhi < -0.87)
            )
            leadingFatJet = ak.firsts(events.FatJet[mask_jets_surviving_HEM15_16_issue]) 
        else:   
            leadingFatJet = ak.firsts(events.FatJet) # for e.g. [0.056304931640625, None, 0.12890625, 0.939453125, 0.0316162109375]
        '''

        leadingFatJet = None

        FatJetParticleNetMD_XbbvsQCD = ak.where(
            (events.FatJet.particleNetMD_Xbb + events.FatJet.particleNetMD_QCD) > 0,
            events.FatJet.particleNetMD_Xbb / (events.FatJet.particleNetMD_Xbb + events.FatJet.particleNetMD_QCD),
            np.full_like(events.FatJet.particleNetMD_Xbb, 0)
        ) 
        FatJetDeepTagMD_ZHbbvsQCD = ak.where(
            events.FatJet.deepTagMD_ZHbbvsQCD >= 0,
            events.FatJet.deepTagMD_ZHbbvsQCD,
            np.full_like(events.FatJet.deepTagMD_ZHbbvsQCD, 0)
        )
        FatJet_ZHbb_plus_Xbb = FatJetDeepTagMD_ZHbbvsQCD + FatJetParticleNetMD_XbbvsQCD
        idx_FatJet_ZHbb_plus_Xbb_max = ak.argmax(FatJet_ZHbb_plus_Xbb, axis=-1, keepdims=True)


        # PNetMD_Hto4b
        idx_FatJet_PNetMD_Hto4b_Haa4bOverQCD_max = None
        if 'particleNetMD_Hto4b_Haa4b' in events.FatJet.fields:
            FatJet_PNetMD_Hto4b_QCD01234b_sum = (
                events.FatJet.particleNetMD_Hto4b_QCD0b + 
                events.FatJet.particleNetMD_Hto4b_QCD1b + 
                events.FatJet.particleNetMD_Hto4b_QCD2b + 
                events.FatJet.particleNetMD_Hto4b_QCD3b + 
                events.FatJet.particleNetMD_Hto4b_QCD4b )  
            FatJet_PNetMD_Hto4b_Htoaa4bOverQCD = ak.where(
                (events.FatJet.particleNetMD_Hto4b_Haa4b + FatJet_PNetMD_Hto4b_QCD01234b_sum) > 0.0,
                (
                    events.FatJet.particleNetMD_Hto4b_Haa4b / 
                    (events.FatJet.particleNetMD_Hto4b_Haa4b + FatJet_PNetMD_Hto4b_QCD01234b_sum)
                ),
                ak.full_like(events.FatJet.particleNetMD_Hto4b_Haa4b, 0) #events.FatJet.particleNetMD_Hto4b_Haa4b
            ) 
            idx_FatJet_PNetMD_Hto4b_Haa4bOverQCD_max = ak.argmax(FatJet_PNetMD_Hto4b_Htoaa4bOverQCD, axis=-1, keepdims=True)
            #leadingBtagFatJet = ak.firsts(events.FatJet[idx_FatJet_PNetMD_Hto4b_Haa4bOverQCD_max])    
            leadingFatJet = ak.firsts(events.FatJet[idx_FatJet_PNetMD_Hto4b_Haa4bOverQCD_max]) 
        else:
  
            #idx_FatJet_PNetMD_XbbvsQCD_max = ak.argmax(FatJetParticleNetMD_XbbvsQCD, axis=-1, keepdims=True)
            #leadingBtagFatJet = ak.firsts(events.FatJet[idx_FatJet_PNetMD_XbbvsQCD_max]) 
            #leadingFatJet = ak.firsts(events.FatJet[idx_FatJet_PNetMD_XbbvsQCD_max]) 
            leadingFatJet = ak.firsts(events.FatJet[idx_FatJet_ZHbb_plus_Xbb_max]) 


        if runMode_OptimizePNetTaggerCut:
            leadingFatJet = ak.firsts(events.FatJet[idx_FatJet_ZHbb_plus_Xbb_max])

        #leadingFatJet = ak.firsts(events.FatJet)        
        leadingFatJet_asSingletons = ak.singletons(leadingFatJet) # for e.g. [[0.056304931640625], [], [0.12890625], [0.939453125], [0.0316162109375]]
        
        leadingFatJetDeepTagMD_ZHbbvsQCD = ak.where(
            leadingFatJet.deepTagMD_ZHbbvsQCD >= 0,
            leadingFatJet.deepTagMD_ZHbbvsQCD,
            np.full_like(leadingFatJet.deepTagMD_ZHbbvsQCD, 0)
        ) 
        leadingFatJetDeepTagMD_ZHccvsQCD = ak.where(
            leadingFatJet.deepTagMD_ZHccvsQCD >= 0,
            leadingFatJet.deepTagMD_ZHccvsQCD,
            np.full_like(leadingFatJet.deepTagMD_ZHccvsQCD, 0)
        )        
        leadingFatJetParticleNetMD_XbbvsQCD = ak.where(
            (leadingFatJet.particleNetMD_Xbb + leadingFatJet.particleNetMD_QCD) > 0,
            leadingFatJet.particleNetMD_Xbb / (leadingFatJet.particleNetMD_Xbb + leadingFatJet.particleNetMD_QCD),
            np.full(len(events), 0)
        )        
        leadingFatJetParticleNetMD_XccvsQCD = ak.where(
            (leadingFatJet.particleNetMD_Xcc + leadingFatJet.particleNetMD_QCD) > 0,
            leadingFatJet.particleNetMD_Xcc / (leadingFatJet.particleNetMD_Xcc + leadingFatJet.particleNetMD_QCD),
            np.full(len(events), 0)
        )         
        leadingFatJetParticleNetMD_XqqvsQCD = ak.where(
            (leadingFatJet.particleNetMD_Xqq + leadingFatJet.particleNetMD_QCD) > 0,
            leadingFatJet.particleNetMD_Xqq / (leadingFatJet.particleNetMD_Xqq + leadingFatJet.particleNetMD_QCD),
            np.full(len(events), 0)
        ) 
        leadingFatJetZHbb_plus_Xbb =  leadingFatJetDeepTagMD_ZHbbvsQCD + leadingFatJetParticleNetMD_XbbvsQCD
        leadingFatJetZHbb_Xbb_avg  = (leadingFatJetDeepTagMD_ZHbbvsQCD + leadingFatJetParticleNetMD_XbbvsQCD) / 2

        # ZHbb = ZHbbvsQCD*(1 - ZHccvsQCD)/(1 - ZHbbvsQCD*ZHccvsQCD)
        leadingFatJetZHbb = leadingFatJetDeepTagMD_ZHbbvsQCD * (1 - leadingFatJetDeepTagMD_ZHccvsQCD)
        leadingFatJetZHbb = leadingFatJetZHbb / (1 - (leadingFatJetDeepTagMD_ZHbbvsQCD * leadingFatJetDeepTagMD_ZHccvsQCD))

        

        # PNetMD_Hto4b
        if 'particleNetMD_Hto4b_Haa4b' in events.FatJet.fields:
            leadingFatJet_PNetMD_Hto4b_QCD01234b_sum = (
                leadingFatJet.particleNetMD_Hto4b_QCD0b + 
                leadingFatJet.particleNetMD_Hto4b_QCD1b + 
                leadingFatJet.particleNetMD_Hto4b_QCD2b + 
                leadingFatJet.particleNetMD_Hto4b_QCD3b + 
                leadingFatJet.particleNetMD_Hto4b_QCD4b )
            leadingFatJet_PNetMD_Hto4b_QCD_avg = (
                leadingFatJet_PNetMD_Hto4b_QCD01234b_sum         + 
                leadingFatJet.particleNetMD_Hto4b_binaryLF_QCDlf + 
                leadingFatJet.particleNetMD_Hto4b_binary_QCD       ) / 3
            leadingFatJet_PNetMD_Hto4b_binary_QCD_avg = (
                leadingFatJet.particleNetMD_Hto4b_binaryLF_QCDlf + 
                leadingFatJet.particleNetMD_Hto4b_binary_QCD       ) / 2

            leadingFatJet_PNetMD_Hto4b_Haa34b_sum = (
                leadingFatJet.particleNetMD_Hto4b_Haa4b          + 
                leadingFatJet.particleNetMD_Hto4b_Haa3b          ) 
            leadingFatJet_PNetMD_Hto4b_Haa4b_avg = (
                leadingFatJet_PNetMD_Hto4b_Haa34b_sum            + 
                leadingFatJet.particleNetMD_Hto4b_binary_Haa4b   + 
                leadingFatJet.particleNetMD_Hto4b_binaryLF_Haa4b   ) / 3
            leadingFatJet_PNetMD_Hto4b_binary_Haa4b_avg = ( 
                leadingFatJet.particleNetMD_Hto4b_binary_Haa4b   + 
                leadingFatJet.particleNetMD_Hto4b_binaryLF_Haa4b   ) / 2

            leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD = ak.where(
                (leadingFatJet.particleNetMD_Hto4b_Haa4b + leadingFatJet_PNetMD_Hto4b_QCD01234b_sum) > 0.0,
                (
                    leadingFatJet.particleNetMD_Hto4b_Haa4b / 
                    (leadingFatJet.particleNetMD_Hto4b_Haa4b + leadingFatJet_PNetMD_Hto4b_QCD01234b_sum)
                ),
                ak.full_like(leadingFatJet.particleNetMD_Hto4b_Haa4b, 0) #leadingFatJet.particleNetMD_Hto4b_Haa4b
            )

            # hLeadingFatJetParticleNet_massA_Hto4b_avg_v013
            leadingFatJet_PNet_massA_Hto4b_avg = calculateAverageOfArrays([
                leadingFatJet.particleNet_massA_Hto4b_v0,
                leadingFatJet.particleNet_massA_Hto4b_v1,
                leadingFatJet.particleNet_massA_Hto4b_v3
            ])
            
            # hLeadingFatJetParticleNet_massH_Hto4b_avg_v0123
            # Scale particleNet_massH_Hto4b_v0 by 1.01 to get better response
            # https://indico.cern.ch/event/1343334/contributions/5655252/attachments/2745224/4781382/2023_11_02_HToAATo4B_Higgs_mass_studies.pdf#page=15
            leadingFatJet_PNet_massH_Hto4b_avg = calculateAverageOfArrays([
                leadingFatJet.particleNet_massH_Hto4b_v0 * 1.01,
                leadingFatJet.particleNet_massH_Hto4b_v1,
                leadingFatJet.particleNet_massH_Hto4b_v2,
                leadingFatJet.particleNet_massH_Hto4b_v3, 
                ])


        # SubJet corresponding to leading FatJet 
        leadingFatJet_subJetIdx_concatenate = ak.concatenate([
                leadingFatJet_asSingletons.subJetIdx1[ leadingFatJet_asSingletons.subJetIdx1 >= 0 ],
                leadingFatJet_asSingletons.subJetIdx2[ leadingFatJet_asSingletons.subJetIdx2 >= 0 ]
            ],
            axis=1 
        )
        leadingFatJet_subJets                    = events.SubJet[leadingFatJet_subJetIdx_concatenate]
        mask_leadingFatJet_subJets_bTagDeepCSV_L = leadingFatJet_subJets.btagDeepB > bTagWPs[self.datasetInfo["era"]]['DeepCSV']['L']
        mask_leadingFatJet_subJets_bTagDeepCSV_M = leadingFatJet_subJets.btagDeepB > bTagWPs[self.datasetInfo["era"]]['DeepCSV']['M']
        leadingFatJet_nSubJets                   = ak.count(leadingFatJet_subJets.btagDeepB, axis=1)
        leadingFatJet_nSubJets_bTag_L            = ak.count(leadingFatJet_subJets[mask_leadingFatJet_subJets_bTagDeepCSV_L].btagDeepB, axis=1)
        leadingFatJet_nSubJets_bTag_M            = ak.count(leadingFatJet_subJets[mask_leadingFatJet_subJets_bTagDeepCSV_M].btagDeepB, axis=1)
        

        # mask satisfying HEM1516 issues conditions
        scaleAK4ToAK8 = 0.4
        mask_HEM1516Issue = ak.fill_none((
            (leadingFatJet.eta > (-3.2  - scaleAK4ToAK8)) & (leadingFatJet.eta < (-1.3  + scaleAK4ToAK8)) & 
            (leadingFatJet.phi > (-1.57 - scaleAK4ToAK8)) & (leadingFatJet.phi < (-0.87 + scaleAK4ToAK8))
        ), False)
        mask_HEM1516Issue_Eta = ak.fill_none((
            (leadingFatJet.eta > (-3.2 - scaleAK4ToAK8)) & (leadingFatJet.eta < (-1.3 + scaleAK4ToAK8))
        ), False)
        mask_HEM1516Issue_Phi = ak.fill_none((
            (leadingFatJet.phi > (-1.57 - scaleAK4ToAK8)) & (leadingFatJet.phi < (-0.87 + scaleAK4ToAK8))
        ), False)
        isRunAffectedBy2018HEM1516Issue = (
            (events.run >= HEM1516Issue2018_AffectedRunRange[0])  & 
            (events.run <= HEM1516Issue2018_AffectedRunRange[1])
        )

        ## SVs matched with leadingFatJet
        mask_SV_matched_leadingFatJet           = leadingFatJet.delta_r(events.SV.p4) < 0.8
        SV_matched_leadingFatJet                = events.SV[mask_SV_matched_leadingFatJet]
        nSV_matched_leadingFatJet               = ak.fill_none( ak.count( SV_matched_leadingFatJet.chi2, axis=1 ), 0)
        idx_SV_matched_leadingFatJet_MaxdxySig  = ak.argmax( SV_matched_leadingFatJet.dxySig, axis=-1, keepdims=True )
        SV_matched_leadingFatJet_MaxdxySig      = ak.firsts( SV_matched_leadingFatJet[ idx_SV_matched_leadingFatJet_MaxdxySig ] )
        SV_matched_leadingFatJet_MaxdxySig      = ak.firsts( SV_matched_leadingFatJet[ idx_SV_matched_leadingFatJet_MaxdxySig ] )
        mass_SV_matched_leadingFatJet_MaxdxySig = ak.fill_none(SV_matched_leadingFatJet_MaxdxySig.mass, 1e-10)
        
        ## match leadingFat jet to genB 
        n_leadingFatJat_matched_genB = np.full(len(events), 0)
        n_leadingFatJat_matched_genB_HToAATo4B = np.full(len(events), 0)
        if self.datasetInfo['isMC'] :
            #mask_leadingFatJat_matched_genB = leadingFatJet.delta_r(vGenBQuarksHardSctred_genBHadronsStatus2_sel) < 0.8
            #n_leadingFatJat_matched_genB = ak.sum(mask_leadingFatJat_matched_genB, axis=1)
            n_leadingFatJat_matched_genB = leadingFatJet.nBHadrons


            if self.datasetInfo['isSignal']:                
                n_leadingFatJat_matched_genB_HToAATo4B = ak.sum(leadingFatJet.delta_r( events.GenPart[ idx_GenB_fromHToAA ] ) < 0.8, axis=1)
                
                mask_FatJet_matched_genB_HToAATo4B = (
                    (events.FatJet.delta_r(LVGenB_0)    < 0.8) &
                    (events.FatJet.delta_r(LVGenBbar_0) < 0.8) &
                    (events.FatJet.delta_r(LVGenB_1)    < 0.8) &
                    (events.FatJet.delta_r(LVGenBbar_1) < 0.8)                     
                )
                mask_events_FatJet_matched_genB_HToAATo4B = ak.any(mask_FatJet_matched_genB_HToAATo4B, axis=1) # events with FatJet matches to 4 GEN B-quarks from HToAATo4B
                idx_FatJet_matched_genB_HToAATo4B = ak.argmax(mask_FatJet_matched_genB_HToAATo4B, axis=-1) # index of FatJet within events that maches to 4 GEN B-quarks from HToAATo4B
                
                       

        ## non-HTo4B FatJet
        nonHto4bFatJet = events.FatJet[(events.FatJet.delta_r(leadingFatJet) > 0.8)]
        leadingNonHto4bFatJet = ak.firsts(nonHto4bFatJet)
        leadingNonHto4bFatJet_asSingletons = ak.singletons(leadingNonHto4bFatJet) # for e.g. [[0.056304931640625], [], [0.12890625], [0.939453125], [0.0316162109375]]
        
        # Calculate W, Z, (W+Z)vsQCD scores from WvsQCD, ZvsQCD and QCD scores
        # Formulas from Andrew on Baylor slack: https://baylorhep.slack.com/archives/C013B0LRAEA/p1706815879028809
        def cal_W_(WQ, Q):
            return WQ * Q / (1 - WQ)
        def cal_Z_(ZQ, Q):
            return ZQ * Q / (1 - ZQ)
        def cal_WZQ_1_(W, Z, Q):
            return (W + Z) / (W + Z + Q)
        def cal_WZQ_2_(WQ, ZQ):
            return (WQ + ZQ + - (2 * WQ * ZQ)) / (1 - (WQ * ZQ))        
        kSmallPositiveNumber = 0.000000001
        leadingNonHto4bFatJet_PNet_WvsQCD        = array_PutLowerBound(leadingNonHto4bFatJet.particleNet_WvsQCD, kSmallPositiveNumber)
        leadingNonHto4bFatJet_PNet_ZvsQCD        = array_PutLowerBound(leadingNonHto4bFatJet.particleNet_ZvsQCD, kSmallPositiveNumber)
        leadingNonHto4bFatJet_PNet_QCD           = array_PutLowerBound(leadingNonHto4bFatJet.particleNet_QCD,    kSmallPositiveNumber)
        leadingNonHto4bFatJet_PNet_W             = cal_W_(leadingNonHto4bFatJet_PNet_WvsQCD, leadingNonHto4bFatJet_PNet_QCD)
        leadingNonHto4bFatJet_PNet_Z             = cal_Z_(leadingNonHto4bFatJet_PNet_ZvsQCD, leadingNonHto4bFatJet_PNet_QCD)
        leadingNonHto4bFatJet_PNet_WZvsQCD       = cal_WZQ_1_(leadingNonHto4bFatJet_PNet_W, leadingNonHto4bFatJet_PNet_Z, leadingNonHto4bFatJet_PNet_QCD)
        leadingNonHto4bFatJet_PNet_WZvsQCD2      = cal_WZQ_2_(leadingNonHto4bFatJet_PNet_WvsQCD, leadingNonHto4bFatJet_PNet_ZvsQCD) 
        leadingNonHto4bFatJet_PNet_VvsQCD_max    = calculateMaxOfTwoArrays(leadingNonHto4bFatJet_PNet_WvsQCD, leadingNonHto4bFatJet_PNet_ZvsQCD)
        leadingNonHto4bFatJet_PNet_V_max         = calculateMaxOfTwoArrays(leadingNonHto4bFatJet_PNet_W, leadingNonHto4bFatJet_PNet_Z)
        #
        leadingNonHto4bFatJet_DeepTag_WvsQCD     = array_PutLowerBound(leadingNonHto4bFatJet.deepTag_WvsQCD, kSmallPositiveNumber)
        leadingNonHto4bFatJet_DeepTag_ZvsQCD     = array_PutLowerBound(leadingNonHto4bFatJet.deepTag_ZvsQCD, kSmallPositiveNumber)
        leadingNonHto4bFatJet_DeepTag_QCD        = array_PutLowerBound(leadingNonHto4bFatJet.deepTag_QCD,    kSmallPositiveNumber)
        leadingNonHto4bFatJet_DeepTag_W          = cal_W_(leadingNonHto4bFatJet_DeepTag_WvsQCD, leadingNonHto4bFatJet_DeepTag_QCD)
        leadingNonHto4bFatJet_DeepTag_Z          = cal_Z_(leadingNonHto4bFatJet_DeepTag_ZvsQCD, leadingNonHto4bFatJet_DeepTag_QCD)
        leadingNonHto4bFatJet_DeepTag_WZvsQCD    = cal_WZQ_1_(leadingNonHto4bFatJet_DeepTag_W, leadingNonHto4bFatJet_DeepTag_Z, leadingNonHto4bFatJet_DeepTag_QCD)
        leadingNonHto4bFatJet_DeepTag_WZvsQCD2   = cal_WZQ_2_(leadingNonHto4bFatJet_DeepTag_WvsQCD, leadingNonHto4bFatJet_DeepTag_ZvsQCD)
        leadingNonHto4bFatJet_DeepTag_VvsQCD_max = calculateMaxOfTwoArrays(leadingNonHto4bFatJet_DeepTag_WvsQCD, leadingNonHto4bFatJet_DeepTag_ZvsQCD)
        leadingNonHto4bFatJet_DeepTag_V_max      = calculateMaxOfTwoArrays(leadingNonHto4bFatJet_DeepTag_W, leadingNonHto4bFatJet_DeepTag_Z)



        if printLevel >= 100:
            printVariable('\n events.FatJet.delta_r(leadingFatJet): ', events.FatJet.delta_r(leadingFatJet))
            printVariable('\n events.FatJet.delta_r(leadingFatJet): ', events.FatJet.delta_r(leadingFatJet) > 0.8)
            printVariable('\n events.FatJet.pt', events.FatJet.pt)
            printVariable('\n leadingFatJet.pt', leadingFatJet.pt)
            printVariable('\n nonHto4bFatJet.pt', nonHto4bFatJet.pt)
            printVariable('\n leadingNonHto4bFatJet.pt', leadingNonHto4bFatJet.pt)
            printVariable('\n leadingNonHto4bFatJet_asSingletons.pt', leadingNonHto4bFatJet_asSingletons.pt)
            printVariable('\n leadingFatJet.delta_phi(leadingNonHto4bFatJet)', leadingFatJet.delta_phi(leadingNonHto4bFatJet))
            printVariable('\n abs(leadingFatJet.delta_phi(leadingNonHto4bFatJet))', abs(leadingFatJet.delta_phi(leadingNonHto4bFatJet)))
            printVariable('\n leadingNonHto4bFatJet.particleNet_WvsQCD', leadingNonHto4bFatJet.particleNet_WvsQCD)
            printVariable('\n leadingNonHto4bFatJet.particleNet_ZvsQCD', leadingNonHto4bFatJet.particleNet_ZvsQCD)
            printVariable('\n leadingNonHto4bFatJet.particleNet_QCD', leadingNonHto4bFatJet.particleNet_QCD)
            printVariable('\n leadingNonHto4bFatJet_PNet_WvsQCD', leadingNonHto4bFatJet_PNet_WvsQCD)
            printVariable('\n leadingNonHto4bFatJet_PNet_ZvsQCD', leadingNonHto4bFatJet_PNet_ZvsQCD)
            printVariable('\n leadingNonHto4bFatJet_PNet_QCD', leadingNonHto4bFatJet_PNet_QCD)
            printVariable('\n leadingNonHto4bFatJet_PNet_W', leadingNonHto4bFatJet_PNet_W)
            printVariable('\n leadingNonHto4bFatJet_PNet_Z', leadingNonHto4bFatJet_PNet_Z)
            printVariable('\n leadingNonHto4bFatJet_PNet_WZvsQCD', leadingNonHto4bFatJet_PNet_WZvsQCD)
            printVariable('\n leadingNonHto4bFatJet_PNet_WZvsQCD2', leadingNonHto4bFatJet_PNet_WZvsQCD2)
            printVariable('\n leadingNonHto4bFatJet_PNet_VvsQCD_max', leadingNonHto4bFatJet_PNet_VvsQCD_max)
            printVariable('\n leadingNonHto4bFatJet_PNet_V_max', leadingNonHto4bFatJet_PNet_V_max)

        
        ## sel leptons
        muonsTight     = selectMuons(
            events.Muon, 
            pT_Thsh    =self.objectSelector.MuonPtThsh, 
            MVAId      =self.objectSelector.MuonMVAId, 
            MiniIsoId  =self.objectSelector.MuonMiniIsoId, 
            MVATTHThsh =self.objectSelector.MuonMVATTHThsh )
        electronsTight = selectElectrons(
            events.Electron, 
            pT_Thsh    =self.objectSelector.ElectronPtThsh, 
            MVAId      =self.objectSelector.ElectronMVAId, 
            MVATTHThsh =self.objectSelector.ElectronMVATTHThsh )
        leptonsTight   = ak.concatenate([muonsTight, electronsTight], axis=1)
        nLeptonsTight  = ak.fill_none(ak.count(leptonsTight.pt, axis=1), 0)
        nLeptons_matched_leadingFatJet    = ak.fill_none(ak.sum(leadingFatJet.delta_r( leptonsTight ) < 0.8, axis=1), 0)
        nLeptons_nonoverlap_leadingFatJet = ak.fill_none(ak.sum(leadingFatJet.delta_r( leptonsTight ) > 0.8, axis=1), 0)
        nLeptons_nonoverlap_selFatJets = ak.fill_none(ak.sum(
            ((leadingFatJet.delta_r( leptonsTight ) > 0.8) &
             (leadingNonHto4bFatJet.delta_r( leptonsTight ) > 0.8)), 
            axis=1), 0)
                    
        

        ## sel AK4 jets
        ak4Jets = selectAK4Jets(Jets=events.Jet, era=self.datasetInfo["era"], pT_Thsh=30)
        mask_ak4Jets_nonoverlaping_leadingFatJet           = ak4Jets.delta_r(leadingFatJet) > 0.8
        ak4Jets_nonoverlaping_leadingFatJet                = ak4Jets[ mask_ak4Jets_nonoverlaping_leadingFatJet ]
        nAk4Jets_nonoverlaping_leadingFatJet               = ak.fill_none(ak.count(ak4Jets_nonoverlaping_leadingFatJet.pt, axis=1), 0)

        ak4JetsCentral_nonoverlaping_leadingFatJet         = ak4Jets_nonoverlaping_leadingFatJet[abs(ak4Jets_nonoverlaping_leadingFatJet.eta) < 2.4]
        nAk4JetsCentral_nonoverlaping_leadingFatJet        = ak.fill_none(ak.count(ak4JetsCentral_nonoverlaping_leadingFatJet.pt, axis=1), 0)

        mask_ak4Jets_bTag_nonoverlaping_leadingFatJet      = ak4Jets_nonoverlaping_leadingFatJet.btagDeepFlavB > self.objectSelector.Ak4JetDeepJetB_Thsh
        ak4Jets_bTag_nonoverlaping_leadingFatJet           = ak4Jets_nonoverlaping_leadingFatJet[mask_ak4Jets_bTag_nonoverlaping_leadingFatJet]
        nAk4Jets_bTag_nonoverlaping_leadingFatJet          = ak.fill_none(ak.count(ak4Jets_bTag_nonoverlaping_leadingFatJet.pt, axis=1), 0)

        ak4JetsCentral_bTag_nonoverlaping_leadingFatJet    = ak4Jets_bTag_nonoverlaping_leadingFatJet[abs(ak4Jets_bTag_nonoverlaping_leadingFatJet.eta) < 2.4]
        nAk4JetsCentral_bTag_nonoverlaping_leadingFatJet   = ak.fill_none(ak.count(ak4JetsCentral_bTag_nonoverlaping_leadingFatJet.pt, axis=1), 0)


        mask_ak4Jets_nonoverlaping_selFatJets           = (ak4Jets.delta_r(leadingFatJet) > 0.8) & (ak4Jets.delta_r(leadingNonHto4bFatJet) > 0.8)
        ak4Jets_nonoverlaping_selFatJets                = ak4Jets[ mask_ak4Jets_nonoverlaping_selFatJets ]
        nAk4Jets_nonoverlaping_selFatJets               = ak.fill_none(ak.count(ak4Jets_nonoverlaping_selFatJets.pt, axis=1), 0)

        ak4JetsCentral_nonoverlaping_selFatJets         = ak4Jets_nonoverlaping_selFatJets[abs(ak4Jets_nonoverlaping_selFatJets.eta) < 2.4]
        nAk4JetsCentral_nonoverlaping_selFatJets        = ak.fill_none(ak.count(ak4JetsCentral_nonoverlaping_selFatJets.pt, axis=1), 0)

        mask_ak4Jets_bTag_nonoverlaping_selFatJets      = ak4Jets_nonoverlaping_selFatJets.btagDeepFlavB > self.objectSelector.Ak4JetDeepJetB_Thsh
        ak4Jets_bTag_nonoverlaping_selFatJets           = ak4Jets_nonoverlaping_selFatJets[mask_ak4Jets_bTag_nonoverlaping_selFatJets]
        nAk4Jets_bTag_nonoverlaping_selFatJets          = ak.fill_none(ak.count(ak4Jets_bTag_nonoverlaping_selFatJets.pt, axis=1), 0)

        ak4JetsCentral_bTag_nonoverlaping_selFatJets    = ak4Jets_bTag_nonoverlaping_selFatJets[abs(ak4Jets_bTag_nonoverlaping_selFatJets.eta) < 2.4]
        nAk4JetsCentral_bTag_nonoverlaping_selFatJets   = ak.fill_none(ak.count(ak4JetsCentral_bTag_nonoverlaping_selFatJets.pt, axis=1), 0)

            



        #####################
        # EVENT SELECTION
        #####################
        
        # reconstruction level cuts for cut-flow table. Order of cuts is IMPORTANT
        cuts_reco = ["dR_LeadingFatJet_GenB_0p8"] + self.sel_names_all["Presel"] #.copy()

       
        # create a PackedSelection object
        # this will help us later in composing the boolean selections easily
        selection = PackedSelection()

        if "run:ls" in self.sel_conditions_all_list:
            # self.datasetInfo['dataLSSelGoldenJSON']
            # using Coffea built-in function: mask_lumi = LumiMask(golden_json_path)(events.run,events.luminosityBlock)
            #selection.add("run:ls", LumiMask(sFilesGoldenJSON[self.datasetInfo["era"]])(events.run,events.luminosityBlock) )
            
            # using selectRunLuminosityBlock function from htoaa_CommonTools
            selection.add("run:ls", selectRunLuminosityBlock(
                dataLSSelGoldenJSON  = self.datasetInfo['dataLSSelGoldenJSON'], 
                runNumber_list       = events.run, 
                luminosityBlock_list = events.luminosityBlock 
                ))

        if "nPV" in self.sel_conditions_all_list:
            # nPVGood >= 1
            selection.add("nPV", events.PV.npvsGood >= 1)

        if "METFilters" in self.sel_conditions_all_list:
            #mask_METFilters = selectMETFilters(events.Flag, self.datasetInfo["era"], self.datasetInfo['isMC'])
            #printVariable('\n mask_METFilters', mask_METFilters)
            selection.add(
                "METFilters", 
                selectMETFilters(events.Flag, self.datasetInfo["era"], self.datasetInfo['isMC'])
            )

        if "leadingFatJetPt" in self.sel_conditions_all_list:
            # >=1 FatJet
            #selection.add("FatJetGet", ak.num(selFatJet) >= self.objectSelector.nFatJetMin)
            selection.add(
                "leadingFatJetPt",
                leadingFatJet.pt > self.objectSelector.FatJetPtThsh
            )
        if "leadingNonHto4bFatJetPt" in self.sel_conditions_all_list:
            selection.add(
                "leadingNonHto4bFatJetPt",
                leadingNonHto4bFatJet.pt > self.objectSelector.FatJetPtThsh
            )

        if "leadingFatJetEta" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetEta",
                abs(leadingFatJet.eta) < self.objectSelector.FatJetEtaThsh
            )
        if "leadingNonHto4bFatJetEta" in self.sel_conditions_all_list:
            selection.add(
                "leadingNonHto4bFatJetEta",
                abs(leadingNonHto4bFatJet.eta) < self.objectSelector.FatJetEtaThsh
            )

        if "JetID"  in self.sel_conditions_all_list:
            selection.add(
                "JetID", 
                leadingFatJet.jetId == self.objectSelector.FatJetJetID
            )
        if "leadingNonHto4bFatJetJetID"  in self.sel_conditions_all_list:
            selection.add(
                "leadingNonHto4bFatJetJetID", 
                leadingNonHto4bFatJet.jetId == self.objectSelector.FatJetJetID
            )

 
        if "leadingFatJetMSoftDrop"  in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetMSoftDrop",
                (leadingFatJet.msoftdrop > self.objectSelector.FatJetMSoftDropThshLow) &
                (leadingFatJet.msoftdrop < self.objectSelector.FatJetMSoftDropThshHigh)
            )
        if "leadingNonHto4bFatMSoftDrop"  in self.sel_conditions_all_list:
            selection.add(
                "leadingNonHto4bFatMSoftDrop",
                (leadingNonHto4bFatJet.msoftdrop > self.objectSelector.FatJetMSoftDropThshLow) &
                (leadingNonHto4bFatJet.msoftdrop < self.objectSelector.FatJetMSoftDropThshHigh)
            )

        if "leadingFatJetParticleNetMD_XbbvsQCD" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetParticleNetMD_XbbvsQCD",
                leadingFatJetParticleNetMD_XbbvsQCD > self.objectSelector.FatJetParticleNetMD_XbbvsQCD_Thsh
            )

        if "leadingFatJetZHbb_Xbb_avg" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetZHbb_Xbb_avg",
                leadingFatJetZHbb_Xbb_avg > self.objectSelector.FatJetZHbb_Xbb_avg_Thsh
            )

        if "leadingFatJetZHbb" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetZHbb",
                leadingFatJetZHbb > self.objectSelector.FatJetZHbb_Thsh
            )       
        
        if "leadingFatJetDeepTagMD_bbvsLight" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetDeepTagMD_bbvsLight",
                leadingFatJet.deepTagMD_bbvsLight > self.objectSelector.FatJetDeepTagMD_bbvsLight_Thsh
            )

        if "leadingFatJetZHbb_plus_Xbb" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetZHbb_plus_Xbb",
                leadingFatJetZHbb_plus_Xbb > self.objectSelector.FatJetZHbb_plus_Xbb_Thsh
            )

        if "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD",
                leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD > self.objectSelector.FatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_Thsh
            )

        if "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP40" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP40",
                leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD > bTagWPs[self.datasetInfo["era"]]['ParticleNetMD_Hto4b_Htoaa4bOverQCD']['WP-40']
            )

        if "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP60" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP60",
                leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD > bTagWPs[self.datasetInfo["era"]]['ParticleNetMD_Hto4b_Htoaa4bOverQCD']['WP-60']
            )

        if "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP80" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP80",
                leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD > bTagWPs[self.datasetInfo["era"]]['ParticleNetMD_Hto4b_Htoaa4bOverQCD']['WP-80']
            )


        if "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP80to40" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP80to40",
                ( (leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD >  bTagWPs[self.datasetInfo["era"]]['ParticleNetMD_Hto4b_Htoaa4bOverQCD']['WP-80']) &
                  (leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD <= bTagWPs[self.datasetInfo["era"]]['ParticleNetMD_Hto4b_Htoaa4bOverQCD']['WP-40']) )
            )

        if "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP95to60" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP95to60",
                ( (leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD >  bTagWPs[self.datasetInfo["era"]]['ParticleNetMD_Hto4b_Htoaa4bOverQCD']['WP-95']) &
                  (leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD <= bTagWPs[self.datasetInfo["era"]]['ParticleNetMD_Hto4b_Htoaa4bOverQCD']['WP-60']) )
            )

        if "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP99to80" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_WP99to80",
                ( (leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD >  bTagWPs[self.datasetInfo["era"]]['ParticleNetMD_Hto4b_Htoaa4bOverQCD']['WP-99']) &
                  (leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD <= bTagWPs[self.datasetInfo["era"]]['ParticleNetMD_Hto4b_Htoaa4bOverQCD']['WP-80']) )
            )

        # mA windows
        for sMassAWindowName, massAWindow in massPseudoscalarA_windows_dict.items():
            sSelName_tmp_ = "leadingFatJetParticleNet_massA_Hto4b_avg_%s" % (sMassAWindowName)
            if sSelName_tmp_ not in self.sel_conditions_all_list: continue

            selection.add(
                sSelName_tmp_,
                ((leadingFatJet_PNet_massA_Hto4b_avg > massAWindow[0]) & 
                 (leadingFatJet_PNet_massA_Hto4b_avg < massAWindow[1] ))
            )

        # mH windows
        for sMassHiggsWindowName, massHiggsWindow in massHiggs_windows_dict. items():
            # msoftdrop
            if 'Msoftdrop' in sMassHiggsWindowName:
                sSelName_tmp_ = "leadingFatJetMsoftdrop_%s" % (sMassHiggsWindowName)
                if sSelName_tmp_ not in self.sel_conditions_all_list: continue

                selection.add(
                    sSelName_tmp_,
                    ((leadingFatJet.msoftdrop > massHiggsWindow[0]) & 
                     (leadingFatJet.msoftdrop < massHiggsWindow[1] ))
                )

            # ParticleNet_massH_Hto4b
            elif 'PNet_massH_Hto4b' in sMassHiggsWindowName:
                sSelName_tmp_ = "leadingFatJetPNet_massH_Hto4b_%s" % (sMassHiggsWindowName)
                if sSelName_tmp_ not in self.sel_conditions_all_list: continue

                selection.add(
                    sSelName_tmp_,
                    ((leadingFatJet_PNet_massH_Hto4b_avg > massHiggsWindow[0]) & 
                     (leadingFatJet_PNet_massH_Hto4b_avg < massHiggsWindow[1] ))
                )

        # leadingNonHto4bFatJet_PNet_WZvsQCD
        if "PNetWZvsQCDLoose" in self.sel_conditions_all_list:
            selection.add(
                "PNetWZvsQCDLoose",
                leadingNonHto4bFatJet_PNet_WZvsQCD > self.objectSelector.NonHto4bFatJetPNet_WZvsQCD_Thsh
            )

        # nLeptonsTight veto
        if "nLeptonsTight" in self.sel_conditions_all_list:
            selection.add(
                "nLeptonsTight",
                ( nLeptonsTight <= self.objectSelector.NLeptonsTight_MaxThsh )
            )



        if "leadingFatJet_nSV" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJet_nSV",
                nSV_matched_leadingFatJet > self.objectSelector.nSV_matched_leadingFatJet_Thsh
            )



        # Trigger selection
        if sTrgSelection in self.sel_conditions_all_list:
            if sTrgSelection not in Triggers_perEra[self.datasetInfo["era"]]:
                logging.critical(f'htoaa_Analysis_VHHadronicMode.py::main():: {sTrgSelection = } not in {Triggers_perEra[self.datasetInfo["era"]] = }.')
                exit(0)  

            mask_Trgs = falses_list
            for HLTName, L1TList in Triggers_perEra[self.datasetInfo["era"]][sTrgSelection].items():
                HLTName_toUse = HLTName.replace('HLT_', '')
                mask_HLT = events.HLT[HLTName_toUse] == True

                mask_L1Ts = falses_list
                for L1TName in L1TList:
                    L1TName_toUse = L1TName.replace('L1_', '')
                    mask_L1T_i = events.L1[L1TName_toUse] == True
                    mask_L1Ts = (mask_L1Ts | mask_L1T_i) # any one of the L1T triggers associated to HLT path should be fired

                mask_Trg_i = (mask_HLT & mask_L1Ts) # HLT path and any of the associated L1T seed should be fired
                mask_Trgs = (mask_Trgs | mask_Trg_i) # Any of the HLT trigger should be fired

            selection.add(
                sTrgSelection,
                mask_Trgs
            )



        if "2018HEM1516Issue" in self.sel_conditions_all_list:
            if not self.datasetInfo['isMC']: # 2018 data
                #mask_HEM1516Issue = mask_HEM1516Issue & (events.run >= 319077)
                # Reject events satisfying HEM1516 issue conditions in 2018 data with run >= 319077
                # Use selection.add('2018HEM1516Issue') for event selection for 2018 data and for event reweight for 2018 MC
                sel_HEM1516Issue = ~ (
                    isRunAffectedBy2018HEM1516Issue & 
                    mask_HEM1516Issue
                )
            else: # 2018 MC
                sel_HEM1516Issue = trues_list
        
            # Reject events satisfying HEM1516 issue conditions in 2018 data with run >= 319077
            # Use selection.add('2018HEM1516Issue') for event selection for 2018 data and for event reweight for 2018 MC
            selection.add(
                "2018HEM1516Issue", 
                sel_HEM1516Issue
            )
            if printLevel >= 10:
                printVariable("\n mask_HEM1516Issue", mask_HEM1516Issue); sys.stdout.flush()
                printVariable("\n ~ mask_HEM1516Issue",~ mask_HEM1516Issue ); sys.stdout.flush()





        if "QCDStitch" in self.sel_conditions_all_list:
            selection.add(
                "QCDStitch",
                #mask_QCD_stitch_eventwise == True
                mask_QCD_stitch_eventwise
            )
            
        
        print(f'{self.sel_names_all["Presel"] = }')
        #sel_SR          = selection.all("nPV", "FatJetGet")
        sel_SR           = selection.all(* self.sel_names_all["Presel"])
        sel_GenHToAATo4B = None

        if self.datasetInfo['isMC'] and self.datasetInfo['isSignal']:
            # max. dR(sel_leadingFatJet, GEN 4B from H->aa)
            dr_LeadingFatJet_GenB = ak.concatenate([leadingFatJet_asSingletons.delta_r(LVGenB_0), leadingFatJet_asSingletons.delta_r(LVGenBbar_0), leadingFatJet_asSingletons.delta_r(LVGenB_1), leadingFatJet_asSingletons.delta_r(LVGenBbar_1)], axis=-1)
            max_dr_LeadingFatJet_GenB = ak.max(dr_LeadingFatJet_GenB, axis=-1)

            if printLevel >= 13:
                printVariable("\n leadingFatJet_asSingletons.delta_r(LVGenB_0)", leadingFatJet_asSingletons.delta_r(LVGenB_0))
                printVariable("\n dr_LeadingFatJet_GenB", dr_LeadingFatJet_GenB)
                printVariable("\n max_dr_LeadingFatJet_GenB", max_dr_LeadingFatJet_GenB)
                
            # GEN level selection
            selection.add("1GenHiggs", ak.num(genHiggs) == 1)
            selection.add("2GenA", ak.num(genACollection) == 2)
            selection.add("2GenAToBBbarPairs", ak.num(genBBar_pairs) == 2)
            selection.add("dR_GenH_GenB_0p8", max_dr_GenH_GenB < 0.8)
            selection.add("dR_LeadingFatJet_GenB_0p8", max_dr_LeadingFatJet_GenB < 0.8)

            # 
            sel_names_GEN = ["1GenHiggs", "2GenA", "2GenAToBBbarPairs", "dR_GenH_GenB_0p8"]
            self.sel_names_all.update( OD([
                ("GenHToAATo4B_1", ["1GenHiggs", "2GenA", "2GenAToBBbarPairs"]),
                ("GenHToAATo4B", [*sel_names_GEN]),
            ]) )
            for idx, cutName in enumerate(cuts_reco):
                if idx == 0:
                    self.sel_names_all.update( OD([
                        ("GenSR_%d" % (idx+1),  [*sel_names_GEN, cutName]),
                    ]) )
                else:
                    self.sel_names_all.update( OD([
                        ("GenSR_%d" % (idx+1),  [*self.sel_names_all["GenSR_%d" % (idx)], cutName]),
                    ]) ) 
            '''
            self.sel_names_all.update( OD([
                ("GenHToAATo4B_1", ["1GenHiggs", "2GenA", "2GenAToBBbarPairs"]),
                ("GenHToAATo4B", [*sel_names_GEN]),

                ("SR_1",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8"]),
                ("SR_2",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt"]),
                ("SR_3",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt", "FatJetEta"]),
                ("SR_4",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt", "FatJetEta", "FatJetBtagDeepB"]),
                ("SR_5",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt", "FatJetEta", "FatJetBtagDeepB", "FatJetMSoftDrop"]),
            ]) )
            '''
            
            sel_GenHToAATo4B = selection.all(* self.sel_names_all["GenHToAATo4B"])

        #print(f"\nsel_SR ({len(sel_SR)}): {sel_SR}   nEventsPass: {ak.sum(sel_SR, axis=0)}")
        #print(f"\nsel_SR_wGenCuts ({len(sel_SR_wGenCuts)}): {sel_SR_wGenCuts}   nEventsPass: {ak.sum(sel_SR_wGenCuts, axis=0)}")
        
        
        # useful debugger for selection efficiency
        if shift_syst is None and printLevel >= 5:
            print(dataset)
            for n in selection.names:
                print(
                    f"- Cut {n} pass {selection.all(n).sum()} of {len(events)} events"
                )
                print(f"selection {n} ({type(selection.all(n))}): {selection.all(n)}")
                #wgt1=np.full(len(events), self.datasetInfo["lumiScale"])
                #print(f"wgt1 ({len(wgt1)}): {wgt1}")
                

        


            
            



        ################
        # EVENT WEIGHTS
        ################
        
        # create a processor Weights object, with the same length as the number of events in the chunk
        weights              = Weights(len(events))
        weights_gen          = Weights(len(events))
        weights_GenHToAATo4B = Weights(len(events))
        weights_woHEM1516Fix = Weights(len(events))
        

        if self.datasetInfo["isMC"]:
            # lumiScale ------------------------------------
            lumiScale_toUse = None
            if self.datasetInfo["MCSamplesStitchOption"] == MCSamplesStitchOptions.PhSpOverlapRewgt and \
               self.datasetInfo['isQCD']:
                mask_PhSp_dict_ = {
                    "QCD_bEnrich": mask_QCD_bEnrich_PhSp,
                    "QCD_bGen": mask_QCD_bGen_PhSp,
                    "QCD_Incl_Remnant": mask_QCD_Incl_Remnant_PhSp,
                }
                lumiScale_toUse = getLumiScaleForPhSpOverlapRewgtMode(
                    hLumiScale      = self.datasetInfo["hMCSamplesStitch"],
                    sample_category = dataset,
                    sample_HT_value = self.datasetInfo['sample_HT_Min'],
                    mask_PhSp_dict  = mask_PhSp_dict_ )
            else:
                lumiScale_toUse = np.full(len(events), self.datasetInfo["lumiScale"])

            # MC wgt for HEM1516Issue --------------------- 
            wgt_HEM1516Issue = None
            if "2018HEM1516Issue" in self.sel_names_all["Presel"]:
                wgt_HEM1516Issue = ak.where(
                    mask_HEM1516Issue, # events w/ jets in HEM15/16 affected phase space 
                    np.full(len(events), (1. - DataFractionAffectedBy2018HEM1516Issue)), 
                    ones_list
                )

            # MC PURewgt ----------------------------------
            wgt_PU = getPURewgts(
                PU_list  = events.Pileup.nTrueInt,
                hPURewgt = self.hPURewgt
            )

            # MC GGF HToAATo4B Higgs pT reweight
            wgt_HiggsPt = None
            if self.datasetInfo['isSignalGGH']:
                wgt_HiggsPt = getHiggsPtRewgtForGGToHToAATo4B(
                    GenHiggsPt_list = ak.firsts(genHiggs.pt)
                )

            # MC QCD_bGen HT reweight ---------------------
            wgt_HT = None
            if self.datasetInfo['isQCD_bGen']:
                wgt_HT = getHTReweight(
                    HT_list            = events.LHE.HT,
                    sFitFunctionFormat = self.datasetInfo['HTRewgt']["fitFunctionFormat"],
                    sFitFunction       = self.datasetInfo['HTRewgt']["fitFunction"],
                    sFitFunctionRange  = self.datasetInfo['HTRewgt']["fitFunctionHTRange"]
                )            

            # MC top pT reweigts for ttbar sample ---------
            if self.datasetInfo['isTTbar']:
                wgt_TopPt = getTopPtRewgt(
                    eventsGenPart = events.GenPart[mask_genTopQuark],
                    isPythiaTuneCP5 = self.datasetInfo['isPythiaTuneCP5']
                )   

            # MC ParticleNetMD_XbbvsQCD SFs      SFs_ParticleNetMD_XbbvsQCD
            if "leadingFatJetParticleNetMD_XbbvsQCD" in self.sel_names_all["Presel"]:
                mask_ParticleNetMD_XbbvsQCD_SFRegion = (
                    (n_leadingFatJat_matched_genB >= 2) &
                    (leadingFatJetParticleNetMD_XbbvsQCD > self.objectSelector.FatJetParticleNetMD_XbbvsQCD_Thsh)
                )
                wgt_ParticleNetMD_XbbvsQCD = ak.fill_none( 
                    ak.where(
                        mask_ParticleNetMD_XbbvsQCD_SFRegion,
                        self.SFs_ParticleNetMD_XbbvsQCD(leadingFatJet.pt),
                        ones_list
                    ), 
                    1
                )



            weights.add(
                "lumiWeight",
                weight = lumiScale_toUse
            )
            weights.add(
                "genWeight",
                weight = np.copysign(np.ones(len(events)), events.genWeight)
            )
            if "2018HEM1516Issue" in self.sel_names_all["Presel"]:
                weights.add(
                    "2018HEM1516IssueWeight",
                    weight = wgt_HEM1516Issue
                )
            weights.add(
                "PUWeight",
                weight = wgt_PU
            )
            if self.datasetInfo['isSignalGGH']:
                weights.add(
                    "GGHPtRewgt",
                    weight = wgt_HiggsPt
                )
            if self.datasetInfo['isQCD_bGen']:
                weights.add(
                    "HTRewgt",
                    weight = wgt_HT
                )            
            if self.datasetInfo['isTTbar']:
                weights.add(
                    "TopPtReWeight",
                    weight = wgt_TopPt
                )         
            if "leadingFatJetParticleNetMD_XbbvsQCD" in self.sel_names_all["Presel"]:
                weights.add(
                    "SF_ParticleNetMD_XbbvsQCD",
                    weight = wgt_ParticleNetMD_XbbvsQCD
                )
            
    
            
 
            ## weights_woHEM1516Fix --------------------
            weights_woHEM1516Fix.add(
                "lumiWeight",
                weight = lumiScale_toUse
            )
            weights_woHEM1516Fix.add(
                "genWeight",
                weight = np.copysign(np.ones(len(events)), events.genWeight)
            )
            #if "2018HEM1516Issue" in self.sel_names_all["Presel"]:
            #    weights.add(
            #        "2018HEM1516IssueWeight",
            #        weight = wgt_HEM1516Issue
            #    )
            weights_woHEM1516Fix.add(
                "PUWeight",
                weight = wgt_PU
            )
            if self.datasetInfo['isSignalGGH']:
                weights_woHEM1516Fix.add(
                    "GGHPtRewgt",
                    weight = wgt_HiggsPt
                )
            if self.datasetInfo['isQCD_bGen']:
                weights_woHEM1516Fix.add(
                    "HTRewgt",
                    weight = wgt_HT
                )  
            if self.datasetInfo['isTTbar']:
                weights_woHEM1516Fix.add(
                    "TopPtReWeight",
                    weight = wgt_TopPt
                )            
            if "leadingFatJetParticleNetMD_XbbvsQCD" in self.sel_names_all["Presel"]:
                weights_woHEM1516Fix.add(
                    "SF_ParticleNetMD_XbbvsQCD",
                    weight = wgt_ParticleNetMD_XbbvsQCD
                )



            ## weights_gen -------------------------------
            weights_gen.add(
                "lumiWeight",
                weight = lumiScale_toUse 
            )
            weights_gen.add(
                "genWeight",
                weight=np.copysign(np.ones(len(events)), events.genWeight)
            )
            if self.datasetInfo['isSignalGGH']:
                weights_gen.add(
                    "GGHPtRewgt",
                    weight = wgt_HiggsPt
                )
            if self.datasetInfo['isQCD_bGen']:
                weights_gen.add(
                    "HTRewgt",
                    weight = wgt_HT
                )

            if printLevel >= 10:
                print(f"\nevents.genWeight ({events.genWeight.fields}) ({len(events.genWeight)}): {events.genWeight.to_list()}")
                genWgt1 = np.copysign(np.ones(len(events)), events.genWeight)
                print(f"genWgt1 ({len(genWgt1)}): {genWgt1}")
                
                print(f"\n\nevents.btagWeight ({events.btagWeight.fields}) ({len(events.btagWeight)}): {events.btagWeight.to_list()}")
                #print(f"events.nPSWeight {events.nPSWeight.to_list()}")
                print(f"\n\nevents.PSWeight ({events.PSWeight.fields}) ({len(events.PSWeight)}): {events.PSWeight.to_list()}")
                print(f"\n\nevents.LHEWeight ({events.LHEWeight.fields}) ({len(events.LHEWeight)}): {events.LHEWeight.to_list()}")
                print(f"\n\nevents.LHEScaleWeight ({events.LHEScaleWeight.fields}) ({len(events.LHEScaleWeight)}): {events.LHEScaleWeight.to_list()}")
                print(f"\n\nevents.LHEReweightingWeight ({events.LHEReweightingWeight.fields}) ({len(events.LHEReweightingWeight)}): {events.LHEReweightingWeight.to_list()}")

                

            if printLevel >= 13:                
                printVariable("\n weights.weight()", weights.weight())
                printVariable("\n weights_gen.weight()", weights_gen.weight())
                printVariable("\n weights_GenHToAATo4B.weight()", weights_GenHToAATo4B.weight())
                #printVariable("\n weights.weight().sum()", weights.weight().sum())
                print(f"weights.weight().sum(): {weights.weight().sum()}")
                printVariable("\n weights.weight()[sel_SR]", weights.weight()[sel_SR])
                print(f"weights.weight()[sel_SR].sum(): {weights.weight()[sel_SR].sum()}")

            if printLevel >= 30:
                print(f"leadingFatJet.tau4[sel_SR]): {leadingFatJet.tau4[sel_SR]}")
                print(f"leadingFatJet.tau3[sel_SR]): {leadingFatJet.tau3[sel_SR]}")
                print(f"divide                     : {np.divide(leadingFatJet.tau4[sel_SR], leadingFatJet.tau3[sel_SR])} \n\n\n")







        ###################
        # FILL HISTOGRAMS
        ###################

        systList = []
        if self.datasetInfo['isMC']:
            if shift_syst is None:
                systList = [
                    "central"
                ]
            else:
                systList = [shift_syst]
        else:
            systList = ["noweight"]

            
        output['cutflow']['all events'] += len(events)
        output['cutflow'][sWeighted+'all events'] += weights.weight().sum() 
        #for n in selection.names:
        #    output['cutflow'][n] += selection.all(n).sum()
        
        for iSelection in self.sel_names_all.keys():
            iName = f"{iSelection}: {self.sel_names_all[iSelection]}"
            sel_i = selection.all(* self.sel_names_all[iSelection])
            output['cutflow'][iName] += sel_i.sum()
            output['cutflow'][sWeighted+iName] +=  weights.weight()[sel_i].sum()

              

        for syst in systList:

            # find the event weight to be used when filling the histograms
            weightSyst = syst
            
            # in the case of 'central', or the jet energy systematics, no weight systematic variation is used (weightSyst=None)
            if syst in ["central", "JERUp", "JERDown", "JESUp", "JESDown"]:
                weightSyst = None

            
            
            if syst == "noweight":
                evtWeight                = np.ones(len(events))
                evtWeight_woHEM1516Fix   = np.ones(len(events))
            else:
                evtWeight                = weights.weight(weightSyst)
                evtWeight_gen            = weights_gen.weight(weightSyst)
                evtWeight_woHEM1516Fix   = weights_woHEM1516Fix.weight(weightSyst)



            ### General or GEN-level histograms ========================================================================================

            # PU
            output['hPV_npvs_beforeSel'].fill(
                dataset=dataset,
                PU=(events.PV.npvs),
                systematic=syst,
                weight=evtWeight
            )            
            output['hPV_npvsGood_beforeSel'].fill(
                dataset=dataset,
                PU=(events.PV.npvsGood),
                systematic=syst,
                weight=evtWeight
            )


            ## isMC --------------------------------------------------
            if self.datasetInfo['isMC'] and runMode_GenLHEPlots:                 
                output['hPileup_nTrueInt'].fill(
                    dataset=dataset,
                    PU=(events.Pileup.nTrueInt),
                    systematic=syst,
                    weight=evtWeight
                )            
                output['hPileup_nPU'].fill(
                    dataset=dataset,
                    PU=(events.Pileup.nPU),
                    systematic=syst,
                    weight=evtWeight
                )      

                output['hGenLHE_HT_all'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenLHE_HTIncoming_all'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HTIncoming),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenLHE_Vpt_all'].fill(
                    dataset=dataset,
                    HT=(events.LHE.Vpt),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenLHE_AlphaS_all'].fill(
                    dataset=dataset,
                    alphaS=(events.LHE.AlphaS),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenLHE_Njets_all'].fill(
                    dataset=dataset,
                    nObject50=(events.LHE.Njets),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenLHE_Nb_all'].fill(
                    dataset=dataset,
                    nObject50=(events.LHE.Nb),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenLHE_Nc_all'].fill(
                    dataset=dataset,
                    nObject50=(events.LHE.Nc),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenLHE_Nuds_all'].fill(
                    dataset=dataset,
                    nObject50=(events.LHE.Nuds),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenLHE_Nglu_all'].fill(
                    dataset=dataset,
                    nObject50=(events.LHE.Nglu),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenLHE_NpNLO_all'].fill(
                    dataset=dataset,
                    nObject200=(events.LHE.NpNLO),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenLHE_NpLO_all'].fill(
                    dataset=dataset,
                    nObject200=(events.LHE.NpLO),
                    systematic=syst,
                    weight=evtWeight_gen
                )


            if (self.datasetInfo['isSignal'] or self.datasetInfo['isHToBB']) and runMode_SignalGenChecks:
                output['hGenHiggsPt_all'].fill(
                    dataset=dataset,
                    Pt2TeV=(ak.firsts(genHiggs.pt)),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenHiggsLog2Pt_all'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(ak.firsts(genHiggs.pt)),
                    systematic=syst,
                    weight=evtWeight_gen
                )                


            ## isMC && isSignal ------------------------------------------------------------------------------------------------------------

            if self.datasetInfo['isSignal'] and runMode_SignalGenChecks: 
                output['hGenHiggsPt_GenHToAATo4B'].fill(
                    dataset=dataset,
                    Pt=(ak.firsts(genHiggs.pt[sel_GenHToAATo4B])),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                output['hGenHiggsPt_sel'].fill(
                    dataset=dataset,
                    Pt=(ak.firsts(genHiggs.pt[sel_SR])),
                    systematic=syst,
                    weight=evtWeight_gen[sel_SR]
                )
                output['hGenHiggsPt_sel_wGenCuts'].fill(
                    dataset=dataset,
                    Pt=(ak.firsts(genHiggs.pt[sel_GenHToAATo4B])),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                
                # m(2b from ATo2B) and m(4b from HToAATo4b) --------------                   
                output['hGenHiggsMass_all_0'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass)),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                    
                output['hMass_GenA_all_0'].fill(
                    dataset=dataset,
                    Mass=(genACollection[:, 0].mass),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hMass_GenA_all_0'].fill(
                    dataset=dataset,
                    Mass=(genACollection[:, 1].mass),
                    systematic=syst,
                    weight=evtWeight_gen
                )


                output['hGenHiggsMass_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                    
                output['hMass_GenA_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[:, 0].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                output['hMass_GenA_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
               
                '''
                output['hMass_GenA_all'].fill(
                    dataset=dataset,
                    Pt=(ak.flatten(genA_Second.mass)),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                '''
                output['hMass_GenAApair_all'].fill(
                    dataset=dataset,
                    Mass=((genA_First + genA_Second).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                output['hMass_GenAToBBbarpair_all'].fill(
                    dataset=dataset,
                    Mass=((events.GenPart[genBBar_pairs['b']] + events.GenPart[genBBar_pairs['bbar']]).mass[:, 0][sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )                
                output['hMass_GenAToBBbarpair_all'].fill(
                    dataset=dataset,
                    Mass=((events.GenPart[genBBar_pairs['b']] + events.GenPart[genBBar_pairs['bbar']]).mass[:, 1][sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                output['hMass_Gen4BFromHToAA_all'].fill(
                    dataset=dataset,
                    Mass=((events.GenPart[genBBar_pairs['b']][:, 0] + events.GenPart[genBBar_pairs['bbar']][:, 0] + events.GenPart[genBBar_pairs['b']][:, 1] + events.GenPart[genBBar_pairs['bbar']][:, 1]).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )


                output['hMass_GenAToBBbarpair_all_1'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_0 + LVGenBbar_0).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )                
                output['hMass_GenAToBBbarpair_all_1'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_1 + LVGenBbar_1).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                output['hMass_Gen4BFromHToAA_all_1'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_0 + LVGenBbar_0 + LVGenB_1 + LVGenBbar_1).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )


                output['hMass_GenA1_vs_GenA2_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[:, 0].mass[sel_GenHToAATo4B]),
                    Mass1=(genACollection[:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                
                output['hMass_GenA1ToBBbar_vs_GenA2ToBBbar_all'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_0 + LVGenBbar_0).mass[sel_GenHToAATo4B]),
                    Mass1=((LVGenB_1 + LVGenBbar_1).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hMass_GenH_vs_GenAHeavy_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    Mass1=(genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )
                
                output['hMass_GenH_vs_GenALight_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    Mass1=(genACollection[idxGenA_sortByMass][:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hMass_GenAHeavy_vs_GenALight_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B]),
                    Mass1=(genACollection[idxGenA_sortByMass][:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hDeltaR_GenH_GenB_max'].fill(
                    dataset=dataset,
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hMassGenH_vs_maxDRGenHGenB_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hMassGenAHeavy_vs_maxDRGenHGenB_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B]),
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )

                output['hMassGenALight_vs_maxDRGenHGenB_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[idxGenA_sortByMass][:, 1].mass[sel_GenHToAATo4B]),
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen[sel_GenHToAATo4B]
                )



            # QCD MC ----------------------------------------------
            if self.datasetInfo['isQCD'] and runMode_QCDGenValidation:
                # all events
                iBin = 0
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list * iBin),
                    systematic=syst
                )
                output['hCutFlowWeighted'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list * iBin),
                    systematic=syst,
                    weight=evtWeight
                )
                
                # genBHadrons_status2 events
                iBin = 1
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_genBHadrons_status2_eventwise] * iBin),
                    systematic=syst
                )
                output['hCutFlowWeighted'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_genBHadrons_status2_eventwise] * iBin),
                    systematic=syst,
                    weight=evtWeight[mask_genBHadrons_status2_eventwise]
                )

                # genBHadrons_status2 events
                iBin = 2
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_genBQuarks_hardSctred_eventwise] * iBin),
                    systematic=syst
                )
                output['hCutFlowWeighted'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_genBQuarks_hardSctred_eventwise] * iBin),
                    systematic=syst,
                    weight=evtWeight[mask_genBQuarks_hardSctred_eventwise]
                )

                # QCD_stitch events
                iBin = 3
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_QCD_stitch_eventwise] * iBin),
                    systematic=syst
                )
                output['hCutFlowWeighted'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[mask_QCD_stitch_eventwise] * iBin),
                    systematic=syst,
                    weight=evtWeight[mask_QCD_stitch_eventwise]
                )
                
                # NEvents in QCD HT samples
                QCDSamplesHTBins_LowEdge = [50, 100, 200, 300, 500, 700, 1000, 1500, 2000]
                idx_QCDSampleHTBin = None
                for idx_ in range(0, len(QCDSamplesHTBins_LowEdge)):
                    if self.datasetInfo['sample_HT_Min'] == QCDSamplesHTBins_LowEdge[idx_]:
                        idx_QCDSampleHTBin = idx_                
                iBin = (idx_QCDSampleHTBin * 5) 
                output['hNEventsQCD'].fill(
                    dataset=dataset,
                    CutFlow50=(ones_list * iBin),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                iBin = (idx_QCDSampleHTBin * 5) + 1
                output['hNEventsQCD'].fill(
                    dataset=dataset,
                    CutFlow50=(ones_list[mask_genBHadrons_status2_eventwise] * iBin),
                    systematic=syst,
                    weight=evtWeight_gen[mask_genBHadrons_status2_eventwise]
                )
                iBin = (idx_QCDSampleHTBin * 5) + 2
                output['hNEventsQCD'].fill(
                    dataset=dataset,
                    CutFlow50=(ones_list[mask_genBQuarks_hardSctred_eventwise] * iBin),
                    systematic=syst,
                    weight=evtWeight_gen[mask_genBQuarks_hardSctred_eventwise]
                )
                iBin = (idx_QCDSampleHTBin * 5) + 3
                output['hNEventsQCD'].fill(
                    dataset=dataset,
                    CutFlow50=(ones_list[mask_genBHadrons_status2_and_noGenBQuarksHardSctred_eventwise] * iBin),
                    systematic=syst,
                    weight=evtWeight_gen[mask_genBHadrons_status2_and_noGenBQuarksHardSctred_eventwise]
                )
                
                iBin = (idx_QCDSampleHTBin * 5)
                output['hNEventsQCDUnweighted'].fill(
                    dataset=dataset,
                    CutFlow50=(ones_list * iBin),
                    systematic=syst
                )
                iBin = (idx_QCDSampleHTBin * 5) + 1
                output['hNEventsQCDUnweighted'].fill(
                    dataset=dataset,
                    CutFlow50=(ones_list[mask_genBHadrons_status2_eventwise] * iBin),
                    systematic=syst
                )
                iBin = (idx_QCDSampleHTBin * 5) + 2
                output['hNEventsQCDUnweighted'].fill(
                    dataset=dataset,
                    CutFlow50=(ones_list[mask_genBQuarks_hardSctred_eventwise] * iBin),
                    systematic=syst
                )
                iBin = (idx_QCDSampleHTBin * 5) + 3
                output['hNEventsQCDUnweighted'].fill(
                    dataset=dataset,
                    CutFlow50=(ones_list[mask_genBHadrons_status2_and_noGenBQuarksHardSctred_eventwise] * iBin),
                    systematic=syst
                )



                output['hGenLHE_HT_SelQCDbHadron'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT[mask_genBHadrons_status2_eventwise]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_genBHadrons_status2_eventwise]
                )
                output['hGenLHE_HT_SelQCDbGen'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT[mask_genBHadrons_status2_and_noGenBQuarksHardSctred_eventwise]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_genBHadrons_status2_and_noGenBQuarksHardSctred_eventwise]
                )
                output['hGenLHE_HT_SelQCDbEnrich'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT[mask_genBQuarks_hardSctred_eventwise]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_genBQuarks_hardSctred_eventwise]
                )
                
                output['hGenLHE_HT_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT[mask_QCD_stitch_CutBQuarkPt_eventwise]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_QCD_stitch_CutBQuarkPt_eventwise]
                )
                output['hGenLHE_HT_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT[mask_QCD_stitch_CutBHadron_eventwise]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_QCD_stitch_CutBHadron_eventwise]
                )
                output['hGenLHE_HT_QCDStitch'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT[mask_QCD_stitch_eventwise]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_QCD_stitch_eventwise]
                )
                output['hGenLHE_HT_QCD_bEnrich_PhSp'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT[mask_QCD_bEnrich_PhSp]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_QCD_bEnrich_PhSp]
                )
                output['hGenLHE_HT_QCD_bGen_PhSp'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT[mask_QCD_bGen_PhSp]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_QCD_bGen_PhSp]
                )
                output['hGenLHE_HT_QCD_Incl_Remnant_PhSp'].fill(
                    dataset=dataset,
                    HT=(events.LHE.HT[mask_QCD_Incl_Remnant_PhSp]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_QCD_Incl_Remnant_PhSp]
                )

                '''
                output['hGenBquark_Status_all'].fill(
                    dataset=dataset,
                    PytPartStatus=(ak.flatten(genBQuarks.status)),
                    systematic=syst
                )
                output['hGenBquark_first_Status_all'].fill(
                    dataset=dataset,
                    PytPartStatus=(ak.to_numpy(genBQuarks_first[mask_genBQuarks].status)),
                    systematic=syst
                )
                output['hGenBquark_first_PdgId_all'].fill(
                    dataset=dataset,
                    PdgId=(ak.to_numpy(abs(genBQuarks_first[mask_genBQuarks].pdgId))),
                    systematic=syst
                )

                for statusFlag_ in GENPART_STATUSFLAGS:
                    output['hGenBquark_first_%s_all' % (statusFlag_)].fill(
                        dataset=dataset,
                        Boolean=(ak.to_numpy(genBQuarks_first[mask_genBQuarks].hasFlags(statusFlag_))),
                        systematic=syst
                    )
                '''
                
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 1)
                output['hGenBquark_leadingPt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[mask_tmp][:, 0]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 2)
                output['hGenBquark_subleadingPt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[mask_tmp][:, 1]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 3)
                output['hGenBquark_thirdLeadingPt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[mask_tmp][:, 2]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 4)
                output['hGenBquark_forthLeadingPt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[mask_tmp][:, 3]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 4)
                output['hGenBquark_forthLeadingPt_UltraLow_all'].fill(
                    dataset=dataset,
                    PtUltraLow=(genBQuarks_pT[mask_tmp][:, 3]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 1)
                output['hGenBquark_leadingPt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 0]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 2)
                output['hGenBquark_subleadingPt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 1]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 3)
                output['hGenBquark_thirdLeadingPt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 2]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 4)
                output['hGenBquark_forthLeadingPt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 3]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 1)
                output['hGenBquark_leadingPt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 0]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 2)
                output['hGenBquark_subleadingPt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 1]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 3)
                output['hGenBquark_thirdLeadingPt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 2]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                mask_tmp = (ak.count(genBQuarks_pT, axis=-1) >= 4)
                output['hGenBquark_forthLeadingPt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_pT[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 3]),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 1)
                output['hLeadingPtGenBquark_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][mask_tmp][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hLeadingPtGenBquark_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][mask_tmp][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 1)
                output['hLeadingPtGenBquarkHardSctred_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][mask_tmp][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hLeadingPtGenBquarkHardSctred_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][mask_tmp][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 1)
                output['hLeadingPtGenBHadron_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][mask_tmp][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hLeadingPtGenBHadron_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][mask_tmp][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 1)
                output['hLeadingPtGenBHadronStatus2_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][mask_tmp][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hLeadingPtGenBHadronStatus2_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][mask_tmp][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 1)
                output['hLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 1)
                output['hLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 1)
                output['hLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 1)
                output['hLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 1)
                output['hLeadingPtGenBquark_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hLeadingPtGenBquark_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 1)
                output['hLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 1)
                output['hLeadingPtGenBHadron_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hLeadingPtGenBHadron_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 1)
                output['hLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 0].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 0].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )


                # subleading pT ----------------------------------------------------------------------------------------------------------------------------
                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBquark_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][mask_tmp][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hSubleadingPtGenBquark_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][mask_tmp][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBquarkHardSctred_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][mask_tmp][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hSubleadingPtGenBquarkHardSctred_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][mask_tmp][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBHadron_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][mask_tmp][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hSubleadingPtGenBHadron_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][mask_tmp][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBHadronStatus2_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][mask_tmp][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hSubleadingPtGenBHadronStatus2_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][mask_tmp][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBquark_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hSubleadingPtGenBquark_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hSubleadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBHadron_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hSubleadingPtGenBHadron_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hSubleadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBquark_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hSubleadingPtGenBquark_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hSubleadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBHadron_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hSubleadingPtGenBHadron_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 2)
                output['hSubleadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 1].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hSubleadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 1].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                # ------------------------------------------------------------------------------------------------------------------------------------------


                # Third-leading pT -------------------------------------------------------------------------------------------------------------------------
                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBquark_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][mask_tmp][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hThirdLeadingPtGenBquark_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][mask_tmp][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBquarkHardSctred_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][mask_tmp][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hThirdLeadingPtGenBquarkHardSctred_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][mask_tmp][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBHadron_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][mask_tmp][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hThirdLeadingPtGenBHadron_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][mask_tmp][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBHadronStatus2_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][mask_tmp][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hThirdLeadingPtGenBHadronStatus2_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][mask_tmp][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hThirdLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hThirdLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hThirdLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hThirdLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBquark_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hThirdLeadingPtGenBquark_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hThirdLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBHadron_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hThirdLeadingPtGenBHadron_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 3)
                output['hThirdLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 2].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hThirdLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 2].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )              
                # ------------------------------------------------------------------------------------------------------------------------------------------


                # Fourth-leading pT ------------------------------------------------------------------------------------------------------------------------
                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBquark_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][mask_tmp][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hFourthLeadingPtGenBquark_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][mask_tmp][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBquarkHardSctred_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][mask_tmp][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hFourthLeadingPtGenBquarkHardSctred_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][mask_tmp][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBHadron_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][mask_tmp][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hFourthLeadingPtGenBHadron_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][mask_tmp][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBHadronStatus2_pt_all'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][mask_tmp][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )
                output['hFourthLeadingPtGenBHadronStatus2_eta_all'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][mask_tmp][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[mask_tmp]
                )


                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hFourthLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hFourthLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hFourthLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )
                output['hFourthLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBQuarkPt_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBquark_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hFourthLeadingPtGenBquark_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )


                mask_tmp = (ak.count(genBQuarks_hardSctred.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hFourthLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )


                mask_tmp = (ak.count(genBHadrons.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBHadron_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hFourthLeadingPtGenBHadron_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons[idx_genBHadrons_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 4)
                output['hFourthLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    PtLow=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 3].pt),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )
                output['hFourthLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron'].fill(
                    dataset=dataset,
                    Eta=(genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 3].eta),
                    systematic=syst,
                    weight=evtWeight_gen[(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)]
                )                
                # ------------------------------------------------------------------------------------------------------------------------------------------






            ### RECO-level histograms ============================================================================================
            
            if self.datasetInfo['isSignal']:
                output['hIdxFatJetMatchedToGenBFromHToAATo4B'].fill(
                    dataset=dataset,
                    nObject=(idx_FatJet_matched_genB_HToAATo4B[mask_events_FatJet_matched_genB_HToAATo4B]),
                    systematic=syst,
                    weight=evtWeight[mask_events_FatJet_matched_genB_HToAATo4B]
                ) 
                if 'particleNetMD_Hto4b_Haa4b' in events.FatJet.fields:
                    output['hIdxFatJetMaxPNetMD_Hto4b_Haa4bOverQCD'].fill(
                        dataset=dataset,
                        nObject=(ak.firsts(idx_FatJet_PNetMD_Hto4b_Haa4bOverQCD_max)[mask_events_FatJet_matched_genB_HToAATo4B]),
                        systematic=syst,
                        weight=evtWeight[mask_events_FatJet_matched_genB_HToAATo4B]
                    ) 
                    output['hIdxFatJetMaxPNetMD_Hto4b_Haa4bOverQCD_1'].fill(
                        dataset=dataset,
                        nObject=(ak.firsts(idx_FatJet_PNetMD_Hto4b_Haa4bOverQCD_max)[(leadingFatJet.nBHadrons >= 4)]),
                        systematic=syst,
                        weight=evtWeight[(leadingFatJet.nBHadrons >= 4)]
                    ) 
                #printVariable('\n mask_events_FatJet_matched_genB_HToAATo4B',mask_events_FatJet_matched_genB_HToAATo4B)
                #printVariable('\n (leadingFatJet.nBHadrons >= 4)', (leadingFatJet.nBHadrons >= 4))
                output['hIdxFatJetMaxZHbb_plus_Xbb'].fill(
                    dataset=dataset,
                    nObject=(ak.firsts(idx_FatJet_ZHbb_plus_Xbb_max)[mask_events_FatJet_matched_genB_HToAATo4B]),
                    systematic=syst,
                    weight=evtWeight[mask_events_FatJet_matched_genB_HToAATo4B]
                ) 
                output['hIdxFatJetMaxZHbb_plus_Xbb_1'].fill(
                    dataset=dataset,
                    nObject=(ak.firsts(idx_FatJet_ZHbb_plus_Xbb_max)[ak.fill_none(leadingFatJet.nBHadrons >= 4, False)]),
                    systematic=syst,
                    weight=evtWeight[ak.fill_none(leadingFatJet.nBHadrons >= 4, False)]
                )
                output['hLeadingBtagFatJetPtOverLeadingFatJetPt_Sig'].fill(
                    dataset=dataset,
                    Ratio=( leadingFatJet[mask_events_FatJet_matched_genB_HToAATo4B].pt / ak.firsts(events.FatJet)[mask_events_FatJet_matched_genB_HToAATo4B].pt),
                    systematic=syst,
                    weight=evtWeight[mask_events_FatJet_matched_genB_HToAATo4B]
                )


            for sel_name in self.sel_names_all.keys(): # loop of list of selections
                

                if sel_name.startswith('Gen'): continue

                sel_SR_toUse = selection.all(* self.sel_names_all[sel_name])
                sel_SR_woSel2018HEM1516_toUse = None
                if "2018HEM1516Issue" in self.sel_names_all[sel_name]:
                    sel_names_all_toUse_ = list(set(self.sel_names_all[sel_name]) - set(["2018HEM1516Issue"])) # all sel_name conditions w/o "2018HEM1516Issue"
                    sel_SR_woSel2018HEM1516_toUse = selection.all(* sel_names_all_toUse_)
                else:
                    sel_SR_woSel2018HEM1516_toUse = sel_SR_toUse


                    


                #printVariable('\n sel_SR', sel_SR)
                for sHExt_0 in self.histosExtensions: # HistogramNameExtensions_QCD = ['_0b', '_1b', '_2b', '_3b', '_4b', '_5bAndMore'], else ['']
                    sHExt = "_%s" % (sel_name)
                    if sHExt_0 != '':
                        sHExt += "_%s" % (sHExt_0)

                    sel_SR_forHExt = None
                    sel_SR_woSel2018HEM1516_forHExt = None

                    if sHExt_0 == '':
                        # No additional GEN-level category
                        sel_SR_forHExt = sel_SR_toUse
                        sel_SR_woSel2018HEM1516_forHExt = sel_SR_woSel2018HEM1516_toUse
                    else:
                        # Split in GEN-level categories
                        nGenBInFatJet = 0
                        if   '0b' in sHExt_0:
                            nGenBInFatJet = 0
                        elif '1b' in sHExt_0:
                            nGenBInFatJet = 1
                        elif '2b' in sHExt_0:
                            nGenBInFatJet = 2
                        elif '3b' in sHExt_0:
                            nGenBInFatJet = 3
                        elif '4b' in sHExt_0:
                            nGenBInFatJet = 4
                        elif '5b' in sHExt_0:
                            nGenBInFatJet = 5

                        if 'AndMore' in sHExt_0:
                            mask_HExt = (n_leadingFatJat_matched_genB >= nGenBInFatJet)
                        else:
                            mask_HExt = (n_leadingFatJat_matched_genB == nGenBInFatJet)

                        mask_HExt = ak.fill_none(mask_HExt, False) # mask for events without FatJet are None. It causes error at the later stage.
                        sel_SR_forHExt = sel_SR_toUse & mask_HExt
                        sel_SR_woSel2018HEM1516_forHExt = sel_SR_woSel2018HEM1516_toUse & mask_HExt

                        #print(f" { n_leadingFatJat_matched_genB[ak.is_none(mask_HExt)] = } "); sys.stdout.flush(); 
                        #print(f" { mask_leadingFatJat_matched_genB[ak.is_none(mask_HExt)] = } "); sys.stdout.flush();
                        #print(f" { vGenBQuarksHardSctred_genBHadronsStatus2_sel[ak.is_none(mask_HExt)] = } "); sys.stdout.flush();
                        #printVariable('\n vGenBQuarksHardSctred_genBHadronsStatus2_sel[ak.is_none(mask_HExt)]', vGenBQuarksHardSctred_genBHadronsStatus2_sel[ak.is_none(mask_HExt)]); sys.stdout.flush();
                        #printVariable('\n leadingFatJet[ak.is_none(mask_HExt)]', leadingFatJet[ak.is_none(mask_HExt)]); sys.stdout.flush();
                    

                        #printVariable('\n mask_HExt (%s)'%sHExt, mask_HExt)
                    #printVariable('\n sel_SR_forHExt (%s)'%sHExt, sel_SR_forHExt)
                    #printVariable('\n evtWeight', evtWeight)
                    #printVariable('\n evtWeight[sel_SR_forHExt]', evtWeight[sel_SR_forHExt])
                    #print(f"{sHExt = }, {len(sel_SR) = }, {ak.sum(sel_SR, axis=0) = }, {ak.sum(mask_HExt, axis=0) = }, {ak.sum(sel_SR_forHExt, axis=0) = }"); sys.stdout.flush();
                    #print(f"{ak.sum(ak.is_none(sel_SR), axis=0) =  }, {ak.sum(ak.is_none(mask_HExt), axis=0) =  }, {ak.sum(ak.is_none(sel_SR_forHExt), axis=0) =  }, "); sys.stdout.flush();

                    sel_SR_forHExt_woGenMatch = sel_SR_forHExt
                    #if runMode_OptimizePNetTaggerCut and self.datasetInfo['isSignal']:
                    #if self.datasetInfo['isSignal']:
                    if (self.datasetInfo['isSignal'] and \
                        (runMode_SignificancsScan2D or \
                         runMode_OptimizePNetTaggerCut or \
                         runMode_SignalGenCuts) ):
                        sel_SR_forHExt = sel_SR_forHExt & (n_leadingFatJat_matched_genB_HToAATo4B >= 4)
                    sel_SR_forHExt = ak.fill_none(sel_SR_forHExt, False) 


                    # all events
                    iBin = 0
                    output['hCutFlow'+sHExt].fill(
                        dataset=dataset,
                        CutFlow=(ones_list * iBin),
                        systematic=syst
                    )
                    output['hCutFlowWeighted'+sHExt].fill(
                        dataset=dataset,
                        CutFlow=(ones_list * iBin),
                        systematic=syst,
                        weight=evtWeight
                    )

                    # events passing SR
                    iBin = 4
                    output['hCutFlow'+sHExt].fill(
                        dataset=dataset,
                        CutFlow=(ones_list[sel_SR_forHExt] * iBin),
                        systematic=syst
                    )
                    output['hCutFlowWeighted'+sHExt].fill(
                        dataset=dataset,
                        CutFlow=(ones_list[sel_SR_forHExt] * iBin),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )

                    output['hPV_npvs'+sHExt].fill(
                        dataset=dataset,
                        PU=(events.PV.npvs[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )            
                    output['hPV_npvsGood'+sHExt].fill(
                        dataset=dataset,
                        PU=(events.PV.npvsGood[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )

                    
                    output['hLeadingFatJetPt'+sHExt].fill(
                        dataset=dataset,
                        #Pt=ak.flatten(selFatJet.pt[sel_SR_forHExt][:, 0]),
                        #Pt=(selFatJet.pt[sel_SR_forHExt][:, 0]),
                        Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )            
                    output['hLeadingFatJetEta'+sHExt].fill(
                        dataset=dataset,
                        Eta=(leadingFatJet.eta[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetPhi'+sHExt].fill(
                        dataset=dataset,
                        Phi=(leadingFatJet.phi[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    if 'hLeadingFatJetEta_vs_Phi'+sHExt in output:
                        output['hLeadingFatJetEta_vs_Phi'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[sel_SR_forHExt]),
                            Phi=(leadingFatJet.phi[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                    output['hLeadingBtagFatJetPtOverLeadingFatJetPt'+sHExt].fill(
                        dataset=dataset,
                        Ratio=( leadingFatJet[sel_SR_forHExt].pt / ak.firsts(events.FatJet)[sel_SR_forHExt].pt),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )



                    # 2018 HEM15/16 issue ----------------------               
                    if self.datasetInfo["era"] == Era_2018 and runMode_2018HEM1516IssueValidation:

                        if printLevel >= 10:
                            printVariable('\n sel_SR_forHExt', sel_SR_forHExt); sys.stdout.flush()
                            printVariable('\n mask_HEM1516Issue_Eta', mask_HEM1516Issue_Eta); sys.stdout.flush()
                            printVariable('\n mask_HEM1516Issue_Phi', mask_HEM1516Issue_Phi); sys.stdout.flush()
                            printVariable('\n sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi', sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi); sys.stdout.flush()
                            printVariable('\n leadingFatJet.pt[ sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi ]', leadingFatJet.pt[ sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi ]); sys.stdout.flush()

                        if not self.datasetInfo['isMC']:
                            mask_DataPreHEM1516Issue  = ~ isRunAffectedBy2018HEM1516Issue
                            mask_DataWithHEM1516Issue = isRunAffectedBy2018HEM1516Issue
                            evtWeight_DataPreHEM1516Issue               = evtWeight
                            evtWeight_DataWithHEM1516Issue              = evtWeight
                            evtWeight_woHEM1516Fix_DataPreHEM1516Issue  = evtWeight_woHEM1516Fix
                            evtWeight_woHEM1516Fix_DataWithHEM1516Issue = evtWeight_woHEM1516Fix
                        else:
                            mask_DataPreHEM1516Issue  = trues_list
                            mask_DataWithHEM1516Issue = trues_list
                            evtWeight_DataPreHEM1516Issue               = evtWeight              * (1. - DataFractionAffectedBy2018HEM1516Issue)
                            evtWeight_DataWithHEM1516Issue              = evtWeight              * (DataFractionAffectedBy2018HEM1516Issue)
                            evtWeight_woHEM1516Fix_DataPreHEM1516Issue  = evtWeight_woHEM1516Fix * (1. - DataFractionAffectedBy2018HEM1516Issue)
                            evtWeight_woHEM1516Fix_DataWithHEM1516Issue = evtWeight_woHEM1516Fix * (DataFractionAffectedBy2018HEM1516Issue)

                        
                        # inclusive, DataPreHEM1516Issue
                        sel_tmp_ = sel_SR_forHExt & mask_DataPreHEM1516Issue
                        output['hLeadingFatJetPt_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataPreHEM1516Issue[ sel_tmp_ ]
                        )            
                        output['hLeadingFatJetEta_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataPreHEM1516Issue[ sel_tmp_ ]
                        )
                        output['hLeadingFatJetPhi_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Phi=(leadingFatJet.phi[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataPreHEM1516Issue[ sel_tmp_ ]
                        )
                        # inclusive, DataWithHEM1516Issue
                        sel_tmp_ = sel_SR_forHExt & mask_DataWithHEM1516Issue
                        output['hLeadingFatJetPt_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataWithHEM1516Issue[ sel_tmp_ ]
                        )            
                        output['hLeadingFatJetEta_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataWithHEM1516Issue[ sel_tmp_ ]
                        )
                        output['hLeadingFatJetPhi_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Phi=(leadingFatJet.phi[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataWithHEM1516Issue[ sel_tmp_ ]
                        )

                        # HEM1516Issue region
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi
                        output['hLeadingFatJetPt_HEM1516IssueEtaPhiCut'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight[ sel_tmp_ ]
                        ) 
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Phi
                        output['hLeadingFatJetEta_HEM1516IssuePhiCut'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta
                        output['hLeadingFatJetPhi_HEM1516IssueEtaCut'+sHExt].fill(
                            dataset=dataset,
                            Phi=(leadingFatJet.phi[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight[ sel_tmp_ ]
                        )
                        # HEM1516Issue region, DataPreHEM1516Issue
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi & mask_DataPreHEM1516Issue
                        output['hLeadingFatJetPt_HEM1516IssueEtaPhiCut_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataPreHEM1516Issue[ sel_tmp_ ]
                        ) 
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Phi & mask_DataPreHEM1516Issue
                        output['hLeadingFatJetEta_HEM1516IssuePhiCut_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataPreHEM1516Issue[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_DataPreHEM1516Issue
                        output['hLeadingFatJetPhi_HEM1516IssueEtaCut_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Phi=(leadingFatJet.phi[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataPreHEM1516Issue[ sel_tmp_ ]
                        )
                        # HEM1516Issue region, DataWithHEM1516Issue
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi & mask_DataWithHEM1516Issue
                        output['hLeadingFatJetPt_HEM1516IssueEtaPhiCut_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataWithHEM1516Issue[ sel_tmp_ ]
                        ) 
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Phi & mask_DataWithHEM1516Issue
                        output['hLeadingFatJetEta_HEM1516IssuePhiCut_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataWithHEM1516Issue[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_DataWithHEM1516Issue
                        output['hLeadingFatJetPhi_HEM1516IssueEtaCut_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Phi=(leadingFatJet.phi[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_DataWithHEM1516Issue[ sel_tmp_ ]
                        )



                        # wo HEM1516 fix ------------------------------------------------------------------------------------
                        sel_tmp_ = sel_SR_woSel2018HEM1516_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi
                        output['hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix[ sel_tmp_ ]
                        )  
                        sel_tmp_ =  sel_SR_woSel2018HEM1516_forHExt & mask_HEM1516Issue_Phi
                        output['hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_woSel2018HEM1516_forHExt & mask_HEM1516Issue_Eta
                        output['hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix'+sHExt].fill(
                            dataset=dataset,
                            Phi=(leadingFatJet.phi[ sel_SR_woSel2018HEM1516_forHExt & mask_HEM1516Issue_Eta ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix[ sel_tmp_ ]
                        )

                        # wo HEM1516 fix, DataPreHEM1516Issue
                        sel_tmp_ = sel_SR_woSel2018HEM1516_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi & mask_DataPreHEM1516Issue
                        output['hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataPreHEM1516Issue[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_woSel2018HEM1516_forHExt & mask_HEM1516Issue_Phi & mask_DataPreHEM1516Issue
                        output['hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataPreHEM1516Issue[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_woSel2018HEM1516_forHExt & mask_HEM1516Issue_Eta & mask_DataPreHEM1516Issue
                        output['hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Phi=(leadingFatJet.phi[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataPreHEM1516Issue[ sel_tmp_ ]
                        )
                        
                        # wo HEM1516 fix, DataWithHEM1516Issue
                        sel_tmp_ = sel_SR_woSel2018HEM1516_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi & mask_DataWithHEM1516Issue
                        output['hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataWithHEM1516Issue[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_woSel2018HEM1516_forHExt & mask_HEM1516Issue_Phi & mask_DataWithHEM1516Issue     
                        output['hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataWithHEM1516Issue[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_woSel2018HEM1516_forHExt & mask_HEM1516Issue_Eta & mask_DataWithHEM1516Issue
                        output['hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Phi=(leadingFatJet.phi[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataWithHEM1516Issue[ sel_tmp_ ]
                        )


                        # w/ HEM1516 fix in data, but w/o HEM1516_MC_Reweights ------------------------------
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi
                        output['hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix[ sel_tmp_ ]
                        )            
                        sel_tmp_ =  sel_SR_forHExt & mask_HEM1516Issue_Phi
                        output['hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta
                        output['hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt'+sHExt].fill(
                            dataset=dataset,
                            Phi=(leadingFatJet.phi[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix[ sel_tmp_ ]
                        )
                        # w/ HEM1516 fix in data (DataPreHEM1516Issue), but w/o HEM1516_MC_Reweights
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi & mask_DataPreHEM1516Issue
                        output['hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataPreHEM1516Issue[ sel_tmp_ ]
                        )            
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Phi & mask_DataPreHEM1516Issue
                        output['hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataPreHEM1516Issue[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_DataPreHEM1516Issue
                        output['hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt_DataPreHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Phi=(leadingFatJet.phi[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataPreHEM1516Issue[ sel_tmp_ ]
                        )
                        # w/ HEM1516 fix in data (DataWithHEM1516Issue), but w/o HEM1516_MC_Reweights
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi & mask_DataWithHEM1516Issue
                        output['hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataWithHEM1516Issue[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Phi & mask_DataWithHEM1516Issue 
                        output['hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataWithHEM1516Issue[ sel_tmp_ ]
                        )
                        sel_tmp_ = sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_DataWithHEM1516Issue
                        output['hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt_DataWithHEM1516Issue'+sHExt].fill(
                            dataset=dataset,
                            Phi=(leadingFatJet.phi[ sel_tmp_ ]),
                            systematic=syst,
                            weight=evtWeight_woHEM1516Fix_DataWithHEM1516Issue[ sel_tmp_ ]
                        )
                    # # 2018 HEM15/16 issue end ------------------
                    

                    output['hLeadingFatJetMass'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetMSoftDrop'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetId'+sHExt].fill(
                        dataset=dataset,
                        nObject=(leadingFatJet.jetId[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetBtagDeepB'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.btagDeepB[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetBtagDDBvLV2'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.btagDDBvLV2[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetBtagDDCvBV2'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.btagDDCvBV2[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetBtagHbb'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.btagHbb[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetDeepTagMD_H4qvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetDeepTagMD_HbbvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )            
                    output['hLeadingFatJetDeepTagMD_ZHbbvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetDeepTagMD_ZHccvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.deepTagMD_ZHccvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetDeepTagMD_ZbbvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.deepTagMD_ZbbvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetDeepTagMD_ZvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.deepTagMD_ZvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetDeepTagMD_bbvsLight'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetDeepTagMD_ccvsLight'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.deepTagMD_ccvsLight[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetDeepTag_H'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.deepTag_H[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetDeepTag_QCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.deepTag_QCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    ) 
                    output['hLeadingFatJetDeepTag_QCDothers'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.deepTag_QCDothers[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )

                    
                    output['hLeadingFatJetN2b1'+sHExt].fill(
                        dataset=dataset,
                        N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetN3b1'+sHExt].fill(
                        dataset=dataset,
                        N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    
                    output['hLeadingFatJetTau1'+sHExt].fill(
                        dataset=dataset,
                        TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetTau2'+sHExt].fill(
                        dataset=dataset,
                        TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetTau3'+sHExt].fill(
                        dataset=dataset,
                        TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetTau4'+sHExt].fill(
                        dataset=dataset,
                        TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )

                    output['hLeadingFatJetTau4by3'+sHExt].fill(
                        dataset=dataset,
                        TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetTau3by2'+sHExt].fill(
                        dataset=dataset,
                        TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetTau2by1'+sHExt].fill(
                        dataset=dataset,
                        TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    
                    
                    output['hLeadingFatJetNConstituents'+sHExt].fill(
                        dataset=dataset,
                        nObject50=(leadingFatJet.nConstituents[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    if self.datasetInfo['isMC']:
                        output['hLeadingFatJetNBHadrons'+sHExt].fill(
                            dataset=dataset,
                            nObject=(leadingFatJet.nBHadrons[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetNCHadrons'+sHExt].fill(
                            dataset=dataset,
                            nObject=(leadingFatJet.nCHadrons[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )     
                        if self.datasetInfo['isSignal']:
                            output['hLeadingFatJetNBHadronsFromHToAA'+sHExt].fill(
                                dataset=dataset,
                                nObject=(n_leadingFatJat_matched_genB_HToAATo4B[sel_SR_forHExt_woGenMatch]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt_woGenMatch]
                            )   
                            output['hLeadingFatJetNBHadrons_Sig'+sHExt].fill(
                                dataset=dataset,
                                nObject=(n_leadingFatJat_matched_genB[sel_SR_forHExt_woGenMatch]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt_woGenMatch]
                            )                           
                        
                    
                    output['hLeadingFatJetParticleNetMD_QCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.particleNetMD_QCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetParticleNetMD_Xbb'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.particleNetMD_Xbb[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetParticleNetMD_Xcc'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.particleNetMD_Xcc[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetParticleNetMD_Xqq'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.particleNetMD_Xqq[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetParticleNetMD_XbbOverQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetParticleNetMD_XccOverQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJetParticleNetMD_XccvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetParticleNetMD_XqqOverQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )


                    output['hLeadingFatJetParticleNet_H4qvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.particleNet_H4qvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetParticleNet_HbbvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.particleNet_HbbvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )

                    output['hLeadingFatJetParticleNet_HccvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.particleNet_HccvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetParticleNet_QCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.particleNet_QCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetParticleNet_mass'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJetZHbb_plus_Xbb'+sHExt].fill(
                        dataset=dataset,
                        MLScore2k=(leadingFatJetZHbb_plus_Xbb[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )

                    ## PNetMD Hto4b
                    if 'particleNetMD_Hto4b_Haa4b' in events.FatJet.fields:
                        output['hLeadingFatJetParticleNetMD_Hto4b_Haa01b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_Haa01b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_Haa2b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_Haa2b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_Haa3b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_Haa3b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_Haa4b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_Haa4b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_QCD0b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_QCD0b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_QCD1b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_QCD1b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_QCD2b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_QCD2b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_QCD3b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_QCD3b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_QCD4b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_QCD4b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_binaryLF_Haa4b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_binaryLF_Haa4b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_binaryLF_QCDlf'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_binaryLF_QCDlf[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_binary_Haa4b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_binary_Haa4b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_binary_QCD'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_binary_QCD[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        output['hLeadingFatJetParticleNetMD_Hto4b_Haa34b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_Haa34b_sum[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_binary_Haa4b_avg'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_binary_Haa4b_avg[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_Haa4b_avg'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_Haa4b_avg[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_QCD01234b'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_QCD01234b_sum[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_binary_QCD_avg'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_binary_QCD_avg[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_QCD_avg'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_QCD_avg[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        output['hLeadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD'+sHExt].fill(
                            dataset=dataset,
                            #MLScore1k=(leadingFatJet.particleNetMD_Hto4b_Haa4b[sel_SR_forHExt] / (leadingFatJet.particleNetMD_Hto4b_Haa4b[sel_SR_forHExt] + leadingFatJet_PNetMD_Hto4b_QCD01234b_sum[sel_SR_forHExt])),
                            MLScore1k=leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD[sel_SR_forHExt],
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        output['hLeadingFatJetParticleNetMD_Hto4b_Htoaa34bOverQCD'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_Haa34b_sum[sel_SR_forHExt] / (leadingFatJet_PNetMD_Hto4b_Haa34b_sum[sel_SR_forHExt] + leadingFatJet_PNetMD_Hto4b_QCD01234b_sum[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_binaryLF_Htoaa4bOverQCD'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_binaryLF_Haa4b[sel_SR_forHExt] / (leadingFatJet.particleNetMD_Hto4b_binaryLF_Haa4b[sel_SR_forHExt] + leadingFatJet.particleNetMD_Hto4b_binaryLF_QCDlf[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_binary_Htoaa4bOverQCD'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_binary_Haa4b[sel_SR_forHExt] / (leadingFatJet.particleNetMD_Hto4b_binary_Haa4b[sel_SR_forHExt] + leadingFatJet.particleNetMD_Hto4b_binary_QCD[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_binary_Htoaa4bOverQCD_avg'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_binary_Haa4b_avg[sel_SR_forHExt] / (leadingFatJet_PNetMD_Hto4b_binary_Haa4b_avg[sel_SR_forHExt] + leadingFatJet_PNetMD_Hto4b_binary_QCD_avg[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_avg'+sHExt].fill(
                            dataset=dataset,
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_Haa4b_avg[sel_SR_forHExt] / (leadingFatJet_PNetMD_Hto4b_Haa4b_avg[sel_SR_forHExt] + leadingFatJet_PNetMD_Hto4b_QCD_avg[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )


                        output['hLeadingFatJetParticleNet_massA_Hto4b_v0'+sHExt].fill(
                            dataset=dataset,
                            Mass1=(leadingFatJet.particleNet_massA_Hto4b_v0[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_v1'+sHExt].fill(
                            dataset=dataset,
                            Mass1=(leadingFatJet.particleNet_massA_Hto4b_v1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_v2'+sHExt].fill(
                            dataset=dataset,
                            Mass1=(leadingFatJet.particleNet_massA_Hto4b_v2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_v3'+sHExt].fill(
                            dataset=dataset,
                            Mass1=(leadingFatJet.particleNet_massA_Hto4b_v3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_v4'+sHExt].fill(
                            dataset=dataset,
                            Mass1=(leadingFatJet.particleNet_massA_Hto4b_v4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        output['hLeadingFatJetParticleNet_massA_Hto4b_avg_v01'+sHExt].fill(
                            dataset=dataset,
                            Mass1=calculateAverageOfArrays([
                                leadingFatJet.particleNet_massA_Hto4b_v0[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v1[sel_SR_forHExt],
                                ]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_avg_v02'+sHExt].fill(
                            dataset=dataset,
                            Mass1=calculateAverageOfArrays([
                                leadingFatJet.particleNet_massA_Hto4b_v0[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v2[sel_SR_forHExt],
                                ]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_avg_v03'+sHExt].fill(
                            dataset=dataset,
                            Mass1=calculateAverageOfArrays([
                                leadingFatJet.particleNet_massA_Hto4b_v0[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v3[sel_SR_forHExt],
                                ]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_avg_v12'+sHExt].fill(
                            dataset=dataset,
                            Mass1=calculateAverageOfArrays([
                                leadingFatJet.particleNet_massA_Hto4b_v1[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v2[sel_SR_forHExt],
                                ]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_avg_v13'+sHExt].fill(
                            dataset=dataset,
                            Mass1=calculateAverageOfArrays([
                                leadingFatJet.particleNet_massA_Hto4b_v1[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v3[sel_SR_forHExt],
                                ]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_avg_v23'+sHExt].fill(
                            dataset=dataset,
                            Mass1=calculateAverageOfArrays([
                                leadingFatJet.particleNet_massA_Hto4b_v2[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v3[sel_SR_forHExt],
                                ]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        output['hLeadingFatJetParticleNet_massA_Hto4b_avg_v012'+sHExt].fill(
                            dataset=dataset,
                            Mass1=calculateAverageOfArrays([
                                leadingFatJet.particleNet_massA_Hto4b_v0[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v1[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v2[sel_SR_forHExt],
                                ]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_avg_v013'+sHExt].fill(
                            dataset=dataset,
                            Mass1=calculateAverageOfArrays([
                                leadingFatJet.particleNet_massA_Hto4b_v0[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v1[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v3[sel_SR_forHExt],
                                ]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_avg_v023'+sHExt].fill(
                            dataset=dataset,
                            Mass1=calculateAverageOfArrays([
                                leadingFatJet.particleNet_massA_Hto4b_v0[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v2[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v3[sel_SR_forHExt],
                                ]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massA_Hto4b_avg_v123'+sHExt].fill(
                            dataset=dataset,
                            Mass1=calculateAverageOfArrays([
                                leadingFatJet.particleNet_massA_Hto4b_v1[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v2[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v3[sel_SR_forHExt],
                                ]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        output['hLeadingFatJetParticleNet_massA_Hto4b_avg_v0123'+sHExt].fill(
                            dataset=dataset,
                            Mass1=calculateAverageOfArrays([
                                leadingFatJet.particleNet_massA_Hto4b_v0[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v1[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v2[sel_SR_forHExt],
                                leadingFatJet.particleNet_massA_Hto4b_v3[sel_SR_forHExt],
                                ]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                                                

                        output['hLeadingFatJetParticleNet_massH_Hto4b_v0'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_massH_Hto4b_v0[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massH_Hto4b_v00'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_massH_Hto4b_v00[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massH_Hto4b_v1'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_massH_Hto4b_v1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massH_Hto4b_v2'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_massH_Hto4b_v2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massH_Hto4b_v3'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_massH_Hto4b_v3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_massH_Hto4b_v4'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_massH_Hto4b_v4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        # Scale particleNet_massH_Hto4b_v0 by 1.01 to get better response
                        # https://indico.cern.ch/event/1343334/contributions/5655252/attachments/2745224/4781382/2023_11_02_HToAATo4B_Higgs_mass_studies.pdf#page=15
                        output['hLeadingFatJetParticleNet_massH_Hto4b_avg_v0123'+sHExt].fill(
                            dataset=dataset,
                            Mass=leadingFatJet_PNet_massH_Hto4b_avg[sel_SR_forHExt],
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                    ## SubJet corresponding to leadingFatJet
                    output['hLeadingFatJet_nSubJets'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(leadingFatJet_nSubJets[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )                    
                    output['hLeadingFatJet_nSubJets_bTag_L'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )                    
                    output['hLeadingFatJet_nSubJets_bTag_M'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )                    


                    ## SV
                    output['hLeadingFatJet_nSV'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJet_mass_SV_MaxdxySig'+sHExt].fill(
                        dataset=dataset,
                        Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingFatJet_logMass_SV_MaxdxySig'+sHExt].fill(
                        dataset=dataset,
                        logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )


                    ## MET
                    output['hMET_pT'+sHExt].fill(
                        dataset=dataset,
                        Pt=(events.MET.pt[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hMET_sumEt'+sHExt].fill(
                        dataset=dataset,
                        Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )


                    ## nLeptons 
                    output['hLeadingFatJet_nLeptons'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nLeptons_matched_leadingFatJet[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hnLeptonsTight'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nLeptonsTight[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hnLeptons_nonoverlap_leadingFatJet'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nLeptons_nonoverlap_leadingFatJet[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hnLeptons_nonoverlap_selFatJets'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nLeptons_nonoverlap_selFatJets[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )


                    ## leadingNonHto4bFatJet                  
                    output['hLeadingNonHto4bFatJetPt'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingNonHto4bFatJet.pt[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )            
                    output['hLeadingNonHto4bFatJetEta'+sHExt].fill(
                        dataset=dataset,
                        Eta=(leadingNonHto4bFatJet.eta[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetPhi'+sHExt].fill(
                        dataset=dataset,
                        Phi=(leadingNonHto4bFatJet.phi[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetMass'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingNonHto4bFatJet.mass[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetMSoftDrop'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingNonHto4bFatJet.msoftdrop[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetId'+sHExt].fill(
                        dataset=dataset,
                        nObject=(leadingNonHto4bFatJet.jetId[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )                                           
                    output['hdPhi_LeadingFJ_LeadingNonHto4bFJ'+sHExt].fill(
                        dataset=dataset,
                        deltaPhi=abs(leadingFatJet.delta_phi(leadingNonHto4bFatJet)[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetDeepTagMD_WvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet.deepTagMD_WvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetDeepTagMD_ZvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet.deepTagMD_ZvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )                    
                    output['hLeadingNonHto4bFatJetDeepTag_WvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet.deepTag_WvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetDeepTag_ZvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet.deepTag_ZvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetDeepTag_VvsQCD_max'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_DeepTag_VvsQCD_max[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetDeepTag_W'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_DeepTag_W[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetDeepTag_Z'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_DeepTag_Z[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetDeepTag_V_max'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_DeepTag_V_max[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetDeepTag_WZvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_DeepTag_WZvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetDeepTag_WZvsQCD2'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_DeepTag_WZvsQCD2[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetParticleNet_WvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet.particleNet_WvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetParticleNet_ZvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet.particleNet_ZvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetPNet_VvsQCD_max'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_PNet_VvsQCD_max[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetPNet_W'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_PNet_W[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetPNet_Z'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_PNet_Z[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetPNet_V_max'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_PNet_V_max[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetPNet_WZvsQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_PNet_WZvsQCD[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    output['hLeadingNonHto4bFatJetPNet_WZvsQCD2'+sHExt].fill(
                        dataset=dataset,
                        MLScore1k=(leadingNonHto4bFatJet_PNet_WZvsQCD2[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )
                    

                    # ak4 jets ----
                    output['hnAK4Jets_NonoverlapLeadingFatJet'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nAk4Jets_nonoverlaping_leadingFatJet[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )    
                    if ak.count(ak.firsts(ak4Jets_nonoverlaping_leadingFatJet).pt[sel_SR_forHExt]) > 0:
                        output['hPtLeadingAK4Jets_NonoverlapLeadingFatJet'+sHExt].fill(
                            dataset=dataset,
                            Pt=(ak.firsts(ak4Jets_nonoverlaping_leadingFatJet).pt[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4Jets_nonoverlaping_leadingFatJet).pt))
                                ]),
                            systematic=syst,
                            weight=evtWeight[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4Jets_nonoverlaping_leadingFatJet).pt))
                                ]
                        )
                    output['hnAK4Jets_bTag_NonoverlapLeadingFatJet'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nAk4Jets_bTag_nonoverlaping_leadingFatJet[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )    
                    if ak.count(ak.firsts(ak4Jets_bTag_nonoverlaping_leadingFatJet).pt[sel_SR_forHExt]) > 0:
                        output['hPtLeadingAK4Jets_bTag_NonoverlapLeadingFatJet'+sHExt].fill(
                            dataset=dataset,
                            Pt=(ak.firsts(ak4Jets_bTag_nonoverlaping_leadingFatJet).pt[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4Jets_bTag_nonoverlaping_leadingFatJet).pt))
                                ]),
                            systematic=syst,
                            weight=evtWeight[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4Jets_bTag_nonoverlaping_leadingFatJet).pt))
                                ]
                        )

                    output['hnAK4JetsCentral_NonoverlapLeadingFatJet'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nAk4JetsCentral_nonoverlaping_leadingFatJet[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )    
                    if ak.count(ak.firsts(ak4JetsCentral_nonoverlaping_leadingFatJet).pt[sel_SR_forHExt]) > 0:
                        output['hPtLeadingAK4JetsCentral_NonoverlapLeadingFatJet'+sHExt].fill(
                            dataset=dataset,
                            Pt=(ak.firsts(ak4JetsCentral_nonoverlaping_leadingFatJet).pt[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4JetsCentral_nonoverlaping_leadingFatJet).pt))
                                ]),
                            systematic=syst,
                            weight=evtWeight[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4JetsCentral_nonoverlaping_leadingFatJet).pt))
                                ]
                        )
                    output['hnAK4JetsCentral_bTag_NonoverlapLeadingFatJet'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nAk4JetsCentral_bTag_nonoverlaping_leadingFatJet[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )    
                    if ak.count(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_leadingFatJet).pt[sel_SR_forHExt]) > 0:
                        output['hPtLeadingAK4JetsCentral_bTag_NonoverlapLeadingFatJet'+sHExt].fill(
                            dataset=dataset,
                            Pt=(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_leadingFatJet).pt[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_leadingFatJet).pt))
                                ]),
                            systematic=syst,
                            weight=evtWeight[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_leadingFatJet).pt))
                                ]
                        )



                    output['hnAK4Jets_NonoverlapSelFatJets'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nAk4Jets_nonoverlaping_selFatJets[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )    
                    if ak.count(ak.firsts(ak4Jets_nonoverlaping_selFatJets).pt[sel_SR_forHExt]) > 0:
                        output['hPtLeadingAK4Jets_NonoverlapSelFatJets'+sHExt].fill(
                            dataset=dataset,
                            Pt=(ak.firsts(ak4Jets_nonoverlaping_selFatJets).pt[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4Jets_nonoverlaping_selFatJets).pt))
                                ]),
                            systematic=syst,
                            weight=evtWeight[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4Jets_nonoverlaping_selFatJets).pt))
                                ]
                        )
                    output['hnAK4Jets_bTag_NonoverlapSelFatJets'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nAk4Jets_bTag_nonoverlaping_selFatJets[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )    
                    if ak.count(ak.firsts(ak4Jets_bTag_nonoverlaping_selFatJets).pt[sel_SR_forHExt]) > 0:
                        output['hPtLeadingAK4Jets_bTag_NonoverlapSelFatJets'+sHExt].fill(
                            dataset=dataset,
                            Pt=(ak.firsts(ak4Jets_bTag_nonoverlaping_selFatJets).pt[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4Jets_bTag_nonoverlaping_selFatJets).pt))
                                ]),
                            systematic=syst,
                            weight=evtWeight[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4Jets_bTag_nonoverlaping_selFatJets).pt))
                                ]
                        )

                    output['hnAK4JetsCentral_NonoverlapSelFatJets'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nAk4JetsCentral_nonoverlaping_selFatJets[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )    
                    if ak.count(ak.firsts(ak4JetsCentral_nonoverlaping_selFatJets).pt[sel_SR_forHExt]) > 0:
                        output['hPtLeadingAK4JetsCentral_NonoverlapSelFatJets'+sHExt].fill(
                            dataset=dataset,
                            Pt=(ak.firsts(ak4JetsCentral_nonoverlaping_selFatJets).pt[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4JetsCentral_nonoverlaping_selFatJets).pt))
                                ]),
                            systematic=syst,
                            weight=evtWeight[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4JetsCentral_nonoverlaping_selFatJets).pt))
                                ]
                        )
                    output['hnAK4JetsCentral_bTag_NonoverlapSelFatJets'+sHExt].fill(
                        dataset=dataset,
                        nObject10=(nAk4JetsCentral_bTag_nonoverlaping_selFatJets[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )    
                    if ak.count(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_selFatJets).pt[sel_SR_forHExt]) > 0:
                        output['hPtLeadingAK4JetsCentral_bTag_NonoverlapSelFatJets'+sHExt].fill(
                            dataset=dataset,
                            Pt=(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_selFatJets).pt[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_selFatJets).pt))
                                ]),
                            systematic=syst,
                            weight=evtWeight[
                                sel_SR_forHExt & 
                                (~ak.is_none(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_selFatJets).pt))
                                ]
                        )



                    ### 2-D distribution ----------------------------------------------------------

                    output['hLeadingFatJetParticleNet_massH_Hto4b_avg_vs_massA_Hto4b_avg'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingFatJet_PNet_massH_Hto4b_avg[sel_SR_forHExt]),
                        Mass2=(leadingFatJet_PNet_massA_Hto4b_avg[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    )                    
                    output['hLeadingFatJetMass_vs_massA_Hto4b_avg'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                        Mass2=(leadingFatJet_PNet_massA_Hto4b_avg[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    ) 
                    output['hLeadingFatJetMSoftDrop_vs_massA_Hto4b_avg'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                        Mass2=(leadingFatJet_PNet_massA_Hto4b_avg[sel_SR_forHExt]),
                        systematic=syst,
                        weight=evtWeight[sel_SR_forHExt]
                    ) 
                    
                    

                    
                    



                    ### 2-D distribution ----------------------------------------------------------

                    if runMode_SignificancsScan2D:

                        ## 2-D hLeadingFatJetDeepTagMD_H4qvsQCD
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetDeepTagMD_HbbvsQCD'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            MLScore1=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetDeepTagMD_ZHbbvsQCD'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            MLScore1=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJetDeepTagMD_bbvsLight'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            MLScore1=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_nSubJets'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_nSV'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_H4qvsQCD_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )


                        ## 2-D hLeadingFatJetDeepTagMD_HbbvsQCD
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetDeepTagMD_ZHbbvsQCD'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            MLScore1=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJetDeepTagMD_bbvsLight'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            MLScore1=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_nSubJets'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_nSV'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_HbbvsQCD_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )



                        ## 2-D hLeadingFatJetDeepTagMD_ZHbbvsQCD
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJetDeepTagMD_bbvsLight'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            MLScore1=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_nSubJets'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_nSV'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_ZHbbvsQCD_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )


                        ## 2-D hLeadingFatJetDeepTagMD_bbvsLight
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_nSubJets'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_nSV'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetDeepTagMD_bbvsLight_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )


                        ## 2-D hLeadingFatJetParticleNetMD_XbbOverQCD
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetParticleNetMD_XqqOverQCD'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            MLScore1=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJetParticleNet_mass'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_nSubJets'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_nSV'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XbbOverQCD_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )


                        ## 2-D hLeadingFatJetParticleNetMD_XqqOverQCD
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJetParticleNet_mass'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_nSubJets'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_nSV'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNetMD_XqqOverQCD_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            MLScore=(leadingFatJetParticleNetMD_XqqvsQCD[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )


                        ## 2-D hLeadingFatJetParticleNet_mass
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            Mass1=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            Mass1=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_nSubJets'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_nSV'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetParticleNet_mass_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.particleNet_mass[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )


                        ## 2-D hLeadingFatJet_nSubJets
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJet_nSubJets_bTag_L'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            nObject10_1=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJet_nSubJets_bTag_M'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            nObject10_1=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJet_nSV'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            nObject10_1=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )


                        ## 2-D hLeadingFatJet_nSubJets_bTag_L
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJet_nSV'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            nObject10_1=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_L_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_L[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )



                        ## 2-D hLeadingFatJet_nSubJets_bTag_M
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJet_nSV'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            nObject10_1=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJet_mass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_LeadingFatJet_logMass_SV_MaxdxySig'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSubJets_bTag_M_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(leadingFatJet_nSubJets_bTag_M[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )



                        ## 2-D hLeadingFatJet_nSV
                        output['hLeadingFatJet_nSV_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_nSV_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nSV_matched_leadingFatJet[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )


                        ## 2-D hLeadingFatJet_mass_SV_MaxdxySig
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_mass_SV_MaxdxySig_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            Mass10=(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )



                        ## 2-D hLeadingFatJet_nSV
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_MET_pT'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJet_logMass_SV_MaxdxySig_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            logMass3=np.log(mass_SV_matched_leadingFatJet_MaxdxySig[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )



                        ## 2-D hMET_pT_
                        output['hMET_pT_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT_vs_MET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.MET.pt[sel_SR_forHExt]),
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )


                        ## 2-D hMET_sumEt
                        output['hMET_sumEt_vs_LeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            Mass=(leadingFatJet.mass[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_sumEt_vs_LeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            Mass=(leadingFatJet.msoftdrop[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_sumEt_vs_LeadingFatJetN2b1'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            N2=(leadingFatJet.n2b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_sumEt_vs_LeadingFatJetN3b1'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            N3=(leadingFatJet.n3b1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_sumEt_vs_LeadingFatJetTau1'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau1[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_sumEt_vs_LeadingFatJetTau2'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau2[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_sumEt_vs_LeadingFatJetTau3'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau3[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_sumEt_vs_LeadingFatJetTau4'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            TauN=(leadingFatJet.tau4[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_sumEt_vs_LeadingFatJetTau4by3'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau4[sel_SR_forHExt], leadingFatJet.tau3[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_sumEt_vs_LeadingFatJetTau3by2'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau3[sel_SR_forHExt], leadingFatJet.tau2[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_sumEt_vs_LeadingFatJetTau2by1'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.MET.sumEt[sel_SR_forHExt]),
                            TauN=(np.divide(leadingFatJet.tau2[sel_SR_forHExt], leadingFatJet.tau1[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        ## FatJet pT vs MLScores
                        output['hLeadingFatJetPt_vs_PNetMD_Hto4b_Haa3b'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_Haa3b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetPt_vs_PNetMD_Hto4b_Haa4b'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_Haa4b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetPt_vs_PNetMD_Hto4b_binaryLF_Haa4b'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_binaryLF_Haa4b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetPt_vs_PNetMD_Hto4b_binary_Haa4b'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore1k=(leadingFatJet.particleNetMD_Hto4b_binary_Haa4b[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetPt_vs_PNetMD_Hto4b_Haa34b'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_Haa34b_sum[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetPt_vs_PNetMD_Hto4b_binary_Haa4b_avg'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_binary_Haa4b_avg[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetPt_vs_PNetMD_Hto4b_Haa4b_avg'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_Haa4b_avg[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetPt_vs_PNetMD_Hto4b_Htoaa4bOverQCD'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_Htoaa4bOverQCD[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetPt_vs_PNetMD_Hto4b_Htoaa34bOverQCD'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_Haa34b_sum[sel_SR_forHExt] / (leadingFatJet_PNetMD_Hto4b_Haa34b_sum[sel_SR_forHExt] + leadingFatJet_PNetMD_Hto4b_QCD01234b_sum[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetPt_vs_PNetMD_Hto4b_binary_Htoaa4bOverQCD_avg'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_binary_Haa4b_avg[sel_SR_forHExt] / (leadingFatJet_PNetMD_Hto4b_binary_Haa4b_avg[sel_SR_forHExt] + leadingFatJet_PNetMD_Hto4b_binary_QCD_avg[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetPt_vs_PNetMD_Hto4b_Htoaa4bOverQCD_avg'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore1k=(leadingFatJet_PNetMD_Hto4b_Haa4b_avg[sel_SR_forHExt] / (leadingFatJet_PNetMD_Hto4b_Haa4b_avg[sel_SR_forHExt] + leadingFatJet_PNetMD_Hto4b_QCD_avg[sel_SR_forHExt])),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetPt_vs_ZHbb_plus_Xbb'+sHExt].fill(
                            dataset=dataset,
                            Pt=(leadingFatJet.pt[sel_SR_forHExt]),
                            MLScore2k=(leadingFatJetZHbb_plus_Xbb[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        

        



















































            '''
            output[''].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            '''

        return output


    def postprocess(self, accumulator):
        #pass
        return accumulator



    
def getLorentzVectorFromAwkArray(awkArray, ptObjectName, etaObjectName, phiObjectName, massObjectName):
    '''
    return ak.Array(
        awkArray,
        with_name='Momentum4D'
    )
    '''
    v1 = ak.zip( {
        'pt':   awkArray[ptObjectName],
        'eta':  awkArray[etaObjectName],
        'phi':  awkArray[phiObjectName],
        'mass': awkArray[massObjectName],
    })
    if printLevel >= 12:
        print(f"getLorentzVectorFromAwkArray(): v1 ({type(v1)}): {v1.to_list()}")

    return ak.Array(v1, with_name='Momentum4D')
                       

def printWithType(sX, X):
    #print(f"{sX} ({type(X)}): {X}")
    print(f"{sX} : {X}")


        
    
    
if __name__ == '__main__':
    print("htoaa_Analysis:: main: {}".format(sys.argv)); sys.stdout.flush()
    print(f"htoaa_Analysis_VHHadronicMode:: here14 {datetime.now() = }")

    if len(sys.argv) != 2:
        print("htoaa_Analysis:: Command-line config file missing.. \t **** ERROR **** \n")

    sConfig = sys.argv[1]

    
    config = GetDictFromJsonFile(sConfig)
    print("Config {}: \n{}".format(sConfig, json.dumps(config, indent=4)))
    print(f"htoaa_Analysis_VHHadronicMode:: here15 {datetime.now() = }")

    nEventsToAnalyze    = config["nEventsToAnalyze"] if "nEventsToAnalyze" in config else nEventToReadInBatch
    sInputFiles         = config["inputFiles"]
    sOutputFile         = config["outputFile"]
    sample_dataset      = config["dataset"] 
    sample_category     = config['sampleCategory']
    isMC                = config["isMC"]
    era                 = config['era']
    downloadIpFiles     = config['downloadIpFiles'] if 'downloadIpFiles' in config else False
    server              = config["server"]
    if isMC:
        sample_crossSection = config["crossSection"]
        sample_nEvents      = config["nEvents"]
        sample_sumEvents    = config["sumEvents"] if config["sumEvents"] > 0 else sample_nEvents
        if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
        
        MCSamplesStitchOption = MCSamplesStitchOptions.PhSpOverlapRewgt if ("MCSamplesStitchOption" in config and \
                                                                            config["MCSamplesStitchOption"] == MCSamplesStitchOptions.PhSpOverlapRewgt.value) \
            else MCSamplesStitchOptions.PhSpOverlapRemove
        
        if MCSamplesStitchOption == MCSamplesStitchOptions.PhSpOverlapRewgt:
            if "MCSamplesStitchInputs" not in config:
                print(frameinfo.filename, frameinfo.lineno, ' ERROR: "MCSamplesStitchInputs" not in config') # https://stackoverflow.com/questions/3056048/filename-and-line-number-of-python-script
                
            MCSamplesStitchInputFileName      = config["MCSamplesStitchInputs"]["inputFile"]
            MCSamplesStitchInputHistogramName = config["MCSamplesStitchInputs"]["histogramName"]
            MCSamplesStitchInputHistogramName = MCSamplesStitchInputHistogramName.replace(
                '$SAMPLECATEGORY', sample_category.split('_')[0]
            )
            print(f"{MCSamplesStitchOption = }, {MCSamplesStitchInputFileName = }, {MCSamplesStitchInputHistogramName = } ")
            if not os.path.exists(MCSamplesStitchInputFileName):
                logging.critical(f'htoaa_Analysis_VHHadronicMode.py::main():: {MCSamplesStitchInputFileName = } does not exists')
                print(f'htoaa_Analysis_VHHadronicMode.py::main() 11:: {MCSamplesStitchInputFileName = } does not exists')
                exit(0)
            print(f"Opening {MCSamplesStitchInputFileName = } "); sys.stdout.flush() 
            with uproot.open(MCSamplesStitchInputFileName) as f_:
                print(f"{f_.keys() = }"); sys.stdout.flush() 
                hMCSamplesStitch = f_[r'%s' % MCSamplesStitchInputHistogramName].to_hist()

    print(f"htoaa_Analysis_VHHadronicMode:: here16 {datetime.now() = }")    
        
        
    #branchesToRead = htoaa_nanoAODBranchesToRead
    #print("branchesToRead: {}".format(branchesToRead))
    sample_dataset = sample_dataset[0] if isinstance(sample_dataset, list) else sample_dataset # dataset is list of datasets w/ same sample name, as they are considered together recently. Those set of datasets are extension of the same samples.

    sInputFiles_toUse = []
    for sInputFile in sInputFiles:
        if "*" in sInputFile:  sInputFiles_toUse.extend( glob.glob( sInputFile ) )
        else:                  sInputFiles_toUse.append( sInputFile )
    sInputFiles = sInputFiles_toUse
    print(f"Initial sInputFiles ({len(sInputFiles)}) (type {type(sInputFiles)}):");
    for sInputFile in sInputFiles:
        print(f"\t{sInputFile}");  sys.stdout.flush()
    print(f"htoaa_Analysis_VHHadronicMode:: here17 {datetime.now() = }")

    for iFile in range(len(sInputFiles)):     
        sInputFile = sInputFiles[iFile]
        sFileLocal = './inputFiles/%s' %(os.path.basename(sInputFile))  
        print(f"  {sInputFile = },  {sFileLocal = }", flush=flushStdout)
        sInputFile, isReadingSuccessful = getNanoAODFile(
            fileName = sInputFile, 
            useLocalFileIfExists = True, 
            downloadFile = downloadIpFiles, 
            fileNameLocal = sFileLocal, 
            nTriesToDownload = 3,
            server = server
            )
        if not isReadingSuccessful:
            logging.critical('htoaa_Analysis_VHHadronicMode:: getNanoAODFile() for input file %s failed. **** CRITICAL ERROR ****. \nAborting...' % (sInputFile)); sys.stdout.flush();
            exit(0)
        
        # Check if input file exists or not
        fileSize = 0
        if os.path.exists(sInputFile):
            try:
                fileSize = os.path.getsize(sInputFile) / (1024 * 1024) # file size in MB
                #print(f"sInputFile: {sInputFile} ({fileSize} MB) ")
            except  FileNotFoundError:
                print(f"sInputFile: {sInputFile} file not found.")
            except OSError: 
                print(f"sInputFile: {sInputFile} OS error occurred.")
        print(f"htoaa_Analysis_VHHadronicMode:: {sInputFile} \t {os.path.exists(sInputFile) = }, {fileSize = } MB");     

        if fileSize > NanoAODFileSize_Min:     
            sInputFiles[iFile] = sInputFile
        else:
            logging.critical('htoaa_Analysis_VHHadronicMode:: Input file %s file size below threshold (%g MB). **** CRITICAL ERROR ****. \nAborting...' % (sInputFile, NanoAODFileSize_Min) ); sys.stdout.flush();
            exit(0)
    
        
    print(f"\nActual  sInputFiles ({len(sInputFiles)}) (type {type(sInputFiles)}):");
    for sInputFile in sInputFiles:
        fileSize = 0
        if os.path.exists(sInputFile):
            try:
                fileSize = os.path.getsize(sInputFile) / (1024 * 1024) # file size in MB
                #print(f"sInputFile: {sInputFile} ({fileSize} MB) ")
            except  FileNotFoundError:
                print(f"sInputFile: {sInputFile} file not found.")
            except OSError: 
                print(f"sInputFile: {sInputFile} OS error occurred.")
        print(f"\t{sInputFile} \t {os.path.exists(sInputFile) = }, {fileSize = } MB");  
    sys.stdout.flush()
    print(f"htoaa_Analysis_VHHadronicMode:: here18 {datetime.now() = }")


    sampleInfo = {
        "era":             era, 
        "isMC":            isMC,
        "sample_category": sample_category,        
        "datasetNameFull": sample_dataset,
    }
    if isMC:
        sampleInfo["sample_crossSection"]   = sample_crossSection
        sampleInfo["sample_sumEvents"]      = sample_sumEvents
        sampleInfo["MCSamplesStitchOption"] = MCSamplesStitchOption
        if MCSamplesStitchOption == MCSamplesStitchOptions.PhSpOverlapRewgt:
            sampleInfo["hMCSamplesStitch"] = hMCSamplesStitch
    print(f"htoaa_Analysis_VHHadronicMode:: here19 {datetime.now() = }", flush=flushStdout)
        
    startTime = time.time()
    tracemalloc.start()


    #client = Client("tls://localhost:8786")
    #executor = processor.DaskExecutor(client=client)
    chunksize = nEventToReadInBatch
    maxchunks = None if nEventsToAnalyze == -1 else int(nEventsToAnalyze/nEventToReadInBatch)
    nWorkers  = 4 if nEventsToAnalyze == -1 else 1
    print(f"nEventsToAnalyze: {nEventsToAnalyze},  nEventToReadInBatch: {nEventToReadInBatch}, chunksize: {chunksize},  maxchunks: {maxchunks},  nWorkers: {nWorkers}", flush=flushStdout)
    run = processor.Runner(
        #executor=executor,
        executor=processor.FuturesExecutor(workers=nWorkers),
        schema=schemas.NanoAODSchema,
        savemetrics=True,
        chunksize=chunksize,  #3 ** 20,  ## Governs the number of times LeptonJetProcessor "process" is called
        maxchunks=maxchunks
    )
    
    output, metrics = run(
        fileset={sample_category: sInputFiles},
        #fileset={"QCD": ["/home/siddhesh/Work/CMS/htoaa/analysis/tmp/20BE2B12-EFF6-8645-AB7F-AFF6A624F816.root"]},
        treename="Events",
        processor_instance=HToAATo4bProcessor(
            datasetInfo=sampleInfo
        )
    )
    print(f"metrics: {metrics}", flush=flushStdout)


    if 'cutflow' in output.keys():
        print("Cutflow::", flush=flushStdout)
        #for key, value in output['cutflow'].items():
        for key in output['cutflow'].keys():
            #print(key, value)
            if key.startswith(sWeighted): continue # to print weighted and unweighted events for cuts on the same line

            print("%10f\t%10d\t%s" % (output['cutflow'][sWeighted+key], output['cutflow'][key], key), flush=flushStdout)
    
    
    if sOutputFile is not None:
        if not sOutputFile.endswith('.root'): sOutputFile += '.root'
        #sOutputFile = sOutputFile.replace('.root', '_wCoffea.root') # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        sample_category_toUse = sample_category
        
        if isMC and \
            MCSamplesStitchOption == MCSamplesStitchOptions.PhSpOverlapRewgt and \
            SplitQCDInGENCats and \
            "QCD" in sample_category:
            sample_category_toUse = "QCD"
        
        sDir1 = 'evt/%s' % (sample_category_toUse)

        
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                #print(f"key: {key},  value ({type(value)}): {value}")
                sHistoName_toUse = key
                sHExt_toUse = ''
                if isMC and \
                    MCSamplesStitchOption == MCSamplesStitchOptions.PhSpOverlapRewgt and \
                    SplitQCDInGENCats and \
                    "QCD" in sample_category:                    
                    for sHExt in HistogramNameExtensions_QCD:
                        if sHExt in key:
                            sHExt_toUse = '_%s' % (sHExt)
                            sHistoName_toUse = sHistoName_toUse.replace(sHExt_toUse, '')
                            break

                sDir1_toUse = '%s%s' % (sDir1, sHExt_toUse)

                #if not (key.startswith('h') or key != 'cutflow'): continue
                if not isinstance(value, hist.Hist): continue
                #print(f"1: key {key}, value ({type(value)})     Hist: {type(hist.Hist)},    isinstance(value, hist.Hist): {isinstance(value, hist.Hist)}") # value: {value}")

                '''
                print(f"value.DEFAULT_DTYPE {value.DEFAULT_DTYPE}")
                print(f"value.fields ({type(value.fields)}): {value.fields}")
                print(f"value.label ({type(value.label)}): {value.label}")
                print(f"value.axes() ({type(value.axes())}): {value.axes()}")
                print(f"value.sparse_axes() ({type(value.sparse_axes())}): {value.sparse_axes()}")
                print(f"value.sparse_dim() {value.sparse_dim()}")
                for ax_label in value.fields:
                    print(f"ax_label ({type(ax_label)}): {ax_label}")
                    print(f"value.axis(ax_label) ({type(value.axis(ax_label))}): {value.axis(ax_label)}")
                    #printWithType('value.axis(ax_label).identifiers()', value.axis(ax_label).identifiers())
                    print(f"value.axis(ax_label).identifiers()  : {value.axis(ax_label).identifiers()}")

                for ax in value.sparse_axes():
                    printWithType('ax', ax)

                h1 = value.to_hist()
                print(f"h1 ({type(h1)}): h1")
                '''


                for _dataset in value.axis('dataset').identifiers():
                    #print(f"_dataset ({type(_dataset)}): {_dataset}")

                    for _syst in value.axis('systematic').identifiers():
                        #print(f"_syst ({type(_syst)}): {_syst}")

                        h1 = value.integrate('dataset',_dataset).integrate('systematic',_syst).to_hist()
                        #print(f"h1 ({type(h1)}): h1")

                        #fOut['%s/%s_%s_%s' % (sDir1, key, _dataset, _syst)] = h1
                        fOut['%s/%s_%s' % (sDir1_toUse, sHistoName_toUse, _syst)] = h1
                
                #fOut['%s%s' % (sDir1, key)] = value
                #fOut['%s%s' % (sDir1, key)] = hist.export1d(value)
                #fOut['%s%s' % (sDir1, key)] = h1

                
        

        #util.save(output, sOutputFile)
            
        
        print("Wrote to sOutputFile {}".format(sOutputFile), flush=flushStdout)
    















    
    current_memory, peak_memory = tracemalloc.get_traced_memory() # https://medium.com/survata-engineering-blog/monitoring-memory-usage-of-a-running-python-program-49f027e3d1ba
    print(f"\n\nMemory usage:: current {current_memory / 10**6}MB;  peak {peak_memory / 10**6}MB", flush=flushStdout)

    endTime = time.time()
    totalTime = endTime - startTime
    totalTime_hr  = int(totalTime/60/60)
    totalTime_min = totalTime - float(totalTime_hr * 60)
    totalTime_min = int(totalTime_min/60)
    totalTime_sec = totalTime - float(totalTime_hr * 60*60) - float(totalTime_min * 60)
    print(f"Total run time: {totalTime_hr}h {totalTime_min}m {totalTime_sec}s = {totalTime}sec ", flush=flushStdout)
    
