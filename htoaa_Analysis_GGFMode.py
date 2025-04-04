#htoaa analysis main code

import os
import sys
from datetime import datetime
#import time
print(f"htoaa_Analysis_GGFMode:: here1 {datetime.now() = }"); sys.stdout.flush()
import subprocess
import json
from urllib.request import urlopen
import glob
from collections import OrderedDict as OD
import time
import tracemalloc
import math
print(f"htoaa_Analysis_GGFMode:: here2 {datetime.now() = }"); sys.stdout.flush()
import numpy as np
from copy import copy, deepcopy
print(f"htoaa_Analysis_GGFMode:: here3 {datetime.now() = }"); sys.stdout.flush()
#import uproot
#import uproot3 as uproot
import uproot as uproot
print(f"htoaa_Analysis_GGFMode:: here4 {datetime.now() = }"); sys.stdout.flush()
#import parse
from parse import *
print(f"htoaa_Analysis_GGFMode:: here4.1 {datetime.now() = }"); sys.stdout.flush()
import logging
print(f"htoaa_Analysis_GGFMode:: here5 {datetime.now() = }"); sys.stdout.flush()

# comment test3
'''
GGF -> H->aa->4b boosted analysis macro

References:
  * Coffea framework used for TTGamma analysis: https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py
* Coffea installation: /home/siddhesh/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea
'''
print(f"htoaa_Analysis_GGFMode:: here6 {datetime.now() = }"); sys.stdout.flush()
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
print(f"htoaa_Analysis_GGFMode:: here7 {datetime.now() = }"); sys.stdout.flush()
from particle import Particle # For PDG particle listing https://github.com/scikit-hep/particle
print(f"htoaa_Analysis_GGFMode:: here8 {datetime.now() = }"); sys.stdout.flush()


from htoaa_Settings import *
print(f"htoaa_Analysis_GGFMode:: here9 {datetime.now() = }"); sys.stdout.flush()
from htoaa_CommonTools import (
    getLorentVector,
    GetDictFromJsonFile, akArray_isin,
    selectRunLuminosityBlock,
    calculate_lumiScale, getLumiScaleForPhSpOverlapRewgtMode, getSampleHTRange, # update_crosssection, 
    getNanoAODFile, setXRootDRedirector,  xrdcpFile,
    selectMETFilters, selectFatJets, getCandidateHiggs, selectAK4Jets, selectMuons, selectElectrons,
    selGenPartsWithStatusFlag,
    getHToAATo4BLundPlaneRewgt, getHiggsPtRewgtForGGToHToAATo4B, 
    getTopPtRewgt, getPURewgts, getHTReweight,
    getPURewgts_variation, get_jetTriggerSF, get_PSWeight, add_pdf_as_weight, get_QCDScaleWeight,
    get_JER_and_JES,
    get_Ak4BtagSF,
    calculateAverageOfArrays, calculateMaxOfTwoArrays, calculateMaxOfArrays,  array_PutLowerBound,
    ak_drop_none,
    fillCoffeaHist, fillCoffeaHist_1,
    printVariable, printVariablePtEtaPhi, printVariablePtEtaPhiM, insertInListBeforeThisElement, stringHasSubstring,
)
print(f"htoaa_Analysis_GGFMode:: here10 {datetime.now() = }"); sys.stdout.flush()
from htoaa_Samples import (
    kData, kQCD_bEnrich, kQCD_bGen, kQCDIncl
)
print(f"htoaa_Analysis_GGFMode:: here11 {datetime.now() = }"); sys.stdout.flush()

from inspect import currentframe, getframeinfo
print(f"htoaa_Analysis_GGFMode:: here12 {datetime.now() = }"); sys.stdout.flush()
frameinfo = getframeinfo(currentframe())
print(f"htoaa_Analysis_GGFMode:: here13 {datetime.now() = }"); sys.stdout.flush()


# use GOldenJSON

 
printLevel = 0
histogramSaveLevel = 1 # 0: hSignal extraction, 1: basic Data-MC validation, 2:..
nEventToReadInBatch = 2*10**4 # 0.5*10**5 # 0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze = -1 # 1000 # 100000 # -1
flushStdout = True
#pd.set_option('display.max_columns', None)  

CrossCheckEvtYieldsWithAndrew = False
LumiSecSelThsh_list = [] # [20, 40, 60, 80, 100]

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
        self.wp_ParticleNetMD_XbbvsQCD = 'VL' # 'L', 'VL' 
        self.wp_PNet_Xto4bv1_Htoaa4bOverQCD = 'WP-40' # 'WP-40' 'WP-80' #'WP-60'

        self.FatJetsPt_Thsh             = 170

        self.FatJetPt_gg0lIncl_MinThsh  = 250
        self.FatJetPt_gg0lIncl_MaxThsh  = 999999
        self.FatJetPt_gg0lLo_MinThsh    = 250
        self.FatJetPt_gg0lLo_MaxThsh    = 400
        self.FatJetPt_gg0lHi_MinThsh    = 400
        self.FatJetPt_gg0lHi_MaxThsh    = 999999        

        self.FatJetEtaThsh = 2.4
        self.FatJetJetID   = int(JetIDs.tightIDPassingLeptonVeto)

        self.FatJetMSoftDropThshLow  = 20 # 50 # 20 # 90
        self.FatJetMSoftDropThshHigh = 9999 #200

        self.FatJetParticleNetMD_Xbb_Thsh       = 0.8
        self.FatJetParticleNetMD_XbbvsQCD_Thsh  = bTagWPs[self.era]['ParticleNetMD_XbbvsQCD'][self.wp_ParticleNetMD_XbbvsQCD]
        self.FatJetDeepTagMD_bbvsLight_Thsh     = 0.98
        self.FatJetPNet_Xto4bv1_Htoaa4bOverQCD_Thsh   = bTagWPs[self.era]['PNet_Xto4bv1_Htoaa4bOverQCD'][self.wp_PNet_Xto4bv1_Htoaa4bOverQCD]
        self.FatJetZHbb_plus_Xbb_Thsh = 0.4
        self.FatJetZHbb_Xbb_avg_Thsh  = 0.4
        self.FatJetZHbb_Thsh          = 0.7
        self.FatJetPNetXto4bv2WorkingPoints = ['40'] #['40', '60', '80'] #['40', '45a', '45b', '50', '60', '65', '70', '80']  # ['40', '50', '60', '65', '70', '80']   ['40', '60', '80']
        
        self.nSV_matched_leadingFatJet_Thsh = 3

        # Vjj veto: nonHto4bFatJet
        self.NonHto4bFatJetPNet_WZvsQCD_Thsh = 0.98 # 0.94
        self.NNonHo4bFatJetPNet_WZvsQCD_MaxThsh = 0

        # Lepton veto
        '''
        # Old selection
        self.MuonPtThsh         = 10
        self.MuonMVAId          =  3 # (1=MvaLoose, 2=MvaMedium, 3=MvaTight, 4=MvaVTight, 5=MvaVVTight)
        self.MuonMiniIsoId      =  3 # (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)
        self.MuonMVATTHThsh     = 0.5
        self.ElectronPtThsh     = 10
        self.ElectronMVAId      = 'mvaFall17V2Iso_WP80' # 'mvaFall17V2Iso_WP80', 'mvaFall17V2Iso_WP90' 'mvaFall17V2Iso_WPL'
        self.ElectronMVATTHThsh = 0.3
        '''
        # New selection
        self.MuonPtThsh         = 26 # TrgMuonPt Thsh: 26,  all: 10  
        self.MuonMiniPFRelIsoId =  0.10 # 0.10: tight-WP 
        self.MuonDxyThsh        = 0.02 # 0.2
        self.MuonDzThsh         = 0.10 # 0.5
        self.ElectronPtThsh     = 35 # TrgElectronPt Thsh: 35,  all: 10
        #self.ElectronMVAId      = 'mvaFall17V2Iso_WP90' # 'mvaFall17V2Iso_WP80', 'mvaFall17V2Iso_WP90' 'mvaFall17V2Iso_WPL'
        self.ElectronDxyThsh    = 0.02 # 0.2
        self.ElectronDzThsh     = 0.10 # 0.5

        self.NLeptonsTight_MaxThsh    = 0
        
        # AK4 b-jet veto
        self.Ak4JetDeepJetB_Thsh = bTagWPs[self.era]['AK4DeepJet']['M']


        # VBF AK4 di-jet veto
        self.VBFDijetMass_MinThsh = 450
        self.VBFDijetEta_MinThsh  = 2.2
         


        # Zvv veto: MET
        self.METPt_ZvvIncl_MinThsh = 200 # 250 # 300
        self.METPt_ZvvIncl_MaxThsh = 999999.0
        self.DPhi_FJHto4b_MET_MinThsh = 1.57
        self.NAK4JetsBtagCentral_MaxThsh  = 0

        self.JetPtThshForHT = 30.0
        self.JetEtaThshForHT = 2.4
        
        self.nFatJetMin = 1
        self.GenHTThsh  = 100.0
        self.LHEHTThsh  = 100.0
        


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




    
    
class HToAATo4bProcessor(processor.ProcessorABC):
    def __init__(self, datasetInfo={}):
        print(f"HToAATo4bProcessor::__init__():: {datasetInfo = }")
         
        global runMode_SignalGenChecks;       runMode_SignalGenChecks  = False; # True
        global runMode_QCDGenValidation;      runMode_QCDGenValidation = False; # True
        global runMode_GenLHEPlots;           runMode_GenLHEPlots      =  False
        global runMode_SignificancsScan2D;    runMode_SignificancsScan2D = False
        global runMode_OptimizePNetTaggerCut; runMode_OptimizePNetTaggerCut = False # False
        global runMode_2018HEM1516IssueValidation; runMode_2018HEM1516IssueValidation = False
        global runMode_SignalGenCuts;         runMode_SignalGenCuts = True; # set False for final round. True for optimization studies.
        
        

        ak.behavior.update(nanoaod.behavior)

        self.datasetInfo = datasetInfo
        self.objectSelector             = ObjectSelection(era=self.datasetInfo["era"])
        datasetName_part1               = self.datasetInfo['datasetNameFull'].split('/')[1]
        self.datasetInfo['datasetName'] = datasetName_part1
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
        #global HLT_AK8PFJet330_name
        #HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4" 
        global sTrgSelection
        sTrgSelection = 'Trg_Combo_AK4AK8Jet_HT'
        
        # sel_names_all = dict of {"selection name" : [list of different cuts]}; for cut-flow table 
        self.sel_names_all = OD([
            ("Presel",                    [
                "nPV",
                "METFilters",
                "candH",
                "leadingFatJetPt",
                sTrgSelection,
                "nLeptonsTight",
                "MetZvvVeto",
                "DijetVBFVeto",
                "BJetVeto"
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
                    #sTrgSelection,
                    #"leadingFatJetBtagDeepB",
                    "leadingFatJetMSoftDrop",
                    #"leadingFatJetZHbb_Xbb_avg",
                    "leadingFatJetZHbb",
                    #"leadingFatJetDeepTagMD_bbvsLight", #"leadingFatJetParticleNetMD_Xbb",
                    #"leadingFatJetParticleNetMD_XbbvsQCD",
                    #"leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD",
                    #"leadingFatJet_nSV"
                ]),
            ])
            self.objectSelector.FatJetPtThsh = 170
            #self.objectSelector.FatJetMSoftDropThshLow  =  50 # 60 # 90
            #self.objectSelector.FatJetMSoftDropThshHigh = 9999 # #160 # 140
            #self.objectSelector.FatJetZHbb_plus_Xbb_Thsh = 0.4

        



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
                addBeforeThisCondition = "leadingFatJetPt"
            )


        
        categories_dict = OD()
        categories_dict["gg0lIncl"] = [ "leadingFatJetPt_gg0lIncl" if s_ == "leadingFatJetPt" else s_     for s_ in self.sel_names_all["Presel"] ]
        categories_dict["gg0lLo"]   = [ "leadingFatJetPt_gg0lLo"   if s_ == "leadingFatJetPt" else s_     for s_ in self.sel_names_all["Presel"] ]
        categories_dict["gg0lHi"]   = [ "leadingFatJetPt_gg0lHi"   if s_ == "leadingFatJetPt" else s_     for s_ in self.sel_names_all["Presel"] ]
                 
        for sCatName, catSels in categories_dict.items():
            #self.sel_names_all["%s" % (sCatName)] = catSels
            
            if CrossCheckEvtYieldsWithAndrew:
                self.sel_names_all["%s_SRWP40" % (sCatName)] = catSels + [
                    "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP40"
                ]    
            '''        
            self.sel_names_all["%s_SRWP60" % (sCatName)] = catSels + [
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP60"
            ]
            self.sel_names_all["%s_SRWP80" % (sCatName)] = catSels + [
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP80"
            ]
            self.sel_names_all["%s_SBWP80to40" % (sCatName)] = catSels + [
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP80to40"
            ]
            self.sel_names_all["%s_SBWP95to60" % (sCatName)] = catSels + [
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP95to60"
            ]
            self.sel_names_all["%s_SBWP99to80" % (sCatName)] = catSels + [
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP99to80"
            ]
            '''

            
            for wp_ in self.objectSelector.FatJetPNetXto4bv2WorkingPoints: 
                self.sel_names_all["%s_Xto4bv2_SRWP%s" % (sCatName, wp_)] = catSels + [ # signal region
                    "leadingFatJetPNet_Xto4bv2_Htoaa4b_SRWP%s" % (wp_)
                ]
                self.sel_names_all["%s_Xto4bv2_SBWP%s" % (sCatName, wp_)] = catSels + [ # side band
                    "leadingFatJetPNet_Xto4bv2_Htoaa4b_SBWP%s" % (wp_)
                ] 
                

        
        #for sCatName, catSels in categories_dict.items():
        #    self.sel_names_all["%s" % (sCatName)] = catSels
        for  selCat_ in []: #["gg0lLo_SRWP40", "gg0lHi_SRWP40"]:
            for LumiSecSelThsh in LumiSecSelThsh_list:
                self.sel_names_all["%s_LSlt%d" % (selCat_, LumiSecSelThsh)] = self.sel_names_all[selCat_] + ["LSlt%d"%(LumiSecSelThsh)]
            

        self.sel_conditions_all_list = set()
        for sel_conditions_ in self.sel_names_all.values():
            self.sel_conditions_all_list.update( sel_conditions_ )
        print(f"{self.sel_conditions_all_list = }")
        
        '''
        # selection region addition each SR conditions successively
        #for iCondition in range(self.sel_names_all["Presel"].index(HLT_AK8PFJet330_name), len(self.sel_names_all["Presel"]) - 1):
        #for iCondition in range(self.sel_names_all["Presel"].index("leadingFatJetPt"), len(self.sel_names_all["Presel"]) - 1):
        for iCondition in range(self.sel_names_all["Presel"].index("leadingFatJetMSoftDrop"), len(self.sel_names_all["Presel"]) - 1):
            conditionName = self.sel_names_all["Presel"][iCondition]
            self.sel_names_all["sel_%s" % conditionName] = self.sel_names_all["Presel"][0 : (iCondition+1)]
        '''
        #nPV
        #for iCondition in range(self.sel_names_all["Presel"].index("nPV"), len(self.sel_names_all["Presel"]) - 1):
        #    conditionName = self.sel_names_all["Presel"][iCondition]
        #    self.sel_names_all["sel_%s" % conditionName] = self.sel_names_all["Presel"][0 : (iCondition+1)]
        


        print(f"self.sel_names_all: {json.dumps(self.sel_names_all, indent=4)}")

        self.evtWeights_list = ['PU', ]
        if self.datasetInfo['isSignalGGH']:
            self.evtWeights_list.extend( ['GGHHiggsPt'] )
        if self.datasetInfo['isQCD_bGen']:
            self.evtWeights_list.extend( ['QCDHT'] )
        if self.datasetInfo['isTTbar']:
            self.evtWeights_list.extend( ['TopPt'] )
        
        




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
                logging.critical(f'htoaa_Analysis_GGFMode.py::main():: {sFilesGoldenJSON[self.datasetInfo["era"]] = } could not read.')
                exit(0) 

            # convert runNumber in str to int
            dataLSSelGoldenJSON = {int(k): v for k, v in dataLSSelGoldenJSON.items()} 
            self.datasetInfo['dataLSSelGoldenJSON'] = dataLSSelGoldenJSON
            #print(f"{dataLSSelGoldenJSON = }")

        else: ## MC

            # lumiScale --------------------------------------------------------------------------------------------------
            if sTrgSelection not in Luminosities_forGGFMode[self.datasetInfo["era"]]:
                logging.critical(f'htoaa_Analysis_GGFMode.py::main():: {sTrgSelection = } not in {Luminosities_forGGFMode[self.datasetInfo["era"]] = }.')
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
            self.SFs_ParticleNetMD_XbbvsQCD = None
            if self.objectSelector.wp_ParticleNetMD_XbbvsQCD in Corrections['ParticleNetMD_XbbvsQCD'][self.datasetInfo["era"]].keys():
                print(f" {Corrections['ParticleNetMD_XbbvsQCD'][self.datasetInfo['era']][self.objectSelector.wp_ParticleNetMD_XbbvsQCD]['SFs'] = } "); sys.stdout.flush()
                print(f" {Corrections['ParticleNetMD_XbbvsQCD'][self.datasetInfo['era']][self.objectSelector.wp_ParticleNetMD_XbbvsQCD]['pT_binEdges'] = } "); sys.stdout.flush()

                self.SFs_ParticleNetMD_XbbvsQCD = dense_lookup(
                    np.array( Corrections['ParticleNetMD_XbbvsQCD'][self.datasetInfo["era"]][self.objectSelector.wp_ParticleNetMD_XbbvsQCD]['SFs'] ), # list of SFs
                    [ np.array(Corrections['ParticleNetMD_XbbvsQCD'][self.datasetInfo["era"]][self.objectSelector.wp_ParticleNetMD_XbbvsQCD]['pT_binEdges']) ] # list of bin edges for all axes of SF histogram
                    )
            print(f"{self.SFs_ParticleNetMD_XbbvsQCD = }")

        
        
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
        mass_axis             = hist.Bin("Mass",                   r"$m$ [GeV]",                 350,       0,     350)
        massCl1_axis          = hist.Bin("MassCl1",                r"$m$ [GeV]",                 350,       0,     350)
        mass_axis1            = hist.Bin("Mass1",                  r"$m$ [GeV]",               20*70,       0,     70)
        mass_axis2            = hist.Bin("Mass2",                  r"$m$ [GeV]",                2*70,       0,     70)
        mass10_axis           = hist.Bin("Mass10",                 r"$m$ [GeV]",                 300,       0,      10)
        mass5_axis            = hist.Bin("Mass5",                  r"$m$ [GeV]",                 150,       0,      5)
        logMass3_axis         = hist.Bin("logMass3",               r"$m$ [GeV]",                 300,       0,       3)
        mlScore_axis          = hist.Bin("MLScore",                r"ML score",                  100,    -1.1,     1.1)
        mlScore_axis1         = hist.Bin("MLScore1",               r"ML score",                  100,    -1.1,     1.1)
        mlScore_axis1k        = hist.Bin("MLScore1k",              r"ML score",                 1100,     0.0,     1.1)
        mlScore_axis2k        = hist.Bin("MLScore2k",              r"ML score",                 2100,     0.0,     2.1)
        mlScore_axis1k_m1To2  = hist.Bin("MLScore1k_m1To2",        r"ML score",                 1100,    -1.1,     2.1)
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
        Weight_axis           = hist.Bin("Weight",                 r"Event weight",              [-10,-3,*np.arange(-2,2,0.05), 3, 10])
        
        sXaxis      = 'xAxis'
        sXaxisLabel = 'xAxisLabel'
        sYaxis      = 'yAxis'
        sYaxisLabel = 'yAxisLabel'

        histos = OD()
        
        # General or GEN-level histograms ---------------------------------------------------------------------------------------------
        if histogramSaveLevel >= 10:
            histos.update(OD([
                # ('histogram_name',  {sXaxis: hist.Bin() axis,  sXaxisLabel: "histogram axis label"})

                ('hPV_npvs_beforeSel',                        {sXaxis: PU_axis,                sXaxisLabel: r"No. of primary vertices - before selection"}),
                ('hPV_npvsGood_beforeSel',                    {sXaxis: PU_axis,                sXaxisLabel: r"No. of good primary vertices - before selection"}),            
            ]))

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
            ]))

        

        # RECO-level histograms --------------------------------------------------------------------------------------------------------------
        if self.datasetInfo['isSignal'] and histogramSaveLevel >= 12:
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

                if histogramSaveLevel >= 1:
                    histos.update(OD([
                        ('hEventWeight_PU'+sHExt,                              {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_PU'}),
                        ('hEventWeight_PUUp'+sHExt,                            {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_PUUp'}),
                        ('hEventWeight_PUDown'+sHExt,                          {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_PUDown'}),
                        ('hEventWeight_PS'+sHExt,                              {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_PS'}),
                        ('hEventWeight_PS_ISRUp'+sHExt,                        {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_PS_ISRUp'}),
                        ('hEventWeight_PS_ISRDown'+sHExt,                      {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_PS_ISRDown'}),
                        ('hEventWeight_PS_FSRUp'+sHExt,                        {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_PS_FSRUp'}),
                        ('hEventWeight_PS_FSRDown'+sHExt,                      {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_PS_FSRDown'}),
                        ('hEventWeight_QCDPdfNom'+sHExt,                       {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDPdfNom'}),
                        ('hEventWeight_QCDPdfUp'+sHExt,                        {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDPdfUp'}),
                        ('hEventWeight_QCDPdfDown'+sHExt,                      {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDPdfDown'}),
                        ('hEventWeight_QCDScale_Nom'+sHExt,                    {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDScale_Nom'}),
                        ('hEventWeight_QCDScale_RenormUp'+sHExt,               {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDScale_RenormUp'}),
                        ('hEventWeight_QCDScale_RenormDown'+sHExt,             {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDScale_RenormDown'}),
                        ('hEventWeight_QCDScale_FactorizationUp'+sHExt,        {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDScale_FactorizationUp'}),
                        ('hEventWeight_QCDScale_FactorizationDown'+sHExt,      {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDScale_FactorizationDown'}),
                        ('hEventWeight_QCDPDFNom'+sHExt,                       {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDPDFNom'}),
                        ('hEventWeight_QCDPDFUp'+sHExt,                        {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDPDFUp'}),
                        ('hEventWeight_QCDPDFDown'+sHExt,                      {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDPDFDown'}),
                        ('hEventWeight_Ak4BtagNom'+sHExt,                      {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_Ak4BtagNom'}),
                        ('hEventWeight_Ak4BtagUp'+sHExt,                       {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_Ak4BtagUp'}),
                        ('hEventWeight_Ak4BtagDown'+sHExt,                     {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_Ak4BtagDown'}),
                        ('hEventWeight_Ak4BtagUpuncorrelated'+sHExt,           {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_Ak4BtagUpuncorrelated'}),
                        ('hEventWeight_Ak4BtagDownuncorrelated'+sHExt,         {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_Ak4BtagDownuncorrelated'}),
                        ('hEventWeight_Ak4BtagUpcorrelated'+sHExt,             {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_Ak4BtagUpcorrelated'}),
                        ('hEventWeight_Ak4BtagDowncorrelated'+sHExt,           {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_Ak4BtagDowncorrelated'}),
                        
                                               
                    ]))
                    if self.datasetInfo['isSignal']:
                        histos.update(OD([
                            ('hEventWeight_LundPlane_Nom'+sHExt,                       {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_LundPlane_Nom'}),
                            ('hEventWeight_LundPlane_Up'+sHExt,                        {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_LundPlane_Up'}),
                            ('hEventWeight_LundPlane_Down'+sHExt,                      {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_LundPlane_Down'}),
                            
                        ]))
                    if self.datasetInfo['isSignalGGH']:
                        histos.update(OD([
                            ('hEventWeight_GGHHiggsPt'+sHExt,                        {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_GGHHiggsPt'}),
                            ('hEventWeight_GGHHiggsPtUp'+sHExt,                      {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_GGHHiggsPtUp'}),
                            ('hEventWeight_GGHHiggsPtDown'+sHExt,                    {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_GGHHiggsPtDown'}),
                            
                        ]))
                    if self.datasetInfo['isQCD_bGen']:
                        histos.update(OD([
                            ('hEventWeight_QCDHT'+sHExt,                        {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_QCDHT'}),
                        ]))
                    if self.datasetInfo['isTTbar']:
                        histos.update(OD([
                            ('hEventWeight_TopPt'+sHExt,                        {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_TopPt'}),
                            ('hEventWeight_TopPtUp'+sHExt,                      {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_TopPtUp'}),
                            ('hEventWeight_TopPtDown'+sHExt,                    {sXaxis: Weight_axis,    sXaxisLabel: 'hEventWeight_TopPtDown'}),
                            
                        ]))



                if histogramSaveLevel >= 0:
                    histos.update(OD([
                        ('hCutFlow'+sHExt,                                  {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),
                        ('hCutFlowWeighted'+sHExt,                          {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),

                        ('hLeadingFatJetParticleNet_massH_Hto4b_avg_vs_massA_Hto4b_avg'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_avg",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg"}),     

                        ('hLeadingFatJetMass_vs_massA_Hto4b_avg'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMass",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg"}),

                        ('hLeadingFatJetMSoftDrop_vs_massA_Hto4b_avg'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMSoftDrop",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg"}),        

                        ## PNetMD Hto4b NanoAOD_v2
                        ('hLeadingFatJetPNet_massH_v2b_vs_massAa'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetPNet_massH",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massAa"}),     

                        ('hLeadingFatJetMass_vs_massAa'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMass",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massAa"}),

                        ('hLeadingFatJetMSoftDrop_vs_massAa'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMSoftDrop",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massAa"}),     

                        ('hLeadingFatJetPNet_massH_v2b_vs_massA34a'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetPNet_massH",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massA34a"}),
                        ('hLeadingFatJetMass_vs_massA34a'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMass",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massA34a"}),
                        ('hLeadingFatJetMSoftDrop_vs_massA34a'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMSoftDrop",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massA34a"}),
                        ('hLeadingFatJetMass_vs_massA34b'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMass",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massA34b"}),
                        ('hLeadingFatJetMass_vs_massA34d'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMass",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massA34d"}),
                        ('hLeadingFatJetMass_vs_massA34ad'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMass",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massA34ad"}),
                        ('hLeadingFatJetMass_vs_massA34d'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMass",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massA34d"}),                        
                        ('hLeadingFatJetMSoftDrop_vs_massA34d'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetMSoftDrop",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massA34d"}),                        
                        ('hLeadingFatJetPNet_massH_v2b_vs_massA34d'+sHExt,     
                        {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetPNet_massH",
                        sYaxis: mass_axis2,       sYaxisLabel: r"hLeadingFatJetPNet_massA34d"}),                        
                         
                    ]))


                if histogramSaveLevel >= 1:
                    histos.update(OD([

                        ('hPV_npvsGood'+sHExt,                           {sXaxis: PU_axis,                sXaxisLabel: r"No. of good primary vertices - signal region"}),

                        ('hLeadingFatJetPt'+sHExt,                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                        ('hLeadingFatJetEta'+sHExt,                         {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                        ('hLeadingFatJetPhi'+sHExt,                         {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
                        ('hLeadingFatJetMass'+sHExt,                        {sXaxis: mass_axis,       sXaxisLabel: r"m (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetMSoftDrop'+sHExt,                   {sXaxis: mass_axis,       sXaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                        ('hLeadingFatJetParticleNet_massH_Hto4b_avg_v0123'+sHExt,     {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_avg_v0123"}), ## selected
                        ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v013'+sHExt,      {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v013"}), ## selected
                        
                        ('hnleadingNonHto4bFatJet_WZvsQCD'+sHExt,             {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of leading non-Hto4b fat jets WZvsQCD-tagged"}),
                        
                        ('hMET_pT'+sHExt,                                   {sXaxis: pt_axis,         sXaxisLabel: r"MET pT [GeV]"}),
                        ('hPuppiMET_pT'+sHExt,                              {sXaxis: pt_axis,         sXaxisLabel: r"PuppiMET pT [GeV]"}),
                        
                        # NanoAODV2
                        ('hLeadingFatJetPNet_X4b_v1_Haa4b_vs_QCD'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetPNet_X4b_v1_Haa4b_vs_QCD"}),
                        ('hLeadingFatJetPNet_X4b_v1_Haa4b_score'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetPNet_X4b_v1_Haa4b_score"}),
                        ('hLeadingFatJetPNet_X4b_v2a_Haa4b_score'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetPNet_X4b_v2a_Haa4b_score"}),
                        ('hLeadingFatJetPNet_X4b_v2b_Haa4b_score'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetPNet_X4b_v2b_Haa4b_score"}),
                        ('hLeadingFatJetPNet_X4b_v2ab_Haa4b_score'+sHExt,   {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetPNet_X4b_v2ab_Haa4b_score"}),
                        ('hLeadingFatJetPNet_X4b_v2a_Haa34b_score'+sHExt,   {sXaxis: mlScore_axis1k,  sXaxisLabel: r"hLeadingFatJetPNet_X4b_v2a_Haa34b_score"}),
                        ('hLeadingFatJetPNet_X4b_v2b_Haa34b_score'+sHExt,   {sXaxis: mlScore_axis1k,  sXaxisLabel: r"hLeadingFatJetPNet_X4b_v2b_Haa34b_score"}),
                        ('hLeadingFatJetPNet_X4b_v2ab_Haa34b_score'+sHExt,   {sXaxis: mlScore_axis1k,  sXaxisLabel: r"hLeadingFatJetPNet_X4b_v2ab_Haa34b_score"}),                        
                        #
                        #('hLeadingFatJetMassH_v1'+sHExt,                    {sXaxis: mass_axis,       sXaxisLabel: r"m (leading FatJet MassH_1) [GeV]"}),
                        #('hLeadingFatJetMassH_v2a'+sHExt,                   {sXaxis: mass_axis,       sXaxisLabel: r"m (leading FatJet MassH_v2a) [GeV]"}),
                        ('hLeadingFatJetMassH_v2b'+sHExt,                   {sXaxis: mass_axis,       sXaxisLabel: r"m (leading FatJet MassH_v2b) [GeV]"}),
                        #('hLeadingFatJetMassH_v2c'+sHExt,                   {sXaxis: mass_axis,       sXaxisLabel: r"m (leading FatJet MassH_v2c) [GeV]"}),
                        #('hLeadingFatJetMassH_v2d'+sHExt,                   {sXaxis: mass_axis,       sXaxisLabel: r"m (leading FatJet MassH_v2d) [GeV]"}),
                        #('hLeadingFatJetMassH_avg'+sHExt,                   {sXaxis: mass_axis,       sXaxisLabel: r"m (leading FatJet MassH_avg) [GeV]"}),
                        #('hLeadingFatJetMassH_std'+sHExt,                   {sXaxis: mass_axis1,      sXaxisLabel: r"m (leading FatJet MassH_std) [GeV]"}),
                        
                        #
                        ('hLeadingFatJetPNet_massAa'+sHExt,         {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_massAa"}), 
                        #('hLeadingFatJetPNet_massAb'+sHExt,         {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_massAb"}), 
                        #('hLeadingFatJetPNet_massAc'+sHExt,         {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_massAc"}), 
                        #('hLeadingFatJetPNet_massAd'+sHExt,         {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_massAd"}), 
                        #('hLeadingFatJetPNet_massA_avg'+sHExt,      {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_massA_avg"}), 
                        #('hLeadingFatJetPNet_massA_std'+sHExt,      {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_massA_std"}), 
                        ('hLeadingFatJetPNet_34massAa'+sHExt,       {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_34massAa"}), 
                        ('hLeadingFatJetPNet_34massAb'+sHExt,       {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_34massAb"}), 
                        #('hLeadingFatJetPNet_34massAc'+sHExt,       {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_34massAc"}), 
                        ('hLeadingFatJetPNet_34massAd'+sHExt,       {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_34massAd"}), 
                        ('hLeadingFatJetPNet_massA1'+sHExt,         {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_massA1"}), 
                        ('hLeadingFatJetPNet_massA2'+sHExt,         {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_massA2"}), 
                        ('hLeadingFatJetPNet_massAA'+sHExt,         {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_massAA"}), 
                        ('hLeadingFatJetPNet_dMassAA'+sHExt,        {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetPNet_dMassAA"}), 
                        ('hLeadingFatJetPNet_dMassAA_relH'+sHExt,   {sXaxis: mass5_axis,       sXaxisLabel: r"hLeadingFatJetPNet_dMassAA_relH"}), 
                        
                        
                        
                    ]))


                if histogramSaveLevel >= 2:
                    histos.update(OD([
                        ('hPV_npvs'+sHExt,                               {sXaxis: PU_axis,                sXaxisLabel: r"No. of primary vertices - signal region"}),

                        ('nSelFatJet'+sHExt,                                {sXaxis: nObject_axis,    sXaxisLabel: 'No. of selected FatJets'}),
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
                        
                        ('hLeadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Htoaa4bOverQCD"}),
                        ('hLeadingFatJetParticleNetMD_Hto4b_Htoaa34bOverQCD'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Htoaa34bOverQCD"}),
                        ('hLeadingFatJetParticleNetMD_Hto4b_binaryLF_Htoaa4bOverQCD'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binaryLF_Htoaa4bOverQCD"}),
                        ('hLeadingFatJetParticleNetMD_Hto4b_binary_Htoaa4bOverQCD'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binary_Htoaa4bOverQCD"}),
                        ('hLeadingFatJetParticleNetMD_Hto4b_binary_Htoaa4bOverQCD_avg'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_binary_Htoaa4bOverQCD_avg"}),
                        ('hLeadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_avg'+sHExt,    {sXaxis: mlScore_axis1k,  sXaxisLabel: r"LeadingFatJetParticleNetMD Hto4b_Htoaa4bOverQCD_avg"}),
                        
                        
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
                        #('hLeadingFatJetParticleNet_massA_Hto4b_avg_v013'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v013"}), ## selected
                        ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v023'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v023"}),
                        ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v123'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v123"}),
                        ('hLeadingFatJetParticleNet_massA_Hto4b_avg_v0123'+sHExt,            {sXaxis: mass_axis1,       sXaxisLabel: r"hLeadingFatJetParticleNet_massA_Hto4b_avg_v0123"}),                                        

                        ('hLeadingFatJetParticleNet_massH_Hto4b_v0'+sHExt,            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v0"}),
                        ('hLeadingFatJetParticleNet_massH_Hto4b_v00'+sHExt,           {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v00"}),
                        ('hLeadingFatJetParticleNet_massH_Hto4b_v1'+sHExt,            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v1"}),
                        ('hLeadingFatJetParticleNet_massH_Hto4b_v2'+sHExt,            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v2"}),
                        ('hLeadingFatJetParticleNet_massH_Hto4b_v3'+sHExt,            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v3"}),
                        ('hLeadingFatJetParticleNet_massH_Hto4b_v4'+sHExt,            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_v4"}),

                        ('hLeadingFatJetParticleNet_massH_Hto4b_avg_v0123'+sHExt,     {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_massH_Hto4b_avg_v0123"}), ## selected
                        
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


                        ## NonHto4bFatJets
                        ('hnNonHto4bFatJet'+sHExt,                            {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of non-Hto4b fat jets "}),
                        #('hnleadingNonHto4bFatJet_WZvsQCD'+sHExt,             {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of leading non-Hto4b fat jets WZvsQCD-tagged"}),
                        

                        ## MET
                        #('hMET_pT'+sHExt,                                   {sXaxis: pt_axis,         sXaxisLabel: r"MET pT [GeV]"}),
                        ('hMET_sumEt'+sHExt,                                {sXaxis: pt4TeV_axis,     sXaxisLabel: r"MET sumEt [GeV]"}),
                        ('hdPhi_MET_leadingFatJet'+sHExt,                   {sXaxis: deltaPhi_axis,   sXaxisLabel: r"deltaPhi(MET, leadingFatJet)"}),
                        #('hPuppiMET_pT'+sHExt,                              {sXaxis: pt_axis,         sXaxisLabel: r"PuppiMET pT [GeV]"}),
                        ('hPuppiMET_sumEt'+sHExt,                           {sXaxis: pt4TeV_axis,     sXaxisLabel: r"PuppiMET sumEt [GeV]"}),
                        ('hdPhi_PuppiMET_leadingFatJet'+sHExt,              {sXaxis: deltaPhi_axis,   sXaxisLabel: r"deltaPhi(PuppiMET, leadingFatJet)"}),
                        ('hPuppiMET_sumEt_minus_FJHto4bPt'+sHExt,           {sXaxis: pt2TeV_axis,     sXaxisLabel: r"PuppiMET sumEt - pT(AK8 jet H->4b) [GeV]"}),
                        ('hMETPhi'+sHExt,                         {sXaxis: phi_axis,        sXaxisLabel: r"\phi (MET)"}),
                        ('hPuppiMETPhi'+sHExt,                         {sXaxis: phi_axis,        sXaxisLabel: r"\phi (MET)"}),
                        
                        ## nLeptons_matched_leadingFatJet
                        ('hLeadingFatJet_nLeptons'+sHExt,                   {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of iso-leptons within leadingFatJet "}),
                        ('hnLeptonsTight'+sHExt,                            {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of tight leptons "}),
                        ('hnLeptons_nonoverlap_leadingFatJet'+sHExt,        {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of tight leptons nonoverlaping leadingFatJet"}),
                        
                        ## AK4 jets
                        ('hnAK4Jets_NonoverlapLeadingFatJet'+sHExt,                     {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of AK4 jets non-overlap FatJet H->4b "}),
                        ('hPtLeadingAK4Jets_NonoverlapLeadingFatJet'+sHExt,             {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading AK4 jets non-overlap FatJet H->4b) [GeV]"}),
                        ('hnAK4Jets_bTag_NonoverlapLeadingFatJet'+sHExt,                {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of bTag AK4 jets non-overlap FatJet H->4b "}),
                        ('hPtLeadingAK4Jets_bTag_NonoverlapLeadingFatJet'+sHExt,        {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading bTag AK4 jets non-overlap FatJet H->4b) [GeV]"}),

                        ('hnAK4JetsCentral_NonoverlapLeadingFatJet'+sHExt,                     {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of central AK4 jets non-overlap FatJet H->4b "}),
                        ('hPtLeadingAK4JetsCentral_NonoverlapLeadingFatJet'+sHExt,             {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading central AK4 jets non-overlap FatJet H->4b) [GeV]"}),
                        ('hnAK4JetsCentral_bTag_NonoverlapLeadingFatJet'+sHExt,                {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of bTag central AK4 jets non-overlap FatJet H->4b "}),
                        ('hPtLeadingAK4JetsCentral_bTag_NonoverlapLeadingFatJet'+sHExt,        {sXaxis: pt_axis,         sXaxisLabel: r"pT(Leading bTag central AK4 jets non-overlap FatJet H->4b) [GeV]"}),

                    ]))

                ### 2-D distribution --------------------------------------------------------------------------------------------------------
                if histogramSaveLevel >= 2:    
                    histos.update(OD([
                        ('hLeadingFatJetEta_vs_Phi'+sHExt,             
                        {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)",
                        sYaxis: phi_axis,        sYaxisLabel: r"\phi (leading FatJet)"}),   

                        ('hMET_pT_vs_dPhi_MET_leadingFatJet'+sHExt,     
                        {sXaxis: pt_axis,         sXaxisLabel: r"MET pT [GeV]",
                        sYaxis: deltaPhi_axis,   sYaxisLabel: r"deltaPhi(MET, leadingFatJet)"}),                                                                              
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


            
        #for statusFlag_ in GENPART_STATUSFLAGS_LIST:
        #    histos['hGenBquark_first_%s_all' % (statusFlag_)] = {sXaxis: boolean_axis,    sXaxisLabel: r"GEN first Bquark %s"  % (statusFlag_)}
        
        self._accumulator = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int)
        })
        
        for histName, histAttributes in histos.items():
            #hXaxis = histAttributes[sXaxis].copy()
            hXaxis = deepcopy(histAttributes[sXaxis])
            hXaxis.label = histAttributes[sXaxisLabel]

            if histName.startswith('hEventWeight_'): 
                # TH1 w/o syst
                self._accumulator.add({
                    histName: hist.Hist(
                        "Counts",
                        dataset_axis,
                        hXaxis, #nObject_axis,
                    )
                })
            elif sYaxis not in histAttributes.keys():
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

    def process11(self, events):
        print(f"process():: *** ", flush=flushStdout)
        output = self.accumulator.identity()
        return

    def process(self, events):
        dataset = events.metadata["dataset"] # dataset label
        print(f"process():: {self.datasetInfo['sample_category'] = }, {dataset = }", flush=flushStdout)

        if printLevel >= 20:
            print(f"nEvents: {len(events)}")
        if printLevel >= 1:
            print(f"\n events.fields ({type(events.fields)}): {events.fields}"); sys.stdout.flush()
             
        if nEventsToAnalyze != -1:
            #print(f"\n (run:ls:event): {ak.zip([events.run, events.luminosityBlock, events.event])}") 
            printVariable('\n (run:ls:event): ', ak.zip([events.run, events.luminosityBlock, events.event])); #sys.stdout.flush()     
        #printVariable('\n (run:ls:event): ', ak.zip([events.run, events.luminosityBlock, events.event])[:20]); #sys.stdout.flush()        

        #if not self.datasetInfo['isMC']:
        if self.datasetInfo['isMC']:
            output = self.accumulator.identity()
            
            #systematics_shift = [None, "JERUp", "JERDown", "JESUp", "JESDown"] # [None, "JERUp", "JERDown", "JESUp", "JESDown"]
            systematics_shift = [None]
            if not self.datasetInfo['systematicsToRun'] == 'no':
                if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['jer', 'full'] ): 
                    systematics_shift.extend( [
                        "JERUp",
                        "JERDown",
                    ] )
                if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['jes', 'full'] ): 
                    systematics_shift.extend( [
                        "JESUp",
                        "JESDown",
                    ] )
                if (stringHasSubstring(self.datasetInfo['systematicsToRun'], ['jeshemissue', 'full'] ) and 
                   (self.datasetInfo["era"] == Era_2018) ):
                    systematics_shift.extend( [
                        "JESHEMIssueUp",
                        "JESHEMIssueDown",
                    ] )                


            #if not self.datasetInfo['systematicsToRun'] == 'no': 
            #    systematics_shift.appen()
            for _syst in systematics_shift:
                output += self.process_shift(events, _syst)
        else:
            print(f" {np.unique(events.run, return_counts=True) = } "); sys.stdout.flush()
            output = self.process_shift(events, None)


        return output


    
    def process_shift(self, events, shift_syst=None):
        
        output = self.accumulator.identity()
        dataset = events.metadata["dataset"] # dataset label
        print(f"process_shift():: {shift_syst = } dataset: {dataset}", flush=flushStdout)

        

        ones_list  = np.ones(len(events))
        trues_list = np.ones(len(events), dtype=bool)
        falses_list = np.full(len(events), False)

        ###########################################
        # UPDATE OBJECTS FOR SYSTEMATICS VARIATIONS
        ###########################################

        #FatJetsToUse = get_JER_and_JES(events, events.FatJet, self.datasetInfo["era"], shift_syst)
        FatJetsToUse = events.FatJet
        JetsToUse    = events.Jet
        METToUse     = events.MET

        
        if 'PNet_X4b_v2a_Haa34b_score' not in FatJetsToUse.fields or CrossCheckEvtYieldsWithAndrew:
            FatJetsToUse['pt_toUse']        = FatJetsToUse.pt 
            JetsToUse['pt_toUse']           = JetsToUse.pt
            METToUse['pt_toUse']            = METToUse.pt
            FatJetsToUse['mass_toUse']      = FatJetsToUse.mass
            FatJetsToUse['msoftdrop_toUse'] = FatJetsToUse.msoftdrop
            JetsToUse['mass_toUse']         = JetsToUse.mass

        else: 
            ## NanoAOD v2: Depending up on systematics to run, set which branches to read 
            # FatJet pt
            if   shift_syst == 'JESUp':            FatJetsToUse['pt_toUse'] = FatJetsToUse.pt_jesTotalUp
            elif shift_syst == 'JESDown':          FatJetsToUse['pt_toUse'] = FatJetsToUse.pt_jesTotalDown
            elif shift_syst == 'JERUp':            FatJetsToUse['pt_toUse'] = FatJetsToUse.pt_jerUp
            elif shift_syst == 'JERDown':          FatJetsToUse['pt_toUse'] = FatJetsToUse.pt_jerDown
            elif shift_syst == 'JESHEMIssueUp':    FatJetsToUse['pt_toUse'] = FatJetsToUse.pt_jesHEMIssueUp
            elif shift_syst == 'JESHEMIssueDown':  FatJetsToUse['pt_toUse'] = FatJetsToUse.pt_jesHEMIssueDown
            else:                                  FatJetsToUse['pt_toUse'] = FatJetsToUse.pt_nom
            # FatJet mass
            if   shift_syst == 'JESUp':            FatJetsToUse['mass_toUse'] = FatJetsToUse.mass_jesTotalUp
            elif shift_syst == 'JESDown':          FatJetsToUse['mass_toUse'] = FatJetsToUse.mass_jesTotalDown
            elif shift_syst == 'JERUp':            FatJetsToUse['mass_toUse'] = FatJetsToUse.mass_jerUp
            elif shift_syst == 'JERDown':          FatJetsToUse['mass_toUse'] = FatJetsToUse.mass_jerDown
            elif shift_syst == 'JESHEMIssueUp':    FatJetsToUse['mass_toUse'] = FatJetsToUse.mass_jesHEMIssueUp
            elif shift_syst == 'JESHEMIssueDown':  FatJetsToUse['mass_toUse'] = FatJetsToUse.mass_jesHEMIssueDown
            else:                                  FatJetsToUse['mass_toUse'] = FatJetsToUse.mass_nom
            # FatJet msoftdrop
            if   shift_syst == 'JESUp':            FatJetsToUse['msoftdrop_toUse'] = FatJetsToUse.msoftdrop_jesTotalUp
            elif shift_syst == 'JESDown':          FatJetsToUse['msoftdrop_toUse'] = FatJetsToUse.msoftdrop_jesTotalDown
            elif shift_syst == 'JERUp':            FatJetsToUse['msoftdrop_toUse'] = FatJetsToUse.msoftdrop_jerUp
            elif shift_syst == 'JERDown':          FatJetsToUse['msoftdrop_toUse'] = FatJetsToUse.msoftdrop_jerDown
            elif shift_syst == 'JESHEMIssueUp':    FatJetsToUse['msoftdrop_toUse'] = FatJetsToUse.msoftdrop_jesHEMIssueUp
            elif shift_syst == 'JESHEMIssueDown':  FatJetsToUse['msoftdrop_toUse'] = FatJetsToUse.msoftdrop_jesHEMIssueDown
            else:                                  FatJetsToUse['msoftdrop_toUse'] = FatJetsToUse.msoftdrop_nom
            
            # AK4 Jet pt
            if   shift_syst == 'AK4JESUp':         JetsToUse['pt_toUse'] = JetsToUse.pt_jesTotalUp
            elif shift_syst == 'AK4JESDown':       JetsToUse['pt_toUse'] = JetsToUse.pt_jesTotalDown
            elif shift_syst == 'AK4JERUp':         JetsToUse['pt_toUse'] = JetsToUse.pt_jerUp
            elif shift_syst == 'AK4JERDown':       JetsToUse['pt_toUse'] = JetsToUse.pt_jerDown
            else:                                  JetsToUse['pt_toUse'] = JetsToUse.pt_nom   
            # AK4 Jet mass
            if   shift_syst == 'AK4JESUp':         JetsToUse['mass_toUse'] = JetsToUse.mass_jesTotalUp
            elif shift_syst == 'AK4JESDown':       JetsToUse['mass_toUse'] = JetsToUse.mass_jesTotalDown
            elif shift_syst == 'AK4JERUp':         JetsToUse['mass_toUse'] = JetsToUse.mass_jerUp
            elif shift_syst == 'AK4JERDown':       JetsToUse['mass_toUse'] = JetsToUse.mass_jerDown
            else:                                  JetsToUse['mass_toUse'] = JetsToUse.mass_nom

            # MET pt
            if   shift_syst == 'METJESUp':         METToUse['pt_toUse'] = METToUse.T1Smear_pt_jesTotalUp
            elif shift_syst == 'METJEDown':        METToUse['pt_toUse'] = METToUse.T1Smear_pt_jesTotalDown
            elif shift_syst == 'METJERUp':         METToUse['pt_toUse'] = METToUse.T1Smear_pt_jerUp
            elif shift_syst == 'METJERDown':       METToUse['pt_toUse'] = METToUse.T1Smear_pt_jerDown
            elif shift_syst == 'METUnclstEnUp':    METToUse['pt_toUse'] = METToUse.T1Smear_pt_unclstEnUp
            elif shift_syst == 'METUnclstEnDown':  METToUse['pt_toUse'] = METToUse.T1Smear_pt_unclstEnDown
            elif not self.datasetInfo['isMC']:     METToUse['pt_toUse'] = METToUse.T1_pt
            else:                                  METToUse['pt_toUse'] = METToUse.T1Smear_pt               
            # MET phi
            if   shift_syst == 'METJESUp':         METToUse['phi_toUse'] = METToUse.T1Smear_phi_jesTotalUp
            elif shift_syst == 'METJEDown':        METToUse['phi_toUse'] = METToUse.T1Smear_phi_jesTotalDown
            elif shift_syst == 'METJERUp':         METToUse['phi_toUse'] = METToUse.T1Smear_phi_jerUp
            elif shift_syst == 'METJERDown':       METToUse['phi_toUse'] = METToUse.T1Smear_phi_jerDown
            elif shift_syst == 'METUnclstEnUp':    METToUse['phi_toUse'] = METToUse.T1Smear_phi_unclstEnUp
            elif shift_syst == 'METUnclstEnDown':  METToUse['phi_toUse'] = METToUse.T1Smear_phi_unclstEnDown
            elif not self.datasetInfo['isMC']:     METToUse['phi_toUse'] = METToUse.T1_phi
            else:                                  METToUse['phi_toUse'] = METToUse.T1Smear_phi
        

        if printLevel >= 100 :
            print(f"{shift_syst = }, {FatJetsToUse.fields = }")
            printVariable('\n events.FatJet.pt \n', events.FatJet.pt)
            printVariable('\n FatJetsToUse.pt\n', FatJetsToUse.pt)
            printVariable('\n FatJetsToUse.pt_toUse\n', FatJetsToUse.pt_toUse)

            print(f"{shift_syst = }, {JetsToUse.fields = }")
            printVariable('\n JetsToUse.pt \n', JetsToUse.pt)
            printVariable('\n JetsToUse.pt_toUse\n', JetsToUse.pt_toUse)


          
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
                #printVariable("\n genB_fromHToAA", genB_fromHToAA)

            
                
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
        mask_FatJetPt = (FatJetsToUse.pt_toUse > self.objectSelector.FatJetPtThsh)
        mask_FatJetEta = (abs(FatJetsToUse.eta) < self.objectSelector.FatJetEtaThsh)
        mask_FatJetBtagDeepB = (FatJetsToUse.btagDeepB > bTagWPs[self.objectSelector.era][self.objectSelector.tagger_btagDeepB][self.objectSelector.wp_btagDeepB])
        mask_FatJetMSoftDrop = (
            (FatJetsToUse.msoftdrop_toUse > self.objectSelector.FatJetMSoftDropThshLow) &
            (FatJetsToUse.msoftdrop_toUse < self.objectSelector.FatJetMSoftDropThshHigh)
        )
        
        selFatJet = FatJetsToUse[(
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
            run_FatJetEta_FatJetPhi = ak.zip({'run': events.run, 'FatJetEta': FatJetsToUse.eta, 'FatJetPhi': FatJetsToUse.phi})   
            mask_jets_surviving_HEM15_16_issue = ~ (
                (run_FatJetEta_FatJetPhi.run >= 319077) &
                (run_FatJetEta_FatJetPhi.FatJetEta > -3.2 ) & (run_FatJetEta_FatJetPhi.FatJetEta < -1.3 ) &
                (run_FatJetEta_FatJetPhi.FatJetPhi > -1.57) & (run_FatJetEta_FatJetPhi.FatJetPhi < -0.87)
            )
            leadingFatJet = ak.firsts(FatJetsToUse[mask_jets_surviving_HEM15_16_issue]) 
        else:   
            leadingFatJet = ak.firsts(FatJetsToUse) # for e.g. [0.056304931640625, None, 0.12890625, 0.939453125, 0.0316162109375]
        '''

        leadingFatJet = None

        FatJetParticleNetMD_XbbvsQCD = ak.where(
            (FatJetsToUse.particleNetMD_Xbb + FatJetsToUse.particleNetMD_QCD) > 0,
            FatJetsToUse.particleNetMD_Xbb / (FatJetsToUse.particleNetMD_Xbb + FatJetsToUse.particleNetMD_QCD),
            np.full_like(FatJetsToUse.particleNetMD_Xbb, 0)
        ) 
        FatJetDeepTagMD_ZHbbvsQCD = ak.where(
            FatJetsToUse.deepTagMD_ZHbbvsQCD >= 0,
            FatJetsToUse.deepTagMD_ZHbbvsQCD,
            np.full_like(FatJetsToUse.deepTagMD_ZHbbvsQCD, 0)
        )
        FatJet_ZHbb_plus_Xbb = FatJetDeepTagMD_ZHbbvsQCD + FatJetParticleNetMD_XbbvsQCD
        idx_FatJet_ZHbb_plus_Xbb_max = ak.argmax(FatJet_ZHbb_plus_Xbb, axis=-1, keepdims=True)


        '''
        # PNetMD_Hto4b
        idx_FatJet_PNetMD_Hto4b_Haa4bOverQCD_max = None
        if 'particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields:
            print(f"\n\nparticleNetMD_Hto4b_Haa4b NanoAODv1\n\n")
            FatJet_PNetMD_Hto4b_QCD01234b_sum = (
                FatJetsToUse.particleNetMD_Hto4b_QCD0b + 
                FatJetsToUse.particleNetMD_Hto4b_QCD1b + 
                FatJetsToUse.particleNetMD_Hto4b_QCD2b + 
                FatJetsToUse.particleNetMD_Hto4b_QCD3b + 
                FatJetsToUse.particleNetMD_Hto4b_QCD4b )  
            FatJet_PNetMD_Hto4b_Htoaa4bOverQCD = ak.where(
                (FatJetsToUse.particleNetMD_Hto4b_Haa4b + FatJet_PNetMD_Hto4b_QCD01234b_sum) > 0.0,
                (
                    FatJetsToUse.particleNetMD_Hto4b_Haa4b / 
                    (FatJetsToUse.particleNetMD_Hto4b_Haa4b + FatJet_PNetMD_Hto4b_QCD01234b_sum)
                ),
                ak.full_like(FatJetsToUse.particleNetMD_Hto4b_Haa4b, 0) #FatJetsToUse.particleNetMD_Hto4b_Haa4b
            ) 
            idx_FatJet_PNetMD_Hto4b_Haa4bOverQCD_max = ak.argmax(FatJet_PNetMD_Hto4b_Htoaa4bOverQCD, axis=-1, keepdims=True)
            #leadingBtagFatJet = ak.firsts(FatJetsToUse[idx_FatJet_PNetMD_Hto4b_Haa4bOverQCD_max])    
            leadingFatJet = ak.firsts(FatJetsToUse[idx_FatJet_PNetMD_Hto4b_Haa4bOverQCD_max]) 
        
        elif 'PNet_X4b_v2a_Haa34b_score' in FatJetsToUse.fields:
            print(f"\n\PNet_X4b_v2a_Haa34b_score NanoAODv2\n\n")
            FatJet_PNet_X4b_v2_Haa34b = FatJetsToUse.PNet_X4b_v2a_Haa34b_score + FatJetsToUse.PNet_X4b_v2b_Haa34b_score
            idx_FatJet_PNet_X4b_v2_Haa34b_max = ak.argmax(FatJet_PNet_X4b_v2_Haa34b, axis=-1, keepdims=True)
            leadingFatJet = ak.firsts(FatJetsToUse[idx_FatJet_PNet_X4b_v2_Haa34b_max]) 
        else:
  
            #idx_FatJet_PNetMD_XbbvsQCD_max = ak.argmax(FatJetParticleNetMD_XbbvsQCD, axis=-1, keepdims=True)
            #leadingBtagFatJet = ak.firsts(FatJetsToUse[idx_FatJet_PNetMD_XbbvsQCD_max]) 
            #leadingFatJet = ak.firsts(FatJetsToUse[idx_FatJet_PNetMD_XbbvsQCD_max]) 
            leadingFatJet = ak.firsts(FatJetsToUse[idx_FatJet_ZHbb_plus_Xbb_max]) 
        '''

        
        #print(f"{self.objectSelector.FatJetsPt_Thsh = }, {self.objectSelector.FatJetEtaThsh = }, {self.objectSelector.FatJetMSoftDropThshLow = }, {self.objectSelector.FatJetJetID = }, {self.objectSelector.FatJetParticleNetMD_XbbvsQCD_Thsh = }, ",flush=True)
        selFatJets = selectFatJets(
            FatJetsToUse, 
            pT_Thsh  = self.objectSelector.FatJetsPt_Thsh, 
            eta_Thsh = self.objectSelector.FatJetEtaThsh, 
            Msd_Thsh = self.objectSelector.FatJetMSoftDropThshLow, 
            JetID    = self.objectSelector.FatJetJetID,
            shift_syst = shift_syst
        )
        leadingFatJet, idx_candHs_PNet_X4b_v2_Haa34b_max = getCandidateHiggs(
            selFatJets, 
            Xbb_Thsh = self.objectSelector.FatJetParticleNetMD_XbbvsQCD_Thsh
            )
        LV_leadingFatJet_wMass      = getLorentVector(leadingFatJet, 'pt_toUse', 'eta', 'phi', 'mass_toUse')   
        LV_leadingFatJet_wMsoftdrop = getLorentVector(leadingFatJet, 'pt_toUse', 'eta', 'phi', 'msoftdrop_toUse') 
        if  'PNet_X4b_v2a_Haa34b_score' in leadingFatJet.fields: # NanoAOD v2
            LV_leadingFatJet_wMpnet = getLorentVector(leadingFatJet, 'pt_toUse', 'eta', 'phi', 'PNet_massH_v2b')
        else:
            LV_leadingFatJet_wMpnet =  LV_leadingFatJet_wMsoftdrop
        
        
        
        #if runMode_OptimizePNetTaggerCut:
        #    leadingFatJet = ak.firsts(FatJetsToUse[idx_FatJet_ZHbb_plus_Xbb_max])

        #leadingFatJet = ak.firsts(FatJetsToUse)        
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
        if 'particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields:
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

            leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD = ak.where(
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
        elif 'PNet_X4b_v2a_Haa34b_score' in FatJetsToUse.fields: # NanoAOD v2 PNet_X4b_v1_Haa4b_vs_QCD
            leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD = leadingFatJet.PNet_X4b_v1_Haa4b_vs_QCD
            leadingFatJet_PNet_Xto4bv2_Htoaa4b        = (leadingFatJet.PNet_X4b_v2a_Haa4b_score + \
                                                         leadingFatJet.PNet_X4b_v2b_Haa4b_score) / 2.0
            leadingFatJet_PNet_Xto34bv2_Htoaa4b       = (leadingFatJet.PNet_X4b_v2a_Haa34b_score + \
                                                         leadingFatJet.PNet_X4b_v2b_Haa34b_score) / 2.0
            


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
        

        ## mask satisfying HEM1516 issues conditions
        '''
        # Iteration 1
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
        '''
        # Iteration 2: Andrew's suggestions https://indico.cern.ch/event/1479951/contributions/6234638/attachments/2968060/5255895/2024_11_15_HToAATo4B_selection_catgories_NanoAODTools.pdf#page=14
        mask_HEM1516Issue = ak.fill_none((
            (leadingFatJet.eta < -1.1) & 
            (np.abs(leadingFatJet.phi + 1.22) < 0.55)
        ), False)
        mask_HEM1516Issue_Eta = ak.fill_none((
            (leadingFatJet.eta < -1.1)
        ), False)
        mask_HEM1516Issue_Phi = ak.fill_none((
            (np.abs(leadingFatJet.phi + 1.22) < 0.55)
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
                    (FatJetsToUse.delta_r(LVGenB_0)    < 0.8) &
                    (FatJetsToUse.delta_r(LVGenBbar_0) < 0.8) &
                    (FatJetsToUse.delta_r(LVGenB_1)    < 0.8) &
                    (FatJetsToUse.delta_r(LVGenBbar_1) < 0.8)                     
                )
                mask_events_FatJet_matched_genB_HToAATo4B = ak.any(mask_FatJet_matched_genB_HToAATo4B, axis=1) # events with FatJet matches to 4 GEN B-quarks from HToAATo4B
                idx_FatJet_matched_genB_HToAATo4B = ak.argmax(mask_FatJet_matched_genB_HToAATo4B, axis=-1) # index of FatJet within events that maches to 4 GEN B-quarks from HToAATo4B
                

    
        ## non-HTo4B FatJet
        nonHto4bFatJet = selFatJets[(selFatJets.delta_r(leadingFatJet) > 0.8)]
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
        '''
        leadingNonHto4bFatJet_DeepTag_WvsQCD     = array_PutLowerBound(leadingNonHto4bFatJet.deepTag_WvsQCD, kSmallPositiveNumber)
        leadingNonHto4bFatJet_DeepTag_ZvsQCD     = array_PutLowerBound(leadingNonHto4bFatJet.deepTag_ZvsQCD, kSmallPositiveNumber)
        leadingNonHto4bFatJet_DeepTag_QCD        = array_PutLowerBound(leadingNonHto4bFatJet.deepTag_QCD,    kSmallPositiveNumber)
        leadingNonHto4bFatJet_DeepTag_W          = cal_W_(leadingNonHto4bFatJet_DeepTag_WvsQCD, leadingNonHto4bFatJet_DeepTag_QCD)
        leadingNonHto4bFatJet_DeepTag_Z          = cal_Z_(leadingNonHto4bFatJet_DeepTag_ZvsQCD, leadingNonHto4bFatJet_DeepTag_QCD)
        leadingNonHto4bFatJet_DeepTag_WZvsQCD    = cal_WZQ_1_(leadingNonHto4bFatJet_DeepTag_W, leadingNonHto4bFatJet_DeepTag_Z, leadingNonHto4bFatJet_DeepTag_QCD)
        leadingNonHto4bFatJet_DeepTag_WZvsQCD2   = cal_WZQ_2_(leadingNonHto4bFatJet_DeepTag_WvsQCD, leadingNonHto4bFatJet_DeepTag_ZvsQCD)
        leadingNonHto4bFatJet_DeepTag_VvsQCD_max = calculateMaxOfTwoArrays(leadingNonHto4bFatJet_DeepTag_WvsQCD, leadingNonHto4bFatJet_DeepTag_ZvsQCD)
        leadingNonHto4bFatJet_DeepTag_V_max      = calculateMaxOfTwoArrays(leadingNonHto4bFatJet_DeepTag_W, leadingNonHto4bFatJet_DeepTag_Z)
        '''
        if 'particleNet_WZvsQCD' in events.FatJet.fields:
            leadingNonHto4bFatJet_PNet_WZvsQCD = leadingNonHto4bFatJet.particleNet_WZvsQCD

        nNonHto4bFatJet  = ak.fill_none(ak.count(nonHto4bFatJet.pt_toUse, axis=1), 0)
        nleadingNonHto4bFatJet_WZvsQCD = ak.fill_none(ak.where(
            (leadingNonHto4bFatJet_PNet_WZvsQCD > self.objectSelector.NonHto4bFatJetPNet_WZvsQCD_Thsh),
            np.full_like(leadingNonHto4bFatJet_PNet_WZvsQCD, 1),
            np.full_like(leadingNonHto4bFatJet_PNet_WZvsQCD, 0),
        ), 0)
        

        ## sel leptons
        '''Old selection
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
        '''
        
        muonsTight     = selectMuons(
            events.Muon, 
            pT_Thsh            = self.objectSelector.MuonPtThsh, 
            MiniPFRelIsoIdThsh = self.objectSelector.MuonMiniPFRelIsoId,
            DxyThsh            = self.objectSelector.MuonDxyThsh,
            DzThsh             = self.objectSelector.MuonDzThsh )
        electronsTight = selectElectrons(
            events.Electron, 
            pT_Thsh            = self.objectSelector.ElectronPtThsh,
            DxyThsh            = self.objectSelector.ElectronDxyThsh,
            DzThsh             = self.objectSelector.ElectronDzThsh )
        leptonsTight   = ak.concatenate([muonsTight, electronsTight], axis=1)
        nLeptonsTight  = ak.fill_none(ak.count(leptonsTight.pt, axis=1), 0)
        #nLeptons_matched_leadingFatJet = ak.fill_none(ak.sum(leadingFatJet.metric_table( leptonsTight, axis=None ) < 0.8, axis=1), 0)
        nLeptons_matched_leadingFatJet    = ak.fill_none(ak.sum(leadingFatJet.delta_r( leptonsTight ) < 0.8, axis=1), 0)
        nLeptons_nonoverlap_leadingFatJet = ak.fill_none(ak.sum(leadingFatJet.delta_r( leptonsTight ) > 0.8, axis=1), 0)
        

        ## sel AK4 jets
        ## 'AK8' Higgs category
        # Don't consider ak4JetHiggsProng yet
        ak4Jets = selectAK4Jets(Jets=JetsToUse, era=self.datasetInfo["era"], pT_Thsh = 30)
        mask_ak4Jets_nonoverlaping_leadingFatJet           = ak4Jets.delta_r(leadingFatJet) > 0.8
        ak4Jets_nonoverlaping_leadingFatJet                = ak4Jets[ mask_ak4Jets_nonoverlaping_leadingFatJet ]        
        nAk4Jets_nonoverlaping_leadingFatJet               = ak.fill_none(ak.count(ak4Jets_nonoverlaping_leadingFatJet.eta, axis=1), 0)

        ak4JetsCentral_nonoverlaping_leadingFatJet         = ak4Jets_nonoverlaping_leadingFatJet[abs(ak4Jets_nonoverlaping_leadingFatJet.eta) < 2.4]
        #nAk4JetsCentral_nonoverlaping_leadingFatJet        = ak.fill_none(ak.count(ak4JetsCentral_nonoverlaping_leadingFatJet.eta, axis=1), 0)

        #mask_ak4Jets_bTag_nonoverlaping_leadingFatJet      = ak4Jets_nonoverlaping_leadingFatJet.btagDeepFlavB > self.objectSelector.Ak4JetDeepJetB_Thsh
        #ak4Jets_bTag_nonoverlaping_leadingFatJet           = ak4Jets_nonoverlaping_leadingFatJet[mask_ak4Jets_bTag_nonoverlaping_leadingFatJet]
        #nAk4Jets_bTag_nonoverlaping_leadingFatJet          = ak.fill_none(ak.count(ak4Jets_bTag_nonoverlaping_leadingFatJet.eta, axis=1), 0)

        #ak4JetsCentral_bTag_nonoverlaping_leadingFatJet    = ak4Jets_bTag_nonoverlaping_leadingFatJet[abs(ak4Jets_bTag_nonoverlaping_leadingFatJet.eta) < 2.4]
        ak4JetsCentral_bTag_nonoverlaping_leadingFatJet    = ak4JetsCentral_nonoverlaping_leadingFatJet[(
            (ak4JetsCentral_nonoverlaping_leadingFatJet.btagDeepFlavB > self.objectSelector.Ak4JetDeepJetB_Thsh)
        )]
        nAk4JetsCentral_bTag_nonoverlaping_leadingFatJet   = ak.fill_none(ak.count(ak4JetsCentral_bTag_nonoverlaping_leadingFatJet.eta, axis=1), 0)

        #mask_ak4Jets_nonbTag_nonoverlaping_leadingFatJet      = ak4Jets_nonoverlaping_leadingFatJet.btagDeepFlavB <= self.objectSelector.Ak4JetDeepJetB_Thsh
        #ak4Jets_nonbTag_nonoverlaping_leadingFatJet           = ak4Jets_nonoverlaping_leadingFatJet[mask_ak4Jets_nonbTag_nonoverlaping_leadingFatJet]
        #nAk4Jets_nonbTag_nonoverlaping_leadingFatJet          = ak.fill_none(ak.count(ak4Jets_nonbTag_nonoverlaping_leadingFatJet.eta, axis=1), 0)

        

        ## VBF jj
        ## 'AK8' Higgs category
        # Don't consider ak4JetHiggsProng yet
        # ak4Jets_nonoverlaping_leadingFatJet
        leading2Ak4Jets_nonoverlaping_leadingFatJet = ak.mask(ak4Jets_nonoverlaping_leadingFatJet, nAk4Jets_nonoverlaping_leadingFatJet >= 2)
        LV_leading2Ak4Jets_nonoverlaping_leadingFatJet = getLorentVector(leading2Ak4Jets_nonoverlaping_leadingFatJet, 'pt_toUse', 'eta', 'phi', 'mass_toUse')
        mass_leadingPair_ak4Jets_nonoverlaping_leadingFatJet = ak.fill_none(
            (LV_leading2Ak4Jets_nonoverlaping_leadingFatJet[:, 0] + LV_leading2Ak4Jets_nonoverlaping_leadingFatJet[:, 1]).mass
        , 0)
        dEta_leadingPair_ak4Jets_nonoverlaping_leadingFatJet = ak.fill_none(
            abs(LV_leading2Ak4Jets_nonoverlaping_leadingFatJet[:, 0].eta - LV_leading2Ak4Jets_nonoverlaping_leadingFatJet[:, 1].eta)
        , 0)
        
        




        #####################
        # EVENT SELECTION
        #####################
        
        # reconstruction level cuts for cut-flow table. Order of cuts is IMPORTANT
        cuts_reco = ["dR_LeadingFatJet_GenB_0p8"] + self.sel_names_all["Presel"] #.copy()

       
        # create a PackedSelection object
        # this will help us later in composing the boolean selections easily
        selection = PackedSelection(dtype='uint32')

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


        if "candH" in self.sel_conditions_all_list:
            selection.add(
                "candH",
                (~ak.is_none(leadingFatJet)) # no. of candH >= 1
            )

        if "leadingFatJetPt" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetPt",
                ((leadingFatJet.pt_toUse > self.objectSelector.FatJetPt_gg0lIncl_MinThsh) &
                 (leadingFatJet.pt_toUse <=  self.objectSelector.FatJetPt_gg0lIncl_MaxThsh))
            )
        if "leadingFatJetPt_gg0lIncl" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetPt_gg0lIncl",
                ((leadingFatJet.pt_toUse > self.objectSelector.FatJetPt_gg0lIncl_MinThsh) &
                 (leadingFatJet.pt_toUse <=  self.objectSelector.FatJetPt_gg0lIncl_MaxThsh))
            )
        if "leadingFatJetPt_gg0lLo" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetPt_gg0lLo",
                ((leadingFatJet.pt_toUse > self.objectSelector.FatJetPt_gg0lLo_MinThsh) &
                 (leadingFatJet.pt_toUse <=  self.objectSelector.FatJetPt_gg0lLo_MaxThsh))
            )
        if "leadingFatJetPt_gg0lHi" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetPt_gg0lHi",
                ((leadingFatJet.pt_toUse > self.objectSelector.FatJetPt_gg0lHi_MinThsh) &
                 (leadingFatJet.pt_toUse <=  self.objectSelector.FatJetPt_gg0lHi_MaxThsh))
            )
            


        if "leadingFatJetEta" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetEta",
                abs(leadingFatJet.eta) < self.objectSelector.FatJetEtaThsh
            )

        if "JetID"  in self.sel_conditions_all_list:
            selection.add(
                "JetID", 
                leadingFatJet.jetId == self.objectSelector.FatJetJetID
            )

 
        if "leadingFatJetMSoftDrop"  in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetMSoftDrop",
                (leadingFatJet.msoftdrop_toUse > self.objectSelector.FatJetMSoftDropThshLow) &
                (leadingFatJet.msoftdrop_toUse < self.objectSelector.FatJetMSoftDropThshHigh)
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

        if "leadingFatJetXbbVsQCD" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetXbbVsQCD",
                leadingFatJet.particleNetMD_XbbvsQCD > self.objectSelector.FatJetParticleNetMD_XbbvsQCD_Thsh
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

        if "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD",
                leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD > self.objectSelector.FatJetPNet_Xto4bv1_Htoaa4bOverQCD_Thsh
            )

        if "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP40" in self.sel_conditions_all_list and \
            (('particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields) or ('PNet_X4b_v1_Haa4b_vs_QCD' in FatJetsToUse.fields) ):
            selection.add(
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP40",
                leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD > bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv1_Htoaa4bOverQCD']['WP-40']
            )

        if "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP60" in self.sel_conditions_all_list and 'particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields:
            selection.add(
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP60",
                leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD > bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv1_Htoaa4bOverQCD']['WP-60']
            )

        if "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP80" in self.sel_conditions_all_list and 'particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields:
            selection.add(
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP80",
                leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD > bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv1_Htoaa4bOverQCD']['WP-80']
            )

        if "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP80to40" in self.sel_conditions_all_list and 'particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields:
            selection.add(
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP80to40",
                ( (leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD >  bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv1_Htoaa4bOverQCD']['WP-80']) &
                  (leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD <= bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv1_Htoaa4bOverQCD']['WP-40']) )
            )

        if "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP95to60" in self.sel_conditions_all_list and 'particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields:
            selection.add(
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP95to60",
                ( (leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD >  bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv1_Htoaa4bOverQCD']['WP-95']) &
                  (leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD <= bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv1_Htoaa4bOverQCD']['WP-60']) )
            )

        if "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP99to80" in self.sel_conditions_all_list and 'particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields:
            selection.add(
                "leadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD_WP99to80",
                ( (leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD >  bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv1_Htoaa4bOverQCD']['WP-99']) &
                  (leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD <= bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv1_Htoaa4bOverQCD']['WP-80']) )
            )

        # PNet Xaa4b v2
        if 'PNet_X4b_v2a_Haa4b_score' in FatJetsToUse.fields:
            for wp_ in self.objectSelector.FatJetPNetXto4bv2WorkingPoints: 
                # leadingFatJetPNet_Xto4bv2_Htoaa4b tagger cut
                if "leadingFatJetPNet_Xto4bv2_Htoaa4b_SRWP%s" % (wp_) in self.sel_conditions_all_list:
                    selection.add(
                        "leadingFatJetPNet_Xto4bv2_Htoaa4b_SRWP%s" % (wp_),
                        leadingFatJet_PNet_Xto4bv2_Htoaa4b > bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv2_Htoaa4b']['SRWP-%s' % wp_]
                    )
                if "leadingFatJetPNet_Xto4bv2_Htoaa4b_SBWP%s" % (wp_) in self.sel_conditions_all_list:
                    selection.add(
                        "leadingFatJetPNet_Xto4bv2_Htoaa4b_SBWP%s" % (wp_),
                        ( (leadingFatJet_PNet_Xto4bv2_Htoaa4b >  bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv2_Htoaa4b']['SBWP-%s' % wp_]) &
                          (leadingFatJet_PNet_Xto4bv2_Htoaa4b <= bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv2_Htoaa4b']['SRWP-%s' % wp_]))
                    )


                

            

           
        if "nLeptonsTight" in self.sel_conditions_all_list:
            selection.add(
                "nLeptonsTight",
                ( nLeptonsTight <= self.objectSelector.NLeptonsTight_MaxThsh )
            )

        if "nonHto4bFatJetVjjVeto" in self.sel_conditions_all_list:
            selection.add(
                "nonHto4bFatJetVjjVeto",
                ( nleadingNonHto4bFatJet_WZvsQCD <= self.objectSelector.NNonHo4bFatJetPNet_WZvsQCD_MaxThsh )
            )

        if "MetZvvVeto" in self.sel_conditions_all_list:
            #selection.add(
            #    "MetZvvVeto",
            #    ~ ( (METToUse.pt_toUse >= self.objectSelector.METPt_ZvvIncl_MinThsh) & 
            #        (METToUse.pt_toUse <  self.objectSelector.METPt_ZvvIncl_MaxThsh) & 
            #        (abs(METToUse.delta_phi(leadingFatJet)) > self.objectSelector.DPhi_FJHto4b_MET_MinThsh) )
            #)
            selection.add(
                "MetZvvVeto",
                ~ ( (METToUse.pt_toUse > self.objectSelector.METPt_ZvvIncl_MinThsh) & 
                    (METToUse.pt_toUse < self.objectSelector.METPt_ZvvIncl_MaxThsh) )
            )



        if "DijetVBFVeto" in self.sel_conditions_all_list:
            selection.add(
                "DijetVBFVeto",
                ~ ( (mass_leadingPair_ak4Jets_nonoverlaping_leadingFatJet > self.objectSelector.VBFDijetMass_MinThsh) & 
                    (dEta_leadingPair_ak4Jets_nonoverlaping_leadingFatJet  > self.objectSelector.VBFDijetEta_MinThsh)   )
            )
        
        if "BJetVeto" in self.sel_conditions_all_list:
            selection.add(
                "BJetVeto",
                nAk4JetsCentral_bTag_nonoverlaping_leadingFatJet == 0
            )

        

        for LumiSecSelThsh in LumiSecSelThsh_list:
            sLSSelCut = "LSlt%d"%(LumiSecSelThsh)
            if sLSSelCut in self.sel_conditions_all_list:
                selection.add(
                    sLSSelCut,
                    events.luminosityBlock < LumiSecSelThsh
                )            




        # Trigger selection
        if sTrgSelection in self.sel_conditions_all_list:
            if sTrgSelection not in Triggers_perEra[self.datasetInfo["era"]]:
                logging.critical(f'htoaa_Analysis_GGFMode.py::main():: {sTrgSelection = } not in {Triggers_perEra[self.datasetInfo["era"]] = }.')
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
            
        
        #print(f'{self.sel_names_all["Presel"] = }')
        #sel_SR          = selection.all("nPV", "FatJetGet")
        sel_SR           = selection.all(* self.sel_names_all["Presel"])
        sel_GenHToAATo4B = None

        if self.datasetInfo['isMC'] and self.datasetInfo['isSignal'] and 1==0:
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
            
            sel_GenHToAATo4B = selection.all(* self.sel_names_all["GenHToAATo4B"])
        
        
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
        weights              = Weights(len(events), storeIndividual=True)
        weights_gen          = Weights(len(events))
        weights_GenHToAATo4B = Weights(len(events))
        weights_woHEM1516Fix = Weights(len(events))
        

        if self.datasetInfo["isMC"]:
            # lumiScale ------------------------------------
            lumiScale_toUse = None
            if self.datasetInfo["MCSamplesStitchOption"] == MCSamplesStitchOptions.PhSpOverlapRewgt and \
               self.datasetInfo['isQCD']:
                mask_PhSp_dict_ = {
                    "QCD_bEnrich":      mask_QCD_bEnrich_PhSp,
                    "QCD_bGen":         mask_QCD_bGen_PhSp,
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
            #wgt_PU = getPURewgts(
            #    PU_list  = events.Pileup.nTrueInt,
            #    hPURewgt = self.hPURewgt
            #)
            wgt_PU, wgt_PUUp, wgt_PUDown = getPURewgts_variation(
                events = events,
                year   = self.datasetInfo["era"]
            )

            # MC TrgEff ------------------------------------
            #wgt_TrgEff, wgt_TrgEffUp, wgt_TrgEffDown  = get_jetTriggerSF(
            #    events = events,
            #    year   = self.datasetInfo["era"],
            #    selection = selection
            #)

            # MC HToAATo4B signal LundPlane reweighting
            wgt_LundPlane_Nom = wgt_LundPlane_Up = wgt_LundPlane_Down = None 
            if self.datasetInfo['isSignal']:
                wgt_LundPlane_Nom, wgt_LundPlane_Up, wgt_LundPlane_Down = getHToAATo4BLundPlaneRewgt(
                    events = events
                )
                #printVariable('LundPlaneWeights: ', ak.zip([wgt_LundPlane_Nom, wgt_LundPlane_Up, wgt_LundPlane_Down]))

            # MC GGF HToAATo4B Higgs pT reweight
            wgt_GGH_HiggsPt = wgt_GGH_HiggsPtUp = wgt_GGH_HiggsPtDown = None
            if self.datasetInfo['isSignalGGH']:
                wgt_GGH_HiggsPt, wgt_GGH_HiggsPtUp, wgt_GGH_HiggsPtDown = getHiggsPtRewgtForGGToHToAATo4B(
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
                wgt_TopPt, wgt_TopPtUp, wgt_TopPtDown = getTopPtRewgt(
                    eventsGenPart = events.GenPart[mask_genTopQuark],
                    isPythiaTuneCP5 = self.datasetInfo['isPythiaTuneCP5']
                )   
                #printVariable('wgt_TopPt: ', ak.zip([wgt_TopPt, wgt_TopPtUp, wgt_TopPtDown]))


            # MC ParticleNetMD_XbbvsQCD SFs      SFs_ParticleNetMD_XbbvsQCD
            #if "leadingFatJetParticleNetMD_XbbvsQCD" in self.sel_names_all["Presel"]:
            if self.SFs_ParticleNetMD_XbbvsQCD != None:
                print(f"{self.SFs_ParticleNetMD_XbbvsQCD = } is not None")
                mask_ParticleNetMD_XbbvsQCD_SFRegion = (
                    (n_leadingFatJat_matched_genB >= 2) &
                    (leadingFatJetParticleNetMD_XbbvsQCD > self.objectSelector.FatJetParticleNetMD_XbbvsQCD_Thsh)
                )
                wgt_ParticleNetMD_XbbvsQCD = ak.fill_none( 
                    ak.where(
                        mask_ParticleNetMD_XbbvsQCD_SFRegion,
                        self.SFs_ParticleNetMD_XbbvsQCD(leadingFatJet.pt_toUse),
                        ones_list
                    ), 
                    1
                )

            # MC Parton Shower weights
            wgt_PS_Nom, wgt_PS_ISRUp, wgt_PS_ISRDown, wgt_PS_FSRUp, wgt_PS_FSRDown  = get_PSWeight(
                events = events,
                dataset = self.datasetInfo['datasetName']
            )

            # MC QCD PDF uncertianty
            wgt_QCDPdfNom, wgt_QCDPdfUp, wgt_QCDPdfDown = add_pdf_as_weight(
                events = events,
                #pdf_weights = events.LHEPdfWeight,
                dataset = self.datasetInfo['datasetName']
            )

            # MC QCD alphaS renormalization and factorization uncertainty
            wgt_QCDScale_Nom, wgt_QCDScale_RenormUp, wgt_QCDScale_RenormDown, wgt_QCDScale_FactorizationUp, wgt_QCDScale_FactorizationDown = get_QCDScaleWeight(
                events = events,
                dataset = self.datasetInfo['datasetName']
            )
                
            # btagSF
            wgt_Ak4Btag_dict = get_Ak4BtagSF(
                jet         = ak4JetsCentral_nonoverlaping_leadingFatJet, 
                btagWPThsh  = self.objectSelector.Ak4JetDeepJetB_Thsh,
                year = self.datasetInfo["era"]
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
                "PU",
                weight     = wgt_PU,
                weightUp   = wgt_PUUp,
                weightDown = wgt_PUDown
            )
            #weights.add(
            #    "TrgEff",
            #    weight     = wgt_TrgEff,
            #    weightUp   = wgt_TrgEffUp,
            #    weightDown = wgt_TrgEffDown
            #)
            
            if self.datasetInfo['isSignal']:
                weights.add(
                    "LPRewgt",
                    weight     = wgt_LundPlane_Nom,
                    weightUp   = wgt_LundPlane_Up,
                    weightDown = wgt_LundPlane_Down
                )
            if self.datasetInfo['isSignalGGH']:
                weights.add(
                    "GGHPtRewgt",
                    weight     = wgt_GGH_HiggsPt,
                    weightUp   = wgt_GGH_HiggsPtUp,
                    weightDown = wgt_GGH_HiggsPtDown
                )
            if self.datasetInfo['isQCD_bGen']:
                weights.add(
                    "HTRewgt",
                    weight = wgt_HT
                )            
            if self.datasetInfo['isTTbar']:
                weights.add(
                    "TopPtReWeight",
                    weight     = wgt_TopPt,
                    weightUp   = wgt_TopPtUp,
                    weightDown = wgt_TopPtDown                    
                )         
            #if "leadingFatJetParticleNetMD_XbbvsQCD" in self.sel_names_all["Presel"]:
            if self.SFs_ParticleNetMD_XbbvsQCD != None:
                weights.add(
                    "SF_ParticleNetMD_XbbvsQCD",
                    weight = wgt_ParticleNetMD_XbbvsQCD
                )
            
            weights.add(
                "ISR",
                weight     = wgt_PS_Nom,
                weightUp   = wgt_PS_ISRUp,
                weightDown = wgt_PS_ISRDown
            )
            weights.add(
                "FSR",
                weight     = wgt_PS_Nom,
                weightUp   = wgt_PS_FSRUp,
                weightDown = wgt_PS_FSRDown
            )
            weights.add(
                "QCDRenorm",
                weight     = wgt_QCDScale_Nom,
                weightUp   = wgt_QCDScale_RenormUp,
                weightDown = wgt_QCDScale_RenormDown
            )
            weights.add(
                "QCDFactr",
                weight     = wgt_QCDScale_Nom,
                weightUp   = wgt_QCDScale_FactorizationUp,
                weightDown = wgt_QCDScale_FactorizationDown
            )            
            weights.add(
                "PDF",
                weight     = wgt_QCDPdfNom,
                weightUp   = wgt_QCDPdfUp,
                weightDown = wgt_QCDPdfDown
            )
            if  kDatasetToAnalyze == DatasetToAnalyze.SingleYear: ## btag 
                weights.add(
                    "Btag",
                    weight     = wgt_Ak4Btag_dict['Nom'],
                    weightUp   = wgt_Ak4Btag_dict['Up'],
                    weightDown = wgt_Ak4Btag_dict['Down']
                )
            elif kDatasetToAnalyze == DatasetToAnalyze.FullRun2:
                weights.add(
                    "BtagUncorr",
                    weight     = wgt_Ak4Btag_dict['Nom'],
                    weightUp   = wgt_Ak4Btag_dict['Upuncorrelated'],
                    weightDown = wgt_Ak4Btag_dict['Downuncorrelated']
                )
                weights.add(
                    "BtagCorr",
                    weight     = ones_list, #wgt_Ak4Btag_dict['Nom'],  #<<<<<< use dummy weights here to avoid application of btag wgt twice
                    weightUp   = wgt_Ak4Btag_dict['Upcorrelated'],
                    weightDown = wgt_Ak4Btag_dict['Downcorrelated']
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
                "PU",
                weight = wgt_PU,
                weightUp = wgt_PUUp,
                weightDown = wgt_PUDown
            )
            #weights_woHEM1516Fix.add(
            #    "TrgEff",
            #    weight     = wgt_TrgEff,
            #    weightUp   = wgt_TrgEffUp,
            #    weightDown = wgt_TrgEffDown
            #)
            if self.datasetInfo['isSignal']:
                weights_woHEM1516Fix.add(
                    "LPRewgt",
                    weight     = wgt_LundPlane_Nom,
                    weightUp   = wgt_LundPlane_Up,
                    weightDown = wgt_LundPlane_Down
                )
            if self.datasetInfo['isSignalGGH']:
                weights_woHEM1516Fix.add(
                    "GGHPtRewgt",
                    weight     = wgt_GGH_HiggsPt,
                    weightUp   = wgt_GGH_HiggsPtUp,
                    weightDown = wgt_GGH_HiggsPtDown
                )
            if self.datasetInfo['isQCD_bGen']:
                weights_woHEM1516Fix.add(
                    "HTRewgt",
                    weight = wgt_HT
                )  
            if self.datasetInfo['isTTbar']:
                weights_woHEM1516Fix.add(
                    "TopPtReWeight",
                    weight     = wgt_TopPt,
                    weightUp   = wgt_TopPtUp,
                    weightDown = wgt_TopPtDown                    
                )         
            #if "leadingFatJetParticleNetMD_XbbvsQCD" in self.sel_names_all["Presel"]:
            if self.SFs_ParticleNetMD_XbbvsQCD != None:
                weights_woHEM1516Fix.add(
                    "SF_ParticleNetMD_XbbvsQCD",
                    weight = wgt_ParticleNetMD_XbbvsQCD
                )
            weights_woHEM1516Fix.add(
                "ISR",
                weight     = wgt_PS_Nom,
                weightUp   = wgt_PS_ISRUp,
                weightDown = wgt_PS_ISRDown
            )
            weights_woHEM1516Fix.add(
                "FSR",
                weight     = wgt_PS_Nom,
                weightUp   = wgt_PS_FSRUp,
                weightDown = wgt_PS_FSRDown
            )
            weights_woHEM1516Fix.add(
                "QCDRenorm",
                weight     = wgt_QCDScale_Nom,
                weightUp   = wgt_QCDScale_RenormUp,
                weightDown = wgt_QCDScale_RenormDown
            )
            weights_woHEM1516Fix.add(
                "QCDFactr",
                weight     = wgt_QCDScale_Nom,
                weightUp   = wgt_QCDScale_FactorizationUp,
                weightDown = wgt_QCDScale_FactorizationDown
            )            
            weights_woHEM1516Fix.add(
                "PDF",
                weight     = wgt_QCDPdfNom,
                weightUp   = wgt_QCDPdfUp,
                weightDown = wgt_QCDPdfDown
            )
            if  kDatasetToAnalyze == DatasetToAnalyze.SingleYear: ## btag 
                weights_woHEM1516Fix.add(
                    "Btag",
                    weight     = wgt_Ak4Btag_dict['Nom'],
                    weightUp   = wgt_Ak4Btag_dict['Up'],
                    weightDown = wgt_Ak4Btag_dict['Down']
                )
            elif kDatasetToAnalyze == DatasetToAnalyze.FullRun2:
                weights_woHEM1516Fix.add(
                    "BtagUncorr",
                    weight     = wgt_Ak4Btag_dict['Nom'],
                    weightUp   = wgt_Ak4Btag_dict['Upuncorrelated'],
                    weightDown = wgt_Ak4Btag_dict['Downuncorrelated']
                )
                weights_woHEM1516Fix.add(
                    "BtagCorr",
                    weight     = ones_list, #wgt_Ak4Btag_dict['Nom'],  #<<<<<< use dummy weights here to avoid application of btag wgt twice
                    weightUp   = wgt_Ak4Btag_dict['Upcorrelated'],
                    weightDown = wgt_Ak4Btag_dict['Downcorrelated']
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
                    weight     = wgt_GGH_HiggsPt,
                    weightUp   = wgt_GGH_HiggsPtUp,
                    weightDown = wgt_GGH_HiggsPtDown
                )
            if self.datasetInfo['isQCD_bGen']:
                weights_gen.add(
                    "HTRewgt",
                    weight = wgt_HT
                )
            if self.datasetInfo['isTTbar']:
                weights_gen.add(
                    "TopPtReWeight",
                    weight     = wgt_TopPt,
                    weightUp   = wgt_TopPtUp,
                    weightDown = wgt_TopPtDown                    
                )         
            








        ###################
        # FILL HISTOGRAMS
        ###################

        systList = []
        if self.datasetInfo['isMC']:
            if shift_syst is None:
                systList = [
                    "Nom",
                ]
                if not self.datasetInfo['systematicsToRun'] == 'no':
                    if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['pu', 'full'] ): 
                        systList.extend( [
                            "PUUp",
                            "PUDown",
                        ] )
                    #if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['trg', 'full'] ):
                    #    systList.extend( [
                    #        "TrgEffUp",
                    #        "TrgEffDown",
                    #    ] )
                    if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['lundplane', 'full'] ) and self.datasetInfo['isSignal']:
                        systList.extend( [
                            "LPRewgtUp",
                            "LPRewgtDown",
                        ] ) 
                    if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['gghptrewgt', 'full'] ) and self.datasetInfo['isSignalGGH']:
                        systList.extend( [
                            "GGHPtRewgtUp",
                            "GGHPtRewgtDown",
                        ] ) 
                    if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['topptrewgt', 'full'] ) and self.datasetInfo['isTTbar']:
                        systList.extend( [
                            "TopPtReWeightUp",
                            "TopPtReWeightDown",
                        ] ) 
                    if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['isr', 'full'] ):
                        systList.extend( [
                            "ISRUp",
                            "ISRDown",
                        ] )
                    if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['fsr', 'full'] ):
                        systList.extend( [
                            "FSRUp",
                            "FSRDown",
                        ] )
                    if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['qcdrenorm', 'full'] ):
                        systList.extend( [
                            "QCDRenormUp",
                            "QCDRenormDown",
                        ] )
                    if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['qcdfactr', 'full'] ):
                        systList.extend( [
                            "QCDFactrUp",
                            "QCDFactrDown",
                        ] )                        
                    if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['pdf', 'full'] ):
                        systList.extend( [
                            "PDFUp",
                            "PDFDown",
                        ] )
                    if stringHasSubstring(self.datasetInfo['systematicsToRun'], ['btag', 'full'] ):
                        if  kDatasetToAnalyze == DatasetToAnalyze.SingleYear:
                            systList.extend( [
                                "BtagUp",
                                "BtagDown",
                            ] ) 
                        elif kDatasetToAnalyze == DatasetToAnalyze.FullRun2:
                            systList.extend( [
                                "BtagUncorrUp",
                                "BtagUncorrDown",
                                "BtagCorrUp",
                                "BtagCorrDown",                                
                            ] )                                           
                    
                    
                
            else:
                systList = [shift_syst]
        else:
            systList = ["noweight"]

            
        if shift_syst is None:
            output['cutflow']['all events'] += len(events)
            output['cutflow'][sWeighted+'all events'] += weights.weight().sum() 
            #for n in selection.names:
            #    output['cutflow'][n] += selection.all(n).sum()

            
            for iSelection in self.sel_names_all.keys():
                iName = f"{iSelection}: {self.sel_names_all[iSelection]}"
                sel_i = selection.all(* self.sel_names_all[iSelection])
                selection.names
                output['cutflow'][iName] += sel_i.sum()
                output['cutflow'][sWeighted+iName] +=  weights.weight()[sel_i].sum()



        if printLevel >= 100:
            print(f"{weights.variations = }")



        for syst in systList:

            # find the event weight to be used when filling the histograms
            weightSyst = syst
            
            # in the case of 'central', or the jet energy systematics, no weight systematic variation is used (weightSyst=None)
            if syst in ["Nom", "JERUp", "JERDown", "JESUp", "JESDown", "JESHEMIssueUp", "JESHEMIssueDown"]:
                weightSyst = None

            
            
            if syst == "noweight":
                evtWeight                = np.ones(len(events))
                evtWeight_woHEM1516Fix   = np.ones(len(events))
            else:
                evtWeight                = weights.weight(weightSyst)
                evtWeight_woHEM1516Fix   = weights_woHEM1516Fix.weight(weightSyst)
                if syst == "Nom":
                    evtWeight_gen            = weights_gen.weight(weightSyst)
                

            #if printLevel >=0:
            #    printVariable('\n evtWeight %s' % (syst), evtWeight)


            ### General or GEN-level histograms ========================================================================================

            # PU
            if histogramSaveLevel >= 10:
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
            if self.datasetInfo['isMC'] and runMode_GenLHEPlots and syst == "Nom":                 
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


            if (self.datasetInfo['isSignal'] or self.datasetInfo['isHToBB']) and runMode_SignalGenChecks and syst == "Nom":
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

            if self.datasetInfo['isSignal'] and runMode_SignalGenChecks and syst == "Nom": 
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
            if self.datasetInfo['isQCD'] and runMode_QCDGenValidation and syst == "Nom":
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
            
            if self.datasetInfo['isSignal'] and histogramSaveLevel >= 12:
                output['hIdxFatJetMatchedToGenBFromHToAATo4B'].fill(
                    dataset=dataset,
                    nObject=(idx_FatJet_matched_genB_HToAATo4B[mask_events_FatJet_matched_genB_HToAATo4B]),
                    systematic=syst,
                    weight=evtWeight[mask_events_FatJet_matched_genB_HToAATo4B]
                ) 
                if 'particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields:
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
                    Ratio=( leadingFatJet[mask_events_FatJet_matched_genB_HToAATo4B].pt_toUse / ak.firsts(FatJetsToUse)[mask_events_FatJet_matched_genB_HToAATo4B].pt_toUse),
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
                #printVariable('\nwgt_LundPlane %s'%(sel_name), ak.zip([wgt_LundPlane_Nom, wgt_LundPlane_Up, wgt_LundPlane_Down ]))
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

                    
                    sel_SR_forHExt_woGenMatch = sel_SR_forHExt
                    #if runMode_OptimizePNetTaggerCut and self.datasetInfo['isSignal']:
                    #if self.datasetInfo['isSignal']:
                    if (self.datasetInfo['isSignal'] and \
                        (runMode_SignificancsScan2D or \
                         runMode_OptimizePNetTaggerCut or \
                         runMode_SignalGenCuts) ):
                        sel_SR_forHExt = sel_SR_forHExt & (n_leadingFatJat_matched_genB_HToAATo4B >= 4)
                    sel_SR_forHExt = ak.fill_none(sel_SR_forHExt, False) 
                    

                    if histogramSaveLevel >= 0:
                        # Cut flow table ------------------------------------------                    
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

                        # 2DAlphabetFit histograms --------------------------------
                        if 'particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields:
                            output['hLeadingFatJetParticleNet_massH_Hto4b_avg_vs_massA_Hto4b_avg'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet_PNet_massH_Hto4b_avg[sel_SR_forHExt]),
                                Mass2=(leadingFatJet_PNet_massA_Hto4b_avg[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )                    
                            output['hLeadingFatJetMass_vs_massA_Hto4b_avg'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.mass_toUse[sel_SR_forHExt]),
                                Mass2=(leadingFatJet_PNet_massA_Hto4b_avg[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            output['hLeadingFatJetMSoftDrop_vs_massA_Hto4b_avg'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.msoftdrop_toUse[sel_SR_forHExt]),
                                Mass2=(leadingFatJet_PNet_massA_Hto4b_avg[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )    
                        ## PNetMD Hto4b NanoAOD_v2
                        if 'PNet_X4b_v2a_Haa4b_score' in FatJetsToUse.fields:
                            # mH vs mAa for different mH versions
                            output['hLeadingFatJetPNet_massH_v2b_vs_massAa'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.PNet_massH_v2b[sel_SR_forHExt]),
                                Mass2=(leadingFatJet.PNet_massAa[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )                    
                            output['hLeadingFatJetMass_vs_massAa'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.mass_toUse[sel_SR_forHExt]),
                                Mass2=(leadingFatJet.PNet_massAa[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            output['hLeadingFatJetMSoftDrop_vs_massAa'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.msoftdrop_toUse[sel_SR_forHExt]),
                                Mass2=(leadingFatJet.PNet_massAa[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )

                            # mH vs mA34a for different mH versions
                            output['hLeadingFatJetPNet_massH_v2b_vs_massA34a'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.PNet_massH_v2b[sel_SR_forHExt]),
                                Mass2=(leadingFatJet.PNet_34massAa[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )                    
                            output['hLeadingFatJetMass_vs_massA34a'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.mass_toUse[sel_SR_forHExt]),
                                Mass2=(leadingFatJet.PNet_34massAa[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            output['hLeadingFatJetMSoftDrop_vs_massA34a'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.msoftdrop_toUse[sel_SR_forHExt]),
                                Mass2=(leadingFatJet.PNet_34massAa[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )

                            # mass vs mA for different mA versions                            
                            output['hLeadingFatJetMass_vs_massA34b'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.mass_toUse[sel_SR_forHExt]),
                                Mass2=(leadingFatJet.PNet_34massAb[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetMass_vs_massA34d'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.mass_toUse[sel_SR_forHExt]),
                                Mass2=(leadingFatJet.PNet_34massAd[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetMass_vs_massA34ad'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.mass_toUse[sel_SR_forHExt]),
                                Mass2=(((leadingFatJet.PNet_34massAa+leadingFatJet.PNet_34massAd)*0.5)[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )

                            # mH vs mA34d for different mH versions                                                
                            output['hLeadingFatJetMass_vs_massA34d'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.mass_toUse[sel_SR_forHExt]),
                                Mass2=(leadingFatJet.PNet_34massAd[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            output['hLeadingFatJetMSoftDrop_vs_massA34d'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.msoftdrop_toUse[sel_SR_forHExt]),
                                Mass2=(leadingFatJet.PNet_34massAd[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetPNet_massH_v2b_vs_massA34d'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.PNet_massH_v2b[sel_SR_forHExt]),
                                Mass2=(leadingFatJet.PNet_34massAd[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            

                    # Event weights histograms ------------------------------------------------------
                    if histogramSaveLevel >= 1 and self.datasetInfo['isMC'] and syst == "Nom":
                        output['hEventWeight_PU'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_PU[sel_SR_forHExt]
                        )
                        output['hEventWeight_PUUp'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_PUUp[sel_SR_forHExt]
                        )
                        output['hEventWeight_PUDown'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_PUDown[sel_SR_forHExt]
                        )
                        if self.datasetInfo['isSignal']:
                            output['hEventWeight_LundPlane_Nom'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_LundPlane_Nom[sel_SR_forHExt]
                            )
                            output['hEventWeight_LundPlane_Up'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_LundPlane_Up[sel_SR_forHExt]
                            )
                            output['hEventWeight_LundPlane_Down'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_LundPlane_Down[sel_SR_forHExt]
                            )                            
                        if self.datasetInfo['isSignalGGH']:
                            output['hEventWeight_GGHHiggsPt'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_GGH_HiggsPt[sel_SR_forHExt]
                            )
                            output['hEventWeight_GGHHiggsPtUp'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_GGH_HiggsPtUp[sel_SR_forHExt]
                            )
                            output['hEventWeight_GGHHiggsPtDown'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_GGH_HiggsPtDown[sel_SR_forHExt]
                            )
                            
                        if self.datasetInfo['isQCD_bGen']:
                            output['hEventWeight_QCDHT'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_HT[sel_SR_forHExt]
                            )
                        if self.datasetInfo['isTTbar']:
                            output['hEventWeight_TopPt'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_TopPt[sel_SR_forHExt]
                            )
                            output['hEventWeight_TopPtUp'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_TopPtUp[sel_SR_forHExt]
                            )
                            output['hEventWeight_TopPtDown'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_TopPtDown[sel_SR_forHExt]
                            )                            
                        output['hEventWeight_PS'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_PS_Nom[sel_SR_forHExt]
                        )
                        output['hEventWeight_PS_ISRUp'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_PS_ISRUp[sel_SR_forHExt]
                        )
                        output['hEventWeight_PS_ISRDown'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_PS_ISRDown[sel_SR_forHExt]
                        )
                        output['hEventWeight_PS_FSRUp'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_PS_FSRUp[sel_SR_forHExt]
                        )
                        output['hEventWeight_PS_FSRDown'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_PS_FSRDown[sel_SR_forHExt]
                        )
                        output['hEventWeight_QCDPdfNom'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_QCDPdfNom[sel_SR_forHExt]
                        )
                        output['hEventWeight_QCDPdfUp'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_QCDPdfUp[sel_SR_forHExt]
                        )
                        output['hEventWeight_QCDPdfDown'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_QCDPdfDown[sel_SR_forHExt]
                        )
                        output['hEventWeight_QCDScale_Nom'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_QCDScale_Nom[sel_SR_forHExt]
                        )
                        output['hEventWeight_QCDScale_RenormUp'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_QCDScale_RenormUp[sel_SR_forHExt]
                        )
                        output['hEventWeight_QCDScale_RenormDown'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_QCDScale_RenormDown[sel_SR_forHExt]
                        )
                        output['hEventWeight_QCDScale_FactorizationUp'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_QCDScale_FactorizationUp[sel_SR_forHExt]
                        )
                        output['hEventWeight_QCDScale_FactorizationDown'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_QCDScale_FactorizationDown[sel_SR_forHExt]
                        )
                        output['hEventWeight_QCDPDFNom'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_QCDPdfNom[sel_SR_forHExt]
                        )
                        output['hEventWeight_QCDPDFUp'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_QCDPdfUp[sel_SR_forHExt]
                        )
                        output['hEventWeight_QCDPDFDown'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_QCDPdfDown[sel_SR_forHExt]
                        )
                        output['hEventWeight_Ak4BtagNom'+sHExt].fill(
                            dataset=dataset,
                            Weight=wgt_Ak4Btag_dict['Nom'][sel_SR_forHExt]
                        )
                        if  kDatasetToAnalyze == DatasetToAnalyze.SingleYear:
                            output['hEventWeight_Ak4BtagUp'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_Ak4Btag_dict['Up'][sel_SR_forHExt]
                            )
                            output['hEventWeight_Ak4BtagDown'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_Ak4Btag_dict['Down'][sel_SR_forHExt]
                            )
                        elif kDatasetToAnalyze == DatasetToAnalyze.FullRun2:    
                            output['hEventWeight_Ak4BtagUpuncorrelated'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_Ak4Btag_dict['Upuncorrelated'][sel_SR_forHExt]
                            )
                            output['hEventWeight_Ak4BtagDownuncorrelated'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_Ak4Btag_dict['Downuncorrelated'][sel_SR_forHExt]
                            )
                            output['hEventWeight_Ak4BtagUpcorrelated'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_Ak4Btag_dict['Upcorrelated'][sel_SR_forHExt]
                            )
                            output['hEventWeight_Ak4BtagDowncorrelated'+sHExt].fill(
                                dataset=dataset,
                                Weight=wgt_Ak4Btag_dict['Downcorrelated'][sel_SR_forHExt]
                            )                        
                    # ------------------------------------------------------    
                        
                        
                        

                    if histogramSaveLevel >= 1: 
                        output['hPV_npvsGood'+sHExt].fill(
                            dataset=dataset,
                            PU=(events.PV.npvsGood[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )                    
                        output['hLeadingFatJetPt'+sHExt].fill(
                            dataset=dataset,
                            #Pt=ak.flatten(selFatJet.pt_toUse[sel_SR_forHExt][:, 0]),
                            #Pt=(selFatJet.pt_toUse[sel_SR_forHExt][:, 0]),
                            Pt=(leadingFatJet.pt_toUse[sel_SR_forHExt]),
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
                        output['hLeadingFatJetMass'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.mass_toUse[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hLeadingFatJetMSoftDrop'+sHExt].fill(
                            dataset=dataset,
                            Mass=(leadingFatJet.msoftdrop_toUse[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        output['hnleadingNonHto4bFatJet_WZvsQCD'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nleadingNonHto4bFatJet_WZvsQCD[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hMET_pT'+sHExt].fill(
                            dataset=dataset,
                            Pt=(METToUse.pt_toUse[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hPuppiMET_pT'+sHExt].fill(
                            dataset=dataset,
                            Pt=(events.PuppiMET.pt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        ## PNetMD Hto4b
                        if 'particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields:
                            # Scale particleNet_massH_Hto4b_v0 by 1.01 to get better response
                            # https://indico.cern.ch/event/1343334/contributions/5655252/attachments/2745224/4781382/2023_11_02_HToAATo4B_Higgs_mass_studies.pdf#page=15
                            output['hLeadingFatJetParticleNet_massH_Hto4b_avg_v0123'+sHExt].fill(
                                dataset=dataset,
                                Mass=leadingFatJet_PNet_massH_Hto4b_avg[sel_SR_forHExt],
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
                        ## PNetMD Hto4b NanoAOD_v2
                        if 'PNet_X4b_v2a_Haa4b_score' in FatJetsToUse.fields:
                            output['hLeadingFatJetPNet_X4b_v1_Haa4b_vs_QCD'+sHExt].fill(
                                dataset=dataset,
                                MLScore1k=leadingFatJet.PNet_X4b_v1_Haa4b_vs_QCD[sel_SR_forHExt],
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetPNet_X4b_v1_Haa4b_score'+sHExt].fill(
                                dataset=dataset,
                                MLScore1k=leadingFatJet.PNet_X4b_v1_Haa4b_score[sel_SR_forHExt],
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetPNet_X4b_v2a_Haa4b_score'+sHExt].fill(
                                dataset=dataset,
                                MLScore1k=leadingFatJet.PNet_X4b_v2a_Haa4b_score[sel_SR_forHExt],
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetPNet_X4b_v2b_Haa4b_score'+sHExt].fill(
                                dataset=dataset,
                                MLScore1k=leadingFatJet.PNet_X4b_v2b_Haa4b_score[sel_SR_forHExt],
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            lFJ_PNet_X4b_v2ab_Haa4ab_score = (leadingFatJet.PNet_X4b_v2a_Haa4b_score + leadingFatJet.PNet_X4b_v2b_Haa4b_score) / 2.0
                            output['hLeadingFatJetPNet_X4b_v2ab_Haa4b_score'+sHExt].fill(
                                dataset=dataset,
                                MLScore1k=lFJ_PNet_X4b_v2ab_Haa4ab_score[sel_SR_forHExt],
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetPNet_X4b_v2a_Haa34b_score'+sHExt].fill(
                                dataset=dataset,
                                MLScore1k=leadingFatJet.PNet_X4b_v2a_Haa34b_score[sel_SR_forHExt],
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetPNet_X4b_v2b_Haa34b_score'+sHExt].fill(
                                dataset=dataset,
                                MLScore1k=leadingFatJet.PNet_X4b_v2b_Haa34b_score[sel_SR_forHExt],
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            lFJ_PNet_X4b_v2ab_Haa34ab_score = (leadingFatJet.PNet_X4b_v2a_Haa34b_score + leadingFatJet.PNet_X4b_v2b_Haa34b_score) / 2.0
                            output['hLeadingFatJetPNet_X4b_v2ab_Haa34b_score'+sHExt].fill(
                                dataset=dataset,
                                MLScore1k=lFJ_PNet_X4b_v2ab_Haa34ab_score[sel_SR_forHExt],
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )

                            # massH                            
                            ''' output['hLeadingFatJetMassH_v1'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.PNet_massH_v1[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]                                
                            )
                            output['hLeadingFatJetMassH_v2a'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.PNet_massH_v2a[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]                                
                            ) '''
                            output['hLeadingFatJetMassH_v2b'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.PNet_massH_v2b[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]                                
                            )
                            ''' output['hLeadingFatJetMassH_v2c'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.PNet_massH_v2c[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]                                
                            )
                            output['hLeadingFatJetMassH_v2d'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.PNet_massH_v2d[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]                                
                            )
                            output['hLeadingFatJetMassH_avg'+sHExt].fill(
                                dataset=dataset,
                                Mass=(leadingFatJet.PNet_massH_avg[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]                                
                            )
                            output['hLeadingFatJetMassH_std'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_massH_std[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]                                
                            ) '''                        
                            
                            # massA
                            output['hLeadingFatJetPNet_massAa'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_massAa[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            ''' output['hLeadingFatJetPNet_massAb'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_massAb[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            output['hLeadingFatJetPNet_massAc'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_massAc[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            output['hLeadingFatJetPNet_massAd'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_massAd[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            output['hLeadingFatJetPNet_massA_avg'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_massA_avg[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            output['hLeadingFatJetPNet_massA_std'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_massA_std[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) '''
                            output['hLeadingFatJetPNet_34massAa'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_34massAa[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            output['hLeadingFatJetPNet_34massAb'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_34massAb[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            ''' output['hLeadingFatJetPNet_34massAc'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_34massAc[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) '''
                            output['hLeadingFatJetPNet_34massAd'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_34massAd[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            ) 
                            output['hLeadingFatJetPNet_massA1'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_massA1[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetPNet_massA2'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_massA2[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetPNet_massAA'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_massAA[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetPNet_dMassAA'+sHExt].fill(
                                dataset=dataset,
                                Mass1=(leadingFatJet.PNet_dMassAA[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                            output['hLeadingFatJetPNet_dMassAA_relH'+sHExt].fill(
                                dataset=dataset,
                                Mass5=(leadingFatJet.PNet_dMassAA_relH[sel_SR_forHExt]),
                                systematic=syst,
                                weight=evtWeight[sel_SR_forHExt]
                            )
                             
                            
                                                       
                            
                            
                            


                    if histogramSaveLevel >= 2:
                        output['hPV_npvs'+sHExt].fill(
                            dataset=dataset,
                            PU=(events.PV.npvs[sel_SR_forHExt]),
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
                        output['hLeadingBtagFatJetPtOverLeadingFatJetPt'+sHExt].fill(
                            dataset=dataset,
                            Ratio=( leadingFatJet[sel_SR_forHExt].pt_toUse / ak.firsts(FatJetsToUse)[sel_SR_forHExt].pt_toUse),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        ## PNetMD Hto4b NanoAOD_v1
                        if 'particleNetMD_Hto4b_Haa4b' in FatJetsToUse.fields:
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

                            output['hLeadingFatJetPNet_Xto4bv1_Htoaa4bOverQCD'+sHExt].fill(
                                dataset=dataset,
                                #MLScore1k=(leadingFatJet.particleNetMD_Hto4b_Haa4b[sel_SR_forHExt] / (leadingFatJet.particleNetMD_Hto4b_Haa4b[sel_SR_forHExt] + leadingFatJet_PNetMD_Hto4b_QCD01234b_sum[sel_SR_forHExt])),
                                MLScore1k=leadingFatJet_PNet_Xto4bv1_Htoaa4bOverQCD[sel_SR_forHExt],
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


                        ## NonHto4bFatJets
                        output['hnNonHto4bFatJet'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nNonHto4bFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )

                        

                        ## MET
                        output['hMET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(METToUse.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )
                        output['hdPhi_MET_leadingFatJet'+sHExt].fill(
                            dataset=dataset,
                            deltaPhi=(abs(METToUse.delta_phi(leadingFatJet))[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )      

                        output['hPuppiMET_sumEt'+sHExt].fill(
                            dataset=dataset,
                            Pt4TeV=(events.PuppiMET.sumEt[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )                    
                        output['hdPhi_PuppiMET_leadingFatJet'+sHExt].fill(
                            dataset=dataset,
                            deltaPhi=(abs(events.PuppiMET.delta_phi(leadingFatJet))[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )  
                        output['hPuppiMET_sumEt_minus_FJHto4bPt'+sHExt].fill(
                            dataset=dataset,
                            Pt2TeV=(events.PuppiMET.sumEt - leadingFatJet.pt_toUse)[sel_SR_forHExt],
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )      
                        output['hMETPhi'+sHExt].fill(
                            dataset=dataset,
                            Phi=(METToUse.phi_toUse[ sel_SR_forHExt ]),
                            systematic=syst,
                            weight=evtWeight[ sel_SR_forHExt ]
                        ) 
                        output['hPuppiMETPhi'+sHExt].fill(
                            dataset=dataset,
                            Phi=(events.PuppiMET.phi[ sel_SR_forHExt ]),
                            systematic=syst,
                            weight=evtWeight[ sel_SR_forHExt ]
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


                        ## AK4 jets
                        output['hnAK4Jets_NonoverlapLeadingFatJet'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nAk4Jets_nonoverlaping_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )    
                        if ak.count(ak.firsts(ak4Jets_nonoverlaping_leadingFatJet).pt_toUse[sel_SR_forHExt]) > 0:
                            output['hPtLeadingAK4Jets_NonoverlapLeadingFatJet'+sHExt].fill(
                                dataset=dataset,
                                Pt=(ak.firsts(ak4Jets_nonoverlaping_leadingFatJet).pt_toUse[
                                    sel_SR_forHExt & 
                                    (~ak.is_none(ak.firsts(ak4Jets_nonoverlaping_leadingFatJet).pt_toUse))
                                    ]),
                                systematic=syst,
                                weight=evtWeight[
                                    sel_SR_forHExt & 
                                    (~ak.is_none(ak.firsts(ak4Jets_nonoverlaping_leadingFatJet).pt_toUse))
                                    ]
                            )
                        output['hnAK4Jets_bTag_NonoverlapLeadingFatJet'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nAk4Jets_bTag_nonoverlaping_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )    
                        if ak.count(ak.firsts(ak4Jets_bTag_nonoverlaping_leadingFatJet).pt_toUse[sel_SR_forHExt]) > 0:
                            output['hPtLeadingAK4Jets_bTag_NonoverlapLeadingFatJet'+sHExt].fill(
                                dataset=dataset,
                                Pt=(ak.firsts(ak4Jets_bTag_nonoverlaping_leadingFatJet).pt_toUse[
                                    sel_SR_forHExt & 
                                    (~ak.is_none(ak.firsts(ak4Jets_bTag_nonoverlaping_leadingFatJet).pt_toUse))
                                    ]),
                                systematic=syst,
                                weight=evtWeight[
                                    sel_SR_forHExt & 
                                    (~ak.is_none(ak.firsts(ak4Jets_bTag_nonoverlaping_leadingFatJet).pt_toUse))
                                    ]
                            )

                        output['hnAK4JetsCentral_NonoverlapLeadingFatJet'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nAk4JetsCentral_nonoverlaping_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )    
                        if ak.count(ak.firsts(ak4JetsCentral_nonoverlaping_leadingFatJet).pt_toUse[sel_SR_forHExt]) > 0:
                            output['hPtLeadingAK4JetsCentral_NonoverlapLeadingFatJet'+sHExt].fill(
                                dataset=dataset,
                                Pt=(ak.firsts(ak4JetsCentral_nonoverlaping_leadingFatJet).pt_toUse[
                                    sel_SR_forHExt & 
                                    (~ak.is_none(ak.firsts(ak4JetsCentral_nonoverlaping_leadingFatJet).pt_toUse))
                                    ]),
                                systematic=syst,
                                weight=evtWeight[
                                    sel_SR_forHExt & 
                                    (~ak.is_none(ak.firsts(ak4JetsCentral_nonoverlaping_leadingFatJet).pt_toUse))
                                    ]
                            )
                        output['hnAK4JetsCentral_bTag_NonoverlapLeadingFatJet'+sHExt].fill(
                            dataset=dataset,
                            nObject10=(nAk4JetsCentral_bTag_nonoverlaping_leadingFatJet[sel_SR_forHExt]),
                            systematic=syst,
                            weight=evtWeight[sel_SR_forHExt]
                        )    
                        if ak.count(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_leadingFatJet).pt_toUse[sel_SR_forHExt]) > 0:
                            output['hPtLeadingAK4JetsCentral_bTag_NonoverlapLeadingFatJet'+sHExt].fill(
                                dataset=dataset,
                                Pt=(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_leadingFatJet).pt_toUse[
                                    sel_SR_forHExt & 
                                    (~ak.is_none(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_leadingFatJet).pt_toUse))
                                    ]),
                                systematic=syst,
                                weight=evtWeight[
                                    sel_SR_forHExt & 
                                    (~ak.is_none(ak.firsts(ak4JetsCentral_bTag_nonoverlaping_leadingFatJet).pt_toUse))
                                    ]
                            )





                    if 'hLeadingFatJetEta_vs_Phi'+sHExt in output:
                        output['hLeadingFatJetEta_vs_Phi'+sHExt].fill(
                            dataset=dataset,
                            Eta=(leadingFatJet.eta[sel_SR_forHExt]),
                            Phi=(leadingFatJet.phi[sel_SR_forHExt]),
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
                            printVariable('\n leadingFatJet.pt_toUse[ sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi ]', leadingFatJet.pt_toUse[ sel_SR_forHExt & mask_HEM1516Issue_Eta & mask_HEM1516Issue_Phi ]); sys.stdout.flush()

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
                            Pt=(leadingFatJet.pt_toUse[ sel_tmp_ ]),
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
                            Pt=(leadingFatJet.pt_toUse[ sel_tmp_ ]),
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
                            Pt=(leadingFatJet.pt_toUse[ sel_tmp_ ]),
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
                            Pt=(leadingFatJet.pt_toUse[ sel_tmp_ ]),
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
                            Pt=(leadingFatJet.pt_toUse[ sel_tmp_ ]),
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
                            Pt=(leadingFatJet.pt_toUse[ sel_tmp_ ]),
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
                            Pt=(leadingFatJet.pt_toUse[ sel_tmp_ ]),
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
                            Pt=(leadingFatJet.pt_toUse[ sel_tmp_ ]),
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
                            Pt=(leadingFatJet.pt_toUse[ sel_tmp_ ]),
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
                            Pt=(leadingFatJet.pt_toUse[ sel_tmp_ ]),
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
                            Pt=(leadingFatJet.pt_toUse[ sel_tmp_ ]),
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
                    






                    ### 2-D distribution ----------------------------------------------------------
                    if histogramSaveLevel >= 2: 
                        output['hMET_pT_vs_dPhi_MET_leadingFatJet'+sHExt].fill(
                            dataset=dataset,
                            Pt=(METToUse.pt_toUse[sel_SR_forHExt]),
                            deltaPhi=(abs(METToUse.delta_phi(leadingFatJet))[sel_SR_forHExt]),
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
    print(f"htoaa_Analysis_GGFMode:: here14 {datetime.now() = }")

    if len(sys.argv) != 2:
        print("htoaa_Analysis:: Command-line config file missing.. \t **** ERROR **** \n")

    sConfig = sys.argv[1]

    
    config = GetDictFromJsonFile(sConfig)
    print("Config {}: \n{}".format(sConfig, json.dumps(config, indent=4)))
    print(f"htoaa_Analysis_GGFMode:: here15 {datetime.now() = }")

    nEventsToAnalyze    = config["nEventsToAnalyze"] if "nEventsToAnalyze" in config else nEventsToAnalyze
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
        systematicsToRun    = config["systematics"].lower() if "systematics" in config else 'no'

        
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
                logging.critical(f'htoaa_Analysis_GGFMode.py::main():: {MCSamplesStitchInputFileName = } does not exists')
                print(f'htoaa_Analysis_GGFMode.py::main() 11:: {MCSamplesStitchInputFileName = } does not exists')
                exit(0)
            print(f"Opening {MCSamplesStitchInputFileName = } "); sys.stdout.flush() 
            with uproot.open(MCSamplesStitchInputFileName) as f_:
                print(f"{f_.keys() = }"); sys.stdout.flush() 
                hMCSamplesStitch = f_[r'%s' % MCSamplesStitchInputHistogramName].to_hist()

    print(f"htoaa_Analysis_GGFMode:: here16 {datetime.now() = }")    
        
        
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
    print(f"htoaa_Analysis_GGFMode:: here17 {datetime.now() = }")

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
            logging.critical('htoaa_Analysis_GGFMode:: getNanoAODFile() for input file %s failed. **** CRITICAL ERROR ****. \nAborting...' % (sInputFile)); sys.stdout.flush();
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
        print(f"htoaa_Analysis_GGFMode:: {sInputFile} \t {os.path.exists(sInputFile) = }, {fileSize = } MB");     

        if fileSize > NanoAODFileSize_Min:     
            sInputFiles[iFile] = sInputFile
        else:
            logging.critical('htoaa_Analysis_GGFMode:: Input file %s file size below threshold (%g MB). **** CRITICAL ERROR ****. \nAborting...' % (sInputFile, NanoAODFileSize_Min) ); sys.stdout.flush();
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
    print(f"htoaa_Analysis_GGFMode:: here18 {datetime.now() = }")


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
        sampleInfo["systematicsToRun"] = systematicsToRun
    print(f"htoaa_Analysis_GGFMode:: here19 {datetime.now() = }", flush=flushStdout)
        
    startTime = time.time()
    tracemalloc.start()


    #client = Client("tls://localhost:8786")
    #executor = processor.DaskExecutor(client=client)
    chunksize = nEventToReadInBatch  if nEventsToAnalyze == -1 else nEventsToAnalyze
    maxchunks = None if nEventsToAnalyze == -1 else int(max(nEventsToAnalyze/nEventToReadInBatch, 1))
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

            #print("%10f\t%10d\t%s" % (output['cutflow'][sWeighted+key], output['cutflow'][key], key), flush=flushStdout)
            print("%10d\t%s" % (output['cutflow'][key], key), flush=flushStdout)
    
    
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

                #print(f"value ({type(value)}): {value}")

                for _dataset in value.axis('dataset').identifiers():
                    #print(f"_dataset ({type(_dataset)}): {_dataset}")

                    #print(f"value.fields : {value.fields}", fl)
                    #print(f"value.fields() ({type(value.fields())}): {value.fields()}")
                    if 'systematic' not in value.fields: # hWeight histogram do not have 'systematics' axes
                        h1 = value.integrate('dataset',_dataset).to_hist()
                        fOut['%s/%s' % (sDir1_toUse, sHistoName_toUse)] = h1
                        continue

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
    
