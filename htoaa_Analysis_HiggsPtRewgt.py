#htoaa analysis main code

import os
import sys
from datetime import datetime
#import time
print(f"htoaa_Analysis_HiggsPtRewgt:: here1 {datetime.now() = }"); sys.stdout.flush()
import subprocess
import json
from urllib.request import urlopen
import glob
from collections import OrderedDict as OD
import time
import tracemalloc
import math
print(f"htoaa_Analysis_HiggsPtRewgt:: here2 {datetime.now() = }"); sys.stdout.flush()
import numpy as np
from copy import copy, deepcopy
print(f"htoaa_Analysis_HiggsPtRewgt:: here3 {datetime.now() = }"); sys.stdout.flush()
#import uproot
#import uproot3 as uproot
import uproot as uproot
print(f"htoaa_Analysis_HiggsPtRewgt:: here4 {datetime.now() = }"); sys.stdout.flush()
#import parse
from parse import *
print(f"htoaa_Analysis_HiggsPtRewgt:: here4.1 {datetime.now() = }"); sys.stdout.flush()
import logging
print(f"htoaa_Analysis_HiggsPtRewgt:: here5 {datetime.now() = }"); sys.stdout.flush()

# comment test3
'''
GGF -> H->aa->4b boosted analysis macro

References:
  * Coffea framework used for TTGamma analysis: https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py
* Coffea installation: /home/siddhesh/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea
'''
print(f"htoaa_Analysis_HiggsPtRewgt:: here6 {datetime.now() = }"); sys.stdout.flush()
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
print(f"htoaa_Analysis_HiggsPtRewgt:: here7 {datetime.now() = }"); sys.stdout.flush()
from particle import Particle # For PDG particle listing https://github.com/scikit-hep/particle
print(f"htoaa_Analysis_HiggsPtRewgt:: here8 {datetime.now() = }"); sys.stdout.flush()


from htoaa_Settings import *
print(f"htoaa_Analysis_HiggsPtRewgt:: here9 {datetime.now() = }"); sys.stdout.flush()
from htoaa_CommonTools import (
    getLorentVector,
    GetDictFromJsonFile, akArray_isin,
    selectRunLuminosityBlock,
    calculate_lumiScale, getLumiScaleForPhSpOverlapRewgtMode, getSampleHTRange, # update_crosssection, 
    getNanoAODFile, setXRootDRedirector,  xrdcpFile,
    selectGenHiggs, selectGenZBoson, selectGenWBoson, selectGenWplusBoson, selectGenWminusBoson, 
    selectGenTop, selectGenAntiTop, selectGenQuarksFromHardScattering,
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
print(f"htoaa_Analysis_HiggsPtRewgt:: here10 {datetime.now() = }"); sys.stdout.flush()
from htoaa_Samples import (
    kData, kQCD_bEnrich, kQCD_bGen, kQCDIncl
)
print(f"htoaa_Analysis_HiggsPtRewgt:: here11 {datetime.now() = }"); sys.stdout.flush()

from inspect import currentframe, getframeinfo
print(f"htoaa_Analysis_HiggsPtRewgt:: here12 {datetime.now() = }"); sys.stdout.flush()
frameinfo = getframeinfo(currentframe())
print(f"htoaa_Analysis_HiggsPtRewgt:: here13 {datetime.now() = }"); sys.stdout.flush()


# use GOldenJSON

 
printLevel = 0
histogramSaveLevel = 2 # 0: hSignal extraction, 1: basic Data-MC validation, 2:..
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
        self.FatJetPNetXto4bv2WorkingPoints = ['40', '60', '80'] #['40', '60', '80'] #['40', '45a', '45b', '50', '60', '65', '70', '80']  # ['40', '50', '60', '65', '70', '80']   ['40', '60', '80']
        
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
         
        global runMode_SignalGenChecks;       runMode_SignalGenChecks  = True; # True
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
        self.datasetInfo['isHToX'         ]  = False
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
            self.datasetInfo['isHToX']           = True if "HTo"     in datasetName_part1 else False
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
            ]),
        ])

        
        



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
                logging.critical(f'htoaa_Analysis_HiggsPtRewgt.py::main():: {sFilesGoldenJSON[self.datasetInfo["era"]] = } could not read.')
                exit(0) 

            # convert runNumber in str to int
            dataLSSelGoldenJSON = {int(k): v for k, v in dataLSSelGoldenJSON.items()} 
            self.datasetInfo['dataLSSelGoldenJSON'] = dataLSSelGoldenJSON
            #print(f"{dataLSSelGoldenJSON = }")

        else: ## MC

            # lumiScale --------------------------------------------------------------------------------------------------
            if sTrgSelection not in Luminosities_forGGFMode[self.datasetInfo["era"]]:
                logging.critical(f'htoaa_Analysis_HiggsPtRewgt.py::main():: {sTrgSelection = } not in {Luminosities_forGGFMode[self.datasetInfo["era"]] = }.')
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
        pt2TeV_axis1          = hist.Bin("Pt2TeV1",                r"$p_{T}$ [GeV]",            2000,       0,    2000)
        log2Pt2TeV_axis       = hist.Bin("Log2Pt2TeV",             r"Log2($p_{T}$) [GeV]",      200,  math.log2(1),    math.log2(2000))
        log2Pt2TeV_axis1      = hist.Bin("Log2Pt2TeV1",            r"Log2($p_{T}$) [GeV]",      200,  math.log2(1),    math.log2(2000))
        log2Ratio_axis        = hist.Bin("Log2Ratio",              r"Log2(Ratio) [GeV]",        120,     -1.5,   1.5)        
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




        if (self.datasetInfo['isSignal'] or self.datasetInfo['isHToX']) and runMode_SignalGenChecks:
            histos.update(OD([
                ('hGenHiggsPt',                           {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]"}),
                ('hGenHiggsLog2Pt',                       {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]"}),
                ('hGenHiggsMass',                         {sXaxis: mass_axis,      sXaxisLabel: r"$Mass(GEN Higgs)$ [GeV]"}),
                ('hGenZPt',                               {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Z)$ [GeV]"}),
                ('hGenZLog2Pt',                           {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Z))$ [GeV]"}),
                ('hGenZMass',                             {sXaxis: mass_axis,      sXaxisLabel: r"$Mass(GEN Z)$ [GeV]"}),
                ('hGenWPt',                               {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN W)$ [GeV]"}),
                ('hGenWLog2Pt',                           {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN W))$ [GeV]"}),
                ('hGenWMass',                             {sXaxis: mass_axis,      sXaxisLabel: r"$Mass(GEN W)$ [GeV]"}),
                ('hGenWplusPt',                           {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Wplus)$ [GeV]"}),
                ('hGenWplusLog2Pt',                       {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Wplus))$ [GeV]"}),
                ('hGenWplusMass',                         {sXaxis: mass_axis,      sXaxisLabel: r"$Mass(GEN Wplus)$ [GeV]"}),
                ('hGenWminusPt',                          {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Wminus)$ [GeV]"}),
                ('hGenWminusLog2Pt',                      {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Wminus))$ [GeV]"}),
                ('hGenWminusMass',                        {sXaxis: mass_axis,      sXaxisLabel: r"$Mass(GEN Wminus)$ [GeV]"}),
                ('hGenTtbarPt',                           {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Ttbar)$ [GeV]"}),
                ('hGenTtbarLog2Pt',                       {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Ttbar))$ [GeV]"}),
                ('hGenQQFromHardScatteringPt',            {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN QQFromHardScattering)$ [GeV]"}),
                ('hGenQQFromHardScatteringLog2Pt',        {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN QQFromHardScattering))$ [GeV]"}),  
                
                
                ('hGenPt_Higgs_vs_Z',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN Z)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_Z',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN Z))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByZ',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(Z)))$"}),  
                
                ('hGenPt_Higgs_vs_Z',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN Z)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_Z',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN Z))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByZ',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(Z)))$"}),  
                
                ('hGenPt_Higgs_vs_W',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN W)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_W',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN W))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByW',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(W)))$"}),  
                
                ('hGenPt_Higgs_vs_Wplus',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN Wplus)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_Wplus',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN Wplus))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByWplus',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(Wplus)))$"}),  
                
                ('hGenPt_Higgs_vs_Wminus',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN Wminus)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_Wminus',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN Wminus))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByWminus',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(Wminus)))$"}),  
                
                ('hGenPt_Higgs_vs_Ttbar',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN Ttbar)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_Ttbar',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN Ttbar))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByTtbar',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(Ttbar)))$"}),  

                ('hGenPt_Higgs_vs_QQFromHardScattering',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN QQFromHardScattering)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_QQFromHardScattering',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN QQFromHardScattering))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByQQFromHardScattering',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(QQFromHardScattering)))$"}),  

                 
                
                             

                ('hnGenHiggs',                            {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of GenHiggs"}),
                ('hnGenZ',                                {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of  GenZ"}),
                ('hnGenW',                                {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of  GenW"}),
                ('hnGenWplus',                            {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of GenWplus"}),
                ('hnGenWminus',                           {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of  GenWminus"}),
                ('hnGenTop',                              {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of  GenTop"}),
                ('hnGenAntiTop',                          {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of  GenAntiTop"}),
                ('hnGenQuarksFromHardScattring',                          {sXaxis: nObject10_axis,  sXaxisLabel: r"No. of GenQuarksFromHardScattring"}),                
                        
                ('hGenHiggsPt_wHiggsPtRewgt',             {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]"}),
                ('hGenHiggsLog2Pt_wHiggsPtRewgt',         {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]"}),
                ('hGenZPt_wHiggsPtRewgt',                               {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Z)$ [GeV]"}),
                ('hGenZLog2Pt_wHiggsPtRewgt',                           {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Z))$ [GeV]"}),
                ('hGenWPt_wHiggsPtRewgt',                               {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN W)$ [GeV]"}),
                ('hGenWLog2Pt_wHiggsPtRewgt',                           {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN W))$ [GeV]"}),
                ('hGenWplusPt_wHiggsPtRewgt',                           {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Wplus)$ [GeV]"}),
                ('hGenWplusLog2Pt_wHiggsPtRewgt',                       {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Wplus))$ [GeV]"}),
                ('hGenWminusPt_wHiggsPtRewgt',                          {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Wminus)$ [GeV]"}),
                ('hGenWminusLog2Pt_wHiggsPtRewgt',                      {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Wminus))$ [GeV]"}),
                ('hGenTtbarPt_wHiggsPtRewgt',                           {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Ttbar)$ [GeV]"}),
                ('hGenTtbarLog2Pt_wHiggsPtRewgt',                       {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Ttbar))$ [GeV]"}),
                ('hGenQQFromHardScatteringPt_wHiggsPtRewgt',            {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN QQFromHardScattering)$ [GeV]"}),
                ('hGenQQFromHardScatteringLog2Pt_wHiggsPtRewgt',        {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN QQFromHardScattering))$ [GeV]"}),                

                ('hGenPt_Higgs_vs_Z_wHiggsPtRewgt',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN Z)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_Z_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN Z))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByZ_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(Z)))$"}),  
                
                ('hGenPt_Higgs_vs_Z_wHiggsPtRewgt',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN Z)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_Z_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN Z))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByZ_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(Z)))$"}),  
                
                ('hGenPt_Higgs_vs_W_wHiggsPtRewgt',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN W)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_W_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN W))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByW_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(W)))$"}),  
                
                ('hGenPt_Higgs_vs_Wplus_wHiggsPtRewgt',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN Wplus)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_Wplus_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN Wplus))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByWplus_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(Wplus)))$"}),  
                
                ('hGenPt_Higgs_vs_Wminus_wHiggsPtRewgt',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN Wminus)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_Wminus_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN Wminus))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByWminus_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(Wminus)))$"}),  
                
                ('hGenPt_Higgs_vs_Ttbar_wHiggsPtRewgt',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN Ttbar)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_Ttbar_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN Ttbar))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByTtbar_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(Ttbar)))$"}),  

                ('hGenPt_Higgs_vs_QQFromHardScattering_wHiggsPtRewgt',                           
                 {sXaxis: pt2TeV_axis,     sXaxisLabel: r"$p_{T}(GEN Higgs)$ [GeV]",
                  sYaxis: pt2TeV_axis1,    sYaxisLabel: r"$p_{T}(GEN QQFromHardScattering)$ [GeV]"}), 
                ('hGenLog2Pt_Higgs_vs_QQFromHardScattering_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Pt2TeV_axis1,sYaxisLabel: r"$log2 p_{T}(GEN QQFromHardScattering))$ [GeV]"}),  
                ('hGenLog2Pt_Higgs_vs_HiggsByQQFromHardScattering_wHiggsPtRewgt',                       
                 {sXaxis: log2Pt2TeV_axis, sXaxisLabel: r"$log2 p_{T}(GEN Higgs))$ [GeV]",
                  sYaxis: log2Ratio_axis,  sYaxisLabel: r"$log2 (2*p_{T}(H)/(p_{T}(H) + p_{T}(QQFromHardScattering)))$"}),  

            ]))


        

        # RECO-level histograms --------------------------------------------------------------------------------------------------------------
       

            
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

        

        ones_list   = np.ones(len(events))
        zeros_list  = np.zeros(len(events))
        trues_list  = np.ones(len(events), dtype=bool)
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
        genHiggses                   = selectGenHiggs(events)
        genZs                       = selectGenZBoson(events)
        genWs                       = selectGenWBoson(events)
        genWpluses                   = selectGenWplusBoson(events)
        genWminuses                  = selectGenWminusBoson(events)
        genTops                     = selectGenTop(events)
        genAntiTops                 = selectGenAntiTop(events)
        genTtbars                   = genTops + genAntiTops
        genQuarksFromHardScattring = selectGenQuarksFromHardScattering(events) # Select quarks coming out of hard scattering

        nGenHiggs                   = ak.fill_none(ak.count(genHiggses.pt, axis=1), 0)
        nGenZ                       = ak.fill_none(ak.count(genZs.pt, axis=1), 0)
        nGenW                       = ak.fill_none(ak.count(genWs.pt, axis=1), 0)
        nGenWplus                   = ak.fill_none(ak.count(genWpluses.pt, axis=1), 0)
        nGenWminus                  = ak.fill_none(ak.count(genWminuses.pt, axis=1), 0)
        nGenTop                     = ak.fill_none(ak.count(genTops.pt, axis=1), 0)
        nGenAntiTop                 = ak.fill_none(ak.count(genAntiTops.pt, axis=1), 0)
        nGenTtbar                   = ak.fill_none(ak.count(genTtbars.pt, axis=1), 0)
        nGenQuarksFromHardScattring = ak.fill_none(ak.count(genQuarksFromHardScattring.pt, axis=1), 0)
        
        # VBF: q1' q2' --> H q1 q2 : Leading two quarks coming out of hard scattering are VBF quarks
        genLeading2QuarksFromHardScattering = ak.mask(genQuarksFromHardScattring, nGenQuarksFromHardScattring >= 2) 
        genQQFromHardScattering = genLeading2QuarksFromHardScattering[:, 0] + genLeading2QuarksFromHardScattering[:, 1]

        # 
        genHiggs = ak.firsts(genHiggses)
        genZ = ak.firsts(genZs)
        genW = ak.firsts(genWs)
        genWplus = ak.firsts(genWpluses)
        genWminus = ak.firsts(genWminuses)
        genTtbar = ak.firsts(genTtbars)
        #genQQFromHardScattering = ak.firsts(genQQFromHardScattering)
        
        
            
            
        

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
        




        #####################
        # EVENT SELECTION
        #####################
        
        # reconstruction level cuts for cut-flow table. Order of cuts is IMPORTANT
        cuts_reco = ["dR_LeadingFatJet_GenB_0p8"] + self.sel_names_all["Presel"] #.copy()

       
        # create a PackedSelection object
        # this will help us later in composing the boolean selections easily
        selection = PackedSelection(dtype='uint32')







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

                   

            # MC PURewgt ----------------------------------
            #wgt_PU, wgt_PUUp, wgt_PUDown = getPURewgts_variation(
            #    events = events,
            #    year   = self.datasetInfo["era"]
            #)


            # MC GGF HToAATo4B Higgs pT reweight
            wgt_GGH_HiggsPt = wgt_GGH_HiggsPtUp = wgt_GGH_HiggsPtDown = None
            if self.datasetInfo['isSignalGGH']:
                wgt_GGH_HiggsPt, wgt_GGH_HiggsPtUp, wgt_GGH_HiggsPtDown = getHiggsPtRewgtForGGToHToAATo4B(
                    GenHiggsPt_list = genHiggs.pt #ak.firsts(genHiggs.pt)
                )


            weights.add(
                "lumiWeight",
                weight = lumiScale_toUse
            )
            weights.add(
                "genWeight",
                weight = np.copysign(np.ones(len(events)), events.genWeight)
            )
            
            #weights.add(
            #    "PU",
            #    weight     = wgt_PU,
            #    weightUp   = wgt_PUUp,
            #    weightDown = wgt_PUDown
            #)

            '''
            if self.datasetInfo['isSignalGGH']:
                weights.add(
                    "GGHPtRewgt",
                    weight     = wgt_GGH_HiggsPt,
                    weightUp   = wgt_GGH_HiggsPtUp,
                    weightDown = wgt_GGH_HiggsPtDown
                )
            '''



            ## weights_gen -------------------------------
            weights_gen.add(
                "lumiWeight",
                weight = lumiScale_toUse 
            )
            weights_gen.add(
                "genWeight",
                weight=np.copysign(np.ones(len(events)), events.genWeight)
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
                if len(self.sel_names_all[iSelection]) == 0: continue
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
                if syst == "Nom":
                    evtWeight_gen            = weights_gen.weight(weightSyst)
                

            #if printLevel >=0:
            #    printVariable('\n evtWeight %s' % (syst), evtWeight)


            ### General or GEN-level histograms ========================================================================================


            

            ## isMC --------------------------------------------------



            if (self.datasetInfo['isSignal'] or self.datasetInfo['isHToX']) and runMode_SignalGenChecks and syst == "Nom":
                output['hGenHiggsPt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenHiggsLog2Pt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenHiggsMass'].fill(
                    dataset=dataset,
                    Mass=(genHiggs.mass),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                mask_ = nGenZ > 0
                output['hGenZPt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genZ.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenZLog2Pt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genZ.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenZMass'].fill(
                    dataset=dataset,
                    Mass=(genZ.mass[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                mask_ = nGenW > 0
                output['hGenWPt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genW.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenWLog2Pt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genW.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenWMass'].fill(
                    dataset=dataset,
                    Mass=(genW.mass[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                mask_ = nGenWplus > 0
                output['hGenWplusPt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genWplus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenWplusLog2Pt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genWplus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenWplusMass'].fill(
                    dataset=dataset,
                    Mass=(genWplus.mass[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                mask_ = nGenWminus > 0
                output['hGenWminusPt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genWminus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenWminusLog2Pt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genWminus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenWminusMass'].fill(
                    dataset=dataset,
                    Mass=(genWminus.mass[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                mask_ = nGenTtbar > 0
                output['hGenTtbarPt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genTtbar.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenTtbarLog2Pt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genTtbar.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                mask_ = nGenQuarksFromHardScattring >= 2
                output['hGenQQFromHardScatteringPt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genQQFromHardScattering.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenQQFromHardScatteringLog2Pt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genQQFromHardScattering.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )

                # 2d plots
                mask_ = nGenZ > 0
                output['hGenPt_Higgs_vs_Z'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genZ.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_Z'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genZ.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByZ'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genZ.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                mask_ = nGenW > 0
                output['hGenPt_Higgs_vs_W'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genW.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_W'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genW.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByW'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genW.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                mask_ = nGenWplus > 0
                output['hGenPt_Higgs_vs_Wplus'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genWplus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_Wplus'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genWplus.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByWplus'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genWplus.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                mask_ = nGenWminus > 0
                output['hGenPt_Higgs_vs_Wminus'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genWminus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_Wminus'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genWminus.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByWminus'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genWminus.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                mask_ = nGenTtbar > 0
                output['hGenPt_Higgs_vs_Ttbar'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genTtbar.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_Ttbar'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genTtbar.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByTtbar'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genTtbar.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                mask_ = nGenQuarksFromHardScattring >= 2
                output['hGenPt_Higgs_vs_QQFromHardScattering'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genQQFromHardScattering.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_QQFromHardScattering'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genQQFromHardScattering.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByQQFromHardScattering'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genQQFromHardScattering.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight_gen[mask_]
                )

                

                
                
                

                



                output['hnGenHiggs'].fill(
                    dataset=dataset,
                    nObject10=(nGenHiggs),
                    systematic=syst,
                    weight=evtWeight_gen
                )   
                output['hnGenZ'].fill(
                    dataset=dataset,
                    nObject10=(nGenZ),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hnGenW'].fill(
                    dataset=dataset,
                    nObject10=(nGenW),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hnGenWplus'].fill(
                    dataset=dataset,
                    nObject10=(nGenWplus),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hnGenWminus'].fill(
                    dataset=dataset,
                    nObject10=(nGenWminus),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hnGenTop'].fill(
                    dataset=dataset,
                    nObject10=(nGenTop),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hnGenAntiTop'].fill(
                    dataset=dataset,
                    nObject10=(nGenAntiTop),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hnGenQuarksFromHardScattring'].fill(
                    dataset=dataset,
                    nObject10=(nGenQuarksFromHardScattring),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                
                   
                                        

                # With HiggsPtRewgt
                output['hGenHiggsPt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt),
                    systematic=syst,
                    weight=evtWeight
                )
                output['hGenHiggsLog2Pt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt),
                    systematic=syst,
                    weight=evtWeight
                ) 
                mask_ = nGenZ > 0
                output['hGenZPt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genZ.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenZLog2Pt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genZ.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                mask_ = nGenW > 0
                output['hGenWPt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genW.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenWLog2Pt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genW.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                mask_ = nGenWplus > 0
                output['hGenWplusPt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genWplus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenWplusLog2Pt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genWplus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                mask_ = nGenWminus > 0
                output['hGenWminusPt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genWminus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenWminusLog2Pt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genWminus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                mask_ = nGenTtbar > 0
                output['hGenTtbarPt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genTtbar.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenTtbarLog2Pt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genTtbar.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                mask_ = nGenQuarksFromHardScattring >= 2
                output['hGenQQFromHardScatteringPt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genQQFromHardScattering.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenQQFromHardScatteringLog2Pt_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genQQFromHardScattering.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )

                # 2d plots
                mask_ = nGenZ > 0
                output['hGenPt_Higgs_vs_Z_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genZ.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_Z_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genZ.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByZ_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genZ.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                mask_ = nGenW > 0
                output['hGenPt_Higgs_vs_W_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genW.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_W_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genW.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByW_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genW.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                mask_ = nGenWplus > 0
                output['hGenPt_Higgs_vs_Wplus_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genWplus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_Wplus_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genWplus.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByWplus_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genWplus.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                mask_ = nGenWminus > 0
                output['hGenPt_Higgs_vs_Wminus_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genWminus.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_Wminus_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genWminus.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByWminus_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genWminus.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                mask_ = nGenTtbar > 0
                output['hGenPt_Higgs_vs_Ttbar_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genTtbar.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_Ttbar_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genTtbar.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByTtbar_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genTtbar.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                mask_ = nGenQuarksFromHardScattring >= 2
                output['hGenPt_Higgs_vs_QQFromHardScattering_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Pt2TeV=(genHiggs.pt[mask_]),
                    Pt2TeV1=(genQQFromHardScattering.pt[mask_]),
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_QQFromHardScattering_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Pt2TeV1=np.log2(genQQFromHardScattering.pt)[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )
                output['hGenLog2Pt_Higgs_vs_HiggsByQQFromHardScattering_wHiggsPtRewgt'].fill(
                    dataset=dataset,
                    Log2Pt2TeV=np.log2(genHiggs.pt)[mask_],
                    Log2Ratio=np.log2( 2*genHiggs.pt / (genHiggs.pt + genQQFromHardScattering.pt) )[mask_],
                    systematic=syst,
                    weight=evtWeight[mask_]
                )


            ## isMC && isSignal ------------------------------------------------------------------------------------------------------------












            ### RECO-level histograms ============================================================================================
            

















































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
    print(f"htoaa_Analysis_HiggsPtRewgt:: here14 {datetime.now() = }")

    if len(sys.argv) != 2:
        print("htoaa_Analysis:: Command-line config file missing.. \t **** ERROR **** \n")

    sConfig = sys.argv[1]

    
    config = GetDictFromJsonFile(sConfig)
    print("Config {}: \n{}".format(sConfig, json.dumps(config, indent=4)))
    print(f"htoaa_Analysis_HiggsPtRewgt:: here15 {datetime.now() = }")

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
                logging.critical(f'htoaa_Analysis_HiggsPtRewgt.py::main():: {MCSamplesStitchInputFileName = } does not exists')
                print(f'htoaa_Analysis_HiggsPtRewgt.py::main() 11:: {MCSamplesStitchInputFileName = } does not exists')
                exit(0)
            print(f"Opening {MCSamplesStitchInputFileName = } "); sys.stdout.flush() 
            with uproot.open(MCSamplesStitchInputFileName) as f_:
                print(f"{f_.keys() = }"); sys.stdout.flush() 
                hMCSamplesStitch = f_[r'%s' % MCSamplesStitchInputHistogramName].to_hist()

    print(f"htoaa_Analysis_HiggsPtRewgt:: here16 {datetime.now() = }")    
        
        
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
    print(f"htoaa_Analysis_HiggsPtRewgt:: here17 {datetime.now() = }")

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
            logging.critical('htoaa_Analysis_HiggsPtRewgt:: getNanoAODFile() for input file %s failed. **** CRITICAL ERROR ****. \nAborting...' % (sInputFile)); sys.stdout.flush();
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
        print(f"htoaa_Analysis_HiggsPtRewgt:: {sInputFile} \t {os.path.exists(sInputFile) = }, {fileSize = } MB");     

        if fileSize > NanoAODFileSize_Min:     
            sInputFiles[iFile] = sInputFile
        else:
            logging.critical('htoaa_Analysis_HiggsPtRewgt:: Input file %s file size below threshold (%g MB). **** CRITICAL ERROR ****. \nAborting...' % (sInputFile, NanoAODFileSize_Min) ); sys.stdout.flush();
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
    print(f"htoaa_Analysis_HiggsPtRewgt:: here18 {datetime.now() = }")


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
    print(f"htoaa_Analysis_HiggsPtRewgt:: here19 {datetime.now() = }", flush=flushStdout)
        
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

            print("%10f\t%10d\t%s" % (output['cutflow'][sWeighted+key], output['cutflow'][key], key), flush=flushStdout)
            #print("%10d\t%s" % (output['cutflow'][key], key), flush=flushStdout)
    
    
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
            if printLevel >= 1: print(f"Opening output file {sOutputFile}")
            for key, value in output.items():
                if printLevel >= 1: print(f"key: {key},  value ({type(value)}): {value}")
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
                    if printLevel >= 1: print(f"_dataset ({type(_dataset)}): {_dataset}")

                    #print(f"value.fields : {value.fields}", fl)
                    #print(f"value.fields() ({type(value.fields())}): {value.fields()}")
                    if 'systematic' not in value.fields: # hWeight histogram do not have 'systematics' axes
                        h1 = value.integrate('dataset',_dataset).to_hist()
                        fOut['%s/%s' % (sDir1_toUse, sHistoName_toUse)] = h1
                        continue

                    for _syst in value.axis('systematic').identifiers():
                        if printLevel >= 1: print(f"_syst ({type(_syst)}): {_syst}")

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
    
