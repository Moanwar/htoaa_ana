#htoaa analysis main code

import os
import sys
from datetime import datetime
#import time
print(f"htoaa_Analysis_GGFMode:: here1 {datetime.now() = }"); sys.stdout.flush()
import subprocess
import json
import csv
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
import hist as histLib
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
    calculate_AbsDeltaPhi, 
    selectRunLuminosityBlock,
    calculate_lumiScale, getLumiScaleForPhSpOverlapRewgtMode, getSampleHTRange, # update_crosssection, 
    getNanoAODFile, setXRootDRedirector,  xrdcpFile,
    selectMETFilters, selectFatJets, getCandidateHiggs, selectAK4Jets, selectMuons, selectElectrons,
    selGenPartsWithStatusFlag,
    getHToAATo4BLundPlaneRewgt, getHiggsPtRewgtForGGToHToAATo4B, 
    getTopPtRewgt, getPURewgts, getHTReweight,
    getPURewgts_variation, get_jetTriggerSF, get_PSWeight, add_pdf_as_weight, get_QCDScaleWeight,
    get_JER_and_JES,
    get_Ak4BtagSF, get_L1TPrefiringWgt,
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


HLTTriggers = {
    Era_2016:[
        "HLT_AK8PFJet360_TrimMass30",
        "HLT_AK8PFJet400_TrimMass30",
        "HLT_AK8PFHT750_TrimMass50",
        "HLT_AK8PFHT800_TrimMass50",
        "HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20",
        "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087",
        "HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087",
        "HLT_AK8DiPFJet300_200_TrimMass30",
        "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50",
        "HLT_AK8PFHT650_TrimR0p1PT0p03Mass50",
        "HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20",
        "HLT_AK8DiPFJet280_200_TrimMass30",
        "HLT_AK8DiPFJet250_200_TrimMass30",
        "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20",
        "HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20",
        "HLT_CaloJet260",
        "HLT_CaloJet500_NoJetID",
        "HLT_Dimuon13_PsiPrime",
        "HLT_Dimuon13_Upsilon",
        "HLT_Dimuon20_Jpsi",
        "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf",
        "HLT_DoubleEle25_CaloIdL_GsfTrkIdVL",
        "HLT_DoubleEle33_CaloIdL",
        "HLT_DoubleEle33_CaloIdL_MW",
        "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW",
        "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",
        "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg",
        "HLT_DoubleTightCombinedIsoPFTau35_Trk1_eta2p1_Reg",
        "HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg",
        "HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1_Reg",
        "HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1",
        "HLT_DoubleTightCombinedIsoPFTau40_Trk1_eta2p1",
        "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
        "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg",
        "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1",
        "HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL",
        "HLT_DoubleMu33NoFiltersNoVtx",
        "HLT_DoubleMu38NoFiltersNoVtx",
        "HLT_DoubleMu23NoFiltersNoVtxDisplaced",
        "HLT_DoubleMu28NoFiltersNoVtxDisplaced",
        "HLT_DoubleMu0",
        "HLT_DoubleMu4_3_Bs",
        "HLT_DoubleMu4_3_Jpsi_Displaced",
        "HLT_DoubleMu4_JpsiTrk_Displaced",
        "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced",
        "HLT_DoubleMu3_Trk_Tau3mu",
        "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
        "HLT_Mu7p5_L2Mu2_Jpsi",
        "HLT_Mu7p5_L2Mu2_Upsilon",
        "HLT_Mu7p5_Track2_Jpsi",
        "HLT_Mu7p5_Track3p5_Jpsi",
        "HLT_Mu7p5_Track7_Jpsi",
        "HLT_Mu7p5_Track2_Upsilon",
        "HLT_Mu7p5_Track3p5_Upsilon",
        "HLT_Mu7p5_Track7_Upsilon",
        "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",
        "HLT_Dimuon0er16_Jpsi_NoVertexing",
        "HLT_Dimuon6_Jpsi_NoVertexing",
        "HLT_Photon150",
        "HLT_Photon90_CaloIdL_HT300",
        "HLT_HT250_CaloMET70",
        "HLT_DoublePhoton60",
        "HLT_DoublePhoton85",
        "HLT_Ele17_Ele8_Gsf",
        "HLT_Ele20_eta2p1_WPLoose_Gsf_LooseIsoPFTau28",
        "HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau29",
        "HLT_Ele22_eta2p1_WPLoose_Gsf",
        "HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1",
        "HLT_Ele23_WPLoose_Gsf",
        "HLT_Ele23_WPLoose_Gsf_WHbbBoost",
        "HLT_Ele24_eta2p1_WPLoose_Gsf",
        "HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20",
        "HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1",
        "HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30",
        "HLT_Ele25_WPTight_Gsf",
        "HLT_Ele25_eta2p1_WPLoose_Gsf",
        "HLT_Ele25_eta2p1_WPTight_Gsf",
        "HLT_Ele27_WPLoose_Gsf",
        "HLT_Ele27_WPLoose_Gsf_WHbbBoost",
        "HLT_Ele27_WPTight_Gsf",
        "HLT_Ele27_WPTight_Gsf_L1JetTauSeeded",
        "HLT_Ele27_eta2p1_WPLoose_Gsf",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1",
        "HLT_Ele27_eta2p1_WPTight_Gsf",
        "HLT_Ele30_WPTight_Gsf",
        "HLT_Ele30_eta2p1_WPLoose_Gsf",
        "HLT_Ele30_eta2p1_WPTight_Gsf",
        "HLT_Ele32_WPTight_Gsf",
        "HLT_Ele32_eta2p1_WPLoose_Gsf",
        "HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1",
        "HLT_Ele32_eta2p1_WPTight_Gsf",
        "HLT_Ele35_WPLoose_Gsf",
        "HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50",
        "HLT_Ele36_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1",
        "HLT_Ele45_WPLoose_Gsf",
        "HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded",
        "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50",
        "HLT_Ele105_CaloIdVT_GsfTrkIdT",
        "HLT_Ele30WP60_SC4_Mass55",
        "HLT_Ele30WP60_Ele8_Mass55",
        "HLT_HT200",
        "HLT_HT275",
        "HLT_HT325",
        "HLT_HT425",
        "HLT_HT575",
        "HLT_HT410to430",
        "HLT_HT430to450",
        "HLT_HT450to470",
        "HLT_HT470to500",
        "HLT_HT500to550",
        "HLT_HT550to650",
        "HLT_HT650",
        "HLT_Mu16_eta2p1_MET30",
        "HLT_IsoMu16_eta2p1_MET30",
        "HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1",
        "HLT_IsoMu17_eta2p1",
        "HLT_IsoMu17_eta2p1_LooseIsoPFTau20",
        "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1",
        "HLT_DoubleIsoMu17_eta2p1",
        "HLT_DoubleIsoMu17_eta2p1_noDzCut",
        "HLT_IsoMu18",
        "HLT_IsoMu19_eta2p1_LooseIsoPFTau20",
        "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1",
        "HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg",
        "HLT_IsoMu19_eta2p1_LooseCombinedIsoPFTau20",
        "HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg",
        "HLT_IsoMu19_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg",
        "HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg",
        "HLT_IsoMu21_eta2p1_TightCombinedIsoPFTau32_Trk1_eta2p1_Reg",
        "HLT_IsoMu20",
        "HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1",
        "HLT_IsoMu21_eta2p1_LooseIsoPFTau50_Trk30_eta2p1_SingleL1",
        "HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg",
        "HLT_IsoMu22",
        "HLT_IsoMu22_eta2p1",
        "HLT_IsoMu24",
        "HLT_IsoMu27",
        "HLT_IsoTkMu18",
        "HLT_IsoTkMu20",
        "HLT_IsoTkMu22",
        "HLT_IsoTkMu22_eta2p1",
        "HLT_IsoTkMu24",
        "HLT_IsoTkMu27",
        "HLT_JetE30_NoBPTX3BX",
        "HLT_JetE30_NoBPTX",
        "HLT_JetE50_NoBPTX3BX",
        "HLT_JetE70_NoBPTX3BX",
        "HLT_L1SingleMu18",
        "HLT_L2Mu10",
        "HLT_L1SingleMuOpen",
        "HLT_L1SingleMuOpen_DT",
        "HLT_L2DoubleMu23_NoVertex",
        "HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10",
        "HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10",
        "HLT_L2Mu10_NoVertex_NoBPTX3BX",
        "HLT_L2Mu10_NoVertex_NoBPTX",
        "HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX",
        "HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120",
        "HLT_PFTau120_eta2p1",
        "HLT_PFTau140_eta2p1",
        "HLT_VLooseIsoPFTau120_Trk50_eta2p1",
        "HLT_VLooseIsoPFTau140_Trk50_eta2p1",
        "HLT_Mu17_Mu8",
        "HLT_Mu17_Mu8_DZ",
        "HLT_Mu17_Mu8_SameSign",
        "HLT_Mu17_Mu8_SameSign_DZ",
        "HLT_Mu20_Mu10",
        "HLT_Mu20_Mu10_DZ",
        "HLT_Mu20_Mu10_SameSign",
        "HLT_Mu20_Mu10_SameSign_DZ",
        "HLT_Mu17_TkMu8_DZ",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
        "HLT_Mu25_TkMu0_dEta18_Onia",
        "HLT_Mu27_TkMu8",
        "HLT_Mu30_TkMu11",
        "HLT_Mu30_eta2p1_PFJet150_PFJet50",
        "HLT_Mu40_TkMu11",
        "HLT_Mu40_eta2p1_PFJet200_PFJet50",
        "HLT_Mu20",
        "HLT_TkMu17",
        "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
        "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
        "HLT_TkMu20",
        "HLT_Mu24_eta2p1",
        "HLT_TkMu24_eta2p1",
        "HLT_Mu27",
        "HLT_TkMu27",
        "HLT_Mu45_eta2p1",
        "HLT_Mu50",
        "HLT_TkMu50",
        "HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL",
        "HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL",
        "HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL",
        "HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL",
        "HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL",
        "HLT_DoubleMu18NoFiltersNoVtx",
        "HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight",
        "HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose",
        "HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose",
        "HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight",
        "HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose",
        "HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose",
        "HLT_Mu28NoFiltersNoVtx_CentralCaloJet40",
        "HLT_PFHT300_PFMET100",
        "HLT_PFHT300_PFMET110",
        "HLT_PFHT550_4JetPt50",
        "HLT_PFHT650_4JetPt50",
        "HLT_PFHT750_4JetPt50",
        "HLT_PFHT750_4JetPt70",
        "HLT_PFHT750_4JetPt80",
        "HLT_PFHT800_4JetPt50",
        "HLT_PFHT850_4JetPt50",
        "HLT_PFJet15_NoCaloMatched",
        "HLT_PFJet25_NoCaloMatched",
        "HLT_DiPFJet15_NoCaloMatched",
        "HLT_DiPFJet25_NoCaloMatched",
        "HLT_DiPFJet15_FBEta3_NoCaloMatched",
        "HLT_DiPFJet25_FBEta3_NoCaloMatched",
        "HLT_DiPFJetAve15_HFJEC",
        "HLT_DiPFJetAve25_HFJEC",
        "HLT_DiPFJetAve35_HFJEC",
        "HLT_AK8PFJet40",
        "HLT_AK8PFJet60",
        "HLT_AK8PFJet80",
        "HLT_AK8PFJet140",
        "HLT_AK8PFJet200",
        "HLT_AK8PFJet260",
        "HLT_AK8PFJet320",
        "HLT_AK8PFJet400",
        "HLT_AK8PFJet450",
        "HLT_AK8PFJet500",
        "HLT_PFJet40",
        "HLT_PFJet60",
        "HLT_PFJet80",
        "HLT_PFJet140",
        "HLT_PFJet200",
        "HLT_PFJet260",
        "HLT_PFJet320",
        "HLT_PFJet400",
        "HLT_PFJet450",
        "HLT_PFJet500",
        "HLT_DiPFJetAve40",
        "HLT_DiPFJetAve60",
        "HLT_DiPFJetAve80",
        "HLT_DiPFJetAve140",
        "HLT_DiPFJetAve200",
        "HLT_DiPFJetAve260",
        "HLT_DiPFJetAve320",
        "HLT_DiPFJetAve400",
        "HLT_DiPFJetAve500",
        "HLT_DiPFJetAve60_HFJEC",
        "HLT_DiPFJetAve80_HFJEC",
        "HLT_DiPFJetAve100_HFJEC",
        "HLT_DiPFJetAve160_HFJEC",
        "HLT_DiPFJetAve220_HFJEC",
        "HLT_DiPFJetAve300_HFJEC",
        "HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140",
        "HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80",
        "HLT_DiCentralPFJet170",
        "HLT_SingleCentralPFJet170_CFMax0p1",
        "HLT_DiCentralPFJet170_CFMax0p1",
        "HLT_DiCentralPFJet220_CFMax0p3",
        "HLT_DiCentralPFJet330_CFMax0p5",
        "HLT_DiCentralPFJet430",
        "HLT_PFHT125",
        "HLT_PFHT200",
        "HLT_PFHT250",
        "HLT_PFHT300",
        "HLT_PFHT350",
        "HLT_PFHT400",
        "HLT_PFHT475",
        "HLT_PFHT600",
        "HLT_PFHT650",
        "HLT_PFHT800",
        "HLT_PFHT900",
        "HLT_PFHT200_PFAlphaT0p51",
        "HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57",
        "HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63",
        "HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55",
        "HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58",
        "HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53",
        "HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54",
        "HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52",
        "HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53",
        "HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51",
        "HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52",
        "HLT_MET60_IsoTrk35_Loose",
        "HLT_MET75_IsoTrk50",
        "HLT_MET90_IsoTrk50",
        "HLT_PFMET120_BTagCSV_p067",
        "HLT_PFMET120_Mu5",
        "HLT_PFMET170_NotCleaned",
        "HLT_PFMET170_NoiseCleaned",
        "HLT_PFMET170_HBHECleaned",
        "HLT_PFMET170_JetIdCleaned",
        "HLT_PFMET170_BeamHaloCleaned",
        "HLT_PFMET170_HBHE_BeamHaloCleaned",
        "HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned",
        "HLT_PFMET90_PFMHT90_IDTight",
        "HLT_PFMET100_PFMHT100_IDTight",
        "HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned",
        "HLT_PFMET110_PFMHT110_IDTight",
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067",
        "HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight",
        "HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200",
        "HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460",
        "HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240",
        "HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500",
        "HLT_QuadPFJet_VBF",
        "HLT_L1_TripleJet_VBF",
        "HLT_QuadJet45_TripleBTagCSV_p087",
        "HLT_QuadJet45_DoubleBTagCSV_p087",
        "HLT_DoubleJet90_Double30_TripleBTagCSV_p087",
        "HLT_DoubleJet90_Double30_DoubleBTagCSV_p087",
        "HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160",
        "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6",
        "HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172",
        "HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6",
        "HLT_DoubleJetsC100_SingleBTagCSV_p026",
        "HLT_DoubleJetsC100_SingleBTagCSV_p014",
        "HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350",
        "HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350",
        "HLT_Photon135_PFMET100",
        "HLT_Photon20_CaloIdVL_IsoL",
        "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40",
        "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF",
        "HLT_Photon250_NoHE",
        "HLT_Photon300_NoHE",
        "HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60",
        "HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15",
        "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40",
        "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF",
        "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40",
        "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF",
        "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40",
        "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF",
        "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40",
        "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF",
        "HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40",
        "HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF",
        "HLT_Mu8_TrkIsoVVL",
        "HLT_Mu17_TrkIsoVVL",
        "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
        "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30",
        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30",
        "HLT_BTagMu_DiJet20_Mu5",
        "HLT_BTagMu_DiJet40_Mu5",
        "HLT_BTagMu_DiJet70_Mu5",
        "HLT_BTagMu_DiJet110_Mu5",
        "HLT_BTagMu_DiJet170_Mu5",
        "HLT_BTagMu_Jet300_Mu5",
        "HLT_BTagMu_AK8Jet300_Mu5",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
        "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL",
        "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL",
        "HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL",
        "HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL",
        "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
        "HLT_Mu12_Photon25_CaloIdL",
        "HLT_Mu12_Photon25_CaloIdL_L1ISO",
        "HLT_Mu12_Photon25_CaloIdL_L1OR",
        "HLT_Mu17_Photon22_CaloIdL_L1ISO",
        "HLT_Mu17_Photon30_CaloIdL_L1ISO",
        "HLT_Mu17_Photon35_CaloIdL_L1ISO",
        "HLT_DiMu9_Ele9_CaloIdL_TrackIdL",
        "HLT_TripleMu_5_3_3",
        "HLT_TripleMu_12_10_5",
        "HLT_Mu3er_PFHT140_PFMET125",
        "HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067",
        "HLT_Mu6_PFHT200_PFMET100",
        "HLT_Mu14er_PFMET100",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Ele17_CaloIdL_GsfTrkIdVL",
        "HLT_Ele17_CaloIdL_TrackIdL_IsoVL",
        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL",
        "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5",
        "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5",
        "HLT_Photon22",
        "HLT_Photon30",
        "HLT_Photon36",
        "HLT_Photon50",
        "HLT_Photon75",
        "HLT_Photon90",
        "HLT_Photon120",
        "HLT_Photon175",
        "HLT_Photon165_HE10",
        "HLT_Photon22_R9Id90_HE10_IsoM",
        "HLT_Photon30_R9Id90_HE10_IsoM",
        "HLT_Photon36_R9Id90_HE10_IsoM",
        "HLT_Photon50_R9Id90_HE10_IsoM",
        "HLT_Photon75_R9Id90_HE10_IsoM",
        "HLT_Photon90_R9Id90_HE10_IsoM",
        "HLT_Photon120_R9Id90_HE10_IsoM",
        "HLT_Photon165_R9Id90_HE10_IsoM",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70",
        "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55",
        "HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55",
        "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55",
        "HLT_Dimuon0_Jpsi_Muon",
        "HLT_Dimuon0_Upsilon_Muon",
        "HLT_QuadMuon0_Dimuon0_Jpsi",
        "HLT_QuadMuon0_Dimuon0_Upsilon",
        "HLT_Rsq0p25_Calo",
        "HLT_RsqMR240_Rsq0p09_MR200_4jet_Calo",
        "HLT_RsqMR240_Rsq0p09_MR200_Calo",
        "HLT_Rsq0p25",
        "HLT_Rsq0p30",
        "HLT_RsqMR240_Rsq0p09_MR200",
        "HLT_RsqMR240_Rsq0p09_MR200_4jet",
        "HLT_RsqMR270_Rsq0p09_MR200",
        "HLT_RsqMR270_Rsq0p09_MR200_4jet",
        "HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200",
        "HLT_Rsq0p02_MR400_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200",
        "HLT_Rsq0p02_MR450_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200",
        "HLT_Rsq0p02_MR500_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200",
        "HLT_Rsq0p02_MR550_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200",
        "HLT_HT200_DisplacedDijet40_DisplacedTrack",
        "HLT_HT250_DisplacedDijet40_DisplacedTrack",
        "HLT_HT350_DisplacedDijet40_DisplacedTrack",
        "HLT_HT350_DisplacedDijet80_DisplacedTrack",
        "HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack",
        "HLT_HT350_DisplacedDijet40_Inclusive",
        "HLT_HT400_DisplacedDijet40_Inclusive",
        "HLT_HT500_DisplacedDijet40_Inclusive",
        "HLT_HT550_DisplacedDijet40_Inclusive",
        "HLT_HT550_DisplacedDijet80_Inclusive",
        "HLT_HT650_DisplacedDijet80_Inclusive",
        "HLT_HT750_DisplacedDijet80_Inclusive",
        "HLT_VBF_DisplacedJet40_DisplacedTrack",
        "HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5",
        "HLT_VBF_DisplacedJet40_TightID_DisplacedTrack",
        "HLT_VBF_DisplacedJet40_Hadronic",
        "HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack",
        "HLT_VBF_DisplacedJet40_TightID_Hadronic",
        "HLT_VBF_DisplacedJet40_VTightID_Hadronic",
        "HLT_VBF_DisplacedJet40_VVTightID_Hadronic",
        "HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack",
        "HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack",
        "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight",
        "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight",
        "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_HT200",
        "HLT_Photon90_CaloIdL_PFHT500",
        "HLT_DoubleMu8_Mass8_PFHT250",
        "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250",
        "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250",
        "HLT_DoubleMu8_Mass8_PFHT300",
        "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300",
        "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300",
        "HLT_Mu10_CentralPFJet30_BTagCSV_p13",
        "HLT_DoubleMu3_PFMET50",
        "HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13",
        "HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400",
        "HLT_Ele15_IsoVVVL_PFHT350_PFMET50",
        "HLT_Ele15_IsoVVVL_PFHT600",
        "HLT_Ele15_IsoVVVL_PFHT350",
        "HLT_Ele15_IsoVVVL_PFHT400_PFMET50",
        "HLT_Ele15_IsoVVVL_PFHT400",
        "HLT_Ele50_IsoVVVL_PFHT400",
        "HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60",
        "HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60",
        "HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400",
        "HLT_Mu15_IsoVVVL_PFHT350_PFMET50",
        "HLT_Mu15_IsoVVVL_PFHT600",
        "HLT_Mu15_IsoVVVL_PFHT350",
        "HLT_Mu15_IsoVVVL_PFHT400_PFMET50",
        "HLT_Mu15_IsoVVVL_PFHT400",
        "HLT_Mu50_IsoVVVL_PFHT400",
        "HLT_Dimuon16_Jpsi",
        "HLT_Dimuon10_Jpsi_Barrel",
        "HLT_Dimuon8_PsiPrime_Barrel",
        "HLT_Dimuon8_Upsilon_Barrel",
        "HLT_Dimuon0_Phi_Barrel",
        "HLT_Mu16_TkMu0_dEta18_Onia",
        "HLT_Mu16_TkMu0_dEta18_Phi",
        "HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx",
        "HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx",
        "HLT_Mu8",
        "HLT_Mu17",
        "HLT_Mu3_PFJet40",
        "HLT_Ele8_CaloIdM_TrackIdM_PFJet30",
        "HLT_Ele12_CaloIdM_TrackIdM_PFJet30",
        "HLT_Ele17_CaloIdM_TrackIdM_PFJet30",
        "HLT_Ele23_CaloIdM_TrackIdM_PFJet30",
        "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140",
        "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165",
        "HLT_PFHT400_SixJet30_DoubleBTagCSV_p056",
        "HLT_PFHT450_SixJet40_BTagCSV_p056",
        "HLT_PFHT400_SixJet30",
        "HLT_PFHT450_SixJet40",
        "HLT_Ele115_CaloIdVT_GsfTrkIdT",
        "HLT_Mu55",
        "HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15",
        "HLT_Photon90_CaloIdL_PFHT600",
        "HLT_PixelTracks_Multiplicity60ForEndOfFill",
        "HLT_PixelTracks_Multiplicity85ForEndOfFill",
        "HLT_PixelTracks_Multiplicity110ForEndOfFill",
        "HLT_PixelTracks_Multiplicity135ForEndOfFill",
        "HLT_PixelTracks_Multiplicity160ForEndOfFill",
        "HLT_FullTracks_Multiplicity80",
        "HLT_FullTracks_Multiplicity100",
        "HLT_FullTracks_Multiplicity130",
        "HLT_FullTracks_Multiplicity150",
        "HLT_ECALHT800",
        "HLT_DiSC30_18_EIso_AND_HE_Mass70",
        "HLT_Photon125",
        "HLT_MET100",
        "HLT_MET150",
        "HLT_MET200",
        "HLT_Ele27_HighEta_Ele20_Mass55",
        "HLT_L1FatEvents",
        "HLT_Physics",
        "HLT_L1FatEvents_part0",
        "HLT_L1FatEvents_part1",
        "HLT_L1FatEvents_part2",
        "HLT_L1FatEvents_part3",
        "HLT_Random",
        "HLT_ZeroBias",
        "HLT_AK4CaloJet30",
        "HLT_AK4CaloJet40",
        "HLT_AK4CaloJet50",
        "HLT_AK4CaloJet80",
        "HLT_AK4CaloJet100",
        "HLT_AK4PFJet30",
        "HLT_AK4PFJet50",
        "HLT_AK4PFJet80",
        "HLT_AK4PFJet100",
        "HLT_HISinglePhoton10",
        "HLT_HISinglePhoton15",
        "HLT_HISinglePhoton20",
        "HLT_HISinglePhoton40",
        "HLT_HISinglePhoton60",
        "HLT_EcalCalibration",
        "HLT_HcalCalibration",
        "HLT_GlobalRunHPDNoise",
        "HLT_L1BptxMinus",
        "HLT_L1BptxPlus",
        "HLT_L1NotBptxOR",
        "HLT_L1BeamGasMinus",
        "HLT_L1BeamGasPlus",
        "HLT_L1BptxXOR",
        "HLT_L1MinimumBiasHF_OR",
        "HLT_L1MinimumBiasHF_AND",
        "HLT_HcalNZS",
        "HLT_HcalPhiSym",
        "HLT_HcalIsolatedbunch",
        "HLT_ZeroBias_FirstCollisionAfterAbortGap",
        "HLT_ZeroBias_FirstCollisionAfterAbortGap_copy",
        "HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS",
        "HLT_ZeroBias_IsolatedBunches",
        "HLT_ZeroBias_FirstCollisionInTrain",
        "HLT_ZeroBias_FirstBXAfterTrain",
        "HLT_Photon500",
        "HLT_Photon600",
        "HLT_Mu300",
        "HLT_Mu350",
        "HLT_MET250",
        "HLT_MET300",
        "HLT_MET600",
        "HLT_MET700",
        "HLT_PFMET300",
        "HLT_PFMET400",
        "HLT_PFMET500",
        "HLT_PFMET600",
        "HLT_Ele250_CaloIdVT_GsfTrkIdT",
        "HLT_Ele300_CaloIdVT_GsfTrkIdT",
        "HLT_HT2000",
        "HLT_HT2500",
        "HLT_IsoTrackHE",
        "HLT_IsoTrackHB"
    ],

    Era_2017: [
        "HLT_AK8PFJet360_TrimMass30",
        "HLT_AK8PFJet380_TrimMass30",
        "HLT_AK8PFJet400_TrimMass30",
        "HLT_AK8PFJet420_TrimMass30",
        "HLT_AK8PFHT750_TrimMass50",
        "HLT_AK8PFHT800_TrimMass50",
        "HLT_AK8PFHT850_TrimMass50",
        "HLT_AK8PFHT900_TrimMass50",
        "HLT_CaloJet500_NoJetID",
        "HLT_CaloJet550_NoJetID",
        "HLT_Trimuon5_3p5_2_Upsilon_Muon",
        "HLT_DoubleEle25_CaloIdL_MW",
        "HLT_DoubleEle27_CaloIdL_MW",
        "HLT_DoubleEle33_CaloIdL_MW",
        "HLT_DoubleEle24_eta2p1_WPTight_Gsf",
        "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350",
        "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350",
        "HLT_Ele27_Ele37_CaloIdL_MW",
        "HLT_Mu27_Ele37_CaloIdL_MW",
        "HLT_Mu37_Ele27_CaloIdL_MW",
        "HLT_Mu37_TkMu27",
        "HLT_DoubleMu4_3_Bs",
        "HLT_DoubleMu4_3_Jpsi_Displaced",
        "HLT_DoubleMu4_JpsiTrk_Displaced",
        "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced",
        "HLT_DoubleMu3_Trk_Tau3mu",
        "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
        "HLT_DoubleMu4_Mass8_DZ_PFHT350",
        "HLT_DoubleMu8_Mass8_PFHT350",
        "HLT_Mu3_PFJet40",
        "HLT_Mu7p5_L2Mu2_Jpsi",
        "HLT_Mu7p5_L2Mu2_Upsilon",
        "HLT_Mu7p5_Track2_Jpsi",
        "HLT_Mu7p5_Track3p5_Jpsi",
        "HLT_Mu7p5_Track7_Jpsi",
        "HLT_Mu7p5_Track2_Upsilon",
        "HLT_Mu7p5_Track3p5_Upsilon",
        "HLT_Mu7p5_Track7_Upsilon",
        "HLT_DoublePhoton33_CaloIdL",
        "HLT_DoublePhoton70",
        "HLT_DoublePhoton85",
        "HLT_Ele20_WPTight_Gsf",
        "HLT_Ele20_WPLoose_Gsf",
        "HLT_Ele20_eta2p1_WPLoose_Gsf",
        "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG",
        "HLT_Ele27_WPTight_Gsf",
        "HLT_Ele32_WPTight_Gsf",
        "HLT_Ele35_WPTight_Gsf",
        "HLT_Ele35_WPTight_Gsf_L1EGMT",
        "HLT_Ele38_WPTight_Gsf",
        "HLT_Ele40_WPTight_Gsf",
        "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
        "HLT_HT450_Beamspot",
        "HLT_HT300_Beamspot",
        "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1",
        "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1",
        "HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1",
        "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1",
        "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1",
        "HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1",
        "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1",
        "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1",
        "HLT_IsoMu20",
        "HLT_IsoMu24",
        "HLT_IsoMu24_eta2p1",
        "HLT_IsoMu27",
        "HLT_IsoMu30",
        "HLT_UncorrectedJetE30_NoBPTX",
        "HLT_UncorrectedJetE30_NoBPTX3BX",
        "HLT_UncorrectedJetE60_NoBPTX3BX",
        "HLT_UncorrectedJetE70_NoBPTX3BX",
        "HLT_L1SingleMu18",
        "HLT_L1SingleMu25",
        "HLT_L2Mu10",
        "HLT_L2Mu10_NoVertex_NoBPTX3BX",
        "HLT_L2Mu10_NoVertex_NoBPTX",
        "HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX",
        "HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX",
        "HLT_L2Mu50",
        "HLT_DoubleL2Mu50",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8",
        "HLT_Mu25_TkMu0_Onia",
        "HLT_Mu30_TkMu0_Onia",
        "HLT_Mu20_TkMu0_Phi",
        "HLT_Mu25_TkMu0_Phi",
        "HLT_Mu20",
        "HLT_Mu27",
        "HLT_Mu50",
        "HLT_Mu55",
        "HLT_OldMu100",
        "HLT_TkMu100",
        "HLT_DiPFJet15_NoCaloMatched",
        "HLT_DiPFJet25_NoCaloMatched",
        "HLT_DiPFJet15_FBEta3_NoCaloMatched",
        "HLT_DiPFJet25_FBEta3_NoCaloMatched",
        "HLT_DiPFJetAve40",
        "HLT_DiPFJetAve60",
        "HLT_DiPFJetAve80",
        "HLT_DiPFJetAve140",
        "HLT_DiPFJetAve200",
        "HLT_DiPFJetAve260",
        "HLT_DiPFJetAve320",
        "HLT_DiPFJetAve400",
        "HLT_DiPFJetAve500",
        "HLT_DiPFJetAve15_HFJEC",
        "HLT_DiPFJetAve25_HFJEC",
        "HLT_DiPFJetAve35_HFJEC",
        "HLT_DiPFJetAve60_HFJEC",
        "HLT_DiPFJetAve80_HFJEC",
        "HLT_DiPFJetAve100_HFJEC",
        "HLT_DiPFJetAve160_HFJEC",
        "HLT_DiPFJetAve220_HFJEC",
        "HLT_DiPFJetAve300_HFJEC",
        "HLT_AK8PFJet40",
        "HLT_AK8PFJet60",
        "HLT_AK8PFJet80",
        "HLT_AK8PFJet140",
        "HLT_AK8PFJet200",
        "HLT_AK8PFJet260",
        "HLT_AK8PFJet320",
        "HLT_AK8PFJet400",
        "HLT_AK8PFJet450",
        "HLT_AK8PFJet500",
        "HLT_AK8PFJet550",
        "HLT_PFJet40",
        "HLT_PFJet60",
        "HLT_PFJet80",
        "HLT_PFJet140",
        "HLT_PFJet200",
        "HLT_PFJet260",
        "HLT_PFJet320",
        "HLT_PFJet400",
        "HLT_PFJet450",
        "HLT_PFJet500",
        "HLT_PFJet550",
        "HLT_PFJetFwd40",
        "HLT_PFJetFwd60",
        "HLT_PFJetFwd80",
        "HLT_PFJetFwd140",
        "HLT_PFJetFwd200",
        "HLT_PFJetFwd260",
        "HLT_PFJetFwd320",
        "HLT_PFJetFwd400",
        "HLT_PFJetFwd450",
        "HLT_PFJetFwd500",
        "HLT_AK8PFJetFwd40",
        "HLT_AK8PFJetFwd60",
        "HLT_AK8PFJetFwd80",
        "HLT_AK8PFJetFwd140",
        "HLT_AK8PFJetFwd200",
        "HLT_AK8PFJetFwd260",
        "HLT_AK8PFJetFwd320",
        "HLT_AK8PFJetFwd400",
        "HLT_AK8PFJetFwd450",
        "HLT_AK8PFJetFwd500",
        "HLT_PFHT180",
        "HLT_PFHT250",
        "HLT_PFHT370",
        "HLT_PFHT430",
        "HLT_PFHT510",
        "HLT_PFHT590",
        "HLT_PFHT680",
        "HLT_PFHT780",
        "HLT_PFHT890",
        "HLT_PFHT1050",
        "HLT_PFHT500_PFMET100_PFMHT100_IDTight",
        "HLT_PFHT500_PFMET110_PFMHT110_IDTight",
        "HLT_PFHT700_PFMET85_PFMHT85_IDTight",
        "HLT_PFHT700_PFMET95_PFMHT95_IDTight",
        "HLT_PFHT800_PFMET75_PFMHT75_IDTight",
        "HLT_PFHT800_PFMET85_PFMHT85_IDTight",
        "HLT_PFMET110_PFMHT110_IDTight",
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_PFMET130_PFMHT130_IDTight",
        "HLT_PFMET140_PFMHT140_IDTight",
        "HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1",
        "HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1",
        "HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1",
        "HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1",
        "HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1",
        "HLT_PFMET120_PFMHT120_IDTight_PFHT60",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
        "HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60",
        "HLT_PFMETTypeOne110_PFMHT110_IDTight",
        "HLT_PFMETTypeOne120_PFMHT120_IDTight",
        "HLT_PFMETTypeOne130_PFMHT130_IDTight",
        "HLT_PFMETTypeOne140_PFMHT140_IDTight",
        "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
        "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight",
        "HLT_L1ETMHadSeeds",
        "HLT_CaloMHT90",
        "HLT_CaloMET80_NotCleaned",
        "HLT_CaloMET90_NotCleaned",
        "HLT_CaloMET100_NotCleaned",
        "HLT_CaloMET110_NotCleaned",
        "HLT_CaloMET250_NotCleaned",
        "HLT_CaloMET70_HBHECleaned",
        "HLT_CaloMET80_HBHECleaned",
        "HLT_CaloMET90_HBHECleaned",
        "HLT_CaloMET100_HBHECleaned",
        "HLT_CaloMET250_HBHECleaned",
        "HLT_CaloMET300_HBHECleaned",
        "HLT_CaloMET350_HBHECleaned",
        "HLT_PFMET200_NotCleaned",
        "HLT_PFMET200_HBHECleaned",
        "HLT_PFMET250_HBHECleaned",
        "HLT_PFMET300_HBHECleaned",
        "HLT_PFMET200_HBHE_BeamHaloCleaned",
        "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned",
        "HLT_MET105_IsoTrk50",
        "HLT_MET120_IsoTrk50",
        "HLT_SingleJet30_Mu12_SinglePFJet40",
        "HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33",
        "HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33",
        "HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33",
        "HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33",
        "HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33",
        "HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33",
        "HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33",
        "HLT_DoublePFJets40_CaloBTagCSV_p33",
        "HLT_DoublePFJets100_CaloBTagCSV_p33",
        "HLT_DoublePFJets200_CaloBTagCSV_p33",
        "HLT_DoublePFJets350_CaloBTagCSV_p33",
        "HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33",
        "HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33",
        "HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33",
        "HLT_Photon300_NoHE",
        "HLT_Mu8_TrkIsoVVL",
        "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
        "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
        "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ",
        "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu17_TrkIsoVVL",
        "HLT_Mu19_TrkIsoVVL",
        "HLT_BTagMu_AK4DiJet20_Mu5",
        "HLT_BTagMu_AK4DiJet40_Mu5",
        "HLT_BTagMu_AK4DiJet70_Mu5",
        "HLT_BTagMu_AK4DiJet110_Mu5",
        "HLT_BTagMu_AK4DiJet170_Mu5",
        "HLT_BTagMu_AK4Jet300_Mu5",
        "HLT_BTagMu_AK8DiJet170_Mu5",
        "HLT_BTagMu_AK8Jet300_Mu5",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu12_DoublePhoton20",
        "HLT_TriplePhoton_20_20_20_CaloIdLV2",
        "HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL",
        "HLT_TriplePhoton_30_30_10_CaloIdLV2",
        "HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL",
        "HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL",
        "HLT_Photon25",
        "HLT_Photon33",
        "HLT_Photon50",
        "HLT_Photon75",
        "HLT_Photon90",
        "HLT_Photon120",
        "HLT_Photon150",
        "HLT_Photon175",
        "HLT_Photon200",
        "HLT_Photon50_R9Id90_HE10_IsoM",
        "HLT_Photon75_R9Id90_HE10_IsoM",
        "HLT_Photon90_R9Id90_HE10_IsoM",
        "HLT_Photon120_R9Id90_HE10_IsoM",
        "HLT_Photon165_R9Id90_HE10_IsoM",
        "HLT_Photon90_CaloIdL_PFHT700",
        "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
        "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95",
        "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55",
        "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55",
        "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55",
        "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55",
        "HLT_Dimuon0_Jpsi_L1_NoOS",
        "HLT_Dimuon0_Jpsi_NoVertexing_NoOS",
        "HLT_Dimuon0_Jpsi",
        "HLT_Dimuon0_Jpsi_NoVertexing",
        "HLT_Dimuon0_Jpsi_L1_4R_0er1p5R",
        "HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R",
        "HLT_Dimuon0_Jpsi3p5_Muon2",
        "HLT_Dimuon0_Upsilon_L1_4p5",
        "HLT_Dimuon0_Upsilon_L1_5",
        "HLT_Dimuon0_Upsilon_L1_4p5NoOS",
        "HLT_Dimuon0_Upsilon_L1_4p5er2p0",
        "HLT_Dimuon0_Upsilon_L1_4p5er2p0M",
        "HLT_Dimuon0_Upsilon_NoVertexing",
        "HLT_Dimuon0_Upsilon_L1_5M",
        "HLT_Dimuon0_LowMass_L1_0er1p5R",
        "HLT_Dimuon0_LowMass_L1_0er1p5",
        "HLT_Dimuon0_LowMass",
        "HLT_Dimuon0_LowMass_L1_4",
        "HLT_Dimuon0_LowMass_L1_4R",
        "HLT_Dimuon0_LowMass_L1_TM530",
        "HLT_Dimuon0_Upsilon_Muon_L1_TM0",
        "HLT_Dimuon0_Upsilon_Muon_NoL1Mass",
        "HLT_TripleMu_5_3_3_Mass3p8to60_DZ",
        "HLT_TripleMu_10_5_5_DZ",
        "HLT_TripleMu_12_10_5",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1",
        "HLT_DoubleMu3_DZ_PFMET50_PFMHT60",
        "HLT_DoubleMu3_DZ_PFMET70_PFMHT70",
        "HLT_DoubleMu3_DZ_PFMET90_PFMHT90",
        "HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass",
        "HLT_DoubleMu4_Jpsi_Displaced",
        "HLT_DoubleMu4_Jpsi_NoVertexing",
        "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
        "HLT_DoubleMu43NoFiltersNoVtx",
        "HLT_DoubleMu48NoFiltersNoVtx",
        "HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL",
        "HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL",
        "HLT_DoubleMu20_7_Mass0to30_L1_DM4",
        "HLT_DoubleMu20_7_Mass0to30_L1_DM4EG",
        "HLT_HT425",
        "HLT_HT430_DisplacedDijet40_DisplacedTrack",
        "HLT_HT430_DisplacedDijet60_DisplacedTrack",
        "HLT_HT430_DisplacedDijet80_DisplacedTrack",
        "HLT_HT400_DisplacedDijet40_DisplacedTrack",
        "HLT_HT650_DisplacedDijet60_Inclusive",
        "HLT_HT550_DisplacedDijet80_Inclusive",
        "HLT_HT550_DisplacedDijet60_Inclusive",
        "HLT_HT650_DisplacedDijet80_Inclusive",
        "HLT_HT750_DisplacedDijet80_Inclusive",
        "HLT_DiJet110_35_Mjj650_PFMET110",
        "HLT_DiJet110_35_Mjj650_PFMET120",
        "HLT_DiJet110_35_Mjj650_PFMET130",
        "HLT_TripleJet110_35_35_Mjj650_PFMET110",
        "HLT_TripleJet110_35_35_Mjj650_PFMET120",
        "HLT_TripleJet110_35_35_Mjj650_PFMET130",
        "HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg",
        "HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg",
        "HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg",
        "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned",
        "HLT_Ele28_eta2p1_WPTight_Gsf_HT150",
        "HLT_Ele28_HighEta_SC20_Mass55",
        "HLT_DoubleMu20_7_Mass0to30_Photon23",
        "HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5",
        "HLT_Ele15_IsoVVVL_PFHT450_PFMET50",
        "HLT_Ele15_IsoVVVL_PFHT450",
        "HLT_Ele50_IsoVVVL_PFHT450",
        "HLT_Ele15_IsoVVVL_PFHT600",
        "HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60",
        "HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60",
        "HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5",
        "HLT_Mu15_IsoVVVL_PFHT450_PFMET50",
        "HLT_Mu15_IsoVVVL_PFHT450",
        "HLT_Mu50_IsoVVVL_PFHT450",
        "HLT_Mu15_IsoVVVL_PFHT600",
        "HLT_Dimuon10_PsiPrime_Barrel_Seagulls",
        "HLT_Dimuon20_Jpsi_Barrel_Seagulls",
        "HLT_Dimuon10_Upsilon_Barrel_Seagulls",
        "HLT_Dimuon12_Upsilon_eta1p5",
        "HLT_Dimuon14_Phi_Barrel_Seagulls",
        "HLT_Dimuon18_PsiPrime",
        "HLT_Dimuon25_Jpsi",
        "HLT_Dimuon18_PsiPrime_noCorrL1",
        "HLT_Dimuon24_Upsilon_noCorrL1",
        "HLT_Dimuon24_Phi_noCorrL1",
        "HLT_Dimuon25_Jpsi_noCorrL1",
        "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
        "HLT_DiMu9_Ele9_CaloIdL_TrackIdL",
        "HLT_DoubleIsoMu20_eta2p1",
        "HLT_DoubleIsoMu24_eta2p1",
        "HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx",
        "HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx",
        "HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx",
        "HLT_Mu8",
        "HLT_Mu17",
        "HLT_Mu19",
        "HLT_Mu17_Photon30_IsoCaloId",
        "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30",
        "HLT_Ele8_CaloIdM_TrackIdM_PFJet30",
        "HLT_Ele17_CaloIdM_TrackIdM_PFJet30",
        "HLT_Ele23_CaloIdM_TrackIdM_PFJet30",
        "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165",
        "HLT_Ele115_CaloIdVT_GsfTrkIdT",
        "HLT_Ele135_CaloIdVT_GsfTrkIdT",
        "HLT_Ele145_CaloIdVT_GsfTrkIdT",
        "HLT_Ele200_CaloIdVT_GsfTrkIdT",
        "HLT_Ele250_CaloIdVT_GsfTrkIdT",
        "HLT_Ele300_CaloIdVT_GsfTrkIdT",
        "HLT_PFHT300PT30_QuadPFJet_75_60_45_40",
        "HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0",
        "HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2",
        "HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2",
        "HLT_PFHT380_SixPFJet32",
        "HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5",
        "HLT_PFHT430_SixPFJet40",
        "HLT_PFHT350",
        "HLT_PFHT350MinPFJet15",
        "HLT_Photon60_R9Id90_CaloIdL_IsoL",
        "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL",
        "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15",
        "HLT_FullTrack_Multiplicity85",
        "HLT_FullTrack_Multiplicity100",
        "HLT_FullTrack_Multiplicity130",
        "HLT_FullTrack_Multiplicity155",
        "HLT_ECALHT800",
        "HLT_DiSC30_18_EIso_AND_HE_Mass70",
        "HLT_Physics",
        "HLT_Physics_part0",
        "HLT_Physics_part1",
        "HLT_Physics_part2",
        "HLT_Physics_part3",
        "HLT_Physics_part4",
        "HLT_Physics_part5",
        "HLT_Physics_part6",
        "HLT_Physics_part7",
        "HLT_Random",
        "HLT_ZeroBias",
        "HLT_ZeroBias_part0",
        "HLT_ZeroBias_part1",
        "HLT_ZeroBias_part2",
        "HLT_ZeroBias_part3",
        "HLT_ZeroBias_part4",
        "HLT_ZeroBias_part5",
        "HLT_ZeroBias_part6",
        "HLT_ZeroBias_part7",
        "HLT_AK4CaloJet30",
        "HLT_AK4CaloJet40",
        "HLT_AK4CaloJet50",
        "HLT_AK4CaloJet80",
        "HLT_AK4CaloJet100",
        "HLT_AK4CaloJet120",
        "HLT_AK4PFJet30",
        "HLT_AK4PFJet50",
        "HLT_AK4PFJet80",
        "HLT_AK4PFJet100",
        "HLT_AK4PFJet120",
        "HLT_HISinglePhoton10_Eta3p1ForPPRef",
        "HLT_HISinglePhoton20_Eta3p1ForPPRef",
        "HLT_HISinglePhoton30_Eta3p1ForPPRef",
        "HLT_HISinglePhoton40_Eta3p1ForPPRef",
        "HLT_HISinglePhoton50_Eta3p1ForPPRef",
        "HLT_HISinglePhoton60_Eta3p1ForPPRef",
        "HLT_Photon20_HoverELoose",
        "HLT_Photon30_HoverELoose",
        "HLT_Photon40_HoverELoose",
        "HLT_Photon50_HoverELoose",
        "HLT_Photon60_HoverELoose",
        "HLT_EcalCalibration",
        "HLT_HcalCalibration",
        "HLT_L1UnpairedBunchBptxMinus",
        "HLT_L1UnpairedBunchBptxPlus",
        "HLT_L1NotBptxOR",
        "HLT_L1MinimumBiasHF_OR",
        "HLT_L1MinimumBiasHF0OR",
        "HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142",
        "HLT_HcalNZS",
        "HLT_HcalPhiSym",
        "HLT_HcalIsolatedbunch",
        "HLT_IsoTrackHB",
        "HLT_IsoTrackHE",
        "HLT_ZeroBias_FirstCollisionAfterAbortGap",
        "HLT_ZeroBias_IsolatedBunches",
        "HLT_ZeroBias_FirstCollisionInTrain",
        "HLT_ZeroBias_LastCollisionInTrain",
        "HLT_ZeroBias_FirstBXAfterTrain",
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1",
        "HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg",
        "HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg",
        "HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg",
        "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg",
        "HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
        "HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
        "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
        "HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr",
        "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr",
        "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1",
        "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_CrossL1",
        "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
        "HLT_Rsq0p35",
        "HLT_Rsq0p40",
        "HLT_RsqMR300_Rsq0p09_MR200",
        "HLT_RsqMR320_Rsq0p09_MR200",
        "HLT_RsqMR300_Rsq0p09_MR200_4jet",
        "HLT_RsqMR320_Rsq0p09_MR200_4jet",
        "HLT_L1_DoubleJet30_Mass_Min400_Mu10",
        "HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1",
        "HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1",
        "HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1",
        "HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50",
        "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3",
        "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3",
        "HLT_PFMET100_PFMHT100_IDTight_PFHT60",
        "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60",
        "HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60",
        "HLT_Mu18_Mu9_SameSign",
        "HLT_Mu18_Mu9_SameSign_DZ",
        "HLT_Mu18_Mu9",
        "HLT_Mu18_Mu9_DZ",
        "HLT_Mu20_Mu10_SameSign",
        "HLT_Mu20_Mu10_SameSign_DZ",
        "HLT_Mu20_Mu10",
        "HLT_Mu20_Mu10_DZ",
        "HLT_Mu23_Mu12_SameSign",
        "HLT_Mu23_Mu12_SameSign_DZ",
        "HLT_Mu23_Mu12",
        "HLT_Mu23_Mu12_DZ",
        "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi",
        "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi",
        "HLT_DoubleMu3_DCA_PFMET50_PFMHT60",
        "HLT_TripleMu_5_3_3_Mass3p8to60_DCA",
        "HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1",
        "HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1",
        "HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1",
        "HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1",
        "HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2",
        "HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2",
        "HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2",
        "HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2",
        "HLT_QuadPFJet98_83_71_15",
        "HLT_QuadPFJet103_88_75_15",
        "HLT_QuadPFJet105_88_76_15",
        "HLT_QuadPFJet111_90_80_15",
        "HLT_AK8PFJet330_PFAK8BTagCSV_p17",
        "HLT_AK8PFJet330_PFAK8BTagCSV_p1",
        "HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55",
        "HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55"
    ],

    Era_2018: [
        "HLT_AK8PFJet360_TrimMass30",
        "HLT_AK8PFJet380_TrimMass30",
        "HLT_AK8PFJet400_TrimMass30",
        "HLT_AK8PFJet420_TrimMass30",
        "HLT_AK8PFHT750_TrimMass50",
        "HLT_AK8PFHT800_TrimMass50",
        "HLT_AK8PFHT850_TrimMass50",
        "HLT_AK8PFHT900_TrimMass50",
        "HLT_CaloJet500_NoJetID",
        "HLT_CaloJet550_NoJetID",
        "HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL",
        "HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon",
        "HLT_Trimuon5_3p5_2_Upsilon_Muon",
        "HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon",
        "HLT_DoubleEle25_CaloIdL_MW",
        "HLT_DoubleEle27_CaloIdL_MW",
        "HLT_DoubleEle33_CaloIdL_MW",
        "HLT_DoubleEle24_eta2p1_WPTight_Gsf",
        "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350",
        "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350",
        "HLT_Ele27_Ele37_CaloIdL_MW",
        "HLT_Mu27_Ele37_CaloIdL_MW",
        "HLT_Mu37_Ele27_CaloIdL_MW",
        "HLT_Mu37_TkMu27",
        "HLT_DoubleMu4_3_Bs",
        "HLT_DoubleMu4_3_Jpsi",
        "HLT_DoubleMu4_JpsiTrk_Displaced",
        "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced",
        "HLT_DoubleMu3_Trk_Tau3mu",
        "HLT_DoubleMu3_TkMu_DsTau3Mu",
        "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
        "HLT_DoubleMu4_Mass3p8_DZ_PFHT350",
        "HLT_Mu3_PFJet40",
        "HLT_Mu7p5_L2Mu2_Jpsi",
        "HLT_Mu7p5_L2Mu2_Upsilon",
        "HLT_Mu7p5_Track2_Jpsi",
        "HLT_Mu7p5_Track3p5_Jpsi",
        "HLT_Mu7p5_Track7_Jpsi",
        "HLT_Mu7p5_Track2_Upsilon",
        "HLT_Mu7p5_Track3p5_Upsilon",
        "HLT_Mu7p5_Track7_Upsilon",
        "HLT_Mu3_L1SingleMu5orSingleMu7",
        "HLT_DoublePhoton33_CaloIdL",
        "HLT_DoublePhoton70",
        "HLT_DoublePhoton85",
        "HLT_Ele20_WPTight_Gsf",
        "HLT_Ele15_WPLoose_Gsf",
        "HLT_Ele17_WPLoose_Gsf",
        "HLT_Ele20_WPLoose_Gsf",
        "HLT_Ele20_eta2p1_WPLoose_Gsf",
        "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG",
        "HLT_Ele27_WPTight_Gsf",
        "HLT_Ele28_WPTight_Gsf",
        "HLT_Ele30_WPTight_Gsf",
        "HLT_Ele32_WPTight_Gsf",
        "HLT_Ele35_WPTight_Gsf",
        "HLT_Ele35_WPTight_Gsf_L1EGMT",
        "HLT_Ele38_WPTight_Gsf",
        "HLT_Ele40_WPTight_Gsf",
        "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1",
        "HLT_HT450_Beamspot",
        "HLT_HT300_Beamspot",
        "HLT_ZeroBias_Beamspot",
        "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1",
        "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1",
        "HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1",
        "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1",
        "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1",
        "HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1",
        "HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1",
        "HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1",
        "HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1",
        "HLT_IsoMu20",
        "HLT_IsoMu24",
        "HLT_IsoMu24_eta2p1",
        "HLT_IsoMu27",
        "HLT_IsoMu30",
        "HLT_UncorrectedJetE30_NoBPTX",
        "HLT_UncorrectedJetE30_NoBPTX3BX",
        "HLT_UncorrectedJetE60_NoBPTX3BX",
        "HLT_UncorrectedJetE70_NoBPTX3BX",
        "HLT_L1SingleMu18",
        "HLT_L1SingleMu25",
        "HLT_L2Mu10",
        "HLT_L2Mu10_NoVertex_NoBPTX3BX",
        "HLT_L2Mu10_NoVertex_NoBPTX",
        "HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX",
        "HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX",
        "HLT_L2Mu50",
        "HLT_L2Mu23NoVtx_2Cha",
        "HLT_L2Mu23NoVtx_2Cha_CosmicSeed",
        "HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4",
        "HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4",
        "HLT_DoubleL2Mu50",
        "HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed",
        "HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched",
        "HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed",
        "HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched",
        "HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4",
        "HLT_DoubleL2Mu23NoVtx_2Cha",
        "HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched",
        "HLT_DoubleL2Mu25NoVtx_2Cha",
        "HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched",
        "HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8",
        "HLT_Mu25_TkMu0_Onia",
        "HLT_Mu30_TkMu0_Psi",
        "HLT_Mu30_TkMu0_Upsilon",
        "HLT_Mu20_TkMu0_Phi",
        "HLT_Mu25_TkMu0_Phi",
        "HLT_Mu12",
        "HLT_Mu15",
        "HLT_Mu20",
        "HLT_Mu27",
        "HLT_Mu50",
        "HLT_Mu55",
        "HLT_OldMu100",
        "HLT_TkMu100",
        "HLT_DiPFJetAve40",
        "HLT_DiPFJetAve60",
        "HLT_DiPFJetAve80",
        "HLT_DiPFJetAve140",
        "HLT_DiPFJetAve200",
        "HLT_DiPFJetAve260",
        "HLT_DiPFJetAve320",
        "HLT_DiPFJetAve400",
        "HLT_DiPFJetAve500",
        "HLT_DiPFJetAve60_HFJEC",
        "HLT_DiPFJetAve80_HFJEC",
        "HLT_DiPFJetAve100_HFJEC",
        "HLT_DiPFJetAve160_HFJEC",
        "HLT_DiPFJetAve220_HFJEC",
        "HLT_DiPFJetAve300_HFJEC",
        "HLT_AK8PFJet15",
        "HLT_AK8PFJet25",
        "HLT_AK8PFJet40",
        "HLT_AK8PFJet60",
        "HLT_AK8PFJet80",
        "HLT_AK8PFJet140",
        "HLT_AK8PFJet200",
        "HLT_AK8PFJet260",
        "HLT_AK8PFJet320",
        "HLT_AK8PFJet400",
        "HLT_AK8PFJet450",
        "HLT_AK8PFJet500",
        "HLT_AK8PFJet550",
        "HLT_PFJet15",
        "HLT_PFJet25",
        "HLT_PFJet40",
        "HLT_PFJet60",
        "HLT_PFJet80",
        "HLT_PFJet140",
        "HLT_PFJet200",
        "HLT_PFJet260",
        "HLT_PFJet320",
        "HLT_PFJet400",
        "HLT_PFJet450",
        "HLT_PFJet500",
        "HLT_PFJet550",
        "HLT_PFJetFwd15",
        "HLT_PFJetFwd25",
        "HLT_PFJetFwd40",
        "HLT_PFJetFwd60",
        "HLT_PFJetFwd80",
        "HLT_PFJetFwd140",
        "HLT_PFJetFwd200",
        "HLT_PFJetFwd260",
        "HLT_PFJetFwd320",
        "HLT_PFJetFwd400",
        "HLT_PFJetFwd450",
        "HLT_PFJetFwd500",
        "HLT_AK8PFJetFwd15",
        "HLT_AK8PFJetFwd25",
        "HLT_AK8PFJetFwd40",
        "HLT_AK8PFJetFwd60",
        "HLT_AK8PFJetFwd80",
        "HLT_AK8PFJetFwd140",
        "HLT_AK8PFJetFwd200",
        "HLT_AK8PFJetFwd260",
        "HLT_AK8PFJetFwd320",
        "HLT_AK8PFJetFwd400",
        "HLT_AK8PFJetFwd450",
        "HLT_AK8PFJetFwd500",
        "HLT_PFHT180",
        "HLT_PFHT250",
        "HLT_PFHT370",
        "HLT_PFHT430",
        "HLT_PFHT510",
        "HLT_PFHT590",
        "HLT_PFHT680",
        "HLT_PFHT780",
        "HLT_PFHT890",
        "HLT_PFHT1050",
        "HLT_PFHT500_PFMET100_PFMHT100_IDTight",
        "HLT_PFHT500_PFMET110_PFMHT110_IDTight",
        "HLT_PFHT700_PFMET85_PFMHT85_IDTight",
        "HLT_PFHT700_PFMET95_PFMHT95_IDTight",
        "HLT_PFHT800_PFMET75_PFMHT75_IDTight",
        "HLT_PFHT800_PFMET85_PFMHT85_IDTight",
        "HLT_PFMET110_PFMHT110_IDTight",
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_PFMET130_PFMHT130_IDTight",
        "HLT_PFMET140_PFMHT140_IDTight",
        "HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1",
        "HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1",
        "HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1",
        "HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1",
        "HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1",
        "HLT_PFMET120_PFMHT120_IDTight_PFHT60",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
        "HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60",
        "HLT_PFMETTypeOne110_PFMHT110_IDTight",
        "HLT_PFMETTypeOne120_PFMHT120_IDTight",
        "HLT_PFMETTypeOne130_PFMHT130_IDTight",
        "HLT_PFMETTypeOne140_PFMHT140_IDTight",
        "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
        "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight",
        "HLT_L1ETMHadSeeds",
        "HLT_CaloMHT90",
        "HLT_CaloMET80_NotCleaned",
        "HLT_CaloMET90_NotCleaned",
        "HLT_CaloMET100_NotCleaned",
        "HLT_CaloMET110_NotCleaned",
        "HLT_CaloMET250_NotCleaned",
        "HLT_CaloMET70_HBHECleaned",
        "HLT_CaloMET80_HBHECleaned",
        "HLT_CaloMET90_HBHECleaned",
        "HLT_CaloMET100_HBHECleaned",
        "HLT_CaloMET250_HBHECleaned",
        "HLT_CaloMET300_HBHECleaned",
        "HLT_CaloMET350_HBHECleaned",
        "HLT_PFMET200_NotCleaned",
        "HLT_PFMET200_HBHECleaned",
        "HLT_PFMET250_HBHECleaned",
        "HLT_PFMET300_HBHECleaned",
        "HLT_PFMET200_HBHE_BeamHaloCleaned",
        "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned",
        "HLT_MET105_IsoTrk50",
        "HLT_MET120_IsoTrk50",
        "HLT_SingleJet30_Mu12_SinglePFJet40",
        "HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71",
        "HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71",
        "HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71",
        "HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71",
        "HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",
        "HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",
        "HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",
        "HLT_DoublePFJets40_CaloBTagDeepCSV_p71",
        "HLT_DoublePFJets100_CaloBTagDeepCSV_p71",
        "HLT_DoublePFJets200_CaloBTagDeepCSV_p71",
        "HLT_DoublePFJets350_CaloBTagDeepCSV_p71",
        "HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",
        "HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",
        "HLT_Photon300_NoHE",
        "HLT_Mu8_TrkIsoVVL",
        "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
        "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
        "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ",
        "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu17_TrkIsoVVL",
        "HLT_Mu19_TrkIsoVVL",
        "HLT_BTagMu_AK4DiJet20_Mu5",
        "HLT_BTagMu_AK4DiJet40_Mu5",
        "HLT_BTagMu_AK4DiJet70_Mu5",
        "HLT_BTagMu_AK4DiJet110_Mu5",
        "HLT_BTagMu_AK4DiJet170_Mu5",
        "HLT_BTagMu_AK4Jet300_Mu5",
        "HLT_BTagMu_AK8DiJet170_Mu5",
        "HLT_BTagMu_AK8Jet170_DoubleMu5",
        "HLT_BTagMu_AK8Jet300_Mu5",
        "HLT_BTagMu_AK4DiJet20_Mu5_noalgo",
        "HLT_BTagMu_AK4DiJet40_Mu5_noalgo",
        "HLT_BTagMu_AK4DiJet70_Mu5_noalgo",
        "HLT_BTagMu_AK4DiJet110_Mu5_noalgo",
        "HLT_BTagMu_AK4DiJet170_Mu5_noalgo",
        "HLT_BTagMu_AK4Jet300_Mu5_noalgo",
        "HLT_BTagMu_AK8DiJet170_Mu5_noalgo",
        "HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo",
        "HLT_BTagMu_AK8Jet300_Mu5_noalgo",
        "HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu12_DoublePhoton20",
        "HLT_TriplePhoton_20_20_20_CaloIdLV2",
        "HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL",
        "HLT_TriplePhoton_30_30_10_CaloIdLV2",
        "HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL",
        "HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL",
        "HLT_Photon20",
        "HLT_Photon33",
        "HLT_Photon50",
        "HLT_Photon75",
        "HLT_Photon90",
        "HLT_Photon120",
        "HLT_Photon150",
        "HLT_Photon175",
        "HLT_Photon200",
        "HLT_Photon100EB_TightID_TightIso",
        "HLT_Photon110EB_TightID_TightIso",
        "HLT_Photon120EB_TightID_TightIso",
        "HLT_Photon100EBHE10",
        "HLT_Photon100EEHE10",
        "HLT_Photon100EE_TightID_TightIso",
        "HLT_Photon50_R9Id90_HE10_IsoM",
        "HLT_Photon75_R9Id90_HE10_IsoM",
        "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3",
        "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3",
        "HLT_Photon90_R9Id90_HE10_IsoM",
        "HLT_Photon120_R9Id90_HE10_IsoM",
        "HLT_Photon165_R9Id90_HE10_IsoM",
        "HLT_Photon90_CaloIdL_PFHT700",
        "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
        "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95",
        "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55",
        "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55",
        "HLT_Photon35_TwoProngs35",
        "HLT_IsoMu24_TwoProngs35",
        "HLT_Dimuon0_Jpsi_L1_NoOS",
        "HLT_Dimuon0_Jpsi_NoVertexing_NoOS",
        "HLT_Dimuon0_Jpsi",
        "HLT_Dimuon0_Jpsi_NoVertexing",
        "HLT_Dimuon0_Jpsi_L1_4R_0er1p5R",
        "HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R",
        "HLT_Dimuon0_Jpsi3p5_Muon2",
        "HLT_Dimuon0_Upsilon_L1_4p5",
        "HLT_Dimuon0_Upsilon_L1_5",
        "HLT_Dimuon0_Upsilon_L1_4p5NoOS",
        "HLT_Dimuon0_Upsilon_L1_4p5er2p0",
        "HLT_Dimuon0_Upsilon_L1_4p5er2p0M",
        "HLT_Dimuon0_Upsilon_NoVertexing",
        "HLT_Dimuon0_Upsilon_L1_5M",
        "HLT_Dimuon0_LowMass_L1_0er1p5R",
        "HLT_Dimuon0_LowMass_L1_0er1p5",
        "HLT_Dimuon0_LowMass",
        "HLT_Dimuon0_LowMass_L1_4",
        "HLT_Dimuon0_LowMass_L1_4R",
        "HLT_Dimuon0_LowMass_L1_TM530",
        "HLT_Dimuon0_Upsilon_Muon_L1_TM0",
        "HLT_Dimuon0_Upsilon_Muon_NoL1Mass",
        "HLT_TripleMu_5_3_3_Mass3p8_DZ",
        "HLT_TripleMu_10_5_5_DZ",
        "HLT_TripleMu_12_10_5",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15",
        "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1",
        "HLT_DoubleMu3_DZ_PFMET50_PFMHT60",
        "HLT_DoubleMu3_DZ_PFMET70_PFMHT70",
        "HLT_DoubleMu3_DZ_PFMET90_PFMHT90",
        "HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass",
        "HLT_DoubleMu4_Jpsi_Displaced",
        "HLT_DoubleMu4_Jpsi_NoVertexing",
        "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
        "HLT_DoubleMu43NoFiltersNoVtx",
        "HLT_DoubleMu48NoFiltersNoVtx",
        "HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL",
        "HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL",
        "HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL",
        "HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL",
        "HLT_DoubleMu33NoFiltersNoVtxDisplaced",
        "HLT_DoubleMu40NoFiltersNoVtxDisplaced",
        "HLT_DoubleMu20_7_Mass0to30_L1_DM4",
        "HLT_DoubleMu20_7_Mass0to30_L1_DM4EG",
        "HLT_HT425",
        "HLT_HT430_DisplacedDijet40_DisplacedTrack",
        "HLT_HT500_DisplacedDijet40_DisplacedTrack",
        "HLT_HT430_DisplacedDijet60_DisplacedTrack",
        "HLT_HT400_DisplacedDijet40_DisplacedTrack",
        "HLT_HT650_DisplacedDijet60_Inclusive",
        "HLT_HT550_DisplacedDijet60_Inclusive",
        "HLT_DiJet110_35_Mjj650_PFMET110",
        "HLT_DiJet110_35_Mjj650_PFMET120",
        "HLT_DiJet110_35_Mjj650_PFMET130",
        "HLT_TripleJet110_35_35_Mjj650_PFMET110",
        "HLT_TripleJet110_35_35_Mjj650_PFMET120",
        "HLT_TripleJet110_35_35_Mjj650_PFMET130",
        "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned",
        "HLT_Ele28_eta2p1_WPTight_Gsf_HT150",
        "HLT_Ele28_HighEta_SC20_Mass55",
        "HLT_DoubleMu20_7_Mass0to30_Photon23",
        "HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5",
        "HLT_Ele15_IsoVVVL_PFHT450_PFMET50",
        "HLT_Ele15_IsoVVVL_PFHT450",
        "HLT_Ele50_IsoVVVL_PFHT450",
        "HLT_Ele15_IsoVVVL_PFHT600",
        "HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60",
        "HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60",
        "HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60",
        "HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5",
        "HLT_Mu15_IsoVVVL_PFHT450_PFMET50",
        "HLT_Mu15_IsoVVVL_PFHT450",
        "HLT_Mu50_IsoVVVL_PFHT450",
        "HLT_Mu15_IsoVVVL_PFHT600",
        "HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight",
        "HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight",
        "HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight",
        "HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight",
        "HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight",
        "HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight",
        "HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight",
        "HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight",
        "HLT_Dimuon10_PsiPrime_Barrel_Seagulls",
        "HLT_Dimuon20_Jpsi_Barrel_Seagulls",
        "HLT_Dimuon12_Upsilon_y1p4",
        "HLT_Dimuon14_Phi_Barrel_Seagulls",
        "HLT_Dimuon18_PsiPrime",
        "HLT_Dimuon25_Jpsi",
        "HLT_Dimuon18_PsiPrime_noCorrL1",
        "HLT_Dimuon24_Upsilon_noCorrL1",
        "HLT_Dimuon24_Phi_noCorrL1",
        "HLT_Dimuon25_Jpsi_noCorrL1",
        "HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8",
        "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
        "HLT_DiMu9_Ele9_CaloIdL_TrackIdL",
        "HLT_DoubleIsoMu20_eta2p1",
        "HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx",
        "HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx",
        "HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx",
        "HLT_Mu8",
        "HLT_Mu17",
        "HLT_Mu19",
        "HLT_Mu17_Photon30_IsoCaloId",
        "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
        "HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30",
        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30",
        "HLT_Ele8_CaloIdM_TrackIdM_PFJet30",
        "HLT_Ele17_CaloIdM_TrackIdM_PFJet30",
        "HLT_Ele23_CaloIdM_TrackIdM_PFJet30",
        "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165",
        "HLT_Ele115_CaloIdVT_GsfTrkIdT",
        "HLT_Ele135_CaloIdVT_GsfTrkIdT",
        "HLT_Ele145_CaloIdVT_GsfTrkIdT",
        "HLT_Ele200_CaloIdVT_GsfTrkIdT",
        "HLT_Ele250_CaloIdVT_GsfTrkIdT",
        "HLT_Ele300_CaloIdVT_GsfTrkIdT",
        "HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5",
        "HLT_PFHT330PT30_QuadPFJet_75_60_45_40",
        "HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94",
        "HLT_PFHT400_SixPFJet32",
        "HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59",
        "HLT_PFHT450_SixPFJet36",
        "HLT_PFHT350",
        "HLT_PFHT350MinPFJet15",
        "HLT_Photon60_R9Id90_CaloIdL_IsoL",
        "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL",
        "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15",
        "HLT_ECALHT800",
        "HLT_DiSC30_18_EIso_AND_HE_Mass70",
        "HLT_Physics",
        "HLT_Physics_part0",
        "HLT_Physics_part1",
        "HLT_Physics_part2",
        "HLT_Physics_part3",
        "HLT_Physics_part4",
        "HLT_Physics_part5",
        "HLT_Physics_part6",
        "HLT_Physics_part7",
        "HLT_Random",
        "HLT_ZeroBias",
        "HLT_ZeroBias_Alignment",
        "HLT_ZeroBias_part0",
        "HLT_ZeroBias_part1",
        "HLT_ZeroBias_part2",
        "HLT_ZeroBias_part3",
        "HLT_ZeroBias_part4",
        "HLT_ZeroBias_part5",
        "HLT_ZeroBias_part6",
        "HLT_ZeroBias_part7",
        "HLT_AK4CaloJet30",
        "HLT_AK4CaloJet40",
        "HLT_AK4CaloJet50",
        "HLT_AK4CaloJet80",
        "HLT_AK4CaloJet100",
        "HLT_AK4CaloJet120",
        "HLT_AK4PFJet30",
        "HLT_AK4PFJet50",
        "HLT_AK4PFJet80",
        "HLT_AK4PFJet100",
        "HLT_AK4PFJet120",
        "HLT_SinglePhoton10_Eta3p1ForPPRef",
        "HLT_SinglePhoton20_Eta3p1ForPPRef",
        "HLT_SinglePhoton30_Eta3p1ForPPRef",
        "HLT_Photon20_HoverELoose",
        "HLT_Photon30_HoverELoose",
        "HLT_EcalCalibration",
        "HLT_HcalCalibration",
        "HLT_L1UnpairedBunchBptxMinus",
        "HLT_L1UnpairedBunchBptxPlus",
        "HLT_L1NotBptxOR",
        "HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142",
        "HLT_CDC_L2cosmic_5_er1p0",
        "HLT_CDC_L2cosmic_5p5_er1p0",
        "HLT_HcalNZS",
        "HLT_HcalPhiSym",
        "HLT_HcalIsolatedbunch",
        "HLT_IsoTrackHB",
        "HLT_IsoTrackHE",
        "HLT_ZeroBias_FirstCollisionAfterAbortGap",
        "HLT_ZeroBias_IsolatedBunches",
        "HLT_ZeroBias_FirstCollisionInTrain",
        "HLT_ZeroBias_LastCollisionInTrain",
        "HLT_ZeroBias_FirstBXAfterTrain",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr",
        "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr",
        "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1",
        "HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1",
        "HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1",
        "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
        "HLT_Rsq0p35",
        "HLT_Rsq0p40",
        "HLT_RsqMR300_Rsq0p09_MR200",
        "HLT_RsqMR320_Rsq0p09_MR200",
        "HLT_RsqMR300_Rsq0p09_MR200_4jet",
        "HLT_RsqMR320_Rsq0p09_MR200_4jet",
        "HLT_IsoMu27_MET90",
        "HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg",
        "HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg",
        "HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg",
        "HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg",
        "HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg",
        "HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1",
        "HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1",
        "HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1",
        "HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50",
        "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3",
        "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3",
        "HLT_PFMET100_PFMHT100_IDTight_PFHT60",
        "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60",
        "HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60",
        "HLT_Mu18_Mu9_SameSign",
        "HLT_Mu18_Mu9_SameSign_DZ",
        "HLT_Mu18_Mu9",
        "HLT_Mu18_Mu9_DZ",
        "HLT_Mu20_Mu10_SameSign",
        "HLT_Mu20_Mu10_SameSign_DZ",
        "HLT_Mu20_Mu10",
        "HLT_Mu20_Mu10_DZ",
        "HLT_Mu23_Mu12_SameSign",
        "HLT_Mu23_Mu12_SameSign_DZ",
        "HLT_Mu23_Mu12",
        "HLT_Mu23_Mu12_DZ",
        "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05",
        "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi",
        "HLT_DoubleMu3_DCA_PFMET50_PFMHT60",
        "HLT_TripleMu_5_3_3_Mass3p8_DCA",
        "HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1",
        "HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1",
        "HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1",
        "HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2",
        "HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2",
        "HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2",
        "HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2",
        "HLT_QuadPFJet98_83_71_15",
        "HLT_QuadPFJet103_88_75_15",
        "HLT_QuadPFJet105_88_76_15",
        "HLT_QuadPFJet111_90_80_15",
        "HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17",
        "HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1",
        "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02",
        "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2",
        "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4",
        "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55",
        "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto",
        "HLT_Mu12_IP6_part0",
        "HLT_Mu12_IP6_part1",
        "HLT_Mu12_IP6_part2",
        "HLT_Mu12_IP6_part3",
        "HLT_Mu12_IP6_part4",
        "HLT_Mu9_IP5_part0",
        "HLT_Mu9_IP5_part1",
        "HLT_Mu9_IP5_part2",
        "HLT_Mu9_IP5_part3",
        "HLT_Mu9_IP5_part4",
        "HLT_Mu7_IP4_part0",
        "HLT_Mu7_IP4_part1",
        "HLT_Mu7_IP4_part2",
        "HLT_Mu7_IP4_part3",
        "HLT_Mu7_IP4_part4",
        "HLT_Mu9_IP4_part0",
        "HLT_Mu9_IP4_part1",
        "HLT_Mu9_IP4_part2",
        "HLT_Mu9_IP4_part3",
        "HLT_Mu9_IP4_part4",
        "HLT_Mu8_IP5_part0",
        "HLT_Mu8_IP5_part1",
        "HLT_Mu8_IP5_part2",
        "HLT_Mu8_IP5_part3",
        "HLT_Mu8_IP5_part4",
        "HLT_Mu8_IP6_part0",
        "HLT_Mu8_IP6_part1",
        "HLT_Mu8_IP6_part2",
        "HLT_Mu8_IP6_part3",
        "HLT_Mu8_IP6_part4",
        "HLT_Mu9_IP6_part0",
        "HLT_Mu9_IP6_part1",
        "HLT_Mu9_IP6_part2",
        "HLT_Mu9_IP6_part3",
        "HLT_Mu9_IP6_part4",
        "HLT_Mu8_IP3_part0",
        "HLT_Mu8_IP3_part1",
        "HLT_Mu8_IP3_part2",
        "HLT_Mu8_IP3_part3",
        "HLT_Mu8_IP3_part4",
        "HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1",
        "HLT_TrkMu6NoFiltersNoVtx",
        "HLT_TrkMu16NoFiltersNoVtx",
        "HLT_DoubleTrkMu_16_6_NoFiltersNoVtx"
    ],
}
HLTTriggers[Era_2016preVFP] = HLTTriggers[Era_2016]
HLTTriggers[Era_2016postVFP]= HLTTriggers[Era_2016]


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
        self.FatJetPNetXto4bv2WorkingPoints = ['60'] #['40', '60', '80'] #['40', '45a', '45b', '50', '60', '65', '70', '80']  # ['40', '50', '60', '65', '70', '80']   ['40', '60', '80']
        
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
        self.MassVBFjj_MinThsh_Loose =   400
        self.DEtaVBFjj_MinThsh_Loose =   2.2         


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
        self.datasetInfo['isSignalZH'     ]  = False
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
            self.datasetInfo['isSignalZH']       = True if "SUSY_ZH_ZToAll_HToAATo4B"   in datasetName_part1 else False
            self.datasetInfo['isSignalTTH']      = True if "SUSY_TTH_TTToAll_HToAATo4B" in datasetName_part1 else False
            self.datasetInfo['isSignal']         = (self.datasetInfo['isSignalGGH']   or \
                                                     self.datasetInfo['isSignalVBFH'] or \
                                                     self.datasetInfo['isSignalWH']   or \
                                                     self.datasetInfo['isSignalZH']   or \
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
            print(f"{self.datasetInfo['isSignalZH'     ] = }")
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
        #sTrgSelection = 'Trg_Combo_AK4AK8Jet_HT'
        sTrgSelection = self.datasetInfo['triggers'] if self.datasetInfo['triggers'] else 'Trg_Combo_AK4AK8Jet_HT'
        print(f"{sTrgSelection = }")
        
        # sel_names_all = dict of {"selection name" : [list of different cuts]}; for cut-flow table 
        self.sel_names_all = OD([
            ("Presel",                    [
                "nPV",
                "METFilters",
                "candH",
                "candHnBHadronsGe4",
                "leadingFatJetPt",
                
            ]),
        ])

        


        if not self.datasetInfo['isMC']: 
            self.sel_names_all["Presel"].insert(0, "run:ls")

               




        
        categories_dict = OD()
        if self.datasetInfo['isSignalGGH']:
            categories_dict["gg0l"] = self.sel_names_all["Presel"] + [
                "nLeptonsTightZero",
            ]
        if self.datasetInfo['isSignalZH']:
            categories_dict["Zvv"]  = self.sel_names_all["Presel"] + [
                "nLeptonsTightZero", 
                "METPt",
                "dPhiLeadingFJHto4bAndMet",
            ]
        if self.datasetInfo['isSignalVBFH']:
            categories_dict["VBF"] = self.sel_names_all["Presel"] + [
                "nLeptonsTightZero",
                "jjvbf_Loose"
            ]
        

        for sCatName, catSels in categories_dict.items():            
            for wp_ in self.objectSelector.FatJetPNetXto4bv2WorkingPoints: 
                '''
                self.sel_names_all["%s_Xto4bv2_SRWP%s" % (sCatName, wp_)] = catSels + [ # signal region
                    "leadingFatJetPNet_Xto4bv2_Htoaa4b_SRWP%s" % (wp_)
                ]
                self.sel_names_all["%s_Xto4bv2_SBWP%s" % (sCatName, wp_)] = catSels + [ # side band
                    "leadingFatJetPNet_Xto4bv2_Htoaa4b_SBWP%s" % (wp_)
                ]
                '''
                self.sel_names_all["%s_SBplusSRWP%s" % (sCatName, wp_)] = catSels + [ # side band + signal region
                    "leadingFatJetPNet_Xto4bv2_Htoaa4b_SBplusSRWP%s" % (wp_)
                ]
                 
                
        ## remove "Presel" category
        self.sel_names_all.pop("Presel", None)
    

        self.sel_conditions_all_list = set()
        for sel_conditions_ in self.sel_names_all.values():
            self.sel_conditions_all_list.update( sel_conditions_ )
        print(f"{self.sel_conditions_all_list = }")
        

        


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
                logging.critical(f'htoaa_Analysis_GGFMode.py::main():: {sFilesGoldenJSON[self.datasetInfo["era"]] = } could not read.')
                exit(0) 

            # convert runNumber in str to int
            dataLSSelGoldenJSON = {int(k): v for k, v in dataLSSelGoldenJSON.items()} 
            self.datasetInfo['dataLSSelGoldenJSON'] = dataLSSelGoldenJSON
            #print(f"{dataLSSelGoldenJSON = }")

        else: ## MC
        
            pass
            '''
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
            
            print(f"{bTagSFEfficiencyDict[self.datasetInfo['era']]['inputFile'] = }")
            '''

        
        
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
        pt_axis               = hist.Bin("Pt",                     r"$p_{T}$ [GeV]",             50,       0,    1000)
        ptLow_axis            = hist.Bin("PtLow",                  r"$p_{T}$ [GeV]",             400,       0,     200)
        ptUltraLow_axis       = hist.Bin("PtUltraLow",             r"$p_{T}$ [GeV]",             200,       0,     0.1)
        pt1to10_axis          = hist.Bin("Pt1to10",                r"$p_{T}$ [GeV]",             100,       0,      10)
        pt2TeV_axis           = hist.Bin("Pt2TeV",                 r"$p_{T}$ [GeV]",            2000,       0,    2000)
        log2Pt2TeV_axis       = hist.Bin("Log2Pt2TeV",             r"Log2($p_{T}$) [GeV]",      200,  math.log2(1),    math.log2(2000))
        eta_axis              = hist.Bin("Eta",                    r"$#eta$",                    100,      -6,       6)
        phi_axis              = hist.Bin("Phi",                    r"$\phi$",                    100,   -3.14,    3.13)
        #mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       200, 0, 600)
        #mass_axis             = hist.Bin("Mass",      r"$m$ [GeV]",       400, 0, 200)
        #mass_axis             = hist.Bin("Mass",                   r"$m$ [GeV]",                 350,       0,     350)
        mass_axis             = hist.Bin("Mass",                   r"$m$ [GeV]",                  48,       0,     240)
        massCl1_axis          = hist.Bin("MassCl1",                r"$m$ [GeV]",                 350,       0,     350)
        mass_axis1            = hist.Bin("Mass1",                  r"$m$ [GeV]",               20*70,       0,      70)
        mass_axis2            = hist.Bin("Mass2",                  r"$m$ [GeV]",                  72,       0,      72)
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
        

        
        for sel_name in self.sel_names_all.keys(): # loop of list of selections
            sHExt = "_%s" % (sel_name)





            histos.update(OD([
                ('hLeadingFatJetPt%s'%(sHExt),                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),                        
                ('hMET_pT%s'%(sHExt),                                   {sXaxis: pt_axis,         sXaxisLabel: r"MET pT [GeV]"}),
                ('hAK4JetPtNonoverlapFatJet%s'%(sHExt),                 {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(AK4 jet nonoverlap candH)$ [GeV]"}),                        
                
            ]))

            for hltPathName_ in HLTTriggers[self.datasetInfo["era"]]:
                histos.update(OD([
                    ('hLeadingFatJetPt%s_%s'%(sHExt, hltPathName_),                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),                        
                    ('hMET_pT%s_%s'%(sHExt, hltPathName_),                                   {sXaxis: pt_axis,         sXaxisLabel: r"MET pT [GeV]"}),
                    ('hAK4JetPtNonoverlapFatJet%s_%s'%(sHExt, hltPathName_),                {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(AK4 jet nonoverlap candH)$ [GeV]"}),                        
                
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
                        #dataset_axis,
                        hXaxis, #nObject_axis,
                        #systematic_axis,
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
        if printLevel >= 0:
            print(f"\n events.fields ({type(events.fields)}): {events.fields}"); sys.stdout.flush()
        if printLevel >= 0:
            hltNames_list = ['HLT_'+f for f in events.HLT.fields ]
            print(f"hltNames_list =")
            print(json.dumps(hltNames_list, indent=4))
             
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








        # Reco-level -----------------------------------------------------------------------------------
        # FatJet selection


        

        ################## 
        # EVENT VARIABLES
        ##################


        leadingFatJet = None

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
        
        if 'PNet_X4b_v2a_Haa34b_score' in FatJetsToUse.fields: # NanoAOD v2 PNet_X4b_v1_Haa4b_vs_QCD
            leadingFatJet_PNet_Xto4bv2_Htoaa4b        = (leadingFatJet.PNet_X4b_v2a_Haa4b_score + \
                                                         leadingFatJet.PNet_X4b_v2b_Haa4b_score) / 2.0
            
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

        if "candHnBHadronsGe4" in self.sel_conditions_all_list:
            selection.add(
                "candHnBHadronsGe4",
                (leadingFatJet.nBHadrons >= 4) # no. of candH >= 1
            )

        

        if "leadingFatJetPt" in self.sel_conditions_all_list:
            selection.add(
                "leadingFatJetPt",
                ((leadingFatJet.pt_toUse > self.objectSelector.FatJetPt_gg0lIncl_MinThsh) &
                 (leadingFatJet.pt_toUse <=  self.objectSelector.FatJetPt_gg0lIncl_MaxThsh))
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
                if "leadingFatJetPNet_Xto4bv2_Htoaa4b_SBplusSRWP%s" % (wp_) in self.sel_conditions_all_list:
                    selection.add(
                        "leadingFatJetPNet_Xto4bv2_Htoaa4b_SBplusSRWP%s" % (wp_),
                        ( (leadingFatJet_PNet_Xto4bv2_Htoaa4b >  bTagWPs[self.datasetInfo["era"]]['PNet_Xto4bv2_Htoaa4b']['SBWP-%s' % wp_]) )
                    )
                


                

            

           
        if "nLeptonsTight" in self.sel_conditions_all_list:
            selection.add(
                "nLeptonsTight",
                ( nLeptonsTight <= self.objectSelector.NLeptonsTight_MaxThsh )
            )
        if "nLeptonsTightZero" in self.sel_conditions_all_list:
            selection.add(
                "nLeptonsTightZero",
                ( nLeptonsTight <= 0 )
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

        if "METPt" in self.sel_conditions_all_list:
            selection.add(
                "METPt",
                ((events.MET.pt >  self.objectSelector.METPt_ZvvIncl_MinThsh) & 
                 (events.MET.pt <= self.objectSelector.METPt_ZvvIncl_MaxThsh) )
            )
        
        if "dPhiLeadingFJHto4bAndMet" in self.sel_conditions_all_list:
            selection.add(
                "dPhiLeadingFJHto4bAndMet",
                ( calculate_AbsDeltaPhi(METToUse.phi_toUse, leadingFatJet.phi) > self.objectSelector.DPhi_FJHto4b_MET_MinThsh )
            )


        if "DijetVBFVeto" in self.sel_conditions_all_list:
            selection.add(
                "DijetVBFVeto",
                ~ ( (mass_leadingPair_ak4Jets_nonoverlaping_leadingFatJet > self.objectSelector.VBFDijetMass_MinThsh) & 
                    (dEta_leadingPair_ak4Jets_nonoverlaping_leadingFatJet  > self.objectSelector.VBFDijetEta_MinThsh)   )
            )
        if "jjvbf_Loose" in self.sel_conditions_all_list:
            selection.add(
                "jjvbf_Loose",
                (
                    (mass_leadingPair_ak4Jets_nonoverlaping_leadingFatJet >  self.objectSelector.MassVBFjj_MinThsh_Loose) &
                    (dEta_leadingPair_ak4Jets_nonoverlaping_leadingFatJet >  self.objectSelector.DEtaVBFjj_MinThsh_Loose) 
                )
            )

        if "BJetVeto" in self.sel_conditions_all_list:
            selection.add(
                "BJetVeto",
                nAk4JetsCentral_bTag_nonoverlaping_leadingFatJet == 0
            )

        




            
            



        ################
        # EVENT WEIGHTS
        ################
        



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
            
            
            for iSelection in self.sel_names_all.keys():
                iName = f"{iSelection}: {self.sel_names_all[iSelection]}"
                sel_i = selection.all(* self.sel_names_all[iSelection])
                selection.names
                output['cutflow'][iName] += sel_i.sum()






               


        for sel_name in self.sel_names_all.keys(): # loop of list of selections
            sHExt = "_%s" % (sel_name)
            sel_SR_toUse = selection.all(* self.sel_names_all[sel_name])





                


                    

                                
            output['hLeadingFatJetPt%s'%(sHExt)].fill(
                #dataset=dataset,
                Pt=(leadingFatJet.pt_toUse[sel_SR_toUse]),
                #systematic=syst
            )                                    
            output['hMET_pT%s'%(sHExt)].fill(
                #dataset=dataset,
                Pt=(METToUse.pt_toUse[sel_SR_toUse]),
                #systematic=syst
            )
            nAk4Jets_nonoverlaping_leadingFatJet = ak.fill_none(ak.count(ak4Jets_nonoverlaping_leadingFatJet.eta, axis=1), 0)
            if ak.sum((sel_SR_toUse & (nAk4Jets_nonoverlaping_leadingFatJet>0))) > 0:
                output['hAK4JetPtNonoverlapFatJet%s'%(sHExt)].fill(
                    #dataset=dataset,
                    Pt=ak.firsts(ak4Jets_nonoverlaping_leadingFatJet.pt_toUse[sel_SR_toUse & (nAk4Jets_nonoverlaping_leadingFatJet>0)]),
                    #systematic=syst
                )                                    
            

            for hltPathName_ in HLTTriggers[self.datasetInfo["era"]]:
                HLTName_toUse = hltPathName_.replace('HLT_', '')
                if HLTName_toUse not in events.HLT.fields: continue
                mask_HLT = events.HLT[HLTName_toUse] == True

                sel_SRplusTrg_toUse = (sel_SR_toUse & mask_HLT)

                output['hLeadingFatJetPt%s_%s'%(sHExt, hltPathName_)].fill(
                    #dataset=dataset,
                    Pt=(leadingFatJet.pt_toUse[sel_SRplusTrg_toUse]),
                    #systematic=syst
                )                                    
                output['hMET_pT%s_%s'%(sHExt, hltPathName_)].fill(
                    #dataset=dataset,
                    Pt=(METToUse.pt_toUse[sel_SRplusTrg_toUse]),
                    #systematic=syst
                )   
                if ak.sum((sel_SRplusTrg_toUse & (nAk4Jets_nonoverlaping_leadingFatJet>0))) > 0:
                    output['hAK4JetPtNonoverlapFatJet%s_%s'%(sHExt, hltPathName_)].fill(
                        #dataset=dataset,
                        Pt=ak.firsts(ak4Jets_nonoverlaping_leadingFatJet.pt_toUse[sel_SRplusTrg_toUse & (nAk4Jets_nonoverlaping_leadingFatJet>0)]),
                        #systematic=syst
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


def calculateEfficiency(N, D): # N: Neumerator, D: Denominator 
    efficiency = np.where(
        (D > 0),
        N / D,
        D
    )     
    # Binomal approach https://indico.cern.ch/event/66256/contributions/2071577/attachments/1017176/1447814/EfficiencyErrors.pdf
    #errorEfficiency = np.sqrt( D * efficiency * (1-efficiency) ) 

    # Baysian approach https://indico.cern.ch/event/66256/contributions/2071577/attachments/1017176/1447814/EfficiencyErrors.pdf
    errorEfficiency = np.sqrt( ((N+1)*(N+2)/((D+2)*(D+3))) - ((N+1)*(N+1)/((D+2)*(D+2))) ) 
    
    # forcefully set error=0 when D=0
    errorEfficiency = np.where(
        (abs(D) < 1e-6),
        np.full_like(D, 0),
        errorEfficiency
    )

    #print(f"{np.stack((efficiency, errorEfficiency, errorEfficiency1), axis=1) = }")


    return efficiency, errorEfficiency
    
    
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
    triggers            = config['triggers'] if 'triggers' in config else ''
    if isMC:
        sample_crossSection = config["crossSection"]
        sample_nEvents      = config["nEvents"]
        sample_sumEvents    = config["sumEvents"] if config["sumEvents"] > 0 else sample_nEvents
        if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
        systematicsToRun    = config["systematics"].lower() if "systematics" in config else 'no'

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
        "triggers":        triggers,
    }
    if isMC:
        sampleInfo["sample_crossSection"]   = sample_crossSection
        sampleInfo["sample_sumEvents"]      = sample_sumEvents
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
        treename="Events",
        processor_instance=HToAATo4bProcessor(
            datasetInfo=sampleInfo
        )
    )
    print(f"metrics: {metrics}", flush=flushStdout)
    #print(f"\noutput: {output}", flush=flushStdout)
    #print(f"\noutput.keys(): {output.keys()}", flush=flushStdout)


    if 'cutflow' in output.keys():
        print("Cutflow::", flush=flushStdout)        
        for key in output['cutflow'].keys():
            #print(f"{key = }", flush=flushStdout)
            if key.startswith(sWeighted): continue # to print weighted and unweighted events for cuts on the same line

            #print("%10f\t%10d\t%s" % (output['cutflow'][sWeighted+key], output['cutflow'][key], key), flush=flushStdout)
            print("%10d\t%s" % (output['cutflow'][key], key), flush=flushStdout)
    
    ##
    datasetName_part1               = sample_dataset.split('/')[1]
    datasetInfo = {}
    datasetInfo['isSignalGGH']      = True if "SUSY_GluGluH_01J_HToAATo4B" in datasetName_part1 else False
    datasetInfo['isSignalVBFH']     = True if "SUSY_VBFH_HToAATo4B"        in datasetName_part1 else False
    datasetInfo['isSignalWH']       = True if "SUSY_WH_WToAll_HToAATo4B"   in datasetName_part1 else False
    datasetInfo['isSignalZH']       = True if "SUSY_ZH_ZToAll_HToAATo4B"   in datasetName_part1 else False
    datasetInfo['isSignalTTH']      = True if "SUSY_TTH_TTToAll_HToAATo4B" in datasetName_part1 else False

    ## Read HLT triggers list saved from TSG twiki(https://twiki.cern.ch/twiki/bin/viewauth/CMS/HLTPathsRunIIList#2017) into csv file
    eraSimple = era
    if '2016' in era: eraSimple = '2016'
    sfHLTTriggersList = 'data/triggerList/HLTTriggers_%s.csv' % (eraSimple)
    HLTTriggersList_TSG_dict = {}
    if os.path.isfile(sfHLTTriggersList):
        reader = csv.reader(open(sfHLTTriggersList, 'r'))

        for row in reader:
            #print(f"{row = }")
            if not row[0].startswith('HLT'): continue

            hltPath = row[0]
            hltPath = hltPath.replace('_v', '')
            lumiActual    = float(row[1])
            lumiEffective = float(row[2])
            dataset       = row[-1]

            HLTTriggersList_TSG_dict[hltPath] = {'lumiActual': lumiActual, 'lumiEffective': lumiEffective, 'dataset': dataset}
            #print(f"{hltPath}: {HLTTriggersList_TSG_dict[hltPath]}")


    
    if sOutputFile is not None:
        if not sOutputFile.endswith('.root'): sOutputFile += '.root'
        #sOutputFile = sOutputFile.replace('.root', '_wCoffea.root') # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        sample_category_toUse = sample_category        
        sDir1 = 'evt/%s' % (sample_category_toUse)
        sDir1_toUse = '%s' % (sDir1) 
        sDir1Effi_toUse = '%s_%s' % (sDir1, 'Effi')     

        with uproot.recreate(sOutputFile) as fOut:
            trgEffi_dict = {}
            trgEffiTotal_dict = {}            
            for key, value in output.items():
                #print(f"key: {key},  value ({type(value)}): {value}")
                sHistoName_toUse = key


                if not isinstance(value, hist.Hist): continue
                if 'HLT' in key: continue # read corresponding trigger ('HLT') histograms togather with efficiency_denominator histograms

                hEffi_Denom = value.to_hist()

                fOut['%s/%s' % (sDir1_toUse, sHistoName_toUse)] = hEffi_Denom

                # read efficiency_numerator histograms and calculate hEfficiency
                for hltPathName_ in HLTTriggers[era]:
                    hEffi_Nume = output['%s_%s'%(sHistoName_toUse, hltPathName_)].to_hist()

                    trgEffi_total, _            = calculateEfficiency(hEffi_Nume.sum(flow=True), hEffi_Denom.sum(flow=True))
                    efficiency, errorEfficiency = calculateEfficiency(hEffi_Nume.values(), hEffi_Denom.values())
                    effi_error_stack = np.stack((efficiency, errorEfficiency), axis=1)

                    hEffi    = histLib.Hist.new.Variable(hEffi_Nume.axes[0].edges, name=hEffi_Nume.axes[0].name, label=hEffi_Nume.axes[0].label).Weight()
                    hEffi[:] = effi_error_stack

                    
                    trgEffiTotal_dict[hltPathName_] = trgEffi_total.item()
                    if hltPathName_ not in trgEffi_dict:
                        trgEffi_dict[hltPathName_] = {sHistoName_toUse: hEffi}
                    else:
                        trgEffi_dict[hltPathName_][sHistoName_toUse] = hEffi


                    '''
                    print(f"\n\n{hltPathName_ = }, {trgEffi_total = }, {trgEffi_total.item() = }")
                    print(f"\n{hEffi_Nume = }, {hEffi_Nume.sum() = }, {hEffi_Nume.sum(flow=True) = }, \n{hEffi_Nume.values() = }, \n{hEffi_Nume.variances() = }", flush=flushStdout)
                    print(f"\n{hEffi_Denom = }, {hEffi_Denom.sum() = }, {hEffi_Denom.sum(flow=True) = }, \n{hEffi_Denom.values() = }, \n{hEffi_Denom.variances() = }", flush=flushStdout)
                    print(f"\n{hEffi = }, {hEffi.sum() = }, \n{hEffi.values() = }, \n{hEffi.variances() = }", flush=flushStdout)
                    #print(f"{list(zip(efficiency, errorEfficiency)) = }")
                    #print(f"{np.stack((efficiency, errorEfficiency), axis=1) = }")
                    '''
                    

                    fOut['%s/%s_%s'             % (sDir1_toUse,     sHistoName_toUse,hltPathName_)] = hEffi_Nume
                    #fOut['%s/hEfficiency_%s_%s' % (sDir1Effi_toUse, sHistoName_toUse,hltPathName_)] = hEffi

            #print(f"{trgEffiTotal_dict = }")
            trgEffiTotal_dict = dict(sorted(trgEffiTotal_dict.items(), key=lambda item: item[1], reverse=True))

            ## First print all triggers, efficiency etc
            print("trgEffiTotal_dict::\n HLTpath: \t trgEffiTotal, \tlumiActual, \tlumiEffective, \tdataset ")
            #print(json.dumps(trgEffiTotal_dict, indent=4))
            for hltPathName_, trgEffiTotal in trgEffiTotal_dict.items():
                if hltPathName_ in HLTTriggersList_TSG_dict:
                    lumiActual    = HLTTriggersList_TSG_dict[hltPathName_]['lumiActual']
                    lumiEffective = HLTTriggersList_TSG_dict[hltPathName_]['lumiEffective']
                    dataset       = HLTTriggersList_TSG_dict[hltPathName_]['dataset']                    
                    trgDetails = '%4.2f, \t%5.2f, \t%5.2f, \t%s' % (trgEffiTotal, lumiActual, lumiEffective, dataset) 
                    print("%-80s: \t %s" % (hltPathName_, trgDetails))


            ## Save/print selected triggers only
            trgSelected_dict = {}
            for hltPathName_, trgEffiTotal in trgEffiTotal_dict.items():
                if trgEffiTotal < 0.10: continue # ignore triggers with efficiency < 10%

                if hltPathName_ in HLTTriggersList_TSG_dict:
                    lumiActual    = HLTTriggersList_TSG_dict[hltPathName_]['lumiActual']
                    lumiEffective = HLTTriggersList_TSG_dict[hltPathName_]['lumiEffective']
                    dataset       = HLTTriggersList_TSG_dict[hltPathName_]['dataset']
                    if datasetInfo['isSignalGGH'] or datasetInfo['isSignalVBFH']:
                        datasetToUse = 'JetHT'
                    if datasetInfo['isSignalZH']:
                        datasetToUse = 'MET'
                    
                    isTrgUnprescale = (lumiEffective > 0.7*lumiActual)                    
                    isDatasetCorrect = (datasetToUse == dataset)

                    if not isTrgUnprescale:  continue
                    #ßif not isDatasetCorrect: continue

                    trgSelected_dict[hltPathName_] = '%4.2f, \t%5.2f, \t%5.2f, \t%s' % (trgEffiTotal, lumiActual, lumiEffective, dataset) #[trgEffiTotal, lumiActual, lumiEffective, dataset ]

                
                for sHistoName_toUse, hEffi in trgEffi_dict[hltPathName_].items():
                    fOut['%s/%s/hEfficiency_%s' % (sDir1Effi_toUse, sHistoName_toUse,hltPathName_)] = hEffi

            print("\n\n\n Selected triggers::\n HLTpath: \t trgEffiTotal, \tlumiActual, \tlumiEffective, \tdataset")
            for hltPathName_, trgDetails in trgSelected_dict.items():
                print("%-80s: \t %s" % (hltPathName_, trgDetails))


                
        

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
    
