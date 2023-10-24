#htoaa triggerStudy main code

import os
import sys
from datetime import datetime
#import time
print(f"htoaa_triggerStudy_GGFMode:: here1 {datetime.now() = }"); sys.stdout.flush()
import subprocess
import json
from urllib.request import urlopen
import glob
from collections import OrderedDict as OD
import time
import tracemalloc
import math
import numpy as np
from copy import copy, deepcopy
#import uproot
#import uproot3 as uproot
import uproot as uproot
#import parse
from parse import *
import logging


'''
Trigger efficiency calculation for 
GGF -> H->aa->4b boosted analysis 

References:
  * Coffea framework used for TTGamma analysis: https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py
* Coffea installation: /home/siddhesh/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea
'''
#import coffea.processor as processor
from coffea import processor, util
from coffea.nanoevents import schemas
from coffea.nanoevents.methods import nanoaod, vector
from coffea.analysis_tools import PackedSelection, Weights
#from coffea.lookup_tools import extractor
from coffea.lookup_tools.dense_lookup import dense_lookup
#from coffea.lumi_tools import LumiMask
#import hist
from coffea import hist
import awkward as ak
#import uproot
#from dask.distributed import Client
from particle import Particle # For PDG particle listing https://github.com/scikit-hep/particle


from htoaa_Settings import *
from htoaa_CommonTools import (
    GetDictFromJsonFile, akArray_isin,
    selectRunLuminosityBlock,
    calculate_lumiScale, getLumiScaleForPhSpOverlapRewgtMode, getSampleHTRange, # update_crosssection, 
    getNanoAODFile, setXRootDRedirector,  xrdcpFile,
    selectMETFilters,
    selGenPartsWithStatusFlag,
    getTopPtRewgt, getPURewgts, getHTReweight,
    fillHist,
    printVariable, printVariablePtEtaPhi,
    insertInListBeforeThisElement,
)
from htoaa_Samples import (
    kData, kQCD_bEnrich, kQCD_bGen, kQCDIncl
)

from inspect import currentframe, getframeinfo
frameinfo = getframeinfo(currentframe())


# use GOldenJSON

 
printLevel = 0
nEventToReadInBatch =  0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze = -1 # 1000 # 100000 # -1
flushStdout = False
#pd.set_option('display.max_columns', None)  

#print("".format())

sWeighted = "Wtd: "




# -----------------------------------------------------------------------------------
class ObjectSelection:
    def __init__(self, era):
        self.era = era
        
        self.tagger_btagDeepB = 'DeepCSV'
        self.wp_btagDeepB = 'M'
        self.wp_ParticleNetMD_XbbvsQCD = 'L'

        self.FatJetPtThsh  = 400 #170
        self.FatJetEtaThsh = 2.4
        self.FatJetJetID   = int(JetIDs.tightIDPassingLeptonVeto)


        self.MuonPtThsh  = 30
        self.MuonEtaThsh = 2.4
        self.MuonMVAId     =  3 # (1=MvaLoose, 2=MvaMedium, 3=MvaTight, 4=MvaVTight, 5=MvaVVTight)
        self.MuonMiniIsoId =  3 # (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)

        self.dRMuonFatJetThsh = 1.5
        


    def selectMuons(self, eventsObj):
        maskSelMuons = (
            (eventsObj.mvaId >= self.MuonMVAId) &
            (eventsObj.miniIsoId >= self.MuonMiniIsoId)
        )
        return eventsObj[maskSelMuons]
    

    def selectElectrons(self, eventsObj):
        maskSelElectrons = (
            (eventsObj[self.ElectronMVAId] > 0)
        )
        return eventsObj[maskSelElectrons]





class HToAATo4bProcessor(processor.ProcessorABC):
    def __init__(self, datasetInfo={}):

         
        global runMode_SignalGenChecks;       runMode_SignalGenChecks  = False; # True
        global runMode_QCDGenValidation;      runMode_QCDGenValidation = False; # True
        global runMode_GenLHEPlots;           runMode_GenLHEPlots      = False
        global runMode_SignificancsScan2D;    runMode_SignificancsScan2D = False
        

        ak.behavior.update(nanoaod.behavior)

        self.datasetInfo = datasetInfo
        self.objectSelector = ObjectSelection(era=self.datasetInfo["era"])

        self.datasetInfo['isSignal'       ]  = False
        self.datasetInfo['isQCD'          ]  = False
        self.datasetInfo['isQCDIncl'      ]  = False
        self.datasetInfo['isQCD_bEnrich'  ]  = False
        self.datasetInfo['isQCD_bGen'     ]  = False
        self.datasetInfo['isTTbar'        ]  = False
        self.datasetInfo['isPythiaTuneCP5']  = False
        
        if self.datasetInfo['isMC']:
            self.datasetInfo['isSignal']         = True if "HToAATo4B"   in self.datasetInfo['sample_category'] else False
            self.datasetInfo['isQCDIncl']        = True if kQCDIncl      in self.datasetInfo['sample_category'] else False
            self.datasetInfo['isQCD_bEnrich']    = True if kQCD_bEnrich  in self.datasetInfo['sample_category'] else False
            self.datasetInfo['isQCD_bGen']       = True if kQCD_bGen     in self.datasetInfo['sample_category'] else False
            self.datasetInfo['isQCD']            = (self.datasetInfo['isQCDIncl']     or \
                                                     self.datasetInfo['isQCD_bEnrich'] or \
                                                     self.datasetInfo['isQCD_bGen'])
            sample_HT_Min, sample_HT_Max = getSampleHTRange( self.datasetInfo["datasetNameFull"] )
            self.datasetInfo['sample_HT_Min']    = sample_HT_Min
            self.datasetInfo['sample_HT_Max']    = sample_HT_Max
            self.datasetInfo['isTTbar']          = True if self.datasetInfo['sample_category'].startswith('TTTo') else False
            self.datasetInfo['isPythiaTuneCP5']  = True if 'TuneCP5' in self.datasetInfo["datasetNameFull"] else False

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


        ## List of all analysis selection condition
        global HLT_Mu_name
        HLT_Mu_name = "HLT_IsoMu27" #"HLT_IsoMu27_v" #"HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4" 
        
        global HLT_AK8PFJet330_name
        HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4"         
        
        # sel_names_all = dict of {"selection name" : [list of different cuts]}; for cut-flow table 
        self.sel_names_all = OD([
            ("SR",                    [
                "nPV",
                "METFilters",
                "leadingMuonPt",
                "leadingMuonEta",
                HLT_Mu_name,
                "dR_Muon_FatJet",
                "leadingFatJetEta",
                "JetID",                 ## Denominator for trigger efficiency calculation
                "L1_SingleJet180",
                HLT_AK8PFJet330_name,    ## Numerator for trigger efficiency calculation
            ]),
        ])
        if not self.datasetInfo['isMC']: 
            self.sel_names_all["SR"].insert(0, "run:ls")

        else:
            if self.datasetInfo['isQCD']: #self.sel_names_all["SR"].append("QCDStitch")
                self.sel_names_all["SR"] = insertInListBeforeThisElement(
                    list1                  = self.sel_names_all["SR"], 
                    sConditionToAdd        = "QCDStitch", 
                    addBeforeThisCondition = "METFilters")                

        if self.datasetInfo["era"] == Era_2018:
            # 2018HEM1516Issue ----------------
            #self.sel_names_all["SR"].append("2018HEM1516Issue")
            # Update self.sel_names_all["SR"] by adding "2018HEM1516Issue" before "leadingFatJetMSoftDrop" in the list    
            '''          
            self.sel_names_all["SR"] = insertInListBeforeThisElement(
                list1                  = self.sel_names_all["SR"], 
                sConditionToAdd        = "2018HEM1516Issue", 
                addBeforeThisCondition = "leadingFatJetMSoftDrop"
            )
            '''
            pass

            

        # selection region addition each SR conditions successively
        for iCondition in range(self.sel_names_all["SR"].index("leadingMuonPt"), len(self.sel_names_all["SR"]) - 1):
            conditionName = self.sel_names_all["SR"][iCondition]
            self.sel_names_all["sel_%s" % conditionName] = self.sel_names_all["SR"][0 : (iCondition+1)]
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
                self.datasetInfo["MCSamplesStitchOption"] == MCSamplesStitchOptions.PhSpOverlapRewgt:
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
        eta_axis              = hist.Bin("Eta",                    r"$#eta$",                    100,      -6,       6)
        phi_axis              = hist.Bin("Phi",                    r"$\phi$",                    100,   -3.14,    3.13)
        #mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       200, 0, 600)
        #mass_axis             = hist.Bin("Mass",      r"$m$ [GeV]",       400, 0, 200)
        mass_axis             = hist.Bin("Mass",                   r"$m$ [GeV]",                 300,       0,     300)
        mass_axis1            = hist.Bin("Mass1",                  r"$m$ [GeV]",                 300,       0,     300)
        mass10_axis           = hist.Bin("Mass10",                 r"$m$ [GeV]",                 300,       0,      10)
        logMass3_axis         = hist.Bin("logMass3",               r"$m$ [GeV]",                 300,       0,       3)
        mlScore_axis          = hist.Bin("MLScore",                r"ML score",                  100,    -1.1,     1.1)
        mlScore_axis1         = hist.Bin("MLScore1",               r"ML score",                  100,    -1.1,     1.1)
        jetN2_axis            = hist.Bin("N2",                     r"N2b1",                      100,       0,       3)
        jetN3_axis            = hist.Bin("N3",                     r"N3b1",                      100,       0,       5)
        jetTau_axis           = hist.Bin("TauN",                   r"TauN",                      100,       0,       1)
        deltaR_axis           = hist.Bin("deltaR",                 r"$delta$ r ",                500,       0,       5)
        #HT_axis               = hist.Bin("HT",                     r"HT",                       3000,       0,    3000)
        HT_axis               = hist.Bin("HT",                     r"HT",                       4000,       0,    4000)
        PytPartStatus_axis    = hist.Bin("PytPartStatus",          r"PytPartStatus",             421,  -210.5,   210.5)
        boolean_axis          = hist.Bin("Boolean",                r"Boolean",                     2,    -0.5,     1.5)
        pdgId_axis            = hist.Bin("PdgId",                  r"PdgId",                     101,    -0.5,   100.5)
        alphaS_axis           = hist.Bin("alphaS",                 r"alphaS",                    101,    0.01,     0.2)
        PU_axis               = hist.Bin("PU",                     r"PU",                         99,     0.0,    99.0)
        
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

        # RECO-level histograms --------------------------------------------------------------------------------------------------------------
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

                    ('nSelMuon'+sHExt,                                {sXaxis: nObject_axis,    sXaxisLabel: 'No. of selected muons'}),
                    ('hLeadingMuonPt'+sHExt,                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading muon)$ [GeV]"}),
                    ('hLeadingMuonEta'+sHExt,                         {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading muon)"}),
                    ('hLeadingMuonPhi'+sHExt,                         {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading muon)"}),

                    #('nSelFatJet'+sHExt,                                {sXaxis: nObject_axis,    sXaxisLabel: 'No. of selected FatJets'}),
                    ('hdR_leadingMuon_leadingFatJet'+sHExt,             {sXaxis: deltaR_axis,     sXaxisLabel: r"$delta$ r(leading muon, leading FatJet) "}),
                    ('hLeadingFatJetPt'+sHExt,                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
                    ('hLeadingFatJetPt_msoftdropGt60_btagHbbGtnp4'+sHExt,                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet, msoftdropGt60_btagHbbGtnp4)$ [GeV]"}),
                    ('hLeadingFatJetEta'+sHExt,                         {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
                    ('hLeadingFatJetPhi'+sHExt,                         {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
                    ('hLeadingFatJetMass'+sHExt,                        {sXaxis: mass_axis,       sXaxisLabel: r"m (leading FatJet) [GeV]"}),
                    ('hLeadingFatJetMSoftDrop'+sHExt,                   {sXaxis: mass_axis,       sXaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
                    ('hLeadingFatJetMSoftDrop_pTGt400_btagHbbGtnp4'+sHExt,                   {sXaxis: mass_axis,       sXaxisLabel: r"m_{soft drop} (leading FatJet, pTGt400_btagHbbGtnp4) [GeV]"}),
                    ('hLeadingFatJetId'+sHExt,                          {sXaxis: nObject_axis,    sXaxisLabel: r"jet Id (leading FatJet)"}),
                    ('hLeadingFatJetParticleNetMD_XbbOverQCD'+sHExt,    {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD Xbb/(Xbb + QCD)"}),
                    ('hLeadingFatJetBtagCSVV2'+sHExt,    {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJet BtagCSVV2"}),
                    ('hLeadingFatJetBtagDDBvLV2'+sHExt,    {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJet BtagDDBvLV2"}),
                    ('hLeadingFatJetBtagDeepB'+sHExt,    {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJet BtagDeepB"}),
                    ('hLeadingFatJetBtagHbb'+sHExt,    {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJet BtagHbb"}),
                    ('hLeadingFatJetBtagHbb_pTGt400_msoftdropGt60'+sHExt,    {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJet BtagHbb pTGt400_msoftdropGt60"}),
                    
                ]))

                ### 2-D distribution --------------------------------------------------------------------------------------------------------
                histos.update(OD([                    
                    ('hLeadingFatJetEta_vs_Phi'+sHExt,             
                     {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)",
                      sYaxis: phi_axis,        sYaxisLabel: r"\phi (leading FatJet)"}),                                        
                    ('hLeadingFatJetPt_vs_Eta'+sHExt,             
                     {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                      sYaxis: eta_axis,        sYaxisLabel: r"\eta (leading FatJet)"}),                    
                    ('hLeadingFatJetPt_vs_Eta_msoftdropGt60_btagHbbGtnp4'+sHExt,             
                     {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet, msoftdropGt60_btagHbbGtnp4)$ [GeV]",
                      sYaxis: eta_axis,        sYaxisLabel: r"\eta (leading FatJet, msoftdropGt60_btagHbbGtnp4)"}),                    
                    ('hLeadingFatJetPt_vs_Phi'+sHExt,             
                     {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                      sYaxis: phi_axis,        sYaxisLabel: r"\phi (leading FatJet)"}),                    
                    ('hLeadingFatJetPt_vs_Mass'+sHExt,             
                     {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                      sYaxis: mass_axis,       sYaxisLabel: r"m (leading FatJet) [GeV]"}),                    
                    ('hLeadingFatJetPt_vs_MSoftDrop'+sHExt,             
                     {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                      sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),                    
                    ('hLeadingFatJetPt_vs_MSoftDrop_btagHbbGtnp4'+sHExt,             
                     {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet, btagHbbGtnp4)$ [GeV]",
                      sYaxis: mass_axis,       sYaxisLabel: r"m_{soft drop} (leading FatJet, btagHbbGtnp4) [GeV]"}),                    
                    ('hLeadingFatJetPt_vs_ParticleNetMD_XbbOverQCD'+sHExt,             
                     {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                      sYaxis: mlScore_axis,    sYaxisLabel: r"LeadingFatJetParticleNetMD Xbb/(Xbb + QCD)"}),      
                    ('hLeadingFatJetPt_vs_BtagHbb'+sHExt,             
                     {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]",
                      sYaxis: mlScore_axis,    sYaxisLabel: r"LeadingFatJet BtagHbb"}),  
                    ('hLeadingFatJetPt_vs_BtagHbb_msoftdropGt60'+sHExt,             
                     {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet, msoftdropGt60)$ [GeV]",
                      sYaxis: mlScore_axis,    sYaxisLabel: r"LeadingFatJet BtagHbb, msoftdropGt60"}),                           
                    ('hLeadingFatJetMSoftDrop_vs_BtagHbb_PtGt400'+sHExt,             
                     {sXaxis: mass_axis,         sXaxisLabel: r"m_{soft drop} (leading FatJet, PtGt400) [GeV]",
                      sYaxis: mlScore_axis,    sYaxisLabel: r"LeadingFatJet BtagHbb, PtGt400"}),                           

                ]))



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
                        systematic_axis,
                        hXaxis, #nObject_axis,                        
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
                        systematic_axis,
                        hXaxis, #nObject_axis,
                        hYaxis,                        
                    )
                })
                


    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        dataset = events.metadata["dataset"] # dataset label
        print(f"process():: {self.datasetInfo['sample_category'] = }, {dataset = }", flush=flushStdout)

        if printLevel >= 1:
            print(f"\n events.fields ({type(events.fields)}): {events.fields}"); sys.stdout.flush()
            print(f"\n events.HLT.fields ({type(events.HLT.fields)}): {events.HLT.fields}", flush=flushStdout)


             
        if nEventsToAnalyze != -1:
            printVariable('\n (run:ls:event): ', ak.zip([events.run, events.luminosityBlock, events.event])); #sys.stdout.flush()           

        if not self.datasetInfo['isMC']:
            print(f" {np.unique(events.run, return_counts=True) = } "); sys.stdout.flush()  



               
            
            output = self.accumulator.identity()
            systematics_shift = [None]
            for _syst in systematics_shift:
                output += self.process_shift(events, _syst)
        else:
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
                      
            genBQuarks = events.GenPart[(
                (abs(events.GenPart.pdgId) == PDGID_BottomQuark )
            )]
            genBQuarks_pT = ak.sort(genBQuarks.pt, axis=-1, ascending=False)
            genBQuarks_first = ak.firsts(genBQuarks)
            mask_genBQuarks = (ak.count(genBQuarks.pdgId, axis=1) >= 1)

            mask_genBQuarks_pTAbvTrsh = ak.any((genBQuarks.pt > 15.0), axis=1)

            idx_genBQuarks_pTsort = ak.argsort(genBQuarks.pt, axis=-1, ascending=False)
            

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
           

            # Check if events has b-quark outgoing from hard subprocess -----------------------------------------------
            mask_genBQuarks_hardSctred = (
                (abs(events.GenPart.pdgId) == PDGID_BottomQuark ) &
                (events.GenPart.status == 23)
            )
            mask_genBQuarks_hardSctred_eventwise = ak.any(mask_genBQuarks_hardSctred, axis=1)

            genBQuarks_hardSctred = events.GenPart[mask_genBQuarks_hardSctred]
            idx_genBQuarks_hardSctred_pTsort = ak.argsort(genBQuarks_hardSctred.pt, axis=-1, ascending=False)
            

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
        # --------------------------------------------------------------------------------------------------
        

        # MC ttbar ----------------------------------------------
        #mask_1 = None
        if self.datasetInfo['isMC'] and self.datasetInfo['isTTbar'] :
            mask_genTopQuark = (
                (abs(events.GenPart.pdgId) == PDGID_TopQuark ) & 
                (events.GenPart.hasFlags("isLastCopy"))
            )   




        ################## 
        # EVENT VARIABLES
        ##################

        ## sel leptons
        muonsTight     = self.objectSelector.selectMuons(events.Muon)
        leadingMuon    = ak.firsts(muonsTight)
        
        # AK8 jet
        leadingFatJet = ak.firsts(events.FatJet)
        leadingFatJet_asSingletons = ak.singletons(leadingFatJet) # for e.g. [[0.056304931640625], [], [0.12890625], [0.939453125], [0.0316162109375]]
        
        leadingFatJetParticleNetMD_XbbvsQCD = ak.where(
            (leadingFatJet.particleNetMD_Xbb + leadingFatJet.particleNetMD_QCD) > 0,
            leadingFatJet.particleNetMD_Xbb / (leadingFatJet.particleNetMD_Xbb + leadingFatJet.particleNetMD_QCD),
            np.full(len(events), 0)
        )        
        
        dR_leadingMuon_leadingFatJet = ak.fill_none(leadingMuon.delta_r(leadingFatJet), 0)

        ## match leadingFat jet to genB 
        n_leadingFatJat_matched_genB = np.full(len(events), 0)
        if self.datasetInfo['isMC'] :
            #mask_leadingFatJat_matched_genB = leadingFatJet.delta_r(vGenBQuarksHardSctred_genBHadronsStatus2_sel) < 0.8
            #n_leadingFatJat_matched_genB = ak.sum(mask_leadingFatJat_matched_genB, axis=1)
            n_leadingFatJat_matched_genB = leadingFatJet.nBHadrons


        if printLevel >= 5:
            printVariablePtEtaPhi('\n leadingMuon',   leadingMuon)
            printVariablePtEtaPhi('\n leadingFatJet', leadingFatJet)
            printVariable('\n dR_leadingMuon_leadingFatJet', dR_leadingMuon_leadingFatJet)



        #####################
        # EVENT SELECTION
        #####################
        
        # create a PackedSelection object
        # this will help us later in composing the boolean selections easily
        selection = PackedSelection()

        if "run:ls" in self.sel_names_all["SR"]:
            # self.datasetInfo['dataLSSelGoldenJSON']
            # using Coffea built-in function: mask_lumi = LumiMask(golden_json_path)(events.run,events.luminosityBlock)
            #selection.add("run:ls", LumiMask(sFilesGoldenJSON[self.datasetInfo["era"]])(events.run,events.luminosityBlock) )
            
            # using selectRunLuminosityBlock function from htoaa_CommonTools
            selection.add("run:ls", selectRunLuminosityBlock(
                dataLSSelGoldenJSON  = self.datasetInfo['dataLSSelGoldenJSON'], 
                runNumber_list       = events.run, 
                luminosityBlock_list = events.luminosityBlock 
                ))

        if "nPV" in self.sel_names_all["SR"]:
            # nPVGood >= 1
            selection.add("nPV", events.PV.npvsGood >= 1)

        if "METFilters" in self.sel_names_all["SR"]:
            #mask_METFilters = selectMETFilters(events.Flag, self.datasetInfo["era"], self.datasetInfo['isMC'])
            #printVariable('\n mask_METFilters', mask_METFilters)
            selection.add(
                "METFilters", 
                selectMETFilters(events.Flag, self.datasetInfo["era"], self.datasetInfo['isMC'])
            )

        if "leadingMuonPt" in self.sel_names_all["SR"]:
            selection.add(
                "leadingMuonPt",
                leadingMuon.pt > self.objectSelector.MuonPtThsh
            )            

        if "leadingMuonEta" in self.sel_names_all["SR"]:
            selection.add(
                "leadingMuonEta",
                abs(leadingMuon.eta) < self.objectSelector.MuonEtaThsh
            )      

        if HLT_Mu_name in self.sel_names_all["SR"]:
            HLT_TRG = HLT_Mu_name.replace('HLT_', '')
            sTmp_ = f"HLT trigger to check: {HLT_TRG}."
            if HLT_TRG in events.HLT.fields:
                selection.add(
                    HLT_Mu_name,
                    events.HLT[HLT_TRG] == True
                )
                sTmp_ += f"    HLT trigger {HLT_TRG} field present in input NanoAOD file."
            else:
                selection.add(
                    HLT_Mu_name,
                    falses_list
                )
            print(sTmp_)

        if "dR_Muon_FatJet" in self.sel_names_all["SR"]:  
            selection.add(
                "dR_Muon_FatJet",
                dR_leadingMuon_leadingFatJet > self.objectSelector.dRMuonFatJetThsh
            )

        if "leadingFatJetEta" in self.sel_names_all["SR"]:
            selection.add(
                "leadingFatJetEta",
                abs(leadingFatJet.eta) < self.objectSelector.FatJetEtaThsh
            )

        if "JetID"  in self.sel_names_all["SR"]:
            selection.add(
                "JetID", 
                leadingFatJet.jetId == self.objectSelector.FatJetJetID
            )


        if "QCDStitch" in self.sel_names_all["SR"]:
            selection.add(
                "QCDStitch",
                #mask_QCD_stitch_eventwise == True
                mask_QCD_stitch_eventwise
            )

        if "L1_SingleJet180" in self.sel_names_all["SR"]:
            selection.add(
                "L1_SingleJet180",
                events.L1.SingleJet180 == True
            )
 
        if HLT_AK8PFJet330_name in self.sel_names_all["SR"]:
            # some files of Run2018A do not have HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 branch
            #HLT_AK8PFJet330_name = None
            if "AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4" in events.HLT.fields:
                #HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4"
                selection.add(
                    HLT_AK8PFJet330_name,
                    events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 == True
                )
#            elif "AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02" in events.HLT.fields:
#                #HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02"
#                selection.add(
#                    HLT_AK8PFJet330_name,
#                    events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02 == True
#                )
            else:
                selection.add(
                    HLT_AK8PFJet330_name,
                    falses_list
                )



        sel_SR           = selection.all(* self.sel_names_all["SR"])



        ################
        # EVENT WEIGHTS
        ################
        
        # create a processor Weights object, with the same length as the number of events in the chunk
        weights              = Weights(len(events))



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

            '''
            # MC wgt for HEM1516Issue --------------------- 
            wgt_HEM1516Issue = None
            if "2018HEM1516Issue" in self.sel_names_all["SR"]:
                wgt_HEM1516Issue = ak.where(
                    mask_HEM1516Issue, # events w/ jets in HEM15/16 affected phase space 
                    np.full(len(events), (1. - DataFractionAffectedBy2018HEM1516Issue)), 
                    ones_list
                )
            '''

            # MC PURewgt ----------------------------------
            wgt_PU = getPURewgts(
                PU_list  = events.Pileup.nTrueInt,
                hPURewgt = self.hPURewgt
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

            '''
            # MC ParticleNetMD_XbbvsQCD SFs      SFs_ParticleNetMD_XbbvsQCD
            if "leadingFatJetParticleNetMD_XbbvsQCD" in self.sel_names_all["SR"]:
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
            '''


            weights.add(
                "lumiWeight",
                weight = lumiScale_toUse
            )
            weights.add(
                "genWeight",
                weight = np.copysign(np.ones(len(events)), events.genWeight)
            )
            '''
            if "2018HEM1516Issue" in self.sel_names_all["SR"]:
                weights.add(
                    "2018HEM1516IssueWeight",
                    weight = wgt_HEM1516Issue
                )
            '''
            weights.add(
                "PUWeight",
                weight = wgt_PU
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
            '''
            if "leadingFatJetParticleNetMD_XbbvsQCD" in self.sel_names_all["SR"]:
                weights.add(
                    "SF_ParticleNetMD_XbbvsQCD",
                    weight = wgt_ParticleNetMD_XbbvsQCD
                )
            '''
            
    



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
            else:
                evtWeight                = weights.weight(weightSyst)


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





            ### RECO-level histograms ============================================================================================
            

            for sel_name in self.sel_names_all.keys(): # loop of list of selections
                

                if sel_name.startswith('Gen'): continue

                sel_SR_toUse = selection.all(* self.sel_names_all[sel_name])
                sel_SR_woSel2018HEM1516_toUse = None
                if "2018HEM1516Issue" in self.sel_names_all[sel_name]:
                    sel_names_all_toUse_ = list(set(self.sel_names_all[sel_name]) - set(["2018HEM1516Issue"])) # all sel_name conditions w/o "2018HEM1516Issue"
                    sel_SR_woSel2018HEM1516_toUse = selection.all(* sel_names_all_toUse_)
                else:
                    sel_SR_woSel2018HEM1516_toUse = sel_SR_toUse

                if len(sel_SR_toUse) == 0: continue 


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



                    if len(sel_SR_forHExt) == 0: continue 


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

                    ## leading Muon
                    sel_tmp_ = sel_SR_forHExt & (~ ak.is_none(leadingMuon.pt))
                    output['hLeadingMuonPt'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingMuon.pt[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    ) 
                    output['hLeadingMuonEta'+sHExt].fill(
                        dataset=dataset,
                        Eta=(leadingMuon.eta[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingMuonPhi'+sHExt].fill(
                        dataset=dataset,
                        Phi=(leadingMuon.phi[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )

                    # "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4"
                    ## leading FatJet
                    sel_tmp_ = sel_SR_forHExt & (~ ak.is_none(leadingFatJet.pt))

                    output['hdR_leadingMuon_leadingFatJet'+sHExt].fill(
                        dataset=dataset,
                        deltaR=(dR_leadingMuon_leadingFatJet[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetPt'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingFatJet.pt[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )   
                    output['hLeadingFatJetPt_msoftdropGt60_btagHbbGtnp4'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingFatJet.pt[( 
                            sel_tmp_ & 
                            (leadingFatJet.msoftdrop > 60) & 
                            (leadingFatJet.btagHbb > -0.4) 
                            )]),
                        systematic=syst,
                        weight=evtWeight[( 
                            sel_tmp_ & 
                            (leadingFatJet.msoftdrop > 60) & 
                            (leadingFatJet.btagHbb > -0.4) 
                            )]
                    )
                    output['hLeadingFatJetEta'+sHExt].fill(
                        dataset=dataset,
                        Eta=(leadingFatJet.eta[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetPhi'+sHExt].fill(
                        dataset=dataset,
                        Phi=(leadingFatJet.phi[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetMass'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingFatJet.mass[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetMSoftDrop'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingFatJet.msoftdrop[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetMSoftDrop_pTGt400_btagHbbGtnp4'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingFatJet.msoftdrop[(
                            sel_tmp_ &
                            (leadingFatJet.pt > 400) &
                            (leadingFatJet.btagHbb > -0.4) 
                            )]),
                        systematic=syst,
                        weight=evtWeight[(
                            sel_tmp_ &
                            (leadingFatJet.pt > 400) &
                            (leadingFatJet.btagHbb > -0.4) 
                            )]
                    )

                    output['hLeadingFatJetId'+sHExt].fill(
                        dataset=dataset,
                        nObject=(leadingFatJet.jetId[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetParticleNetMD_XbbOverQCD'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetBtagCSVV2'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.btagCSVV2[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetBtagDDBvLV2'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.btagDDBvLV2[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetBtagDeepB'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.btagDeepB[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetBtagHbb'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.btagHbb[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetBtagHbb_pTGt400_msoftdropGt60'+sHExt].fill(
                        dataset=dataset,
                        MLScore=(leadingFatJet.btagHbb[(
                            sel_tmp_ &
                            (leadingFatJet.pt > 400) &
                            (leadingFatJet.msoftdrop > 60) 
                            )]),
                        systematic=syst,
                        weight=evtWeight[(
                            sel_tmp_ &
                            (leadingFatJet.pt > 400) &
                            (leadingFatJet.msoftdrop > 60) 
                            )]
                    )
                    
                    output['hLeadingFatJetEta_vs_Phi'+sHExt].fill(
                        dataset=dataset,
                        Eta=(leadingFatJet.eta[sel_tmp_]),
                        Phi=(leadingFatJet.phi[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetPt_vs_Eta'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingFatJet.pt[sel_tmp_]),
                        Eta=(leadingFatJet.eta[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetPt_vs_Eta_msoftdropGt60_btagHbbGtnp4'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingFatJet.pt[(
                            sel_tmp_ & 
                            (leadingFatJet.msoftdrop > 60) & 
                            (leadingFatJet.btagHbb > -0.4) 
                            )]),
                        Eta=(leadingFatJet.eta[(
                            sel_tmp_ & 
                            (leadingFatJet.msoftdrop > 60) & 
                            (leadingFatJet.btagHbb > -0.4) 
                            )]),
                        systematic=syst,
                        weight=evtWeight[(
                            sel_tmp_ & 
                            (leadingFatJet.msoftdrop > 60) & 
                            (leadingFatJet.btagHbb > -0.4) 
                            )]
                    )                    
                    output['hLeadingFatJetPt_vs_Phi'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingFatJet.pt[sel_tmp_]),
                        Phi=(leadingFatJet.phi[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )                    
                    output['hLeadingFatJetPt_vs_Mass'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingFatJet.pt[sel_tmp_]),
                        Mass=(leadingFatJet.mass[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )                    
                    output['hLeadingFatJetPt_vs_MSoftDrop'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingFatJet.pt[sel_tmp_]),
                        Mass=(leadingFatJet.msoftdrop[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetPt_vs_MSoftDrop_btagHbbGtnp4'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingFatJet.pt[(
                            sel_tmp_ &
                            (leadingFatJet.btagHbb > -0.4)
                            )]),
                        Mass=(leadingFatJet.msoftdrop[(
                            sel_tmp_ &
                            (leadingFatJet.btagHbb > -0.4)
                            )]),
                        systematic=syst,
                        weight=evtWeight[(
                            sel_tmp_ &
                            (leadingFatJet.btagHbb > -0.4)
                            )]
                    )                                        
                    output['hLeadingFatJetPt_vs_ParticleNetMD_XbbOverQCD'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingFatJet.pt[sel_tmp_]),
                        MLScore=(leadingFatJetParticleNetMD_XbbvsQCD[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetPt_vs_BtagHbb'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingFatJet.pt[sel_tmp_]),
                        MLScore=(leadingFatJet.btagHbb[sel_tmp_]),
                        systematic=syst,
                        weight=evtWeight[sel_tmp_]
                    )
                    output['hLeadingFatJetPt_vs_BtagHbb_msoftdropGt60'+sHExt].fill(
                        dataset=dataset,
                        Pt=(leadingFatJet.pt[(
                            sel_tmp_ &
                            (leadingFatJet.msoftdrop > 60)
                            )]),
                        MLScore=(leadingFatJet.btagHbb[(
                            sel_tmp_ &
                            (leadingFatJet.msoftdrop > 60)
                            )]),
                        systematic=syst,
                        weight=evtWeight[(
                            sel_tmp_ &
                            (leadingFatJet.msoftdrop > 60)
                            )]
                    )
                    output['hLeadingFatJetMSoftDrop_vs_BtagHbb_PtGt400'+sHExt].fill(
                        dataset=dataset,
                        Mass=(leadingFatJet.msoftdrop[(
                            sel_tmp_ & 
                            (leadingFatJet.pt > 400)
                            )]),
                        MLScore=(leadingFatJet.btagHbb[(
                            sel_tmp_ & 
                            (leadingFatJet.pt > 400)
                            )]),
                        systematic=syst,
                        weight=evtWeight[(
                            sel_tmp_ & 
                            (leadingFatJet.pt > 400)
                            )]
                    )
                    






        return output




    def postprocess(self, accumulator):
        #pass
        return accumulator












if __name__ == '__main__':
    print("htoaa_Analysis:: main: {}".format(sys.argv)); sys.stdout.flush()
    print(f"htoaa_triggerStudy_GGFMode:: here14 {datetime.now() = }")

    if len(sys.argv) != 2:
        print("htoaa_Analysis:: Command-line config file missing.. \t **** ERROR **** \n")

    sConfig = sys.argv[1]

    
    config = GetDictFromJsonFile(sConfig)
    print("Config {}: \n{}".format(sConfig, json.dumps(config, indent=4)))
    print(f"htoaa_triggerStudy_GGFMode:: here15 {datetime.now() = }")

    lumiScale = 1
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
        luminosity          = Luminosities_forGGFMode[era]['HLT_IsoMu27'][0]  # Luminosities_Inclusive[era][0]
        sample_crossSection = config["crossSection"]
        sample_nEvents      = config["nEvents"]
        sample_sumEvents    = config["sumEvents"] if config["sumEvents"] > 0 else sample_nEvents
        if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
        lumiScale = calculate_lumiScale(luminosity=luminosity, crossSection=sample_crossSection, sumEvents=sample_sumEvents)

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
                logging.critical(f'htoaa_triggerStudy_GGFMode.py::main():: {MCSamplesStitchInputFileName = } does not exists')
                print(f'htoaa_triggerStudy_GGFMode.py::main() 11:: {MCSamplesStitchInputFileName = } does not exists')
                exit(0)
            print(f"Opening {MCSamplesStitchInputFileName = } "); sys.stdout.flush() 
            with uproot.open(MCSamplesStitchInputFileName) as f_:
                print(f"{f_.keys() = }"); sys.stdout.flush() 
                hMCSamplesStitch = f_[r'%s' % MCSamplesStitchInputHistogramName].to_hist()

        print(f"isMC: {isMC}, luminosity: {luminosity}, lumiScale: {lumiScale}")
    print(f"htoaa_triggerStudy_GGFMode:: here16 {datetime.now() = }")    
        
        
    #branchesToRead = htoaa_nanoAODBranchesToRead
    #print("branchesToRead: {}".format(branchesToRead))
    sample_dataset = sample_dataset[0] if isinstance(sample_dataset, list) else sample_dataset # dataset is list of datasets w/ same sample name, as they are considered together recently. Those set of datasets are extension of the same samples.

    print(f"isMC: {isMC}, lumiScale: {lumiScale}")
    sInputFiles_toUse = []
    for sInputFile in sInputFiles:
        if "*" in sInputFile:  sInputFiles_toUse.extend( glob.glob( sInputFile ) )
        else:                  sInputFiles_toUse.append( sInputFile )
    sInputFiles = sInputFiles_toUse
    print(f"Initial sInputFiles ({len(sInputFiles)}) (type {type(sInputFiles)}):");
    for sInputFile in sInputFiles:
        print(f"\t{sInputFile}");  sys.stdout.flush()
    print(f"htoaa_triggerStudy_GGFMode:: here17 {datetime.now() = }")

    for iFile in range(len(sInputFiles)):     
        sInputFile = sInputFiles[iFile]
        sFileLocal = './inputFiles/%s' %(os.path.basename(sInputFile))   
        #cp_command = 'eos cp' if server in ['lxplus'] else 'xrdcp'
        sInputFile, isReadingSuccessful = getNanoAODFile(
            fileName = sInputFile, 
            useLocalFileIfExists = True, 
            downloadFile = downloadIpFiles, 
            fileNameLocal = './inputFiles/%s' %(os.path.basename(sInputFile)), 
            nTriesToDownload = 3,
            server = server
            )
        if not isReadingSuccessful:
            logging.critical('htoaa_triggerStudy_GGFMode:: getNanoAODFile() for input file %s failed. **** CRITICAL ERROR ****. \nAborting...' % (sInputFile)); sys.stdout.flush();
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
        print(f"htoaa_triggerStudy_GGFMode:: {sInputFile} \t {os.path.exists(sInputFile) = }, {fileSize = } MB");     

        if fileSize > NanoAODFileSize_Min:     
            sInputFiles[iFile] = sInputFile
        else:
            logging.critical('htoaa_triggerStudy_GGFMode:: Input file %s file size below threshold (%g MB). **** CRITICAL ERROR ****. \nAborting...' % (sInputFile, NanoAODFileSize_Min) ); sys.stdout.flush();
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
    print(f"htoaa_triggerStudy_GGFMode:: here18 {datetime.now() = }")


    sampleInfo = {
        "era":             era, 
        "isMC":            isMC,
        "sample_category": sample_category,        
        "datasetNameFull": sample_dataset,
        "lumiScale":       lumiScale,
    }
    if isMC:
        sampleInfo["MCSamplesStitchOption"] = MCSamplesStitchOption
        if MCSamplesStitchOption == MCSamplesStitchOptions.PhSpOverlapRewgt:
            sampleInfo["hMCSamplesStitch"] = hMCSamplesStitch
    print(f"htoaa_triggerStudy_GGFMode:: here19 {datetime.now() = }")
        
    startTime = time.time()
    tracemalloc.start()


    #client = Client("tls://localhost:8786")
    #executor = processor.DaskExecutor(client=client)
    chunksize = nEventToReadInBatch
    maxchunks = None if nEventsToAnalyze == -1 else int(nEventsToAnalyze/nEventToReadInBatch)
    nWorkers  = 4 if nEventsToAnalyze == -1 else 1
    print(f"nEventsToAnalyze: {nEventsToAnalyze},  nEventToReadInBatch: {nEventToReadInBatch}, chunksize: {chunksize},  maxchunks: {maxchunks},  nWorkers: {nWorkers}")
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
        for key in output['cutflow'].keys():
            if key.startswith(sWeighted): continue # to print weighted and unweighted events for cuts on the same line

            print("%10f\t%10d\t%s" % (output['cutflow'][sWeighted+key], output['cutflow'][key], key), flush=flushStdout)
    
    
    if sOutputFile is not None:
        if not sOutputFile.endswith('.root'): sOutputFile += '.root'
        sample_category_toUse = sample_category
        
        if isMC and \
            MCSamplesStitchOption == MCSamplesStitchOptions.PhSpOverlapRewgt and \
            "QCD" in sample_category:
            sample_category_toUse = "QCD"
        
        sDir1 = 'evt/%s' % (sample_category_toUse)

        
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                sHistoName_toUse = key
                sHExt_toUse = ''
                if isMC and \
                    MCSamplesStitchOption == MCSamplesStitchOptions.PhSpOverlapRewgt and \
                    "QCD" in sample_category:                    
                    for sHExt in HistogramNameExtensions_QCD:
                        if sHExt in key:
                            sHExt_toUse = '_%s' % (sHExt)
                            sHistoName_toUse = sHistoName_toUse.replace(sHExt_toUse, '')
                            break

                sDir1_toUse = '%s%s' % (sDir1, sHExt_toUse)

                if not isinstance(value, hist.Hist): continue

                for _dataset in value.axis('dataset').identifiers():

                    for _syst in value.axis('systematic').identifiers():

                        h1 = value.integrate('dataset',_dataset).integrate('systematic',_syst).to_hist()

                        fOut['%s/%s_%s' % (sDir1_toUse, sHistoName_toUse, _syst)] = h1
            
        
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
    
