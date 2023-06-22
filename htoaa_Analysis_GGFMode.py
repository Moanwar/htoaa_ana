#htoaa analysis main code

import os
import sys
from datetime import datetime
print(f"htoaa_Analysis_GGFMode:: here1 {datetime.now() = }")
import subprocess
import json
from urllib.request import urlopen
import glob
from collections import OrderedDict as OD
import time
import tracemalloc
import math
print(f"htoaa_Analysis_GGFMode:: here2 {datetime.now() = }")
import numpy as np
from copy import copy, deepcopy
print(f"htoaa_Analysis_GGFMode:: here3 {datetime.now() = }")
#import uproot
#import uproot3 as uproot
import uproot as uproot
print(f"htoaa_Analysis_GGFMode:: here4 {datetime.now() = }")
#import parse
from parse import *
import logging
print(f"htoaa_Analysis_GGFMode:: here5 {datetime.now() = }")

# comment test3
'''
GGF -> H->aa->4b boosted analysis macro

References:
  * Coffea framework used for TTGamma analysis: https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py
* Coffea installation: /home/siddhesh/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea
'''
print(f"htoaa_Analysis_GGFMode:: here6 {datetime.now() = }")
#import coffea.processor as processor
from coffea import processor, util
from coffea.nanoevents import schemas
from coffea.nanoevents.methods import nanoaod, vector
from coffea.analysis_tools import PackedSelection, Weights
#from coffea.lumi_tools import LumiMask
#import hist
from coffea import hist
import awkward as ak
#import uproot
#from dask.distributed import Client
print(f"htoaa_Analysis_GGFMode:: here7 {datetime.now() = }")
from particle import Particle # For PDG particle listing https://github.com/scikit-hep/particle
print(f"htoaa_Analysis_GGFMode:: here8 {datetime.now() = }")


from htoaa_Settings import *
print(f"htoaa_Analysis_GGFMode:: here9 {datetime.now() = }")
from htoaa_CommonTools import (
    GetDictFromJsonFile, selectRunLuminosityBlock,
    calculate_lumiScale, getLumiScaleForPhSpOverlapRewgtMode, update_crosssection, getSampleHTRange,
    getNanoAODFile, setXRootDRedirector,  xrdcpFile,
    getHTReweight
)
print(f"htoaa_Analysis_GGFMode:: here10 {datetime.now() = }")
from htoaa_Samples import (
    kData, kQCD_bEnrich, kQCD_bGen, kQCDIncl
)
print(f"htoaa_Analysis_GGFMode:: here11 {datetime.now() = }")

from inspect import currentframe, getframeinfo
print(f"htoaa_Analysis_GGFMode:: here12 {datetime.now() = }")
frameinfo = getframeinfo(currentframe())
print(f"htoaa_Analysis_GGFMode:: here13 {datetime.now() = }")


# use GOldenJSON

 
printLevel = 0
nEventToReadInBatch = 0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze = -1 # 1000 # 100000 # -1
#pd.set_option('display.max_columns', None)

#print("".format())

sWeighted = "Wtd: "




    
# -----------------------------------------------------------------------------------
def get_GenPartDaughters(awkArray, index_GenPart):
    if printLevel >= 9:
        print(f"\n get_GenPartDaughters:: awkArray: {awkArray},   index_GenPart: {index_GenPart}"); sys.stdout.flush()
    
    return False

# -----------------------------------------------------------------------------------
def printVariable(sName, var):
    printInDetail=True
    #if nEventsToAnalyze == -1: printInDetail = False
    if str(type(var)) in ['numpy.ndarray', "<class 'numpy.ndarray'>"]: printInDetail = False # as gave error
    #print(f"printInDetail: {printInDetail} {sName} ({type(var)}) ({len(var)}): {var}")
    if not printInDetail:
        print(f"{sName} ({type(var)}) ({len(var)}): {var}")
    else:
        print(f"{sName} ({type(var)}) ({len(var)}): {var.to_list()}")

# -----------------------------------------------------------------------------------
class ObjectSelection:
    def __init__(self, era):
        self.era = era
        
        self.tagger_btagDeepB = 'DeepCSV'
        self.wp_btagDeepB = 'M'

        self.FatJetPtThsh  = 400 #170
        self.FatJetEtaThsh = 2.4

        self.FatJetMSoftDropThshLow  = 90
        self.FatJetMSoftDropThshHigh = 200

        self.FatJetParticleNetMD_Xbb_Thsh = 0.8
        
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
        if printLevel >= 15:
            print(f"\n maskGenHiggs:  {maskGenHiggs.to_list()} ")
            print(f"\n events.GenPart[maskGenHiggs]:  {events.GenPart[maskGenHiggs].to_list()} ")
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

        ak.behavior.update(nanoaod.behavior)

        self.datasetInfo = datasetInfo
        self.objectSelector = ObjectSelection(era=self.datasetInfo["era"])
        dataLSSelGoldenJSON = None
        if kData in self.datasetInfo:
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
            self.datasetInfo[kData]['dataLSSelGoldenJSON'] = dataLSSelGoldenJSON
            #print(f"{dataLSSelGoldenJSON = }")
                
        
        # set self.pdgId_BHadrons for 'QCD_bGenFilter' sample requirement
        self.pdgId_BHadrons = []
        bHadrons_ = Particle.findall(lambda p: p.pdgid.has_bottom) # Find all bottom hadrons
        print(f"List of B-hadrons for QCD B-GEN-filter ({len(bHadrons_)}):")
        print("%s %-20s %15s %15s" %(" "*4, "pdg name", "pdgId", "Mass in MeV"))
        for bHadron in bHadrons_:
            self.pdgId_BHadrons.append(bHadron.pdgid.abspid)
            print("%s %-20s %15d %15s" % (" "*4, str(bHadron), bHadron.pdgid.abspid, str(bHadron.mass)))
        print(f"self.pdgId_BHadrons ({len(self.pdgId_BHadrons)}): {self.pdgId_BHadrons}")
        self.pdgId_BHadrons = list(set(self.pdgId_BHadrons))
        print(f" after duplicate removal --> \nself.pdgId_BHadrons ({len(self.pdgId_BHadrons)}): {self.pdgId_BHadrons}")
        
        
        #dataset_axis = hist.axis.StrCategory(name="dataset", label="", categories=[], growth=True)
        #muon_axis = hist.axis.Regular(name="massT", label="Transverse Mass [GeV]", bins=50, start=15, stop=250)
        dataset_axis    = hist.Cat("dataset", "Dataset")
        systematic_axis = hist.Cat("systematic", "Systematic Uncertatinty")

        cutFlow_axis          = hist.Bin("CutFlow",                r"Cuts",                       21,    -0.5,    20.5)
        cutFlow50_axis        = hist.Bin("CutFlow50",              r"Cuts",                       51,    -0.5,    50.5)
        nObject_axis          = hist.Bin("nObject",                r"No. of object",              21,    -0.5,    20.5)
        nObject50_axis        = hist.Bin("nObject50",              r"No. of object",              51,    -0.5,    50.5)
        nObject200_axis       = hist.Bin("nObject200",             r"No. of object",             201,    -0.5,   200.5)
        pt_axis               = hist.Bin("Pt",                     r"$p_{T}$ [GeV]",             200,       0,    1000)
        ptLow_axis            = hist.Bin("PtLow",                  r"$p_{T}$ [GeV]",             400,       0,     200)
        ptUltraLow_axis         = hist.Bin("PtUltraLow",           r"$p_{T}$ [GeV]",             200,       0,     0.1)
        pt1to10_axis          = hist.Bin("Pt1to10",                r"$p_{T}$ [GeV]",             100,       0,      10)
        eta_axis              = hist.Bin("Eta",                    r"$#eta$",                    100,      -6,       6)
        phi_axis              = hist.Bin("Phi",                    r"$\phi$",                    100,   -3.14,    3.13)
        #mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       200, 0, 600)
        #mass_axis             = hist.Bin("Mass",      r"$m$ [GeV]",       400, 0, 200)
        mass_axis             = hist.Bin("Mass",                   r"$m$ [GeV]",                 300,       0,     300)
        mass_axis1            = hist.Bin("Mass1",                  r"$m$ [GeV]",                 300,       0,     300)
        mlScore_axis          = hist.Bin("MLScore",                r"ML score",                  100,    -1.1,     1.1)
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
        
        sXaxis      = 'xAxis'
        sXaxisLabel = 'xAxisLabel'
        sYaxis      = 'yAxis'
        sYaxisLabel = 'yAxisLabel'
        histos = OD([
            # ('histogram_name',  {sXaxis: hist.Bin() axis,  sXaxisLabel: "histogram axis label"})
            ('hCutFlow',                                  {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),
            ('hCutFlowWeighted',                          {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),
            ('hNEventsQCD',                               {sXaxis: cutFlow50_axis,  sXaxisLabel: 'Cuts'}),
            ('hNEventsQCDUnweighted',                     {sXaxis: cutFlow50_axis,  sXaxisLabel: 'Cuts'}),

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
            
            ('nSelFatJet',                                {sXaxis: nObject_axis,    sXaxisLabel: 'No. of selected FatJets'}),
            ('hLeadingFatJetPt',                          {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(leading FatJet)$ [GeV]"}),
            ('hLeadingFatJetEta',                         {sXaxis: eta_axis,        sXaxisLabel: r"\eta (leading FatJet)"}),
            ('hLeadingFatJetPhi',                         {sXaxis: phi_axis,        sXaxisLabel: r"\phi (leading FatJet)"}),
            ('hLeadingFatJetMass',                        {sXaxis: mass_axis,       sXaxisLabel: r"m (leading FatJet) [GeV]"}),
            ('hLeadingFatJetMSoftDrop',                   {sXaxis: mass_axis,       sXaxisLabel: r"m_{soft drop} (leading FatJet) [GeV]"}),
            
            ('hLeadingFatJetBtagDeepB',                   {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagDeepB"}),
            ('hLeadingFatJetBtagDDBvLV2',                 {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagDDBvLV2"}),

            ('hLeadingFatJetBtagDDCvBV2',                 {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagDDCvBV2"}),
            ('hLeadingFatJetBtagHbb',                     {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetBtagHbb"}),
            ('hLeadingFatJetDeepTagMD_H4qvsQCD',          {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_H4qvsQCD"}),
            ('hLeadingFatJetDeepTagMD_HbbvsQCD',          {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_HbbvsQCD"}),
            ('hLeadingFatJetDeepTagMD_ZHbbvsQCD',         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHbbvsQCD"}),
            ('hLeadingFatJetDeepTagMD_ZHccvsQCD',         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZHccvsQCD"}),
            
            ('hLeadingFatJetDeepTagMD_ZbbvsQCD',          {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetdeepTagMD_ZbbvsQCD"}),
            ('hLeadingFatJetDeepTagMD_ZvsQCD',            {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ZvsQCD"}),
            ('hLeadingFatJetDeepTagMD_bbvsLight',         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_bbvsLight"}),
            ('hLeadingFatJetDeepTagMD_ccvsLight',         {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTagMD_ccvsLight"}),
            ('hLeadingFatJetDeepTag_H',                   {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTag_H"}),
            ('hLeadingFatJetDeepTag_QCD',                 {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTag_QCD"}),
            ('hLeadingFatJetDeepTag_QCDothers',           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetDeepTag_QCDothers"}),
            
            ('hLeadingFatJetN2b1',                        {sXaxis: jetN2_axis,      sXaxisLabel: r"LeadingFatJetn2b1"}),
            ('hLeadingFatJetN3b1',                        {sXaxis: jetN3_axis,      sXaxisLabel: r"LeadingFatJetn3b1"}),
            ('hLeadingFatJetTau1',                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau1"}),
            ('hLeadingFatJetTau2',                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau2"}),
            ('hLeadingFatJetTau3',                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau3"}),
            ('hLeadingFatJetTau4',                        {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau4"}),

            ('hLeadingFatJetTau4by3',                     {sXaxis: jetTau_axis,     sXaxisLabel: r"LeadingFatJetTau4by3"}),
            ('hLeadingFatJetTau3by2',                     {sXaxis: jetTau_axis,     sXaxisLabel: r"hLeadingFatJetTau3by2"}),
            ('hLeadingFatJetTau2by1',                     {sXaxis: jetTau_axis,     sXaxisLabel: r"hLeadingFatJetTau2by1"}),
            
            ('hLeadingFatJetNBHadrons',                   {sXaxis: nObject_axis,    sXaxisLabel: r"LeadingFatJetNBHadrons"}),
            ('hLeadingFatJetNCHadrons',                   {sXaxis: nObject_axis,    sXaxisLabel: r"LeadingFatJetNCHadrons"}),            
            ('hLeadingFatJetNConstituents',               {sXaxis: nObject_axis,    sXaxisLabel: r"LeadingFatJetNConstituents"}),
            
            ('hLeadingFatJetParticleNetMD_QCD',           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_QCD"}),
            ('hLeadingFatJetParticleNetMD_Xbb',           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_Xbb"}),
            ('hLeadingFatJetParticleNetMD_Xcc',           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_Xcc"}),

            ('hLeadingFatJetParticleNetMD_Xqq',           {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNetMD_Xqq"}),
            ('hLeadingFatJetParticleNet_H4qvsQCD',        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_H4qvsQCD"}),
            ('hLeadingFatJetParticleNet_HbbvsQCD',        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_HbbvsQCD"}),
            ('hLeadingFatJetParticleNet_HccvsQCD',        {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_HccvsQCD"}),
            ('hLeadingFatJetParticleNet_QCD',             {sXaxis: mlScore_axis,    sXaxisLabel: r"LeadingFatJetParticleNet_QCD"}),
            
            ('hLeadingFatJetParticleNet_mass',            {sXaxis: mass_axis,       sXaxisLabel: r"LeadingFatJetParticleNet_mass"}),

            ('hGenHiggsPt_all',                           {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(GEN Higgs (pdgId: 25, status=62))$ [GeV]"}),
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

            
        ])

        for statusFlag_ in GENPART_STATUSFLAGS:
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

        if printLevel >= 20:
            print(f"nEvents: {len(events)}")
        if printLevel >= 0:
            print(f"\n events.fields ({type(events.fields)}): {events.fields}")
            #print(f"\n events.GenPart.fields: {events.GenPart.fields}")
            #print(f"\n events.HLT.fields: {events.HLT.fields}")
            #printVariable('\n events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4', events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4)
            #print(f"\n events.L1.fields: {events.L1.fields}")
            #printVariable('\n events.L1.SingleJet180', events.L1.SingleJet180)
            #print(f"\n events.FatJet.fields: {events.FatJet.fields}")
            #print(f"\n events.LHE.fields: {events.LHE.fields}")
            #print(f"\n events.LHE.HT: {events.LHE.HT.to_list()}")
            #print(f"{events.LHE.fields = } ")
            #print(f"{events.LHE = } ")
            #printVariable('events.LHE', events.LHE)
            print(f"{events.run.fields = }")
            printVariable('events.run', events.run)
            print(f"{events.luminosityBlock.fields = }")
            printVariable('events.luminosityBlock', events.luminosityBlock)
             
        if nEventsToAnalyze != -1:
            print(f"\n (run:ls:event): {ak.zip([events.run, events.luminosityBlock, events.event])}")

        #print(f"htoaa_Analysis_GGFMode.py::process():: {self.datasetInfo = }"); sys.stdout.flush()

        self.datasetInfo[dataset]['isSignal'     ]  = False
        self.datasetInfo[dataset]['isQCD'        ]  = False
        self.datasetInfo[dataset]['isQCDIncl'    ]  = False
        self.datasetInfo[dataset]['isQCD_bEnrich']  = False
        self.datasetInfo[dataset]['isQCD_bGen'   ]  = False
        
        if self.datasetInfo[dataset]['isMC']:
            self.datasetInfo[dataset]['isSignal']      = True if "HToAATo4B"   in dataset else False
            self.datasetInfo[dataset]['isQCDIncl']     = True if kQCDIncl      in dataset else False
            self.datasetInfo[dataset]['isQCD_bEnrich'] = True if kQCD_bEnrich  in dataset else False
            self.datasetInfo[dataset]['isQCD_bGen']    = True if kQCD_bGen     in dataset else False
            self.datasetInfo[dataset]['isQCD']         = (self.datasetInfo[dataset]['isQCDIncl']     or \
                                                          self.datasetInfo[dataset]['isQCD_bEnrich'] or \
                                                          self.datasetInfo[dataset]['isQCD_bGen'])
            sample_HT_Min, sample_HT_Max = getSampleHTRange( self.datasetInfo[dataset]["datasetNameFull"] )
            self.datasetInfo[dataset]['sample_HT_Min'] = sample_HT_Min
            self.datasetInfo[dataset]['sample_HT_Max'] = sample_HT_Max

            if self.datasetInfo[dataset]['isQCD_bGen']:
                # 'Corrections' variable defined in htoaa_Settings.py
                fitFunctionFormat_  = Corrections["HTRewgt"]["QCD_bGen"][self.datasetInfo["era"]]["FitFunctionFormat"] 
                fitFunctionHTRange_ = ""
                for sHTBin in Corrections["HTRewgt"]["QCD_bGen"][self.datasetInfo["era"]]:
                    if "HT%dto" % (self.datasetInfo[dataset]['sample_HT_Min']) in sHTBin:
                        fitFunctionHTRange_ = sHTBin
                        fitFunction_  = Corrections["HTRewgt"]["QCD_bGen"][self.datasetInfo["era"]][sHTBin]

                        
                self.datasetInfo[dataset]['HTRewgt'] = {
                    "fitFunctionFormat":  fitFunctionFormat_,
                    "fitFunction":        fitFunction_,
                    "fitFunctionHTRange": fitFunctionHTRange_,                    
                }

               
            
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
        #print(f"dataset: {dataset}")

        

        ones_list  = np.ones(len(events))
        trues_list = np.ones(len(events), dtype=bool)


          
        ##################
        # OBJECT SELECTION
        ##################


        # Gen-level selection ---------------------------------------------------------------------
        genHiggs = None
        genHT    = None
        if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isSignal']: 
            genHiggs  = self.objectSelector.selectGenHiggs(events)        
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

            # LorentVector of GenB quarks from HToAATo4b
            nEvents_11 = ak.num(events.GenPart[genBBar_pairs['b']][:, 0].pt, axis=0)
            mass_bQuark = 4.18
            #print(f"\n np.full(nEvents_11, mass_bQuark): {np.full(nEvents_11, mass_bQuark)}")

            # https://coffeateam.github.io/coffea/modules/coffea.nanoevents.methods.vector.html
            # bQuark from 1st A
            LVGenB_0 = ak.zip(
                {
                    "pt"  : events.GenPart[genBBar_pairs['b']][:, 0].pt,
                    "eta" : events.GenPart[genBBar_pairs['b']][:, 0].eta,
                    "phi" : events.GenPart[genBBar_pairs['b']][:, 0].phi,
                    "mass": np.full(nEvents_11, mass_bQuark),
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
                    "mass": np.full(nEvents_11, mass_bQuark),
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
                    "mass": np.full(nEvents_11, mass_bQuark),
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
                    "mass": np.full(nEvents_11, mass_bQuark),
                },
                with_name="PtEtaPhiMLorentzVector",
                behavior=vector.behavior,
            )            

            dr_GenH_GenB = ak.concatenate([genHiggs.delta_r(LVGenB_0), genHiggs.delta_r(LVGenBbar_0), genHiggs.delta_r(LVGenB_1), genHiggs.delta_r(LVGenBbar_1)], axis=-1)
            max_dr_GenH_GenB = ak.max(dr_GenH_GenB, axis=-1)            

                
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
        if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isQCD'] :
            mask_genLHEHTLt100 = (events.LHE.HT < 100)

            if printLevel >= 2:
                printVariable('\n events.LHE.HT', events.LHE.HT); sys.stdout.flush()
                printVariable('\n mask_genLHEHTLt100', mask_genLHEHTLt100); sys.stdout.flush()
            
            
            genBQuarks = events.GenPart[(
                (abs(events.GenPart.pdgId) == 5)
            )]
            genBQuarks_pT = ak.sort(genBQuarks.pt, axis=-1, ascending=False)
            genBQuarks_first = ak.firsts(genBQuarks)
            #genBQuarks_first = ak.firsts(genBQuarks, axis=-1)
            mask_genBQuarks = (ak.count(genBQuarks.pdgId, axis=1) >= 1)

            mask_genBQuarks_pTAbvTrsh = ak.any((genBQuarks.pt > 15.0), axis=1)

            idx_genBQuarks_pTsort = ak.argsort(genBQuarks.pt, axis=-1, ascending=False)
            
            
            if printLevel >= 2:
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
           
                    
            if printLevel >= 2:
                #genBHadrons_status2 = events.GenPart[mask_genBHadrons_status2]
                printVariable('\n genBHadrons_status2', genBHadrons_status2); sys.stdout.flush()
                printVariable(' genBHadrons_status2.pdgId', genBHadrons_status2.pdgId); sys.stdout.flush()
                printVariable(' ak.any(mask_genBHadrons_status2, axis=1)', ak.any(mask_genBHadrons_status2, axis=1)); sys.stdout.flush()


            # Check if events has b-quark outgoing from hard subprocess -----------------------------------------------
            mask_genBQuarks_hardSctred = (
                (abs(events.GenPart.pdgId) == 5) &
                (events.GenPart.status == 23)
            )
            mask_genBQuarks_hardSctred_eventwise = ak.any(mask_genBQuarks_hardSctred, axis=1)

            genBQuarks_hardSctred = events.GenPart[mask_genBQuarks_hardSctred]
            idx_genBQuarks_hardSctred_pTsort = ak.argsort(genBQuarks_hardSctred.pt, axis=-1, ascending=False)
            
            if printLevel >= 2:
                #genBQuarks_hardSctred = events.GenPart[mask_genBQuarks_hardSctred]
                printVariable('\n genBQuarks_hardSctred', genBQuarks_hardSctred); sys.stdout.flush()
                printVariable(' genBQuarks_hardSctred.hasFlags("isHardProcess")', genBQuarks_hardSctred.hasFlags("isHardProcess")); sys.stdout.flush()
                printVariable(' genBQuarks_hardSctred.hasFlags("fromHardProcess")', genBQuarks_hardSctred.hasFlags("fromHardProcess")); sys.stdout.flush()
                printVariable(' ak.any(mask_genBQuarks_hardSctred, axis=1)', ak.any(mask_genBQuarks_hardSctred, axis=1)); sys.stdout.flush()

                printVariable(' genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort].pt', genBQuarks_hardSctred[idx_genBQuarks_hardSctred_pTsort].pt); sys.stdout.flush()

            # QCD stitch cut-based: conditions -----------------------------------------------------------------------------------
            # option 1: GEN b-quark pT > 15 GeV for QCD BEnrich and QCD bGen samples.
            if self.datasetInfo[dataset]['isQCD_bEnrich'] or self.datasetInfo[dataset]['isQCD_bGen']:
                mask_QCD_stitch_CutBQuarkPt_eventwise = mask_genBQuarks_pTAbvTrsh
            elif self.datasetInfo[dataset]['isQCDIncl']:
                mask_QCD_stitch_CutBQuarkPt_eventwise = (
                    (mask_genBQuarks_pTAbvTrsh == False) |
                    (mask_genLHEHTLt100 == True)
                )
                
            if printLevel >= 2:
                printVariable('\n mask_QCD_stitch_CutBQuarkPt_eventwise', mask_QCD_stitch_CutBQuarkPt_eventwise); sys.stdout.flush()


            # Option 2: B-Hadron for QCD bGEN. b-quark from hard subprocess for QCD bEnrich
            if self.datasetInfo[dataset]['isQCD_bEnrich']:
                mask_QCD_stitch_CutBHadron_eventwise = trues_list
            elif self.datasetInfo[dataset]['isQCD_bGen']:
                mask_QCD_stitch_CutBHadron_eventwise = (
                    (mask_genBQuarks_hardSctred_eventwise == False)
                )
            elif self.datasetInfo[dataset]['isQCDIncl']:
                mask_QCD_stitch_CutBHadron_eventwise = (
                    (
                        (mask_genBQuarks_hardSctred_eventwise == False) &
                        (mask_genBHadrons_status2_eventwise == False)
                    ) |
                    (
                        (mask_genLHEHTLt100 == True)
                    )
                )
                
            if printLevel >= 2:
                printVariable('\n mask_QCD_stitch_CutBHadron_eventwise', mask_QCD_stitch_CutBHadron_eventwise); sys.stdout.flush()

            if printLevel >= 10:
                mask_tmp = (ak.count(genBHadrons_status2.pt, axis=-1) >= 4)
                #mask_tmp1 = (genBHadrons_status2[idx_genBHadrons_status2_pTsort].pt < 0.01)
                
                #printVariable('\n (genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp)].pt)', (genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp)].pt)); sys.stdout.flush()
                printVariable('\n (genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp)][:, 3].pt)', (genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp)][:, 2].pt)); sys.stdout.flush()

                genBHadrons_status2_fourth_pt = genBHadrons_status2[idx_genBHadrons_status2_pTsort][(mask_tmp)][:, 2].pt
                mask_tmp1_ = (genBHadrons_status2_fourth_pt < 0.01)
                printVariable('\n genBHadrons_status2_fourth_pt[mask_tmp1_]', genBHadrons_status2_fourth_pt[mask_tmp1_]); sys.stdout.flush()

            if printLevel >= 3:
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

            if printLevel >= 1:
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
                
            if self.datasetInfo[dataset]["MCSamplesStitchOption"] == MCSamplesStitchOptions.PhSpOverlapRewgt:
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
            # --------------------------------------------------------------------------------------------------

                    
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
        
        leadingFatJet = None
        if ( not self.datasetInfo[dataset]['isMC']) and (self.datasetInfo["era"] == Era_2018):
            '''
            HEM15/16 issue in jets with -3.2<eta<-1.3 and -1.57<phi< -0.87 in 2018 data (runs>=319077, i.e. last certified run of 2018B, and all of 2018C+D) https://twiki.cern.ch/twiki/bin/view/CMS/JetMET#
            '''
            run_FatJetEta_FatJetPhi = ak.zip({'run': events.run, 'FatJetEta': events.FatJet.eta, 'FatJetPhi': events.FatJet.phi})   
            mask_jets_surviving_HEM15_16_issue = ~ (
                (run_FatJetEta_FatJetPhi.run >= 319077) &
                (run_FatJetEta_FatJetPhi.FatJetEta > -3.2 ) & (run_FatJetEta_FatJetPhi.FatJetEta < -1.3 ) &
                (run_FatJetEta_FatJetPhi.FatJetPhi > -1.57) & (run_FatJetEta_FatJetPhi.FatJetPhi < -0.87)
            )
            leadingFatJet = ak.firsts(events.FatJet[mask_jets_surviving_HEM15_16_issue]) 
        else:   
            leadingFatJet = ak.firsts(events.FatJet) # for e.g. [0.056304931640625, None, 0.12890625, 0.939453125, 0.0316162109375]


        leadingFatJet_asSingletons = ak.singletons(leadingFatJet) # for e.g. [[0.056304931640625], [], [0.12890625], [0.939453125], [0.0316162109375]]
        
        if printLevel >= 13:
            #printVariable("\n ", )
            printVariable("\n ak.firsts(events.FatJet)", ak.firsts(events.FatJet))
            printVariable("\n ak.singletons( ak.firsts(events.FatJet) )", ak.singletons( ak.firsts(events.FatJet) ))

            printVariable("\n leadingFatJet.btagDeepB", leadingFatJet.btagDeepB)
            printVariable("\n ak.singletons(leadingFatJet.btagDeepB)", ak.singletons(leadingFatJet.btagDeepB))
            printVariable("\n leadingFatJet.btagDeepB > bTagWPs[self.objectSelector.era][self.objectSelector.tagger_btagDeepB][self.objectSelector.wp_btagDeepB]", leadingFatJet.btagDeepB > bTagWPs[self.objectSelector.era][self.objectSelector.tagger_btagDeepB][self.objectSelector.wp_btagDeepB])

        #if printLevel >= 2:
        #    print(f" : {}")
            
        #####################
        # EVENT SELECTION
        #####################

        HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4" 
        
        # sel_names_all = dict of {"selection name" : [list of different cuts]}; for cut-flow table
        sel_names_all = OD([
            ("SR",                    [
                "nPV",
                "leadingFatJetPt",
                "leadingFatJetEta",
                #"leadingFatJetBtagDeepB",
                "leadingFatJetMSoftDrop",
                "leadingFatJetParticleNetMD_Xbb",
                "L1_SingleJet180",
                HLT_AK8PFJet330_name
            ]),
        ])
        if not self.datasetInfo[dataset]['isMC']: 
            sel_names_all["SR"].insert(0, "run:ls")
        else:
            if self.datasetInfo[dataset]['isQCD']: sel_names_all["SR"].append("QCDStitch")
        # reconstruction level cuts for cut-flow table. Order of cuts is IMPORTANT
        cuts_reco = ["dR_LeadingFatJet_GenB_0p8"] + sel_names_all["SR"] #.copy()

       
        # create a PackedSelection object
        # this will help us later in composing the boolean selections easily
        selection = PackedSelection()

        if "run:ls" in sel_names_all["SR"]:
            # self.datasetInfo[dataset]['dataLSSelGoldenJSON']
            # using Coffea built-in function: mask_lumi = LumiMask(golden_json_path)(events.run,events.luminosityBlock)
            #selection.add("run:ls", LumiMask(sFilesGoldenJSON[self.datasetInfo["era"]])(events.run,events.luminosityBlock) )
            
            # using selectRunLuminosityBlock function from htoaa_CommonTools
            selection.add("run:ls", selectRunLuminosityBlock(
                dataLSSelGoldenJSON  = self.datasetInfo[dataset]['dataLSSelGoldenJSON'], 
                runNumber_list       = events.run, 
                luminosityBlock_list = events.luminosityBlock 
                ))

        if "nPV" in sel_names_all["SR"]:
            # nPVGood >= 1
            selection.add("nPV", events.PV.npvsGood >= 1)

        

        if "leadingFatJetPt" in sel_names_all["SR"]:
            # >=1 FatJet
            #selection.add("FatJetGet", ak.num(selFatJet) >= self.objectSelector.nFatJetMin)
            selection.add(
                "leadingFatJetPt",
                leadingFatJet.pt > self.objectSelector.FatJetPtThsh
            )


        if "leadingFatJetEta" in sel_names_all["SR"]:
            selection.add(
                "leadingFatJetEta",
                abs(leadingFatJet.eta) < self.objectSelector.FatJetEtaThsh
            )


        if "leadingFatJetBtagDeepB" in sel_names_all["SR"]:
            selection.add(
                "leadingFatJetBtagDeepB",
                leadingFatJet.btagDeepB > bTagWPs[self.objectSelector.era][self.objectSelector.tagger_btagDeepB][self.objectSelector.wp_btagDeepB]
            )

 
        if"leadingFatJetMSoftDrop"  in sel_names_all["SR"]:
           selection.add(
                "leadingFatJetMSoftDrop",
                (leadingFatJet.msoftdrop > self.objectSelector.FatJetMSoftDropThshLow) &
                (leadingFatJet.msoftdrop < self.objectSelector.FatJetMSoftDropThshHigh)
            )

 
        if "leadingFatJetParticleNetMD_Xbb" in sel_names_all["SR"]:
            selection.add(
                "leadingFatJetParticleNetMD_Xbb",
                leadingFatJet.particleNetMD_Xbb > self.objectSelector.FatJetParticleNetMD_Xbb_Thsh
            )
            
            
        if "L1_SingleJet180" in sel_names_all["SR"]:
            selection.add(
                "L1_SingleJet180",
                events.L1.SingleJet180 == True
            )

 
        if HLT_AK8PFJet330_name in sel_names_all["SR"]:
            # some files of Run2018A do not have HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 branch
            #HLT_AK8PFJet330_name = None
            if "AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4" in events.HLT.fields:
                #HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4"
                selection.add(
                    HLT_AK8PFJet330_name,
                    events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 == True
                )
            elif "AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02" in events.HLT.fields:
                #HLT_AK8PFJet330_name = "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02"
                selection.add(
                    HLT_AK8PFJet330_name,
                    events.HLT.AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02 == True
                )
        

        if "QCDStitch" in sel_names_all["SR"]:
            selection.add(
                "QCDStitch",
                #mask_QCD_stitch_eventwise == True
                mask_QCD_stitch_eventwise
            )
            
            

            
        #sel_SR          = selection.all("nPV", "FatJetGet")
        sel_SR           = selection.all(* sel_names_all["SR"])
        sel_GenHToAATo4B = None

        if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isSignal']:
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
            sel_names_all.update( OD([
                ("GenHToAATo4B_1", ["1GenHiggs", "2GenA", "2GenAToBBbarPairs"]),
                ("GenHToAATo4B", [*sel_names_GEN]),
            ]) )
            for idx, cutName in enumerate(cuts_reco):
                if idx == 0:
                    sel_names_all.update( OD([
                        ("SR_%d" % (idx+1),  [*sel_names_GEN, cutName]),
                    ]) )
                else:
                    sel_names_all.update( OD([
                        ("SR_%d" % (idx+1),  [*sel_names_all["SR_%d" % (idx)], cutName]),
                    ]) ) 
            '''
            sel_names_all.update( OD([
                ("GenHToAATo4B_1", ["1GenHiggs", "2GenA", "2GenAToBBbarPairs"]),
                ("GenHToAATo4B", [*sel_names_GEN]),

                ("SR_1",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8"]),
                ("SR_2",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt"]),
                ("SR_3",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt", "FatJetEta"]),
                ("SR_4",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt", "FatJetEta", "FatJetBtagDeepB"]),
                ("SR_5",  [*sel_names_GEN, "nPV", "dR_LeadingFatJet_GenB_0p8", "FatJetPt", "FatJetEta", "FatJetBtagDeepB", "FatJetMSoftDrop"]),
            ]) )
            '''
            
            sel_GenHToAATo4B = selection.all(* sel_names_all["GenHToAATo4B"])

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
                #wgt1=np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
                #print(f"wgt1 ({len(wgt1)}): {wgt1}")
                

        if printLevel >= 15:
            print(f"selFatJet.fields: {selFatJet.fields}")
            print(f"sel_SR ({len(sel_SR)}): {sel_SR}")
            print(f"selFatJet.pt[sel_SR].to_list(): {selFatJet.pt[sel_SR].to_list()} ")

        if printLevel >= 30:
            #print(f" ({type()}) ({len()}): {} \n")
            print(f"leadingFatJet ({type(leadingFatJet)}) ({len(leadingFatJet)}): {leadingFatJet} \n")
            print(f"leadingFatJet_asSingletons ({type(leadingFatJet_asSingletons)}) ({len(leadingFatJet_asSingletons)}): {leadingFatJet_asSingletons} \n")
            print(f"sel_SR ({type(sel_SR)}) ({len(sel_SR)}): {sel_SR} \n")
            print(f"sel_GenHToAATo4B ({type(sel_GenHToAATo4B)}) ({len(sel_GenHToAATo4B)}): {sel_GenHToAATo4B} \n")
            print(f"leadingFatJet ({type(leadingFatJet)}) ({len(leadingFatJet)}): {leadingFatJet} \n")
            print(f"leadingFatJet.pt[sel_SR] ({type(leadingFatJet.pt[sel_SR])}) ({len(leadingFatJet.pt[sel_SR])}): {leadingFatJet.pt[sel_SR]} \n")
            #print(f" ({type()}) ({len()}): {} \n")
                    
        
        if printLevel >= 30:
            #printVariable("", )
            printVariable("\n events.FatJet.pt", events.FatJet.pt)
            printVariable("\n mask_FatJetPt", mask_FatJetPt)
            printVariable("\n events.FatJet[mask_FatJetPt].pt", events.FatJet[mask_FatJetPt].pt)            
            printVariable("\n ak.num(events.FatJet[mask_FatJetPt])", ak.num(events.FatJet[mask_FatJetPt]))
            printVariable("\n ak.num(events.FatJet[mask_FatJetPt]) >= self.objectSelector.nFatJetMin", ak.num(events.FatJet[mask_FatJetPt]) >= self.objectSelector.nFatJetMin)
            print(f'\n selection.all("FatJetPt") ({type(selection.all("FatJetPt"))}) ({len(selection.all("FatJetPt"))}): {selection.all("FatJetPt")}')

            printVariable("\n\n events.FatJet.eta", events.FatJet.eta)
            printVariable("\n mask_FatJetEta", mask_FatJetEta)
            printVariable("\n events.FatJet[mask_FatJetEta].eta", events.FatJet[mask_FatJetEta].eta)
            printVariable("\n ak.num(events.FatJet[mask_FatJetEta]) >= self.objectSelector.nFatJetMin", ak.num(events.FatJet[mask_FatJetEta]) >= self.objectSelector.nFatJetMin)
            #printVariable('\n selection.all("FatJetEta")', selection.all("FatJetEta"))            
            print(f'\n selection.all("FatJetEta") ({type(selection.all("FatJetEta"))}) ({len(selection.all("FatJetEta"))}): {selection.all("FatJetEta")}')
            

            printVariable("\n events.FatJet.btagDeepB", events.FatJet.btagDeepB)
            printVariable("\n mask_FatJetBtagDeepB", mask_FatJetBtagDeepB)
            printVariable("\n events.FatJet[mask_FatJetBtagDeepB].btagDeepB", events.FatJet[mask_FatJetBtagDeepB].btagDeepB)
            printVariable("\n ak.num(events.FatJet[mask_FatJetBtagDeepB])", ak.num(events.FatJet[mask_FatJetBtagDeepB]))
            printVariable("\n ak.num(events.FatJet[mask_FatJetBtagDeepB]) >= self.objectSelector.nFatJetMin ", ak.num(events.FatJet[mask_FatJetBtagDeepB]) >= self.objectSelector.nFatJetMin )
            print(f'\n selection.all("FatJetBtagDeepB") ({type(selection.all("FatJetBtagDeepB"))}) ({len(selection.all("FatJetBtagDeepB"))}): {selection.all("FatJetBtagDeepB")}')

            printVariable("\n events.FatJet.msoftdrop", events.FatJet.msoftdrop)
            printVariable("\n mask_FatJetMSoftDrop", mask_FatJetMSoftDrop)
            printVariable("\n events.FatJet[mask_FatJetMSoftDrop].msoftdrop", events.FatJet[mask_FatJetMSoftDrop].msoftdrop)
            printVariable("\n ak.num(events.FatJet[mask_FatJetMSoftDrop])", ak.num(events.FatJet[mask_FatJetMSoftDrop]))
            printVariable("\n ak.num(events.FatJet[mask_FatJetMSoftDrop]) >= self.objectSelector.nFatJetMin ", ak.num(events.FatJet[mask_FatJetMSoftDrop]) >= self.objectSelector.nFatJetMin )
            print(f'\n selection.all("FatJetMSoftDrop") ({type(selection.all("FatJetMSoftDrop"))}) ({len(selection.all("FatJetMSoftDrop"))}): {selection.all("FatJetMSoftDrop")}')


            print(f"\n\nsel_names_all: {sel_names_all}")
            print(f'\n selection.all(* sel_names_all["SR"]) ({type(selection.all(* sel_names_all["SR"]))}) ({len(selection.all(* sel_names_all["SR"]))}): {selection.all(* sel_names_all["SR"])}')
            
            print(f'\n selection.all(* sel_names_all["GenHToAATo4B"]) ({type(selection.all(* sel_names_all["GenHToAATo4B"]))}) ({len(selection.all(* sel_names_all["GenHToAATo4B"]))}): {selection.all(* sel_names_all["GenHToAATo4B"])}')
            
            printVariable("\n events.FatJet.pt[sel_SR]", events.FatJet.pt[sel_SR])
            printVariable("\n selFatJet.pt", selFatJet.pt)
            printVariable("\n selFatJet.pt[sel_SR]", selFatJet.pt[sel_SR])

            printVariable("\n selFatJet.pt[mask_FatJetPt]", selFatJet.pt[mask_FatJetPt])

            printVariable("\n mask_FatJetPt", mask_FatJetPt)
            printVariable("\n mask_FatJetEta", mask_FatJetEta)
            printVariable("\n mask_FatJetBtagDeepB", mask_FatJetBtagDeepB)
            printVariable("\n mask_FatJetMSoftDrop", mask_FatJetMSoftDrop)
            printVariable("\n mask_all &:", (
                mask_FatJetPt &
                mask_FatJetEta &
                mask_FatJetBtagDeepB &
                mask_FatJetMSoftDrop
            ))
            #printVariable("\n ", )

            
            
            
            '''
            printVariable("events.FatJet.btagDeepB", events.FatJet.btagDeepB)
            printVariable("", )
            printVariable("", )

            printVariable("mask_FatJetBtagDeepB", mask_FatJetBtagDeepB)
            printVariable("", )
            '''

                    
        


            
            


        
        ################
        # EVENT WEIGHTS
        ################

        # create a processor Weights object, with the same length as the number of events in the chunk
        weights     = Weights(len(events))
        weights_gen = Weights(len(events))
        weights_GenHToAATo4B = Weights(len(events))
        

        if self.datasetInfo[dataset]["isMC"]:
            # lumiScale --
            lumiScale_toUse = None
            if self.datasetInfo[dataset]["MCSamplesStitchOption"] == MCSamplesStitchOptions.PhSpOverlapRewgt and \
               self.datasetInfo[dataset]['isQCD']:
                mask_PhSp_dict_ = {
                    "QCD_bEnrich": mask_QCD_bEnrich_PhSp,
                    "QCD_bGen": mask_QCD_bGen_PhSp,
                    "QCD_Incl_Remnant": mask_QCD_Incl_Remnant_PhSp,
                }
                lumiScale_toUse = getLumiScaleForPhSpOverlapRewgtMode(
                    hLumiScale      = self.datasetInfo[dataset]["hMCSamplesStitch"],
                    sample_category = dataset,
                    sample_HT_value = self.datasetInfo[dataset]['sample_HT_Min'],
                    mask_PhSp_dict  = mask_PhSp_dict_ )
            else:
                lumiScale_toUse = np.full(len(events), self.datasetInfo[dataset]["lumiScale"])

            
            weights.add(
                "lumiWeight",
                weight = lumiScale_toUse
            )
            weights.add(
                "genWeight",
                weight = np.copysign(np.ones(len(events)), events.genWeight)
            )
            weights.add(
                "btagWeight",
                weight = (events.btagWeight.DeepCSVB)
            )
            
 
            weights_gen.add(
                "lumiWeight",
                weight = lumiScale_toUse 
            )
            weights_gen.add(
                "genWeight",
                weight=np.copysign(np.ones(len(events)), events.genWeight)
            )
            
            
            if printLevel >= 5:
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



            if self.datasetInfo[dataset]['isQCD_bGen']:
                wgt_HT = getHTReweight(
                    HT_list            = events.LHE.HT,
                    sFitFunctionFormat = self.datasetInfo[dataset]['HTRewgt']["fitFunctionFormat"],
                    sFitFunction       = self.datasetInfo[dataset]['HTRewgt']["fitFunction"],
                    sFitFunctionRange  = self.datasetInfo[dataset]['HTRewgt']["fitFunctionHTRange"]
                )
                weights.add(
                    "HTRewgt",
                    weight=wgt_HT
                )

                weights_gen.add(
                    "HTRewgt",
                    weight=wgt_HT
                )
               
                if printLevel >= 30:
                    printVariable("\n events.LHE.HT", events.LHE.HT)
                    printVariable("wgt_HT", wgt_HT)


                
        ###################
        # FILL HISTOGRAMS
        ###################

        systList = []
        if self.datasetInfo[dataset]['isMC']:
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

        for iSelection in sel_names_all.keys():
            iName = f"{iSelection}: {sel_names_all[iSelection]}"
            sel_i = selection.all(* sel_names_all[iSelection])
            output['cutflow'][iName] += sel_i.sum()
            output['cutflow'][sWeighted+iName] +=  weights.weight()[sel_i].sum()
            

        for syst in systList:

            # find the event weight to be used when filling the histograms
            weightSyst = syst
            
            # in the case of 'central', or the jet energy systematics, no weight systematic variation is used (weightSyst=None)
            if syst in ["central", "JERUp", "JERDown", "JESUp", "JESDown"]:
                weightSyst = None

            
            
            if syst == "noweight":
                evtWeight = np.ones(len(events))
            else:
                evtWeight     = weights.weight(weightSyst)
                evtWeight_gen = weights_gen.weight(weightSyst)


            if printLevel >= 15:
                print(f"ak.num(selFatJet[sel_SR]) ({len(ak.num(selFatJet[sel_SR]))}): {ak.num(selFatJet[sel_SR])}")
                print(f"ak.to_numpy(ak.num(selFatJet[sel_SR])) ({len(ak.to_numpy(ak.num(selFatJet[sel_SR])))}): {ak.to_numpy(ak.num(selFatJet[sel_SR]))}")
                print(f"selFatJet.pt[sel_SR] ({len(selFatJet.pt[sel_SR])}): {selFatJet.pt[sel_SR]}")
                print(f"ak.flatten(selFatJet.pt[sel_SR]) ({len(ak.flatten(selFatJet.pt[sel_SR]))}): {ak.flatten(selFatJet.pt[sel_SR])}")
                print(f"selFatJet.pt[sel_SR][:, 0] ({len(selFatJet.pt[sel_SR][:, 0])}): {selFatJet.pt[sel_SR][:, 0]}")
                #print(f"ak.flatten(selFatJet.pt[sel_SR][:, 0]) : {ak.flatten(selFatJet.pt[sel_SR][:, 0])}")
                #print(f"ak.flatten(selFatJet.pt[sel_SR][:, 0]) ({len(ak.flatten(selFatJet.pt[sel_SR][:, 0]))}): {ak.flatten(selFatJet.pt[sel_SR][:, 0])}")

                print(f"evtWeight ({len(evtWeight)}): {evtWeight}")
                print(f"evtWeight[sel_SR] ({len(evtWeight[sel_SR])}): {evtWeight[sel_SR]}")

            if printLevel >= 30:
                printVariable("\n weights.weight()", weights.weight())
                printVariable("\n weights.weight()[sel_SR]", weights.weight()[sel_SR])
                printVariable("\n selFatJet.pt", selFatJet.pt)
                printVariable("\n selFatJet.pt[sel_SR]", selFatJet.pt[sel_SR])
                printVariable("\n leadingFatJet.pt", leadingFatJet.pt)
                printVariable("\n leadingFatJet.pt[sel_SR]", leadingFatJet.pt[sel_SR])

            if printLevel >= 30:
                printVariable("\n ones_list", ones_list)
                printVariable("\n ones_list[sel_SR]", ones_list[sel_SR])

                printVariable("\n evtWeight", evtWeight)
                printVariable("\n evtWeight[sel_SR]", evtWeight[sel_SR])

                k = 3
                print(f"k: {k}")
                printVariable("\n ones_list * k", ones_list * k)
                printVariable("\n ones_list[sel_SR] * k", ones_list[sel_SR] * k)
                
                printVariable("\n evtWeight * k", evtWeight * k)
                printVariable("\n evtWeight[sel_SR]", evtWeight[sel_SR] * k)
                

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

            # QCD MC ----------------------------------------------
            if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isQCD'] :
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
                    if self.datasetInfo[dataset]['sample_HT_Min'] == QCDSamplesHTBins_LowEdge[idx_]:
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


                

            # events passing SR
            iBin = 4
            output['hCutFlow'].fill(
                dataset=dataset,
                CutFlow=(ones_list[sel_SR] * iBin),
                systematic=syst
            )
            output['hCutFlowWeighted'].fill(
                dataset=dataset,
                CutFlow=(ones_list[sel_SR] * iBin),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )

            
            
            '''
            output['nSelFatJet'].fill(
                dataset=dataset,
                nObject=ak.to_numpy(ak.num(selFatJet[sel_SR])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            '''
            
            output['hLeadingFatJetPt'].fill(
                dataset=dataset,
                #Pt=ak.flatten(selFatJet.pt[sel_SR][:, 0]),
                #Pt=(selFatJet.pt[sel_SR][:, 0]),
                Pt=(leadingFatJet.pt[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            
            output['hLeadingFatJetEta'].fill(
                dataset=dataset,
                Eta=(leadingFatJet.eta[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetPhi'].fill(
                dataset=dataset,
                Phi=(leadingFatJet.phi[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetMass'].fill(
                dataset=dataset,
                Mass=(leadingFatJet.mass[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetMSoftDrop'].fill(
                dataset=dataset,
                Mass=(leadingFatJet.msoftdrop[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetBtagDeepB'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.btagDeepB[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetBtagDDBvLV2'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.btagDDBvLV2[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetBtagDDCvBV2'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.btagDDCvBV2[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetBtagHbb'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.btagHbb[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_H4qvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_H4qvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_HbbvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_HbbvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )            
            output['hLeadingFatJetDeepTagMD_ZHbbvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_ZHbbvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_ZHccvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_ZHccvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_ZbbvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_ZbbvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_ZvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_ZvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_bbvsLight'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_bbvsLight[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTagMD_ccvsLight'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTagMD_ccvsLight[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTag_H'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTag_H[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTag_QCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTag_QCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetDeepTag_QCDothers'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.deepTag_QCDothers[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )

            
            output['hLeadingFatJetN2b1'].fill(
                dataset=dataset,
                N2=(leadingFatJet.n2b1[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetN3b1'].fill(
                dataset=dataset,
                N3=(leadingFatJet.n3b1[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            
            output['hLeadingFatJetTau1'].fill(
                dataset=dataset,
                TauN=(leadingFatJet.tau1[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetTau2'].fill(
                dataset=dataset,
                TauN=(leadingFatJet.tau2[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetTau3'].fill(
                dataset=dataset,
                TauN=(leadingFatJet.tau3[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetTau4'].fill(
                dataset=dataset,
                TauN=(leadingFatJet.tau4[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )

            output['hLeadingFatJetTau4by3'].fill(
                dataset=dataset,
                TauN=(np.divide(leadingFatJet.tau4[sel_SR], leadingFatJet.tau3[sel_SR])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetTau3by2'].fill(
                dataset=dataset,
                TauN=(np.divide(leadingFatJet.tau3[sel_SR], leadingFatJet.tau2[sel_SR])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetTau2by1'].fill(
                dataset=dataset,
                TauN=(np.divide(leadingFatJet.tau2[sel_SR], leadingFatJet.tau1[sel_SR])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            
            
            output['hLeadingFatJetNConstituents'].fill(
                dataset=dataset,
                nObject=(leadingFatJet.nConstituents[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            if self.datasetInfo[dataset]['isMC']:
                output['hLeadingFatJetNBHadrons'].fill(
                    dataset=dataset,
                    nObject=(leadingFatJet.nBHadrons[sel_SR]),
                    systematic=syst,
                    weight=evtWeight[sel_SR]
                )
                output['hLeadingFatJetNCHadrons'].fill(
                    dataset=dataset,
                    nObject=(leadingFatJet.nCHadrons[sel_SR]),
                    systematic=syst,
                    weight=evtWeight[sel_SR]
                )            
                
            
            output['hLeadingFatJetParticleNetMD_QCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNetMD_QCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNetMD_Xbb'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNetMD_Xbb[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNetMD_Xcc'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNetMD_Xcc[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNetMD_Xqq'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNetMD_Xqq[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNet_H4qvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNet_H4qvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNet_HbbvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNet_HbbvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )

            output['hLeadingFatJetParticleNet_HccvsQCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNet_HccvsQCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNet_QCD'].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.particleNet_QCD[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output['hLeadingFatJetParticleNet_mass'].fill(
                dataset=dataset,
                Mass=(leadingFatJet.particleNet_mass[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )


            if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isSignal']: 
                output['hGenHiggsPt_all'].fill(
                    dataset=dataset,
                    #Pt=(ak.flatten(genHiggs.pt)),
                    Pt=(ak.firsts(genHiggs.pt)),
                    systematic=syst,
                    weight=evtWeight_gen
                )
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



            if self.datasetInfo[dataset]['isMC']:
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
                 
                
            # QCD MC ----------------------------------------------
            if self.datasetInfo[dataset]['isMC'] and self.datasetInfo[dataset]['isQCD']:
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






















                
                
                if printLevel >= 100:                    
                    printVariable("\n genBQuarks_pT", genBQuarks_pT)
                    printVariable("\n ak.count(genBQuarks_pT, axis=-1)", ak.count(genBQuarks_pT, axis=-1))
                    printVariable("\n (ak.count(genBQuarks_pT, axis=-1) >= 2)", (ak.count(genBQuarks_pT, axis=-1) >= 2))
                    printVariable("\n evtWeight_gen", evtWeight_gen)
                    printVariable("\n evtWeight_gen[(ak.count(genBQuarks_pT, axis=-1) >= 2)]", evtWeight_gen[(ak.count(genBQuarks_pT, axis=-1) >= 2)])
                    printVariable("\n genBQuarks_pT[(ak.count(genBQuarks_pT, axis=-1) >= 2)]", genBQuarks_pT[(ak.count(genBQuarks_pT, axis=-1) >= 2)])
                    printVariable("\n genBQuarks_pT[(ak.count(genBQuarks_pT, axis=-1) >= 2)][:, 1]", genBQuarks_pT[(ak.count(genBQuarks_pT, axis=-1) >= 2)][:, 1])



                if printLevel >= 20:
                    mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 1)
                    printVariable("\n genBQuarks[idx_genBQuarks_pTsort][mask_tmp][:, 0].pt", genBQuarks[idx_genBQuarks_pTsort][mask_tmp][:, 0].pt)
                    printVariable("\n mask_tmp", mask_tmp)

                    mask_tmp = (ak.count(genBQuarks.pt, axis=-1) >= 1)
                    printVariable("\n genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 0].pt", genBQuarks[idx_genBQuarks_pTsort][(mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)][:, 0].pt)
                    printVariable("\n (mask_tmp & mask_QCD_stitch_CutBHadron_eventwise)", (mask_tmp & mask_QCD_stitch_CutBHadron_eventwise))

                    
                    
                '''
                Notes: Check status of first occuriing b-quark in event, and check how it is produced.

                '''
                    


                
            '''
            output[''].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output[''].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            output[''].fill(
                dataset=dataset,
                MLScore=(leadingFatJet.[sel_SR]),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            '''

            

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

    lumiScale = 1
    nEventsToAnalyze    = config["nEventsToAnalyze"] if "nEventsToAnalyze" in config else nEventToReadInBatch
    sInputFiles         = config["inputFiles"]
    sOutputFile         = config["outputFile"]
    sample_dataset      = config["dataset"] 
    sample_category     = config['sampleCategory']
    isMC                = config["isMC"]
    era                 = config['era']
    downloadIpFiles     = config['downloadIpFiles'] if 'downloadIpFiles' in config else False
    if isMC:
        luminosity          = Luminosities[era][0]
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

    print(f"isMC: {isMC}, lumiScale: {lumiScale}")
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
        sInputFiles[iFile] = getNanoAODFile(
            fileName = sInputFile, 
            useLocalFileIfExists = True, 
            downloadFile = True, 
            fileNameLocal = './inputFiles/%s' %(os.path.basename(sInputFile)), 
            nTriesToDownload = 3)
        
    print(f"\nActual  sInputFiles ({len(sInputFiles)}) (type {type(sInputFiles)}):");
    for sInputFile in sInputFiles:
        print(f"\t{sInputFile} \t {os.path.exists(sInputFile) = }");  
    sys.stdout.flush()
    print(f"htoaa_Analysis_GGFMode:: here18 {datetime.now() = }")


    sampleInfo = {
        "isMC": isMC,
        "datasetNameFull": sample_dataset,
        "lumiScale": lumiScale
    }
    if isMC:
        sampleInfo["MCSamplesStitchOption"] = MCSamplesStitchOption
        if MCSamplesStitchOption == MCSamplesStitchOptions.PhSpOverlapRewgt:
            sampleInfo["hMCSamplesStitch"] = hMCSamplesStitch
    print(f"htoaa_Analysis_GGFMode:: here19 {datetime.now() = }")
        
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
            datasetInfo={
                "era":           era, 
                sample_category: sampleInfo,
            }
        )
    )

    print(f"metrics: {metrics}")


    if 'cutflow' in output.keys():
        print("Cutflow::")
        #for key, value in output['cutflow'].items():
        for key in output['cutflow'].keys():
            #print(key, value)
            if key.startswith(sWeighted): continue

            print("%10f\t%10d\t%s" % (output['cutflow'][sWeighted+key], output['cutflow'][key], key))


    
    if sOutputFile is not None:
        if not sOutputFile.endswith('.root'): sOutputFile += '.root'
        #sOutputFile = sOutputFile.replace('.root', '_wCoffea.root') # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        sDir1 = 'evt/%s' % (sample_category)

        
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                #print(f"key: {key},  value ({type(value)}): {value}")
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
                        fOut['%s/%s_%s' % (sDir1, key, _syst)] = h1
                
                #fOut['%s%s' % (sDir1, key)] = value
                #fOut['%s%s' % (sDir1, key)] = hist.export1d(value)
                #fOut['%s%s' % (sDir1, key)] = h1

                
        

        #util.save(output, sOutputFile)
            
        
        print("Wrote to sOutputFile {}".format(sOutputFile))
   















    
    current_memory, peak_memory = tracemalloc.get_traced_memory() # https://medium.com/survata-engineering-blog/monitoring-memory-usage-of-a-running-python-program-49f027e3d1ba
    print(f"\n\nMemory usage:: current {current_memory / 10**6}MB;  peak {peak_memory / 10**6}MB")

    endTime = time.time()
    totalTime = endTime - startTime
    totalTime_hr  = int(totalTime/60/60)
    totalTime_min = totalTime - float(totalTime_hr * 60)
    totalTime_min = int(totalTime_min/60)
    totalTime_sec = totalTime - float(totalTime_hr * 60*60) - float(totalTime_min * 60)
    print(f"Total run time: {totalTime_hr}h {totalTime_min}m {totalTime_sec}s = {totalTime}sec ")
    
