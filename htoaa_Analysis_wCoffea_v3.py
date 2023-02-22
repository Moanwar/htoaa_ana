#htoaa analysis main code

import os
import sys
import json
import glob
from collections import OrderedDict as OD
import time
import tracemalloc
import math
import numpy as np
from copy import deepcopy

'''
H->aa->4b boosted analysis macro

References:
  * Coffea framework used for TTGamma analysis: https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py
* Coffea installation: /home/siddhesh/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea
'''

#import coffea.processor as processor
from coffea import processor, util
from coffea.nanoevents import schemas
from coffea.nanoevents.methods import nanoaod, vector
from coffea.analysis_tools import PackedSelection, Weights
#import hist
from coffea import hist
import awkward as ak
import uproot
#from dask.distributed import Client


from htoaa_Settings import *
from htoaa_CommonTools import GetDictFromJsonFile, calculate_lumiScale
'''
from htoaa_Settings import *
from htoaa_NanoAODBranches import htoaa_nanoAODBranchesToRead
from htoaa_CommonTools import GetDictFromJsonFile, DfColLabel_convert_bytes_to_string, calculate_lumiScale
from htoaa_CommonTools import cut_ObjectMultiplicity, cut_ObjectPt, cut_ObjectEta, cut_ObjectPt_1
'''



 
printLevel = 0
nEventToReadInBatch =  0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze =  -1 # 1000 # 100000 # -1
#pd.set_option('display.max_columns', None)

#print("".format())



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

        self.FatJetPtThsh  = 170
        self.FatJetEtaThsh = 2.4

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
        #self.isMC = isMC
        self.objectSelector = ObjectSelection(era=self.datasetInfo["era"])
        
        #dataset_axis = hist.axis.StrCategory(name="dataset", label="", categories=[], growth=True)
        #muon_axis = hist.axis.Regular(name="massT", label="Transverse Mass [GeV]", bins=50, start=15, stop=250)
        dataset_axis    = hist.Cat("dataset", "Dataset")
        systematic_axis = hist.Cat("systematic", "Systematic Uncertatinty")

        nObject_axis  = hist.Bin("nObject",   r"No. of object",   21, -0.5, 20.5)
        pt_axis       = hist.Bin("Pt",        r"$p_{T}$ [GeV]",   200, 0, 1000)
        eta_axis      = hist.Bin("Eta",       r"$#eta$",          100, -5, 5)
        phi_axis      = hist.Bin("Phi",       r"$\phi$",          100, -3.14, 3.13)
        #mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       200, 0, 600)
        #mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       400, 0, 200)
        mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       300, 0, 300)
        mass_axis1    = hist.Bin("Mass1",     r"$m$ [GeV]",       300, 0, 300)
        mlScore_axis  = hist.Bin("MLScore",   r"ML score",        100, -1.1, 1.1)
        jetN2_axis    = hist.Bin("N2",        r"N2b1",            100, 0, 3)
        jetN3_axis    = hist.Bin("N3",        r"N3b1",            100, 0, 5)
        jetTau_axis   = hist.Bin("TauN",      r"TauN",            100, 0, 1)
        deltaR_axis   = hist.Bin("deltaR",     r"$delta$ r ",     500, 0, 5)

        sXaxis      = 'xAxis'
        sXaxisLabel = 'xAxisLabel'
        sYaxis      = 'yAxis'
        sYaxisLabel = 'yAxisLabel'
        histos = OD([
            # ('histogram_name',  {sXaxis: hist.Bin() axis,  sXaxisLabel: "histogram axis label"})
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
            ('hGenHiggsPt_sel',                           {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(GEN Higgs (pdgId: 25, status=62))$ [GeV]"}),
            ('hGenHiggsPt_sel_wGenCuts',                  {sXaxis: pt_axis,         sXaxisLabel: r"$p_{T}(GEN Higgs (pdgId: 25, status=62))$ [GeV]"}),

            ('hGenHiggsMass_all',                        {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN H) [GeV]"}),
            ('hMass_GenA_all',                        {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN A) [GeV]"}),
            ('hMass_GenAApair_all',                        {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN HToAA) [GeV]"}),
            ('hMass_GenAToBBbarpair_all',                        {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN AToBB) [GeV]"}),
            ('hMass_Gen4BFromHToAA_all',                        {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN HTOAATo4B) [GeV]"}),
            ('hMass_GenAToBBbarpair_all_1',                        {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN AToBB) [GeV]"}),
            ('hMass_Gen4BFromHToAA_all_1',                        {sXaxis: mass_axis,       sXaxisLabel: r"m (GEN HTOAATo4B) [GeV]"}),
            ('hDeltaR_GenH_GenB_max',                        {sXaxis: deltaR_axis,       sXaxisLabel: r"$Delta$r (GEN H, GEN B)_{max}"}),

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

        if printLevel >= 5:
            print(f"nEvents: {len(events)}")
        if printLevel >= 5:
            print(f"events.fields: {events.fields}")
        
        if self.datasetInfo[dataset]['isMC']:
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

        
        if printLevel >= 6:
            print(f"events.GenPart: {events.GenPart.fields}") # ['eta', 'mass', 'phi', 'pt', 'genPartIdxMother', 'pdgId', 'status', 'statusFlags', 'genPartIdxMotherG', 'distinctParentIdxG', 'childrenIdxG', 'distinctChildrenIdxG']
            print(f"events.GenPart.nGenPart: {ak.count(events.GenPart.pdgId, axis=1).to_list()}")
            '''
            print(f"\nevents.GenPart.pdgId: {events.GenPart.pdgId.to_list()}")
            print(f"\nevents.GenPart.pt: {events.GenPart.pt.to_list()}")
            print(f"\nevents.GenPart.eta: {events.GenPart.eta.to_list()}")
            print(f"\nevents.GenPart.phi: {events.GenPart.phi.to_list()}")
            print(f"\nevents.GenPart.mass: {events.GenPart.mass.to_list()}")
            print(f"\nevents.GenPart.genPartIdxMother: {events.GenPart.genPartIdxMother.to_list()}")
            print(f"\nevents.GenPart.status: {events.GenPart.status.to_list()}")
            print(f"\nevents.GenPart.statusFlags: {events.GenPart.statusFlags.to_list()}")
            print(f"\nevents.GenPart.genPartIdxMotherG: {events.GenPart.genPartIdxMotherG.to_list()}")
            print(f"\nevents.GenPart.distinctParentIdxG: {events.GenPart.distinctParentIdxG.to_list()}")
            print(f"\nevents.GenPart.childrenIdxG: {events.GenPart.childrenIdxG.to_list()}")
            print(f"\nevents.GenPart.distinctChildrenIdxG: {events.GenPart.distinctChildrenIdxG.to_list()}"); sys.stdout.flush()
            #print(f"events.GenPart.: {events.GenPart.}")
            '''
            print(f"\nevents.GenPart: {events.GenPart.to_list()} ")
            print(f"\nevents.GenPart.eta: {events.GenPart.eta.to_list()} ")

        if printLevel >= 9:
            mask_GenHiggs_1 = (
                (events.GenPart.pdgId  == 25) #& # pdgId:: 25: H0
                #(events.GenPart.status == 62)   # statu 62: outgoing subprocess particle with primordial kT included https://pythia.org/latest-manual/ParticleProperties.html
            )
            #print(f"\n mask_GenHiggs_1:  {mask_GenHiggs.to_list()} ")
            print(f"\n events.GenPart[mask_GenHiggs_1]:  {events.GenPart[mask_GenHiggs_1].to_list()} ")
            print(f"\n events.GenPart[mask_GenHiggs_1].mass:  {events.GenPart[mask_GenHiggs_1].mass.to_list()} "); sys.stdout.flush()
            print(f"\n events.GenPart[mask_GenHiggs_1].status:  {events.GenPart[mask_GenHiggs_1].status.to_list()} "); sys.stdout.flush()

            mask_GenHiggs = (
                (events.GenPart.pdgId  == 25) & # pdgId:: 25: H0
                (events.GenPart.status == 62)   # statu 62: outgoing subprocess particle with primordial kT included https://pythia.org/latest-manual/ParticleProperties.html
            )
            print(f"\n mask_GenHiggs:  {mask_GenHiggs.to_list()} ")
            print(f"\n events.GenPart[mask_GenHiggs]:  {events.GenPart[mask_GenHiggs].to_list()} ")
            print(f"\n events.GenPart[mask_GenHiggs].mass:  {events.GenPart[mask_GenHiggs].mass.to_list()} "); sys.stdout.flush()

            GenHiggsCollection = events.GenPart[mask_GenHiggs]
            print(f"GenHiggsCollection: {GenHiggsCollection.to_list()}")

        if printLevel >= 9:
            mask_GenA = (
                (events.GenPart.pdgId == 36)
            )

            print(f"\n mask_GenA: {mask_GenA.to_list()}")
            print(f"\n events.GenPart[mask_GenA]: {events.GenPart[mask_GenA].to_list()}")
            print(f"\n events.GenPart[mask_GenA].mass: {events.GenPart[mask_GenA].mass.to_list()}")

            print(f"\n events.GenPart[0]: {events.GenPart[0]} ")
            print(f"events.GenPart[:, 0]: {events.GenPart[:, 0]} ")
            print(f"events.GenPart[:, 0].eta: {events.GenPart[:, 0].eta} ")


            #  events.GenPart[mask_GenHiggs]:  [[{'eta': 1.03515625, 'mass': 125.0, 'phi': -1.37890625, 'pt': 181.5, 'genPartIdxMother': 10, 'pdgId': 25, 'status': 62, 'statusFlags': 10497, 'genPartIdxMotherG': 10, 'distinctParentIdxG': 0, 'childrenIdxG': [12, 13], 'distinctChildrenIdxG': [12, 13]}]]
            #  events.GenPart[mask_GenA]: [[{'eta': 1.0625, 'mass': 20.0, 'phi': -1.49609375, 'pt': 175.5, 'genPartIdxMother': 11, 'pdgId': 36, 'status': 22, 'statusFlags': 14721, 'genPartIdxMotherG': 11, 'distinctParentIdxG': 11, 'childrenIdxG': [14, 15], 'distinctChildrenIdxG': [14, 15, 20, 21]}, {'eta': 0.023193359375, 'mass': 20.0, 'phi': -0.1494140625, 'pt': 21.625, 'genPartIdxMother': 11, 'pdgId': 36, 'status': 22, 'statusFlags': 14721, 'genPartIdxMotherG': 11, 'distinctParentIdxG': 11, 'childrenIdxG': [16, 17], 'distinctChildrenIdxG': [16, 17]}]]


            print(f"\n mask_GenHiggs: {mask_GenHiggs.to_list()} ")
            print(f"\n ak.local_index(mask_GenHiggs, axis=-1): {ak.local_index(mask_GenHiggs, axis=-1).to_list()} ")
            print(f"ak.argmax(mask_GenHiggs, axis=-1): {ak.argmax(mask_GenHiggs, axis=-1).to_list()} ")

            index_GenHiggs = ak.argmax(mask_GenHiggs, axis=-1) # ak.argmax(mask_GenHiggs, axis=-1): [11, 10]
            #index_GenHiggs = ak.argmax(mask_GenHiggs, axis=-1, keepdims=True) # ak.argmax(mask_GenHiggs, axis=-1): [11, 10]
            
            print(f"index_GenHiggs: {index_GenHiggs.to_list()}")
            print(f"\n events.GenPart.genPartIdxMotherG: {events.GenPart.genPartIdxMotherG.to_list()}")
            
            '''
            # did not work
            mask_GenA_1 = (
                (events.GenPart.genPartIdxMotherG == index_GenHiggs)
            )
            print(f"mask_GenA_1: {mask_GenA_1.to_list()}")
            print(f"events.GenPart[mask_GenA_1]: {events.GenPart[mask_GenA_1].to_list()}")

            mask_GenA_2 = np.vectorize(get_GenPartDaughters)(events.GenPart.genPartIdxMotherG, index_GenHiggs)
            '''


            selGenA = events.GenPart[mask_GenA]
            selGenA_first  = selGenA[:, 0]
            selGenA_second = selGenA[:, 1]

            '''
            mask_GenB = (
                abs(events.GenPart.pdgId) == 5                
            )

            '''


            selGenB    = events.GenPart[( (events.GenPart.pdgId ==  5) & (events.GenPart.status == 23) )] # status = 23: outgoing particles of the hardest subprocess
            selGenBbar = events.GenPart[( (events.GenPart.pdgId == -5) & (events.GenPart.status == 23) )] # status = 71: partons in preparation of hadronization process, 71 : copied partons to collect into contiguous colour singlet
            selGenBBbar_pairs = ak.cartesian((selGenB, selGenBbar) )

            selGenBBbar_pair_B, selGenBBbar_pair_Bbar = ak.unzip(selGenBBbar_pairs)

            selGenBBbar_LV = selGenBBbar_pair_B + selGenBBbar_pair_Bbar

            dr_GenAFirst_GenBBbar  = selGenA_first.delta_r(selGenBBbar_LV)
            dr_GenASecond_GenBBbar = selGenA_second.delta_r(selGenBBbar_LV)

            
            
            print(f"\n selGenA: {selGenA.to_list()}")
            print(f"selGenA_first: {selGenA_first.to_list()}")
            print(f"selGenA_second: {selGenA_second.to_list()}")

            print(f"\n events.GenPart[(abs(events.GenPart.pdgId) ==  5)]: {events.GenPart[(abs(events.GenPart.pdgId) ==  5)].to_list()} ")
            print(f"\n events.GenPart[(abs(events.GenPart.pdgId) ==  5)].status: {events.GenPart[(abs(events.GenPart.pdgId) ==  5)].status.to_list()} ")
            print(f"\n selGenB: {selGenB.to_list()} ")
            print(f"\n selGenBbar: {selGenBbar.to_list()} ")
            print(f"\n selGenBBbar_pairs: {selGenBBbar_pairs.to_list()} ")

            print(f"\n selGenBBbar_pair_B: {selGenBBbar_pair_B.to_list()} ")
            print(f"\n selGenBBbar_pair_Bbar: {selGenBBbar_pair_Bbar.to_list()} ")

            print(f"\n selGenBBbar_LV: {selGenBBbar_LV.to_list()} ")
            print(f"\n selGenBBbar_LV.mass: {selGenBBbar_LV.mass.to_list()} ")

            print(f"\n dr_GenAFirst_GenBBbar: {dr_GenAFirst_GenBBbar.to_list()} ")
            print(f"\n dr_GenASecond_GenBBbar: {dr_GenASecond_GenBBbar.to_list()} ")

            
            selGenBQuarks = events.GenPart[( (abs(events.GenPart.pdgId) ==  5) & (events.GenPart.status == 23) )] # status = 23: outgoing particles of the hardest subprocess.  status = 71: partons in preparation of hadronization process, 71 : copied partons to collect into contiguous colour singlet
            genBBar_pairs_all = ak.argcombinations(selGenBQuarks, 2, fields=['b', 'bbar'])
            genBBar_pairs = genBBar_pairs_all[(
                ((selGenBQuarks[genBBar_pairs_all['b']].pdgId) == (-1*selGenBQuarks[genBBar_pairs_all['bbar']].pdgId)  ) &
                (selGenBQuarks[genBBar_pairs_all['b']].genPartIdxMother == selGenBQuarks[genBBar_pairs_all['bbar']].genPartIdxMother)
            )]
            #selGenBBbar_pairs_1 = selGenBQuarks[genBBar_pairs['b']]

            print(f"\n\n\n MethodII]: \n selGenBQuarks: {selGenBQuarks.to_list()}")
            #print(f"\n genBBar_pairs_all: {genBBar_pairs_all.to_list()} ")
            print(f"\n genBBar_pairs: {genBBar_pairs.to_list()} ")
            print(f"\n selGenBQuarks[genBBar_pairs['b']: {selGenBQuarks[genBBar_pairs['b']].to_list()} ")
            print(f"\n selGenBQuarks[genBBar_pairs['bbar']: {selGenBQuarks[genBBar_pairs['bbar']].to_list()} ")
            print(f"\n (selGenBQuarks[genBBar_pairs['b']] + selGenBQuarks[genBBar_pairs['bbar']]): { (selGenBQuarks[genBBar_pairs['b']] + selGenBQuarks[genBBar_pairs['bbar']]).to_list() } ")
            print(f"\n (selGenBQuarks[genBBar_pairs['b']] + selGenBQuarks[genBBar_pairs['bbar']]).mass: { (selGenBQuarks[genBBar_pairs['b']] + selGenBQuarks[genBBar_pairs['bbar']]).mass.to_list() } ")
            print(f"\n (selGenBQuarks[:, 0] + selGenBQuarks[:, 1] + selGenBQuarks[:, 2] + selGenBQuarks[:, 3]).mass: {(selGenBQuarks[:, 0] + selGenBQuarks[:, 1] + selGenBQuarks[:, 2] + selGenBQuarks[:, 3]).mass.to_list()} ")


            genBBar_pairs_all_1 = ak.argcombinations(events.GenPart, 2, fields=['b', 'bbar'])
            genBBar_pairs_1 = genBBar_pairs_all_1[(
                (abs(events.GenPart[genBBar_pairs_all_1['b']].pdgId) == 5) &
                (abs(events.GenPart[genBBar_pairs_all_1['bbar']].pdgId) == 5) &
                ((events.GenPart[genBBar_pairs_all_1['b']].pdgId) == (-1*events.GenPart[genBBar_pairs_all_1['bbar']].pdgId)  ) &
                (events.GenPart[genBBar_pairs_all_1['b']].genPartIdxMother == events.GenPart[genBBar_pairs_all_1['bbar']].genPartIdxMother) &
                (events.GenPart[ events.GenPart[genBBar_pairs_all_1['b']].genPartIdxMother ].pdgId == 36) &
                (events.GenPart[ events.GenPart[genBBar_pairs_all_1['bbar']].genPartIdxMother ].pdgId == 36) 
            )]

           
            #print(f"\n\n Check here:: genBBar_pairs_all_1: {genBBar_pairs_all_1.to_list()} ")
            print(f"\n genBBar_pairs_1: {genBBar_pairs_1.to_list()} ")



            idxGenA_sortByMass = ak.argsort(selGenA.mass, axis=-1, ascending=False) 
            print(f"\n\n\n selGenA: {selGenA.to_list()}")
            print(f"\n selGenA.mass: {selGenA.mass.to_list()}")
            print(f"\n idxGenA_sortByMass: {idxGenA_sortByMass.to_list()}")
            print(f"\n selGenA[idxGenA_sortByMass]. {selGenA[idxGenA_sortByMass].to_list()}")
            print(f"\n selGenA[idxGenA_sortByMass].mass: {selGenA[idxGenA_sortByMass].mass.to_list()}")
            print(f"\n idxGenA_sortByMass[:, 0]: {idxGenA_sortByMass[:, 0].to_list()}")
            print(f"\n selGenA[idxGenA_sortByMass[:, 0]].mass: {selGenA[idxGenA_sortByMass[:, 0]].mass.to_list()}") # wrong
            print(f"\n selGenA[idxGenA_sortByMass][:, 0].mass: {selGenA[idxGenA_sortByMass][:, 0].mass.to_list()}") # correct
            
            
            
            #print(f": {} ")
            
            

            



            '''
# find distance between leading jet and all electrons in each event
dr = events.Jet[:, 0].delta_r(events.Electron)
dr
 <Array [[], [3.13], [3.45, ... 0.0858], [], []] type='40 * var * float32'>
￼
[8]:
# find minimum distance
ak.min(dr, axis=1)
￼
[8]:
<Array [None, 3.13, 2.18, ... None, None] type='40 * ?float32'>
￼
[9]:
# a convenience method for this operation on all jets is available
events.Jet.nearest(events.Electron)           
<ElectronArray [[None, None, None, ... [None, None]] type='40 * var * ?electron'>
            


For generated particles, the parent index is similarly mapped:

[14]:
events.GenPart.parent.pdgId
￼
[14]:
<Array [[None, None, 1, 1, ... 111, 111, 111]] type='40 * var * ?int32[parameter...'>
￼
In addition, using the parent index, a helper method computes the inverse mapping, namely, children. As such, one can find particle siblings with:

[15]:
events.GenPart.parent.children.pdgId
# notice this is a doubly-jagged array
￼
[15]:
<Array [[None, None, [23, 21, ... [22, 22]]] type='40 * var * option[var * ?int3...'>
￼
Since often one wants to shortcut repeated particles in a decay sequence, a helper method distinctParent is also available. Here we use it to find the parent particle ID for all prompt electrons:

[16]:
events.GenPart[
    (abs(events.GenPart.pdgId) == 11)
    & events.GenPart.hasFlags(['isPrompt', 'isLastCopy'])
].distinctParent.pdgId
￼
[16]:
<Array [[], [23, 23], [23, ... [23, 23], []] type='40 * var * ?int32[parameters=...'>
￼

            
            '''
            
        ##################
        # OBJECT SELECTION
        ##################

        # FatJet selection
        selFatJet = self.objectSelector.selectFatJets(events)

        genHiggs = None
        genHT    = None
        if self.datasetInfo[dataset]['isMC']:
            genHiggs  = self.objectSelector.selectGenHiggs(events)        
            genHT     = self.objectSelector.GenHT(events)

            # m(bbar from A) and m(4b from HToAA) ----------
            genHiggsCollection = events.GenPart[(
                (events.GenPart.pdgId  == 25) & # pdgId:: 25: H0
                (events.GenPart.status == 62)   # statu 62: outgoing subprocess particle with primordial kT included https://pythia.org/latest-manual/ParticleProperties.html
            )]

            genACollection = events.GenPart[(
                (events.GenPart.pdgId == 36)
            )]
            genA_First  = genACollection[:, 0]
            genA_Second = genACollection[:, 1]
            
            idxGenA_sortByMass = ak.argsort(genACollection.mass, axis=-1, ascending=False)
            # genACollection[idxGenA_sortByMass[0]] : Leading mass GenA
            # genACollection[idxGenA_sortByMass[1]] : Subleading mass GenA
            
            '''
            genBQuarkCollection = events.GenPart[(
                (abs(events.GenPart.pdgId) ==  5) &
                (events.GenPart.status == 23)
            )] # status = 23: outgoing particles of the hardest subprocess.  status = 71: partons in preparation of hadronization process, 71 : copied partons to collect into contiguous colour singlet
            genBBar_pairs_all = ak.argcombinations(genBQuarkCollection, 2, fields=['b', 'bbar'])
            genBBar_pairs = genBBar_pairs_all[(
                ((genBQuarkCollection[genBBar_pairs_all['b']].pdgId) == (-1*genBQuarkCollection[genBBar_pairs_all['bbar']].pdgId)  ) &
                (genBQuarkCollection[genBBar_pairs_all['b']].genPartIdxMother == genBQuarkCollection[genBBar_pairs_all['bbar']].genPartIdxMother)
            )]
            '''

            
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

            if printLevel >= 9:
                print(f" \n ak.num(genHiggs): {ak.num(genHiggs).to_list()} ")
                print(f" \n ak.num(genHiggs) == 1: {(ak.num(genHiggs) == 1).to_list()} ")
                print(f"\n ak.num(genBBar_pairs['b']: {ak.num(genBBar_pairs['b']).to_list()} ")
                print(f"\n ak.num(genBBar_pairs: {ak.num(genBBar_pairs).to_list()} ")
                print(f"\n genBBar_pairs: {genBBar_pairs.to_list()}")
            

            if printLevel >= 9:
                print(f"\n events.GenPart[genBBar_pairs['b']]: {events.GenPart[genBBar_pairs['b']].to_list()} ")
                print(f"\n events.GenPart[genBBar_pairs['b']][:, 0].pt: {events.GenPart[genBBar_pairs['b']][:, 0].pt.to_list()} ")
                print(f"\n ak.num(events.GenPart[genBBar_pairs['b']][:, 0].pt, axis=0): {ak.num(events.GenPart[genBBar_pairs['b']][:, 0].pt, axis=0)} ")

                print(f"\n LVGenB_0: {LVGenB_0.to_list()} ")
                print(f"\n LVGenBbar_0: {LVGenBbar_0.to_list()} ")
                add1 = LVGenB_0.add(LVGenBbar_0)
                print(f"\n after addition \n LVGenB_0: {LVGenB_0.to_list()} ")
                print(f"\n LVGenBbar_0: {LVGenBbar_0.to_list()} ")
                print(f"\n add1: {add1.to_list()} ")
                print(f"\n add1.mass: {add1.mass.to_list()} ")
                print(f"\n (LVGenB_0 + LVGenBbar_0).mass: {(LVGenB_0 + LVGenBbar_0).mass.to_list()} ")
                
                

            if printLevel >= 9:
                print(f"\n LVGenB_0: {LVGenB_0.to_list()} ")
                print(f"\n LVGenBbar_0: {LVGenBbar_0.to_list()} ")
                print(f"\n LVGenB_1: {LVGenB_1.to_list()} ")
                print(f"\n LVGenBbar_1: {LVGenBbar_1.to_list()} ")
                
                # genBQuarks_sel1 = ak.concatenate([LVGenB_0, LVGenBbar_0, LVGenB_1, LVGenBbar_1], axis=1) # not working
                genBQuarks_sel1 = ak.zip([LVGenB_0, LVGenBbar_0, LVGenB_1, LVGenBbar_1])
                print(f"\n genBQuarks_sel1: {genBQuarks_sel1.to_list()} ")
                
                #gen_Bquarks
                print(f"\n genHiggs.delta_r(LVGenB_0): {genHiggs.delta_r(LVGenB_0).to_list()}")
                print(f"\n genHiggs.delta_r(LVGenBbar_0): {genHiggs.delta_r(LVGenBbar_0).to_list()}")
                print(f"\n genHiggs.delta_r(LVGenB_1): {genHiggs.delta_r(LVGenB_1).to_list()}")
                print(f"\n genHiggs.delta_r(LVGenBbar_1): {genHiggs.delta_r(LVGenBbar_1).to_list()}")
                #print(f"\n genHiggs.delta_r(genBQuarks_sel1): {genHiggs.delta_r(genBQuarks_sel1).to_list()}")
                #print(f"\n genHiggs.delta_r(LVGenB_0): {genHiggs.delta_r(events.GenPart[]).to_list()}")
                print(f"\n ak.concatenate([genHiggs.delta_r(LVGenB_0), genHiggs.delta_r(LVGenBbar_0)], axis=-1): {ak.concatenate([genHiggs.delta_r(LVGenB_0), genHiggs.delta_r(LVGenBbar_0), genHiggs.delta_r(LVGenB_1), genHiggs.delta_r(LVGenBbar_1)], axis=-1).to_list()}")
                
                print(f"\n dr_GenH_GenB: {dr_GenH_GenB.to_list()}")                
                print(f"\n max_dr_GenH_GenB: {max_dr_GenH_GenB.to_list()} ")
                


            
        #####################
        # EVENT SELECTION
        #####################

        # create a PackedSelection object
        # this will help us later in composing the boolean selections easily
        selection = PackedSelection()

        if printLevel >= 12:
            print(f"events.PV.npvsGood.to_list(): {events.PV.npvsGood.to_list()}")
        # nPVGood >= 1
        selection.add("nPV", events.PV.npvsGood >= 1)
        
        # >=1 FatJet
        selection.add("FatJetGet", ak.num(selFatJet) >= self.objectSelector.nFatJetMin)

        if self.datasetInfo[dataset]['isMC']:
            #selection.add("GenHT", genHT >= self.objectSelector.GenHTThsh)
            selection.add("LHEHT", events.LHE.HT  >= self.objectSelector.LHEHTThsh)

            if printLevel >= 15:
                print(f"\nevents.LHE.fields: {events.LHE.fields}")
                print(f"\nevents.LHE.HT: {events.LHE.HT.to_list()}")
                print(f"\ngenHT: {genHT.to_list()} ")


            selection.add("1GenHiggs", ak.num(genHiggs) == 1)
            selection.add("2GenA", ak.num(genACollection) == 2)
            selection.add("2GenAToBBbarPairs", ak.num(genBBar_pairs) == 2)

            
        sel_names_all = OD([
            ("SR",           ["nPV", "FatJetGet"]),
            ("SR_wGenCuts",  ["LHEHT"]),
            ("GenHToAATo4B", ["1GenHiggs", "2GenA", "2GenAToBBbarPairs"]),
        ])
        #sel_SR          = selection.all("nPV", "FatJetGet")
        sel_SR           = selection.all(* sel_names_all["SR"])
        sel_SR_wGenCuts  = selection.all(* sel_names_all["SR_wGenCuts"])
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



        #if printLevel >= 5:
        #    print(f" {events.}")
            
            
        ##################
        # EVENT VARIABLES
        ##################

        leadingFatJet = ak.firsts(selFatJet)


        
        ################
        # EVENT WEIGHTS
        ################

        # create a processor Weights object, with the same length as the number of events in the chunk
        weights     = Weights(len(events))
        weights_gen = Weights(len(events))
        weights_GenHToAATo4B = Weights(len(events))


        if self.datasetInfo[dataset]["isMC"]:
            weights.add(
                "lumiWeight",
                weight=np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
            )
            weights.add(
                "genWeight",
                weight=np.copysign(np.ones(len(events)), events.genWeight)
            )
            weights.add(
                "btagWeight",
                weight=(events.btagWeight.DeepCSVB)
            )
            
 
            weights_gen.add(
                "lumiWeight",
                weight=np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
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
        for n in selection.names:
            output['cutflow'][n] += selection.all(n).sum()

        for iSelection in sel_names_all.keys():
            iName = f"{iSelection}: {sel_names_all[iSelection]}"
            output['cutflow'][iName] += selection.all(* sel_names_all[iSelection]).sum()
            

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

                
                
            output['nSelFatJet'].fill(
                dataset=dataset,
                nObject=ak.to_numpy(ak.num(selFatJet[sel_SR])),
                systematic=syst,
                weight=evtWeight[sel_SR]
            )
            
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
            output['hLeadingFatJetNConstituents'].fill(
                dataset=dataset,
                nObject=(leadingFatJet.nConstituents[sel_SR]),
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


            if self.datasetInfo[dataset]["isMC"]:
                output['hGenHiggsPt_all'].fill(
                    dataset=dataset,
                    Pt=(ak.flatten(genHiggs.pt)),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hGenHiggsPt_sel'].fill(
                    dataset=dataset,
                    Pt=(ak.flatten(genHiggs.pt[sel_SR])),
                    systematic=syst,
                    weight=evtWeight_gen[sel_SR]
                )
                output['hGenHiggsPt_sel_wGenCuts'].fill(
                    dataset=dataset,
                    Pt=(ak.flatten(genHiggs.pt[sel_SR_wGenCuts])),
                    systematic=syst,
                    weight=evtWeight_gen[sel_SR_wGenCuts]
                )

                
                # m(2b from ATo2B) and m(4b from HToAATo4b) --------------
                if printLevel >= 9:
                    print(f"\n Fill histograms: \ngenHiggs.mass: {genHiggs.mass.to_list()} ")
                    print(f"\n ak.flatten(genHiggs.mass): {ak.flatten(genHiggs.mass.to_list())} ")
                    print(f"\n output['hGenHiggsMass_all']: {output['hGenHiggsMass_all']} ")
                    
                output['hGenHiggsMass_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                if printLevel >= 9:
                    print(f"\n genACollection.mass ({len(genACollection.mass)}): {genACollection.mass}")
                    print(f"\n ak.flatten(genACollection.mass) ({len(ak.flatten(genACollection.mass))}): {ak.flatten(genACollection.mass)}")
                    print(f"\n evtWeight_gen ({len(evtWeight_gen)}): {evtWeight_gen} ")
                    print(f"\n genACollection[:, 0].mass ({len(genACollection[:, 0].mass)}): {genACollection[:, 0].mass}")
                    #print(f"\n ak.flatten(genACollection[:, 0].mass) ({len(ak.flatten(genACollection[:, 0].mass))}): {ak.flatten(genACollection[:, 0].mass)}")
                    '''
                    print(f"\n (genBQuarkCollection[genBBar_pairs['b']] + genBQuarkCollection[genBBar_pairs['bbar']]).mass ({len((genBQuarkCollection[genBBar_pairs['b']] + genBQuarkCollection[genBBar_pairs['bbar']]).mass)}): {(genBQuarkCollection[genBBar_pairs['b']] + genBQuarkCollection[genBBar_pairs['bbar']]).mass.to_list()} ")
                    print(f"\n (genBQuarkCollection[genBBar_pairs['b']] + genBQuarkCollection[genBBar_pairs['bbar']]).mass[:, 0] ({len((genBQuarkCollection[genBBar_pairs['b']] + genBQuarkCollection[genBBar_pairs['bbar']]).mass[:, 0])}): {(genBQuarkCollection[genBBar_pairs['b']] + genBQuarkCollection[genBBar_pairs['bbar']]).mass[:, 0].to_list()} ")
                    print(f"\n genBQuarkCollection[genBBar_pairs['b']]: {genBQuarkCollection[genBBar_pairs['b']]}")
                    print(f"\n genBQuarkCollection[genBBar_pairs['b']][:, 0]: {genBQuarkCollection[genBBar_pairs['b']][:, 0]}")
                    print(f"\n (genBQuarkCollection[genBBar_pairs['b']][:, 0] + genBQuarkCollection[genBBar_pairs['bbar']][:, 0]).mass ({len((genBQuarkCollection[genBBar_pairs['b']][:, 0] + genBQuarkCollection[genBBar_pairs['bbar']][:, 0]).mass)}): {(genBQuarkCollection[genBBar_pairs['b']][:, 0] + genBQuarkCollection[genBBar_pairs['bbar']][:, 0]).mass.to_list()} ")
                    '''
                    print(f"\n\n CheckAgain:: \n genBBar_pairs: {genBBar_pairs}")
                    print(f"\n events.GenPart[genBBar_pairs['b']].pdgId: {events.GenPart[genBBar_pairs['b']].pdgId}")
                    print(f"\n events.GenPart[genBBar_pairs['bbar']].pdgId: {events.GenPart[genBBar_pairs['bbar']].pdgId}")
                    print(f"\n ((events.GenPart[genBBar_pairs['b']][:, 0] + events.GenPart[genBBar_pairs['bbar']][:, 0] + events.GenPart[genBBar_pairs['b']][:, 1] + events.GenPart[genBBar_pairs['bbar']][:, 1]).mass): {((events.GenPart[genBBar_pairs['b']][:, 0] + events.GenPart[genBBar_pairs['bbar']][:, 0] + events.GenPart[genBBar_pairs['b']][:, 1] + events.GenPart[genBBar_pairs['bbar']][:, 1]).mass)}")
                    print(f"\n Fill histogram ends....\n\n")
                    
                output['hMass_GenA_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[:, 0].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hMass_GenA_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
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
                    weight=evtWeight_gen
                )
                output['hMass_GenAToBBbarpair_all'].fill(
                    dataset=dataset,
                    Mass=((events.GenPart[genBBar_pairs['b']] + events.GenPart[genBBar_pairs['bbar']]).mass[:, 0][sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )                
                output['hMass_GenAToBBbarpair_all'].fill(
                    dataset=dataset,
                    Mass=((events.GenPart[genBBar_pairs['b']] + events.GenPart[genBBar_pairs['bbar']]).mass[:, 1][sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hMass_Gen4BFromHToAA_all'].fill(
                    dataset=dataset,
                    Mass=((events.GenPart[genBBar_pairs['b']][:, 0] + events.GenPart[genBBar_pairs['bbar']][:, 0] + events.GenPart[genBBar_pairs['b']][:, 1] + events.GenPart[genBBar_pairs['bbar']][:, 1]).mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                if printLevel >= 9:
                    print(f"\n\n (events.GenPart[genBBar_pairs['b']] + events.GenPart[genBBar_pairs['bbar']]).mass[:, 0]: {(events.GenPart[genBBar_pairs['b']] + events.GenPart[genBBar_pairs['bbar']]).mass[:, 0].to_list()}")
                    print(f"\n\n (events.GenPart[genBBar_pairs['b']] + events.GenPart[genBBar_pairs['bbar']]).mass[:, 0][sel_GenHToAATo4B]: {(events.GenPart[genBBar_pairs['b']] + events.GenPart[genBBar_pairs['bbar']]).mass[:, 0][sel_GenHToAATo4B].to_list()}")



                output['hMass_GenAToBBbarpair_all_1'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_0 + LVGenBbar_0).mass),
                    systematic=syst,
                    weight=evtWeight_gen
                )                
                output['hMass_GenAToBBbarpair_all_1'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_1 + LVGenBbar_1).mass),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                output['hMass_Gen4BFromHToAA_all_1'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_0 + LVGenBbar_0 + LVGenB_1 + LVGenBbar_1).mass),
                    systematic=syst,
                    weight=evtWeight_gen
                )


                output['hMass_GenA1_vs_GenA2_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[:, 0].mass[sel_GenHToAATo4B]),
                    Mass1=(genACollection[:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                
                output['hMass_GenA1ToBBbar_vs_GenA2ToBBbar_all'].fill(
                    dataset=dataset,
                    Mass=((LVGenB_0 + LVGenBbar_0).mass),
                    Mass1=((LVGenB_1 + LVGenBbar_1).mass),
                    systematic=syst,
                    weight=evtWeight_gen
                )

                if printLevel >= 9:
                    print(f"\n\n genHiggs.mass[sel_GenHToAATo4B]: {genHiggs.mass[sel_GenHToAATo4B].to_list()}")
                    print(f"\n ak.flatten(genHiggs.mass[sel_GenHToAATo4B]) ({len(ak.flatten(genHiggs.mass[sel_GenHToAATo4B]))}): {ak.flatten(genHiggs.mass[sel_GenHToAATo4B]).to_list()}")
                    print(f"\n genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B] ({len(genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B])}): {genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B].to_list()}")

                output['hMass_GenH_vs_GenAHeavy_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    Mass1=(genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                
                output['hMass_GenH_vs_GenALight_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    Mass1=(genACollection[idxGenA_sortByMass][:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )

                output['hMass_GenAHeavy_vs_GenALight_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B]),
                    Mass1=(genACollection[idxGenA_sortByMass][:, 1].mass[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )

                output['hDeltaR_GenH_GenB_max'].fill(
                    dataset=dataset,
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )

                output['hMassGenH_vs_maxDRGenHGenB_all'].fill(
                    dataset=dataset,
                    Mass=(ak.flatten(genHiggs.mass[sel_GenHToAATo4B])),
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )

                output['hMassGenAHeavy_vs_maxDRGenHGenB_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[idxGenA_sortByMass][:, 0].mass[sel_GenHToAATo4B]),
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )

                output['hMassGenALight_vs_maxDRGenHGenB_all'].fill(
                    dataset=dataset,
                    Mass=(genACollection[idxGenA_sortByMass][:, 1].mass[sel_GenHToAATo4B]),
                    deltaR=(max_dr_GenH_GenB[sel_GenHToAATo4B]),
                    systematic=syst,
                    weight=evtWeight_gen
                )
                

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


    if len(sys.argv) != 2:
        print("htoaa_Analysis:: Command-line config file missing.. \t **** ERROR **** \n")

    sConfig = sys.argv[1]

    
    config = GetDictFromJsonFile(sConfig)
    print("Config {}: \n{}".format(sConfig, json.dumps(config, indent=4)))
    
    sInputFiles         = config["inputFiles"]
    sOutputFile         = config["outputFile"]
    sample_category     = config['sampleCategory']
    isMC                = config["isMC"]
    era                 = config['era']
    luminosity          = Luminosities[era][0]
    sample_crossSection = config["crossSection"]
    sample_nEvents      = config["nEvents"]
    sample_sumEvents    = config["sumEvents"] if config["sumEvents"] != -1 else sample_nEvents
    if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
    lumiScale = calculate_lumiScale(luminosity=luminosity, crossSection=sample_crossSection, sumEvents=sample_sumEvents)
    #branchesToRead = htoaa_nanoAODBranchesToRead
    #print("branchesToRead: {}".format(branchesToRead))

    print(f"isMC: {isMC}, lumiScale: {lumiScale}")
    sInputFiles_toUse = []
    for sInputFile in sInputFiles:
        if "*" in sInputFile:  sInputFiles_toUse.extend( glob.glob( sInputFile ) )
        else:                  sInputFiles_toUse.append( sInputFile )
    sInputFiles = sInputFiles_toUse
    for iFile in range(len(sInputFiles)):
        if sInputFiles[iFile].startswith("/store/"): # LFN: Logical File Name
            sInputFiles[iFile] = xrootd_redirectorName + sInputFiles[iFile]
    print(f"sInputFiles ({len(sInputFiles)}) (type {type(sInputFiles)}): {sInputFiles}"); sys.stdout.flush()


    startTime = time.time()
    
    tracemalloc.start()
    



    #client = Client("tls://localhost:8786")
    #executor = processor.DaskExecutor(client=client)
    chunksize = nEventToReadInBatch
    maxchunks = None if nEventsToAnalyze == -1 else int(nEventsToAnalyze/nEventToReadInBatch)
    print(f"nEventsToAnalyze: {nEventsToAnalyze},  nEventToReadInBatch: {nEventToReadInBatch}, chunksize: {chunksize},  maxchunks: {maxchunks}")
    run = processor.Runner(
        #executor=executor,
        executor=processor.FuturesExecutor(workers=4),
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
                "era": era, 
                sample_category: {"isMC": isMC, "lumiScale": lumiScale}
            }
        )
    )

    print(f"metrics: {metrics}")


    if 'cutflow' in output.keys():
        print("Cutflow::")
        for key, value in output['cutflow'].items():
            print(key, value)    


    
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
    
