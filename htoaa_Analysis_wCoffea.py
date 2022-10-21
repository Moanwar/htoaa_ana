#htoaa analysis main code

import os
import sys
import json
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
from coffea.nanoevents.methods import nanoaod
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


class ObjectSelection:
    def __init__(self, era):
        self.era = era
        

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
        tagger_btagDeepB = 'DeepCSV'
        wp_btagDeepB = 'M'

        FatJetPtThsh  = 170
        FatJetEtaThsh = 2.4

        maskSelFatJets = (
            (events.FatJet.pt > FatJetPtThsh) &
            (abs(events.FatJet.eta) < FatJetEtaThsh) &
            (events.FatJet.btagDeepB > bTagWPs[self.era][tagger_btagDeepB][wp_btagDeepB])
        )
        if printLevel >= 15:
            #print(f"era: {self.era}, bTagWPs[self.era]: {bTagWPs[self.era]}")
            print(f"selectFatJets()::maskSelFatJets {len(maskSelFatJets)}: {maskSelFatJets.to_list()}")
        return events.FatJet[maskSelFatJets]


    
    
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
        eta_axis      = hist.Bin("Eta",       r"$\eta$",          100, -5, 5)
        phi_axis      = hist.Bin("Phi",       r"$\phi$",          100, -3.14, 3.13)
        mass_axis     = hist.Bin("Mass",      r"$m$ [GeV]",       200, 0, 600)
        mlScore_axis  = hist.Bin("MLScore",   r"ML score",        100, -1.1, 1.1)
        jetN2_axis    = hist.Bin("N2",        r"N2b1",            100, 0, 3)
        jetN3_axis    = hist.Bin("N3",        r"N3b1",            100, 0, 5)
        jetTau_axis   = hist.Bin("TauN",      r"TauN",            100, 0, 1)

        sXaxis      = 'xAxis'
        sXaxisLabel = 'xAxisLabel'
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
        ])

        
        

        self._accumulator = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int)
        })

        for histName, histAttributes in histos.items():
            #hXaxis = histAttributes[sXaxis].copy()
            hXaxis = deepcopy(histAttributes[sXaxis])
            hXaxis.label = histAttributes[sXaxisLabel]
            
            self._accumulator.add({
                histName: hist.Hist(
                    "Counts",
                    dataset_axis,
                    hXaxis, #nObject_axis,
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
            'hLeadingFatJetEta': hist.Hist(
                "Counts",
                dataset_axis,
                eta_axis,
                systematic_axis,
            ),
            'hLeadingFatJetPhi': hist.Hist(
                "Counts",
                dataset_axis,
                phi_axis,
                systematic_axis,
            ),
             'hLeadingFatJetMass': hist.Hist(
                "Counts",
                dataset_axis,
                mass_axis,
                systematic_axis,
            ),
             'hLeadingFatJetMSoftDrop': hist.Hist(
                "Counts",
                dataset_axis,
                mass_axis,
                systematic_axis,
            ),
             'hLeadingFatJetBtagDeepB': hist.Hist(
                "Counts",
                dataset_axis,
                mlScore_axis,
                systematic_axis,
            ),
             'hLeadingFatJetBtagDDBvLV2': hist.Hist(
                "Counts",
                dataset_axis,
                mlScore_axis,
                systematic_axis,
            ),
             'hLeadingFatJetBtagDDCvBV2': hist.Hist(
                "Counts",
                dataset_axis,
                mlScore_axis,
                systematic_axis,
            ),
             'hLeadingFatJetBtagHbb': hist.Hist(
                "Counts",
                dataset_axis,
                mlScore_axis,
                systematic_axis,
            ),
             'hLeadingFatJetDeepTagMD_H4qvsQCD': hist.Hist(
                "Counts",
                dataset_axis,
                mlScore_axis,
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
        nFatJetMin = 1

        output = self.accumulator.identity()
        dataset = events.metadata["dataset"] # dataset label
        #print(f"dataset: {dataset}")


        
        ##################
        # OBJECT SELECTION
        ##################

        # FatJet selection
        selFatJet = self.objectSelector.selectFatJets(events)
        

        #####################
        # EVENT SELECTION
        #####################

        # create a PackedSelection object
        # this will help us later in composing the boolean selections easily
        selection = PackedSelection()

        if printLevel >= 5:
            print(f"events.PV.npvsGood.to_list(): {events.PV.npvsGood.to_list()}")
        # nPVGood >= 1
        selection.add("nPV", events.PV.npvsGood >= 1)
        
        # >=1 FatJet
        selection.add("FatJetGe1", ak.num(selFatJet) >= nFatJetMin)

        sel_SR = selection.all("nPV", "FatJetGe1")

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

        ##################
        # EVENT VARIABLES
        ##################

        leadingFatJet = ak.firsts(selFatJet)


        
        ################
        # EVENT WEIGHTS
        ################

        # create a processor Weights object, with the same length as the number of events in the chunk
        weights = Weights(len(events))


        if self.datasetInfo[dataset]["isMC"]:
            weights.add(
                "lumiWeight",
                weight=np.full(len(events), self.datasetInfo[dataset]["lumiScale"])
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

                

            weights.add(
                "genWeight",
                weight=np.copysign(np.ones(len(events)), events.genWeight)
            )

            weights.add(
                "btagWeight",
                weight=(events.btagWeight.DeepCSVB)
            )
                
        
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

        for syst in systList:

            # find the event weight to be used when filling the histograms
            weightSyst = syst
            
            # in the case of 'central', or the jet energy systematics, no weight systematic variation is used (weightSyst=None)
            if syst in ["central", "JERUp", "JERDown", "JESUp", "JESDown"]:
                weightSyst = None

            if syst == "noweight":
                evtWeight = np.ones(len(events))
            else:
                evtWeight = weights.weight(weightSyst)


            if printLevel >= 5:
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


            
            
        '''
       
        nFatJet = 1
        FatJetPtThsh = 220
        FatJetEtaThsh = 2.4
        
        dataset = events.metadata["dataset"]
        print(f"dataset: {dataset}")
    
        # Keep track of muons and electrons by tagging them 0/1.
        #electrons = ak.with_field(events.Electron, 11, 'flavor')
        #muons     = ak.with_field(events.Muon, 13, 'flavor')
        
        if printLevel >= 5:
            print(f"events.fields: {events.fields}")
            print(f"events.FatJet.fields: {events.FatJet.fields}")
            print(f"events.FatJet.pt: {events.FatJet.pt}")


        self.output['cutflow']['All events'] += len(events)


        self.output['hnFatJet_level0'].fill( ak.to_numpy(ak.num(events.FatJet)) )
        
        #events = events[(
        #    (ak.num(events.FatJet) >= 1) )]
        #self.output['cutflow']['nFatJet >= 1'] += len(events)

        if printLevel >= 5:
            print(f"events.FatJet.pt _0 ({len(events)}): {events.FatJet.pt.to_list()}")

        #events = events[(
        #    (events.FatJet.pt[:, 0] > FatJetPtThsh) )]
        #self.output['cutflow']['LeadingFatJetPt > %s' % (FatJetPtThsh)] += len(events)

        if printLevel >= 5:
            print(f"\n\nevents.FatJet.pt _1 ({len(events)}): {events.FatJet.pt.to_list()}")
            cut1  = (events.FatJet.pt > FatJetPtThsh)
            print(f"\n\ncut1 ({len(cut1)}): {cut1.to_list()}")
            #events1 = events.FatJet[cut1]
            events1 = events[ak.any(cut1, axis=1)]
            print(f"\n\nevents1 ({len(events1)}): {events1.to_list()}")
            print(f"\n\nevents1.FatJet.pt _1 ({len(events1)}): {events1.FatJet.pt.to_list()}")

        events = events[(
            (events.FatJet.pt > FatJetPtThsh) )]
        self.output['cutflow']['LeadingFatJetPt _1 > %s' % (FatJetPtThsh)] += len(events)
        
        if printLevel >= 5:
            print(f"events.FatJet.pt _2 ({len(events)}): {events.FatJet.pt.to_list()}")

        
        events = events[(
            (abs(events.FatJet.eta) < 2.4) )]
        self.output['cutflow']['FatJetAbsEta < 2.4'] += len(events)

        if printLevel >= 5:
            print(f"nEvents FatJetEta cut ({len(events)}) ")
            print(f"events: {events.to_list()}")
            
        self.output['hnFatJet_level1'].fill( ak.to_numpy(ak.num(events.FatJet)) )
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
    
    sInputFiles  = config["inputFiles"]
    sOutputFile  = config["outputFile"]
    sample_category = config['sampleCategory']
    isMC = config["isMC"]
    era = config['era']
    luminosity = Luminosities[era][0]
    sample_crossSection = config["crossSection"]
    sample_nEvents = config["nEvents"]
    sample_sumEvents = config["sumEvents"] if config["sumEvents"] != -1 else sample_nEvents
    if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
    lumiScale = calculate_lumiScale(luminosity=luminosity, crossSection=sample_crossSection, sumEvents=sample_sumEvents)
    #branchesToRead = htoaa_nanoAODBranchesToRead
    #print("branchesToRead: {}".format(branchesToRead))

    print(f"isMC: {isMC}, lumiScale: {lumiScale}")
    for iFile in range(len(sInputFiles)):
        if sInputFiles[iFile].startswith("/store/"): # LFN: Logical File Name
            sInputFiles[iFile] = xrootd_redirectorName + sInputFiles[iFile]
    print(f"sInputFiles ({len(sInputFiles)}): {sInputFiles}"); sys.stdout.flush()

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
        sOutputFile = sOutputFile.replace('.root', '_wCoffea.root') # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        sDir1 = 'evt/%s' % (sample_category)

        
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                print(f"key: {key},  value ({type(value)}): {value}")
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
    
