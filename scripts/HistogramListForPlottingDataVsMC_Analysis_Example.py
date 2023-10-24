import os
import numpy as np
from collections import OrderedDict as OD

sXRange = "xAxisRange"; sYRange = "yAxisRange";
sXLabel = 'xAxisLabel'; sYLabel = 'yAxisLabel';
sNRebinX = 'nRebinX';  sNRebinY = 'nRebinY'; 


ExpData_dict = {
    'Data ABCD': ['SingleMuon_Run2018A', 'SingleMuon_Run2018B', 'SingleMuon_Run2018C', 'SingleMuon_Run2018D'],
    #'Data A': ['SingleMuon_Run2018A'],
    #'Data B': ['SingleMuon_Run2018B'],
    #'Data C': ['SingleMuon_Run2018C'],
    #'Data D': ['SingleMuon_Run2018D']
}
MCBkg_list = [
    'QCD_0bCat', 'QCD_1bCat', 'QCD_2bCat', 'QCD_3bCat', 'QCD_4bCat', 'QCD_5bAndMoreCat',  
    'TTToHadronic_powheg', 'TTToSemiLeptonic_powheg', 'TTTo2L2Nu_powheg', "SingleTop", 
    'ZJetsToQQ_HT', "DYJets_M-10to50_Incl_NLO", "DYJets_M-50_Incl_NLO", 
    'WJetsToQQ_HT', 'WJetsToLNu_Incl_NLO', 
]
MCSig_list = [] #['SUSY_GluGluH_01J_HToAATo4B_M-20_HPtAbv150', 'SUSY_GluGluH_01J_HToAATo4B_M-30_HPtAbv150', ]
sLableSig = [] #['HToAATo4B_M-20', 'HToAATo4B_M-30']
systematics_list = ['central']
systematics_forData = 'noweight'
selectionTags = ['sel_JetID', 'sel_L1_SingleJet180', 'SR'] # ['sel_JetID', 'sel_L1_SingleJet180', 'SR'] # ['sel_leadingFatJetMSoftDrop', 'sel_leadingFatJetParticleNetMD_XbbvsQCD', 'SR'] #['SR', 'sel_leadingFatJetMSoftDrop', 'sel_leadingFatJetParticleNetMD_XbbvsQCD', 'sel_2018HEM1516Issue']

HLT_toUse = 'HLT_IsoMu27'

scale_MCSig = 50 #1000
yRatioLimit = [0.4, 1.6]

logYMinScaleFactor = 1 # 10 # 100 # 1 # scale yMin by factor logYMinScaleFactor to not concentrate lowest stats background processes



histograms_dict = OD([
    # format:
    # ("HistogramName", {sXLabel: 'xAxisLabel', sYLabel: 'yAxisLabel', sXRange: [xMin, xMax], sNRebinX: NoOfBinsToRebin }),

    #("hLeadingMuonPt", {sXLabel: 'hLeadingMuonPt', sYLabel: 'Events', sXRange: [0, 500], sNRebinX: 4 }),
    #("hLeadingMuonEta", {sXLabel: 'hLeadingMuonEta', sYLabel: 'Events', sXRange: [-3.5, 3.5], sNRebinX: 2 }),
    #("hLeadingMuonPhi", {sXLabel: 'hLeadingMuonPhi', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),

    ("hLeadingFatJetPt", {sXLabel: 'hLeadingFatJetPt', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    ("hLeadingFatJetEta", {sXLabel: 'hLeadingFatJetEta', sYLabel: 'Events', sXRange: [-3.5, 3.5], sNRebinX: 2 }),
    ("hLeadingFatJetPhi", {sXLabel: 'hLeadingFatJetPhi', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),
    #("hLeadingFatJetMass", {sXLabel: 'hLeadingFatJetMass', sYLabel: 'Events', sXRange: [50, 300], sNRebinX: 5}),
    ("hLeadingFatJetMSoftDrop", {sXLabel: 'hLeadingFatJetMSoftDrop', sYLabel: 'Events', sXRange: [50, 300], sNRebinX: 5 }),
    #("hLeadingFatJetParticleNetMD_XbbOverQCD", {sXLabel: 'hLeadingFatJetParticleNetMD_XbbOverQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }), 
])


