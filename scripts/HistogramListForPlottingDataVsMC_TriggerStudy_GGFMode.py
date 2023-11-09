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
    'ZJetsToQQ_HT', 
    #"DYJets_M-10to50_Incl_LO", "DYJets_M-50_Incl_LO", #"DYJets_M-10to50_Incl_NLO", "DYJets_M-50_Incl_NLO", 
    "DYJets_M-10to50_Incl_LO", "DYJets_HT_LO", 
    'WJetsToQQ_HT', 'WJetsToLNu_HT_LO', #'WJetsToLNu_Incl_NLO', 
]
MCSig_list = [] #['SUSY_GluGluH_01J_HToAATo4B_M-20_HPtAbv150', 'SUSY_GluGluH_01J_HToAATo4B_M-30_HPtAbv150', ]
sLableSig = [] #['HToAATo4B_M-20', 'HToAATo4B_M-30']
systematics_list = ['central']
systematics_forData = 'noweight'
selectionTags = ['SR', 'SR_TrgAK8330_M30_BDBnp4', 'SR_Trg2AK4116_DCSVp71', 'SR_TrgAK8400_M30', 'SR_TrgAK8500', 'SR_TrgComb2', 'SR_TrgComb4' ] # ['sel_JetID', 'sel_lFJPNetXbbPlusDZHbb', 'sel_L1_SingleJet180', 'SR'] # ['sel_JetID', 'sel_L1_SingleJet180', 'SR'] #['sel_HLT_IsoMu27'] # ['sel_JetID', 'sel_L1_SingleJet180', 'SR'] # ['sel_JetID', 'sel_L1_SingleJet180', 'SR'] # ['sel_leadingFatJetMSoftDrop', 'sel_leadingFatJetParticleNetMD_XbbvsQCD', 'SR'] #['SR', 'sel_leadingFatJetMSoftDrop', 'sel_leadingFatJetParticleNetMD_XbbvsQCD', 'sel_2018HEM1516Issue']

HLT_toUse = 'HLT_IsoMu24'

scale_MCSig = 50 #1000 
yRatioLimit = [0.4, 1.6]

logYMinScaleFactor = 1 # 10 # 100 # 1 # scale yMin by factor logYMinScaleFactor to not concentrate lowest stats background processes



histograms_dict = OD([
    ("hLeadingMuonPt", {sXLabel: 'hLeadingMuonPt', sYLabel: 'Events', sXRange: [0, 500], sNRebinX: 4 }),
    ("hLeadingMuonEta", {sXLabel: 'hLeadingMuonEta', sYLabel: 'Events', sXRange: [-3.5, 3.5], sNRebinX: 2 }),
    ("hLeadingMuonPhi", {sXLabel: 'hLeadingMuonPhi', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),

    ("hdR_leadingMuon_leadingFatJet", {sXLabel: 'hdR_leadingMuon_leadingFatJet', sYLabel: 'Events',  sNRebinX: 2 }),
    
    ("hLeadingFatJetPt", {sXLabel: 'hLeadingFatJetPt', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    ("hLeadingFatJetPt_msoftdropGt60", {sXLabel: 'hLeadingFatJetPt_msoftdropGt60', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    ("hLeadingFatJetPt_msoftdropGt60_btagHbbGtnp1", {sXLabel: 'hLeadingFatJetPt_msoftdropGt60_btagHbbGtnp1', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    ("hLeadingFatJetEta", {sXLabel: 'hLeadingFatJetEta', sYLabel: 'Events', sXRange: [-3.5, 3.5], sNRebinX: 2 }),
    ("hLeadingFatJetPhi", {sXLabel: 'hLeadingFatJetPhi', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),
    ("hLeadingFatJetMass", {sXLabel: 'hLeadingFatJetMass', sYLabel: 'Events', sXRange: [0, 300], sNRebinX: 5}),
    ("hLeadingFatJetMSoftDrop", {sXLabel: 'hLeadingFatJetMSoftDrop', sYLabel: 'Events', sXRange: [0, 300], sNRebinX: 5 }),
    ("hLeadingFatJetMass_pTGt400_btagHbbGtnp1", {sXLabel: 'hLeadingFatJetMass_pTGt400_btagHbbGtnp1', sYLabel: 'Events', sXRange: [00, 300], sNRebinX: 5}),
    ("hLeadingFatJetMSoftDrop_pTGt400", {sXLabel: 'hLeadingFatJetMSoftDrop_pTGt400', sYLabel: 'Events', sXRange: [00, 300], sNRebinX: 5}),
    ("hLeadingFatJetMSoftDrop_pTGt400_btagHbbGtnp1", {sXLabel: 'hLeadingFatJetMSoftDrop_pTGt400_btagHbbGtnp1', sYLabel: 'Events', sXRange: [00, 300], sNRebinX: 5}),
    
    ("hLeadingFatJetParticleNetMD_XbbOverQCD", {sXLabel: 'hLeadingFatJetParticleNetMD_XbbOverQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }), 
    ("hLeadingFatJetBtagCSVV2", {sXLabel: 'hLeadingFatJetBtagCSVV2', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }), 
    ("hLeadingFatJetBtagDDBvLV2", {sXLabel: 'hLeadingFatJetBtagDDBvLV2', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }), 
    ("hLeadingFatJetBtagDeepB", {sXLabel: 'hLeadingFatJetBtagDeepB', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }), 
    ("hLeadingFatJetBtagHbb", {sXLabel: 'hLeadingFatJetBtagHbb', sYLabel: 'Events', sXRange: [-1, 1], sNRebinX: 2 }), 
    ("hLeadingFatJetBtagHbb_pTGt400_msoftdropGt60", {sXLabel: 'hLeadingFatJetBtagHbb_pTGt400_msoftdropGt60', sYLabel: 'Events', sXRange: [-1, 1], sNRebinX: 2 }), 
    ("hLeadingFatJetParticleNetMD_XbbOverQCD_pTGt400_msoftdropGt60", {sXLabel: 'hLeadingFatJetParticleNetMD_XbbOverQCD_pTGt400_msoftdropGt60', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }), 
    
    
])

#histograms_dict = OD([
#    ("hdR_leadingMuon_leadingFatJet", {sXLabel: 'hdR_leadingMuon_leadingFatJet', sYLabel: 'Events',  sNRebinX: 2 }),
#    
#])