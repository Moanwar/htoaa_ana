import os
import numpy as np
from collections import OrderedDict as OD

#import import_ipynb
#import nbformat
#from importnb import imports
#with imports("ipynb"):
#    from PlotHistos1D_DataVsMC import era

sXRange = "xAxisRange"; sYRange = "yAxisRange";
sXLabel = 'xAxisLabel'; sYLabel = 'yAxisLabel';
sNRebinX = 'nRebinX';  sNRebinY = 'nRebinY'; 
sXRebinning = 'xRebinning'; sYRebinning = 'yRebinning'; 

'''
with open('PlotHistos1D_DataVsMC.ipynb') as f_:
    nbPlotHistos1D_DataVsMC = nbformat.read(f_, as_version=4)

print(f"HistogramListForPlottingDataVsMC_Analysis_GGFMode.py:: {nbPlotHistos1D_DataVsMC['cells'] = }", flush=True)
print(f"HistogramListForPlottingDataVsMC_Analysis_GGFMode.py:: {nbPlotHistos1D_DataVsMC['cells'][0]['source'] = }", flush=True)
'''

ExpData_dict = {
    'Data': ['JetHT_Run2018A', 'JetHT_Run2018B', 'JetHT_Run2018C', 'JetHT_Run2018D'],
    #'Data A': ['JetHT_Run2018A'],
    #'Data B': ['JetHT_Run2018B'],
    #'Data C': ['JetHT_Run2018C'],
    #'Data D': ['JetHT_Run2018D']
}
#ExpDatasetName = 'JetHT'
#ExpDatasetNames = ['JetHT', 'BTagCSV']
ExpData_dict = {}
#MCBkg_list = [
#    'QCD_0bCat', 'QCD_1bCat', 'QCD_2bCat', 'QCD_3bCat', 'QCD_4bCat', 'QCD_5bAndMoreCat',  
#    'TTToHadronic_powheg', 'TTToSemiLeptonic_powheg', 'TTTo2L2Nu_powheg', "SingleTop", 
#    'ZJetsToQQ_HT', "DYJets_M-50_Incl_NLO", 
#    'WJetsToQQ_HT', 'WJetsToLNu_HT_LO'
#    ]
#MCBkg_list = ['QCD_0bCat', 'QCD_1bCat', 'QCD_2bCat', 'QCD_3bCat', 'QCD_4bCat', 'QCD_5bAndMoreCat',  ]
MCBkg_list_1 = [
    "QCD_bEnr", "QCD_BGen", "QCD_Incl", 
    "TT0l", "TT1l", "TT2l", 
    "STop_t", "STbar_t", "ST_s_0l", "ST_s_1l", "STop_tW_Incl", "STbar_tW_Incl", #"STop_tW_12l", "STbar_tW_12l", 
    #"ttZ", "ttW", "tZq", 
    "Zqq", "Zvv", "Zll", "Wqq", "Wlv", 
    "ZZ", "WZ", "WW", #"ZZZ", "WZZ", "WWZ", "WWW", 
    #"ggH", #"VBFH", "VBFH_dipoleRecoilOn", "VBFWH_dipoleRecoilOn", 
    #"WHbbqq", "WHbblv", 
    #"ZH", #"ggZH", 
    #"ttH", #"", "", "", 
    #"", "", "", "", "", "",     
]
MCBkg_dict = {
    'QCD': ["QCD_bEnr", "QCD_BGen", "QCD_Incl"],
    #r't$\bar{t}$+X': ["TT0l", "TT1l", "TT2l"],
    r't$\bar{t}$+X (0l)': ["TT0l"],
    r't$\bar{t}$+X (1l)': ["TT1l"],
    r't$\bar{t}$+X (2l)': ["TT2l"],
    'Single top': ["STop_t", "STbar_t", "ST_s_0l", "ST_s_1l", "STop_tW_Incl", "STbar_tW_Incl"], #"STop_tW_12l", "STbar_tW_12l"],
    'V+X': ["Zqq", "Zvv", "Zll", "Wqq", "Wlv"],
    'Diboson': ["ZZ", "WZ", "WW"],
    #'': [],
}
MCSig_list = [
    #'SUSY_GluGluH_01J_HToAATo4B_M-15_HPtAbv150', 
    #'SUSY_GluGluH_01J_HToAATo4B_M-20_HPtAbv150', 
    #'SUSY_GluGluH_01J_HToAATo4B_M-25_HPtAbv150', 
    #'SUSY_GluGluH_01J_HToAATo4B_M-30_HPtAbv150', 
    #'SUSY_GluGluH_01J_HToAATo4B_M-50_HPtAbv150', 
    #'SUSY_GluGluH_01J_HToAATo4B_M-55_HPtAbv150', 
    #'ggHtoaato4b_mA_20'
    'ggHtoaato4b_mA_15',
    'ggHtoaato4b_mA_30',
    'ggHtoaato4b_mA_55',
    
]
MCSig_list = []
sLableSig = [
    #'HToAATo4B_M-15', 
    #'HToAATo4B_M-20', 
    #'HToAATo4B_M-25', 
    #'HToAATo4B_M-30', 
    #'HToAATo4B_M-50',
    #'HToAATo4B_M-55',
    #'ggHtoaato4b_mA_20'
    r'ggH signal $m_a=15$  GeV', #'ggHtoaato4b_mA_15',
    r'ggH signal $m_a=30$  GeV',
    r'ggH signal $m_a=55$  GeV',
]
sLableSig = []
systematics_list = ['Nom'] #['central'] ['Nom'] 
systematics_forData = 'noweight'
#selectionTags = ['SRWP40', 'SBWP80to40'] # ['SRWP40_mA30Window'] ['SRWP40'] ['SRWP40_mA55Window']  ['sel_leadingFatJetMSoftDrop', 'sel_leadingFatJetParticleNetMD_XbbvsQCD', 'SR'] #['SR', 'sel_leadingFatJetMSoftDrop', 'sel_leadingFatJetParticleNetMD_XbbvsQCD', 'sel_2018HEM1516Issue']

#HLT_toUse = 'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4'
HLT_toUse = 'Trg_Combo_AK4AK8Jet_HT'

scale_MCSig = 5 #20 # 1# #15 # 2 #5 #50 #1000
yRatioLimit = [0.4, 1.6]
ySignfLimit = [1e-1, 5e2]

logYMinScaleFactor = 10 # 100 # 1 # scale yMin by factor logYMinScaleFactor to not concentrate lowest stats background processes



histograms_dict = OD([
    #("hCutFlowPerCat", {sXLabel: r'Event cut flow', sYLabel: 'Events (unweighted)', sXRange: [-0.5, 20.5], sNRebinX: 1 }),
    #("hCutFlowPerCatWeighted", {sXLabel: r'Event cut flow', sYLabel: 'Events', sXRange: [-0.5, 20.5], sNRebinX: 1 }),

    ("hLeadingFatJetPt", {sXLabel: r'$p_{T}$(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    ("hLeadingFatJetEta", {sXLabel: r'$\eta$(Higgs candidate AK8 jet)', sYLabel: 'Events', sXRange: [-3.5, 3.5], sNRebinX: 2 }),
    ("hLeadingFatJetPhi", {sXLabel: r'$\phi$(Higgs candidate AK8 jet)', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),
    ("hLeadingFatJetMass", {sXLabel: r'Mass(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events', sXRange: [0, 250], sNRebinX: 5}),
    ("hLeadingFatJetMSoftDrop", {sXLabel: r'Mass$_{Soft\, drop}$(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events', sXRange: [0, 250], sNRebinX: 5 }),
    

    ("hMET_pT", {sXLabel: 'hMET_pT', sYLabel: 'Events', sXRange: [0, 1000], sNRebinX: 5 }),
    ("hPuppiMET_pT", {sXLabel: 'hPuppiMET_pT', sYLabel: 'Events', sXRange: [0, 1000], sNRebinX: 5 }),
    ("hHTtrig", {sXLabel: 'hHTtrig', sYLabel: 'Events', sXRange: [0, 2000], sNRebinX: 15 }),
    
    ("hLeadingFatJetParticleNetMD_XbbOverQCD", {sXLabel: 'hLeadingFatJetParticleNetMD_XbbOverQCD', sYLabel: 'Events', sXRange: [0.75, 1], sNRebinX: 1 }),
    ("hLeadingFatJetPNet_X4b_v2ab_Haa34b_score", {sXLabel: r'$X\to 3,4b$ tagger score', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 50 }),
    
    ("hLeadingFatJetMassH_v2b", {sXLabel: r'Mass$_{PNet\, X\to 4b}$(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events', sXRange: [0, 250], sNRebinX: 5}),
    
    ("hLeadingFatJetPNet_massAa", {sXLabel: 'hLeadingFatJetPNet_massAa', sYLabel: 'Events', sXRange: [5, 70], sNRebinX: 20}),
    ("hLeadingFatJetPNet_34massAa", {sXLabel: r'Mass$_{version\, a}$(a)  [GeV]', sYLabel: 'Events', sXRange: [5, 70], sNRebinX: 20}),
    ("hLeadingFatJetPNet_34massAb", {sXLabel: 'hLeadingFatJetPNet_34massAb', sYLabel: 'Events', sXRange: [5, 70], sNRebinX: 20}),
    ("hLeadingFatJetPNet_34massAd", {sXLabel: r'Mass$_{version\, d}$(a)  [GeV]', sYLabel: 'Events', sXRange: [5, 70], sNRebinX: 20}),
    
])

'''
histograms_dict = OD([
    #("hLeadingFatJetPt", {sXLabel: r'$p_{T}$(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    ("hHTtrig", {sXLabel: 'hHTtrig', sYLabel: 'Events', sXRange: [0, 2000], sNRebinX: 50 }),
    
])    
'''