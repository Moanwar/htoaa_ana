import os
import numpy as np
from collections import OrderedDict as OD

sXRange = "xAxisRange"; sYRange = "yAxisRange";
sXLabel = 'xAxisLabel'; sYLabel = 'yAxisLabel';
sNRebinX = 'nRebinX';  sNRebinY = 'nRebinY'; 
sXRebinning = 'xRebinning'; sYRebinning = 'yRebinning'; 

ExpData_dict = {
    #'Data ABCD': ['JetHT_Run2018A', 'JetHT_Run2018B', 'JetHT_Run2018C', 'JetHT_Run2018D'],
    'Data': ['MET_Run2018A', 'MET_Run2018B', 'MET_Run2018C', 'MET_Run2018D'],
    #'Data A': ['JetHT_Run2018A'],
    #'Data B': ['JetHT_Run2018B'],
    #'Data C': ['JetHT_Run2018C'],
    #'Data D': ['JetHT_Run2018D']
}
#MCBkg_list = [
#    'QCD_0bCat', 'QCD_1bCat', 'QCD_2bCat', 'QCD_3bCat', 'QCD_4bCat', 'QCD_5bAndMoreCat',  
#    'TTToHadronic_powheg', 'TTToSemiLeptonic_powheg', 'TTTo2L2Nu_powheg', "SingleTop", 
#    'ZJetsToQQ_HT', "ZJetsToNuNuQ_HT", "DYJets_M-50_Incl_NLO", 
#    'WJetsToQQ_HT', 'WJetsToLNu_HT_LO',
#    'ZZ','WZ','WW', 'ZZZ','WZZ','WWZ','WWW',
#    'WH_HToBB_WToLNu',
#    'ZH_HToBB',
#    'ttH'
#]
MCBkg_list_1 = [
    "QCD_bEnr", "QCD_BGen", "QCD_Incl", 
    "TT0l", "TT1l", "TT2l", 
    "STop_t", "STbar_t", "ST_s_0l", "ST_s_1l", "STop_tW_Incl", "STbar_tW_Incl", #"STop_tW_12l", "STbar_tW_12l", 
    "ttZ", "ttW", "tZq", 
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
    r'$t\bar{t}$': ["TT0l", "TT1l", "TT2l"],
    't': ["STop_t", "STbar_t", "ST_s_0l", "ST_s_1l", "STop_tW_Incl", "STbar_tW_Incl"], # "STop_tW_12l", "STbar_tW_12l"],
    'V': ["Zqq", "Zvv", "Zll", "Wqq", "Wlv"],
    'VV': ["ZZ", "WZ", "WW"],
    #'': [],
}
MCBkg_dict = {
    'QCD': ["QCD_bEnr", "QCD_BGen", "QCD_Incl"],
    r't$\bar{t}$+X': ["TT0l", "TT1l", "TT2l"],
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
    #'SUSY_WH_WToAll_HToAATo4B_M-20_HPtAbv150',
    #'SUSY_ZH_ZToAll_HToAATo4B_M-20_HPtAbv150'
    'ZHtoaato4b_mA_15',
    'ZHtoaato4b_mA_30',
    'ZHtoaato4b_mA_55',
    
    ]
sLableSig = [
    #'HToAATo4B_M-15', 
    #'HToAATo4B_M-20', 
    #'HToAATo4B_M-25', 
    #'HToAATo4B_M-30', 
    #'HToAATo4B_M-50',
    #'HToAATo4B_M-55',
    #'WH_HToAATo4B_M-20',
    r'ZH signal $m_a=15$  GeV', #'ggHtoaato4b_mA_15',
    r'ZH signal $m_a=30$  GeV',
    r'ZH signal $m_a=55$  GeV',
     ]
systematics_list = ['Nom'] #['central'] ['Nom'] 
systematics_forData = 'noweight'
#selectionTags = ['SRWP80'] #['Presel'] #['SRWP60'] # ['SRWP80'] # ['SRWP40_mA30Window'] #['SRWP40'] # ['SRWP40_mA55Window']  ['sel_leadingFatJetMSoftDrop', 'sel_leadingFatJetParticleNetMD_XbbvsQCD', 'SR'] #['SR', 'sel_leadingFatJetMSoftDrop', 'sel_leadingFatJetParticleNetMD_XbbvsQCD', 'sel_2018HEM1516Issue']

#HLT_toUse = 'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4'
HLT_toUse = 'Trg_Combo_MET' # 'Trg_Combo_AK4AK8Jet_HT'

scale_MCSig = 5 #10 #50 #1000
yRatioLimit = [0., 2.0] #[0.4, 1.6]
ySignfLimit = [1e-1, 5e2]

logYMinScaleFactor = 10 # 100 # 1 # scale yMin by factor logYMinScaleFactor to not concentrate lowest stats background processes


histograms_dict = OD([
    #("hLeadingFatJetMass", {sXLabel: 'Leading FatJet mass [GeV]', sYLabel: 'Events', sXRange: [50, 250], sYRange: [1e-2, 1e8]})
    #("hLeadingFatJetMass", {sXLabel: 'Leading FatJet mass [GeV]', sYLabel: 'Events', sXRange: [50, 250]}),
    
    #("", {sXLabel: '', sYLabel: 'Events', sXRange: []}),
    
    #("hCutFlow", {sXLabel: 'hCutFlow', sYLabel: 'Events'}),
    #("hCutFlowWeighted", {sXLabel: 'hCutFlowWeighted', sYLabel: 'Events'}),

    #("hPV_npvs_beforeSel", {sXLabel: 'No. of primary vertices before selection', sYLabel: 'Events', sXRange: [0, 100] }),
    #("hPV_npvsGood_beforeSel", {sXLabel: 'No. of good primary vertices before selection', sYLabel: 'Events', sXRange: [0, 100] }),
    #("hPV_npvs_SR", {sXLabel: 'No. of primary vertices in SR', sYLabel: 'Events', sXRange: [0, 100] }),
    #("hPV_npvsGood_SR", {sXLabel: 'No. of good primary vertices in SR', sYLabel: 'Events', sXRange: [0, 100] }),

    ("hLeadingFatJetPt", {sXLabel: r'$p_{T}$(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    ("hLeadingFatJetEta", {sXLabel: r'$\eta$(Higgs candidate AK8 jet)', sYLabel: 'Events', sXRange: [-3.5, 3.5], sNRebinX: 2 }),
    ("hLeadingFatJetPhi", {sXLabel: r'$\phi$(Higgs candidate AK8 jet)', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),

    #("hLeadingFatJetPt_HEM1516IssueEtaPhiCut", {sXLabel: 'hLeadingFatJetPt_HEM1516IssueEtaPhiCut', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    #("hLeadingFatJetEta_HEM1516IssuePhiCut", {sXLabel: 'hLeadingFatJetEta_HEM1516IssuePhiCut', sYLabel: 'Events', sXRange: [-3.5, 3.5] }),
    #("hLeadingFatJetPhi_HEM1516IssueEtaCut", {sXLabel: 'hLeadingFatJetPhi_HEM1516IssueEtaCut', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),
    #("hLeadingFatJetPt_HEM1516IssueEtaPhiCut_DataPreHEM1516Issue", {sXLabel: 'hLeadingFatJetPt_HEM1516IssueEtaPhiCut_DataPreHEM1516Issue', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    #("hLeadingFatJetEta_HEM1516IssuePhiCut_DataPreHEM1516Issue", {sXLabel: 'hLeadingFatJetEta_HEM1516IssuePhiCut_DataPreHEM1516Issue', sYLabel: 'Events', sXRange: [-3.5, 3.5] }),
    #("hLeadingFatJetPhi_HEM1516IssueEtaCut_DataPreHEM1516Issue", {sXLabel: 'hLeadingFatJetPhi_HEM1516IssueEtaCut_DataPreHEM1516Issue', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),
    #("hLeadingFatJetPt_HEM1516IssueEtaPhiCut_DataWithHEM1516Issue", {sXLabel: 'hLeadingFatJetPt_HEM1516IssueEtaPhiCut_DataWithHEM1516Issue', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    #("hLeadingFatJetEta_HEM1516IssuePhiCut_DataWithHEM1516Issue", {sXLabel: 'hLeadingFatJetEta_HEM1516IssuePhiCut_DataWithHEM1516Issue', sYLabel: 'Events', sXRange: [-3.5, 3.5] }),
    #("hLeadingFatJetPhi_HEM1516IssueEtaCut_DataWithHEM1516Issue", {sXLabel: 'hLeadingFatJetPhi_HEM1516IssueEtaCut_DataWithHEM1516Issue', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),

    #("hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix", {sXLabel: 'hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4  }),
    #("hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix", {sXLabel: 'hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix', sYLabel: 'Events', sXRange: [-3.5, 3.5] }),
    #("hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix", {sXLabel: 'hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),
    #("hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix_DataPreHEM1516Issue", {sXLabel: 'hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix_DataPreHEM1516Issue', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4  }),
    #("hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix_DataPreHEM1516Issue", {sXLabel: 'hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix_DataPreHEM1516Issue', sYLabel: 'Events', sXRange: [-3.5, 3.5] }),
    #("hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix_DataPreHEM1516Issue", {sXLabel: 'hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix_DataPreHEM1516Issue', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),
    #("hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix_DataWithHEM1516Issue", {sXLabel: 'hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516Fix_DataWithHEM1516Issue', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    #("hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix_DataWithHEM1516Issue", {sXLabel: 'hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516Fix_DataWithHEM1516Issue', sYLabel: 'Events', sXRange: [-3.5, 3.5] }),
    #("hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix_DataWithHEM1516Issue", {sXLabel: 'hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516Fix_DataWithHEM1516Issue', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),

    # 2018 HEM15/16 issue validation
    #("hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt", {sXLabel: 'hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    #("hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt", {sXLabel: 'hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt', sYLabel: 'Events', sXRange: [-3.5, 3.5] }),
    #("hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt", {sXLabel: 'hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),
    #("hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt_DataPreHEM1516Issue", {sXLabel: 'hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt_DataPreHEM1516Issue', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    #("hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt_DataPreHEM1516Issue", {sXLabel: 'hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt_DataPreHEM1516Issue', sYLabel: 'Events', sXRange: [-3.5, 3.5] }),
    #("hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt_DataPreHEM1516Issue", {sXLabel: 'hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt_DataPreHEM1516Issue', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),
    #("hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt_DataWithHEM1516Issue", {sXLabel: 'hLeadingFatJetPt_HEM1516IssueEtaPhiCut_woHEM1516MCRewgt_DataWithHEM1516Issue', sYLabel: 'Events', sXRange: [180, 1000], sNRebinX: 4 }),
    #("hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt_DataWithHEM1516Issue", {sXLabel: 'hLeadingFatJetEta_HEM1516IssuePhiCut_woHEM1516MCRewgt_DataWithHEM1516Issue', sYLabel: 'Events', sXRange: [-3.5, 3.5] }),
    #("hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt_DataWithHEM1516Issue", {sXLabel: 'hLeadingFatJetPhi_HEM1516IssueEtaCut_woHEM1516MCRewgt_DataWithHEM1516Issue', sYLabel: 'Events', sXRange: [-3.14, 3.14], sNRebinX: 2 }),

    ("hLeadingFatJetMass", {sXLabel: r'Mass(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events', sXRange: [0, 250], sNRebinX: 5}),
    ("hLeadingFatJetMSoftDrop", {sXLabel: r'Mass$_{Soft\, drop}$(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events', sXRange: [0, 250], sNRebinX: 5 }),
    #("hLeadingFatJetBtagDeepB", {sXLabel: 'hLeadingFatJetBtagDeepB', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetBtagDDBvLV2", {sXLabel: 'hLeadingFatJetBtagDDBvLV2', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetBtagDDCvBV2", {sXLabel: 'hLeadingFatJetBtagDDCvBV2', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    
    #("hLeadingFatJetBtagHbb", {sXLabel: 'hLeadingFatJetBtagHbb', sYLabel: 'Events', sXRange: [-1, 1], sNRebinX: 2 }),
    #("hLeadingFatJetDeepTagMD_H4qvsQCD", {sXLabel: 'hLeadingFatJetDeepTagMD_H4qvsQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetDeepTagMD_HbbvsQCD", {sXLabel: 'hLeadingFatJetDeepTagMD_HbbvsQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetDeepTagMD_ZHbbvsQCD", {sXLabel: 'hLeadingFatJetDeepTagMD_ZHbbvsQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetDeepTagMD_ZHccvsQCD", {sXLabel: 'hLeadingFatJetDeepTagMD_ZHccvsQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    
    #("hLeadingFatJetDeepTagMD_ZbbvsQCD", {sXLabel: 'hLeadingFatJetDeepTagMD_ZbbvsQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetDeepTagMD_ZvsQCD", {sXLabel: 'hLeadingFatJetDeepTagMD_ZvsQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetDeepTagMD_bbvsLight", {sXLabel: 'hLeadingFatJetDeepTagMD_bbvsLight', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetDeepTagMD_ccvsLight", {sXLabel: 'hLeadingFatJetDeepTagMD_ccvsLight', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetDeepTag_H", {sXLabel: 'hLeadingFatJetDeepTag_H', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    
    #("hLeadingFatJetDeepTag_QCD", {sXLabel: 'hLeadingFatJetDeepTag_QCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetDeepTag_QCDothers", {sXLabel: 'hLeadingFatJetDeepTag_QCDothers', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetN2b1", {sXLabel: 'hLeadingFatJetN2b1', sYLabel: 'Events', sXRange: [0, 0.6], sNRebinX: 2 }),
    #("hLeadingFatJetN3b1", {sXLabel: 'hLeadingFatJetN3b1', sYLabel: 'Events', sXRange: [0.5, 3.5], sNRebinX: 2 }),
    #("hLeadingFatJetTau1", {sXLabel: 'hLeadingFatJetTau1', sYLabel: 'Events', sXRange: [0, 0.6], sNRebinX: 2 }),    
    #("hLeadingFatJetTau2", {sXLabel: 'hLeadingFatJetTau2', sYLabel: 'Events', sXRange: [0, 0.5], sNRebinX: 2 }),
    #("hLeadingFatJetTau3", {sXLabel: 'hLeadingFatJetTau3', sYLabel: 'Events', sXRange: [0, 0.3], sNRebinX: 2 }),
    #("hLeadingFatJetTau4", {sXLabel: 'hLeadingFatJetTau4', sYLabel: 'Events', sXRange: [0, 0.4], sNRebinX: 2 }),
    #("hLeadingFatJetTau4by3", {sXLabel: 'hLeadingFatJetTau4by3', sYLabel: 'Events', sXRange: [0.2, 1], sNRebinX: 2 }),
    #("hLeadingFatJetTau3by2", {sXLabel: 'hLeadingFatJetTau3by2', sYLabel: 'Events', sXRange: [0.2, 1], sNRebinX: 2 }),    
    #("hLeadingFatJetTau2by1", {sXLabel: 'hLeadingFatJetTau2by1', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),

    #("hLeadingFatJetNConstituents", {sXLabel: 'hLeadingFatJetNConstituents', sYLabel: 'Events'}),
    #("hLeadingFatJetNBHadrons", {sXLabel: 'hLeadingFatJetNBHadrons', sYLabel: 'Events', sXRange: [-0.5, 10.5]}),
    #("hLeadingFatJetNCHadrons", {sXLabel: 'hLeadingFatJetNCHadrons', sYLabel: 'Events', sXRange: [-0.5, 10.5]}),

    #("hLeadingFatJetParticleNetMD_QCD", {sXLabel: 'hLeadingFatJetParticleNetMD_QCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetParticleNetMD_Xbb", {sXLabel: 'hLeadingFatJetParticleNetMD_Xbb', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetParticleNetMD_Xcc", {sXLabel: 'hLeadingFatJetParticleNetMD_Xcc', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),       
    #("hLeadingFatJetParticleNetMD_Xqq", {sXLabel: 'hLeadingFatJetParticleNetMD_Xqq', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2}),

    #("hLeadingFatJetParticleNetMD_XbbOverQCD", {sXLabel: 'hLeadingFatJetParticleNetMD_XbbOverQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }), 
    #("hLeadingFatJetParticleNetMD_XccOverQCD", {sXLabel: 'hLeadingFatJetParticleNetMD_XccOverQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }), 
    #("hLeadingFatJetParticleNetMD_XqqOverQCD", {sXLabel: 'hLeadingFatJetParticleNetMD_XqqOverQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }), 

    #("hLeadingFatJetParticleNet_H4qvsQCD", {sXLabel: 'hLeadingFatJetParticleNet_H4qvsQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetParticleNet_HbbvsQCD", {sXLabel: 'hLeadingFatJetParticleNet_HbbvsQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetParticleNet_HccvsQCD", {sXLabel: 'hLeadingFatJetParticleNet_HccvsQCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    #("hLeadingFatJetParticleNet_QCD", {sXLabel: 'hLeadingFatJetParticleNet_QCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 2 }),
    
    #("hLeadingFatJetParticleNet_mass", {sXLabel: 'hLeadingFatJetParticleNet_mass', sYLabel: 'Events', sXRange: [0, 300], sNRebinX: 5 }),

    #("hLeadingFatJet_nSubJets", {sXLabel: 'hLeadingFatJet_nSubJets', sYLabel: 'Events', sXRange: [-0.5, 6.5] }),
    #("hLeadingFatJet_nSubJets_bTag_L", {sXLabel: 'hLeadingFatJet_nSubJets_bTag_L', sYLabel: 'Events', sXRange: [-0.5, 10.5] }),
    #("hLeadingFatJet_nSubJets_bTag_M", {sXLabel: 'hLeadingFatJet_nSubJets_bTag_M', sYLabel: 'Events', sXRange: [-0.5, 10.5] }),
    #("hLeadingFatJet_nSV", {sXLabel: 'hLeadingFatJet_nSV', sYLabel: 'Events', sXRange: [-0.5, 10.5] }),

    ("hMET_pT", {sXLabel: r'$E_{T}^{miss}$', sYLabel: 'Events', sXRange: [0, 1000], sNRebinX: 5 }),
    ("hPuppiMET_pT", {sXLabel: 'hPuppiMET_pT', sYLabel: 'Events', sXRange: [0, 1000], sNRebinX: 5 }),
    #("hMET_sumEt", {sXLabel: 'hMET_sumEt', sYLabel: 'Events', sXRange: [1000, 4000], sNRebinX: 5 }),
    ("hMETPhi", {sXLabel: r'$\phi$($\vec{E}_{T}^{miss}$)', sYLabel: 'Events', sXRange: [-3.5, 3.5], sNRebinX: 2 }),
    ("hdPhi_MET_leadingFatJet", {sXLabel: r'$\Delta\phi$(H candidiate AK8 jet and $\vec{E}_{T}^{miss}$)', sYLabel: 'Events', sXRange: [2.2, 3.14], sNRebinX: 10 }),

    #("hLeadingFatJet_nLeptons", {sXLabel: 'hLeadingFatJet_nLeptons', sYLabel: 'Events', sXRange: [-0.5, 6.5] }),

    #("hLeadingFatJetParticleNet_massA_Hto4b_avg_v013", {sXLabel: 'hLeadingFatJetParticleNet_massA_Hto4b_avg_v013', sYLabel: 'Events', sXRange: [0, 70], sNRebinX: 5}),    
    #("hLeadingFatJetParticleNet_massH_Hto4b_avg_v0123", {sXLabel: 'hLeadingFatJetParticleNet_massH_Hto4b_avg_v0123', sYLabel: 'Events', sXRange: [50, 300], sNRebinX: 5}),

    #("hLeadingFatJetPNet_X4b_v1_Haa4b_vs_QCD", {sXLabel: 'hLeadingFatJetPNet_X4b_v1_Haa4b_vs_QCD', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 50 }),
    ("hLeadingFatJetPNet_X4b_v1_Haa4b_score", {sXLabel: 'hLeadingFatJetPNet_X4b_v1_Haa4b_score', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 50 }),
    ("hLeadingFatJetPNet_X4b_v2a_Haa4b_score", {sXLabel: 'hLeadingFatJetPNet_X4b_v2a_Haa4b_score', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 50 }),
    ("hLeadingFatJetPNet_X4b_v2b_Haa4b_score", {sXLabel: 'hLeadingFatJetPNet_X4b_v2b_Haa4b_score', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 50 }),
    ("hLeadingFatJetPNet_X4b_v2ab_Haa4b_score", {sXLabel: r'$X\to 4b$ tagger score', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 50 }),
    ("hLeadingFatJetPNet_X4b_v2a_Haa34b_score", {sXLabel: 'hLeadingFatJetPNet_X4b_v2a_Haa34b_score', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 50 }),
    ("hLeadingFatJetPNet_X4b_v2b_Haa34b_score", {sXLabel: 'hLeadingFatJetPNet_X4b_v2b_Haa34b_score', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 50 }),
    ("hLeadingFatJetPNet_X4b_v2ab_Haa34b_score", {sXLabel: r'$X\to 3,4b$ tagger score', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 50 }),
    
    ("hLeadingFatJetMassH_v2b", {sXLabel: r'Mass$_{PNet\, X\to 4b}$(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events', sXRange: [0, 250], sNRebinX: 5}),
    
    ("hLeadingFatJetPNet_massAa", {sXLabel: 'hLeadingFatJetPNet_massAa', sYLabel: 'Events', sXRange: [5, 70], sNRebinX: 20}),
    ("hLeadingFatJetPNet_34massAa", {sXLabel: r'Mass$_{version\, a}$(a)  [GeV]', sYLabel: 'Events', sXRange: [5, 70], sNRebinX: 20}),
    ("hLeadingFatJetPNet_34massAb", {sXLabel: 'hLeadingFatJetPNet_34massAb', sYLabel: 'Events', sXRange: [5, 70], sNRebinX: 20}),
    ("hLeadingFatJetPNet_34massAd", {sXLabel: r'Mass$_{version\, d}$(a)  [GeV]', sYLabel: 'Events', sXRange: [5, 70], sNRebinX: 20}),
    #("hLeadingFatJetPNet_34massAad", {sXLabel: 'hLeadingFatJetPNet_34massAad', sYLabel: 'Events', sXRange: [5, 70], sNRebinX: 20}),
    #("hLeadingFatJetPNet_massA1", {sXLabel: 'hLeadingFatJetPNet_massA1', sYLabel: 'Events', sXRange: [0, 80], sNRebinX: 20}),
    #("hLeadingFatJetPNet_massA2", {sXLabel: 'hLeadingFatJetPNet_massA2', sYLabel: 'Events', sXRange: [0, 80], sNRebinX: 20}),
    #("hLeadingFatJetPNet_massAA", {sXLabel: 'hLeadingFatJetPNet_massAA', sYLabel: 'Events', sXRange: [0, 80], sNRebinX: 20}),
    #("hLeadingFatJetPNet_dMassAA", {sXLabel: 'hLeadingFatJetPNet_dMassAA', sYLabel: 'Events', sXRange: [0, 80], sNRebinX: 20}),
    #("hLeadingFatJetPNet_dMassAA_relH", {sXLabel: 'hLeadingFatJetPNet_dMassAA_relH', sYLabel: 'Events', sXRange: [0, 1], sNRebinX: 1}),
    
    #("hLeadingFatJetMass_H34bCat", {sXLabel: 'm(Leading AK8 jet (mass) + nearest AK4 jet) [GeV]', sYLabel: 'Events', sXRange: [0, 250], sNRebinX: 5}),
    #("hLeadingFatJetMSoftDrop_H34bCat", {sXLabel: 'm(Leading AK8 jet (m-soft-drop) + nearest AK4 jet) [GeV]', sYLabel: 'Events', sXRange: [0, 250], sNRebinX: 5}),
    #("hLeadingFatJetMassH_H34bCat", {sXLabel: 'm(Leading AK8 jet (m-PNet) + nearest AK4 jet) [GeV]', sYLabel: 'Events', sXRange: [0, 250], sNRebinX: 5}),
    
    
    

    #("", {sXLabel: '', sYLabel: 'Events'}),

])

'''
histograms_dict = OD([
    ("hLeadingFatJetPNet_34massAa", {sXLabel: r'Mass$_{version\, a}$(a)  [GeV]', sYLabel: 'Events', sXRange: [5, 70], sNRebinX: 20}),
    

])
'''