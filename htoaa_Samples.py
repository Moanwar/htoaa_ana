#import OrderedDict as OD
from collections import OrderedDict as OD

kData         = "Data" # dict key for Datasets
kQCD_bEnrich  = "QCD_bEnr" # "QCD_bEnrich"
kQCD_bGen     = "QCD_BGen" # "QCD_bGen"
kQCDIncl      = "QCD_Incl"
kZJets        = "ZJets"
kWJets        = "WJets"

QCDInclMode = 2 # 1: run "QCD_Incl", 2: run "QCD_Incl_PSWeight", 0: run both "QCD_Incl" and "QCD_Incl_PSWeight". Use 2 as "QCD_Incl_PSWeight" (recommended)

Samples2018 = OD([

    #(kData, [
    #    "JetHT_Run2018A",
    #    "JetHT_Run2018B",
    #    "JetHT_Run2018C",
    #    "JetHT_Run2018D"
    #]),

    ## JetHT data
    ("JetHT_Run2018A", ["JetHT_Run2018A"]),
    ("JetHT_Run2018B", ["JetHT_Run2018B"]),
    ("JetHT_Run2018C", ["JetHT_Run2018C"]),
    ("JetHT_Run2018D", ["JetHT_Run2018D"]),

    ## MET data
    ("MET_Run2018A", ["MET_Run2018A"]),
    ("MET_Run2018B", ["MET_Run2018B"]),
    ("MET_Run2018C", ["MET_Run2018C"]),
    ("MET_Run2018D", ["MET_Run2018D"]),

    ## SingleMuon data
    ("SingleMuon_Run2018A", ["SingleMuon_Run2018A"]),
    ("SingleMuon_Run2018B", ["SingleMuon_Run2018B"]),
    ("SingleMuon_Run2018C", ["SingleMuon_Run2018C"]),
    ("SingleMuon_Run2018D", ["SingleMuon_Run2018D"]),

    ## EGamma data
    ("EGamma_Run2018A", ["EGamma_Run2018A"]),
    ("EGamma_Run2018B", ["EGamma_Run2018B"]),
    ("EGamma_Run2018C", ["EGamma_Run2018C"]),
    ("EGamma_Run2018D", ["EGamma_Run2018D"]),
    

    
    (kQCD_bEnrich,[
        "QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT200to300_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT300to500_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8"
    ]),

    (kQCD_bGen, [
        "QCD_HT100to200_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT200to300_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT300to500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT700to1000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT1000to1500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT1500to2000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT2000toInf_BGenFilter_TuneCP5_13TeV-madgraph-pythia8"
    ]),

    ("QCD_Incl", [
        "QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8",
        "QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),

    ("QCD_Incl_PSWeight", [
        "QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
        "QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
        "QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
        "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
        "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
        "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
        "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
        "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
        "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
     ]),


    ("TT0l", [    #("TTToHadronic_powheg", [
        "TTToHadronic_TuneCP5_13TeV-powheg-pythia8",
    ]),
    ("TT1l", [   #("TTToSemiLeptonic_powheg", [
        "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8",
    ]),
    ("TT2l", [  #("TTTo2L2Nu_powheg", [
        "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
    ]),


#    ("TTJets_Incl_NLO", [
#        "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8"
#    ]),

#    ("TTJets_Incl_LO", [
#        "TTJets_TuneCP5_13TeV-madgraphMLM-pythia8",
#    ]),

#    ("TTJets_HT_LO", [
#        "TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
#        "TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8",
#        "TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8",
#        "TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
#    ]),

#    ("TTJets_Lep_LO", [
#        "TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8",
#        "TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8",
#        "TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8",
#    ]),    

    ("STop_t", [
        "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8",
    ]),
    ("STbar_t", [
        "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8",
    ]),    
    ("ST_s_0l", [
        "ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8",
    ]),    
    ("ST_s_1l", [
        "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8",
    ]),    
    ("STop_tW_Incl", [
        "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",
    ]),    
    ("STbar_tW_Incl", [
        "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",
    ]),    
    ("STop_tW_12l", [
        "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
    ]),
    ("STbar_tW_12l", [
        "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
    ]),

    ("ttZ", [
        "ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8",
    ]),
    ("ttW", [
        "ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8",
    ]),
    ("tZq", [
        "tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8",
    ]),
    
    

    ("Zqq", [  #("ZJetsToQQ_HT", [
        "ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),

    ("Zvv", [    #("ZJetsToNuNu_HT", [
        "ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",        
    ]),

    ("Zll", [
        "DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8",
        "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    ]),

    #("DYJets_Incl_NLO", [
    #    "DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    #    "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    #]),
    #("DYJets_M-10to50_Incl_NLO", [
    #    "DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    #]),
    #("DYJets_M-50_Incl_NLO", [
    #    "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    #]),

    #("DYJets_Incl_LO", [
    #    "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8",
    #    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
    #]),
    #("DYJets_M-10to50_Incl_LO", [
    #    "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8",
    #]),
    #("DYJets_M-50_Incl_LO", [
    #    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
    #]),

    #("DYJets_M-50_HT_LO", [
    #    "DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
    #    "DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
    #    "DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
    #    "DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
    #    "DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
    #    "DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
    #    "DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
    #    "DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
    #]),        
    
    ("Wqq", [  #("WJetsToQQ_HT", [
        "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),

    #('Wlv', [   #('WJetsToLNu_Incl_NLO', [
    #    "WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    #]),

    #('WJetsToLNu_Incl_LO', [
    #    "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
    #]),

    ('Wlv', [  #('WJetsToLNu_HT_LO', [
        "WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",
    ]),
 
    #('W1JetsToLNu_LO', [
    #    "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
    #]),   
 
    #('W2JetsToLNu_LO', [
    #    "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
    #]),   
 
    #('W3JetsToLNu_LO', [
    #    "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
    #]),   
 
    #('W4JetsToLNu_LO', [
    #    "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
    #]),


    ('ZZ', [
        "ZZ_TuneCP5_13TeV-pythia8",
    ]),
    ('WZ', [
        "WZ_TuneCP5_13TeV-pythia8",
    ]),
    ('WW', [
        "WW_TuneCP5_13TeV-pythia8",
    ]),
    ('ZZZ', [
        "ZZZ_TuneCP5_13TeV-amcatnlo-pythia8",
    ]),
    ('WZZ', [
        "WZZ_TuneCP5_13TeV-amcatnlo-pythia8",
    ]),
    ('WWZ', [
        "WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8",
    ]),
    ('WWW', [
        "WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8",
    ]),        




    #('GluGluHToBB_Incl', [
    #    "GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8",
    #]),
    ('ggH', [   #('GluGluHToBB_Pt-200ToInf', [
        "GluGluHToBB_Pt-200ToInf_M-125_TuneCP5_MINLO_13TeV-powheg-pythia8",
    ]),
    ('VBFH', [   #('VBFHToBB_powheg', [
        "VBFHToBB_M-125_TuneCP5_13TeV-powheg-pythia8",
    ]),
    #('VBFHToBB_herwig', [
    #    "VBFHToBB_M-125_TuneCH3_13TeV-powheg-herwig",
    #]),
    ('VBFH_dipoleRecoilOn', [
        "VBFHToBB_M-125_dipoleRecoilOn_TuneCP5_13TeV-powheg-pythia8",
    ]),
    ('VBFWH_dipoleRecoilOn', [
        "VBFWH_HToBB_WToLNu_M-125_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8",
    ]),

    ('WHbbqq', [
        "WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8",
        "WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8"
    ]),
    ('WHbblv', [
        "WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8",
        "WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8"
    ]),
    ('ZH', [
        "ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8",
        "ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8",
        "ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8",
        "ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8"
    ]),
    ('ggZH', [
        "ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8",
        "ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8",
        "ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8",
        "ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8",
    ]),
    ('ttH', [
        "ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8",
        "ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8"
    ]),


    ("ggHtoaato4b_mA_12", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ggHtoaato4b_mA_15", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ggHtoaato4b_mA_20", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ggHtoaato4b_mA_25", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ggHtoaato4b_mA_30", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ggHtoaato4b_mA_35", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ggHtoaato4b_mA_40", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ggHtoaato4b_mA_45", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ggHtoaato4b_mA_50", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ggHtoaato4b_mA_55", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ggHtoaato4b_mA_60", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8"]),

    ("VBFHtoaato4b_mA_12", ["SUSY_VBFH_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8"]),
    ("VBFHtoaato4b_mA_15", ["SUSY_VBFH_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8"]),
    ("VBFHtoaato4b_mA_20", ["SUSY_VBFH_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"]),
    ("VBFHtoaato4b_mA_25", ["SUSY_VBFH_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8"]),
    ("VBFHtoaato4b_mA_30", ["SUSY_VBFH_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8"]),
    ("VBFHtoaato4b_mA_35", ["SUSY_VBFH_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8"]),
    ("VBFHtoaato4b_mA_40", ["SUSY_VBFH_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8"]),
    ("VBFHtoaato4b_mA_45", ["SUSY_VBFH_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8"]),
    ("VBFHtoaato4b_mA_50", ["SUSY_VBFH_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8"]),
    ("VBFHtoaato4b_mA_55", ["SUSY_VBFH_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8"]),
    ("VBFHtoaato4b_mA_60", ["SUSY_VBFH_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8"]),

    ("WHtoaato4b_mA_12", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8"]),
    ("WHtoaato4b_mA_15", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8"]),
    ("WHtoaato4b_mA_20", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"]),
    ("WHtoaato4b_mA_25", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8"]),
    ("WHtoaato4b_mA_30", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8"]),
    ("WHtoaato4b_mA_35", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8"]),
    ("WHtoaato4b_mA_40", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8"]),
    ("WHtoaato4b_mA_45", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8"]),
    ("WHtoaato4b_mA_50", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8"]),
    ("WHtoaato4b_mA_55", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8"]),
    ("WHtoaato4b_mA_60", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8"]),

    ("ZHtoaato4b_mA_12", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ZHtoaato4b_mA_15", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ZHtoaato4b_mA_20", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ZHtoaato4b_mA_25", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ZHtoaato4b_mA_30", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ZHtoaato4b_mA_35", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ZHtoaato4b_mA_40", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ZHtoaato4b_mA_45", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ZHtoaato4b_mA_50", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ZHtoaato4b_mA_55", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ZHtoaato4b_mA_60", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8"]),

    ("ttHtoaato4b_mA_12", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ttHtoaato4b_mA_15", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ttHtoaato4b_mA_20", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ttHtoaato4b_mA_25", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ttHtoaato4b_mA_30", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ttHtoaato4b_mA_35", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ttHtoaato4b_mA_40", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ttHtoaato4b_mA_45", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ttHtoaato4b_mA_50", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ttHtoaato4b_mA_55", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8"]),
    ("ttHtoaato4b_mA_60", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8"]),


    #("ggHtoaato4b_Incl_mA_12", ["SUSY_GluGluH_01J_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ggHtoaato4b_Incl_mA_15", ["SUSY_GluGluH_01J_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ggHtoaato4b_Incl_mA_20", ["SUSY_GluGluH_01J_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ggHtoaato4b_Incl_mA_25", ["SUSY_GluGluH_01J_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ggHtoaato4b_Incl_mA_30", ["SUSY_GluGluH_01J_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ggHtoaato4b_Incl_mA_35", ["SUSY_GluGluH_01J_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ggHtoaato4b_Incl_mA_40", ["SUSY_GluGluH_01J_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ggHtoaato4b_Incl_mA_45", ["SUSY_GluGluH_01J_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ggHtoaato4b_Incl_mA_50", ["SUSY_GluGluH_01J_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ggHtoaato4b_Incl_mA_55", ["SUSY_GluGluH_01J_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ggHtoaato4b_Incl_mA_60", ["SUSY_GluGluH_01J_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"]),

    #("VBFHtoaato4b_Incl_mA_12", ["SUSY_VBFH_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8"]),
    #("VBFHtoaato4b_Incl_mA_15", ["SUSY_VBFH_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8"]),
    #("VBFHtoaato4b_Incl_mA_20", ["SUSY_VBFH_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8"]),
    #("VBFHtoaato4b_Incl_mA_25", ["SUSY_VBFH_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8"]),
    #("VBFHtoaato4b_Incl_mA_30", ["SUSY_VBFH_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8"]),
    #("VBFHtoaato4b_Incl_mA_35", ["SUSY_VBFH_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8"]),
    #("VBFHtoaato4b_Incl_mA_40", ["SUSY_VBFH_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8"]),
    #("VBFHtoaato4b_Incl_mA_45", ["SUSY_VBFH_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8"]),
    #("VBFHtoaato4b_Incl_mA_50", ["SUSY_VBFH_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8"]),
    #("VBFHtoaato4b_Incl_mA_55", ["SUSY_VBFH_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8"]),
    #("VBFHtoaato4b_Incl_mA_60", ["SUSY_VBFH_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"]),

    #("WHtoaato4b_Incl_mA_12", ["SUSY_WH_WToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8"]),
    #("WHtoaato4b_Incl_mA_15", ["SUSY_WH_WToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8"]),
    #("WHtoaato4b_Incl_mA_20", ["SUSY_WH_WToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8"]),
    #("WHtoaato4b_Incl_mA_25", ["SUSY_WH_WToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8"]),
    #("WHtoaato4b_Incl_mA_30", ["SUSY_WH_WToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8"]),
    #("WHtoaato4b_Incl_mA_35", ["SUSY_WH_WToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8"]),
    #("WHtoaato4b_Incl_mA_40", ["SUSY_WH_WToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8"]),
    #("WHtoaato4b_Incl_mA_45", ["SUSY_WH_WToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8"]),
    #("WHtoaato4b_Incl_mA_50", ["SUSY_WH_WToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8"]),
    #("WHtoaato4b_Incl_mA_55", ["SUSY_WH_WToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8"]),
    #("WHtoaato4b_Incl_mA_60", ["SUSY_WH_WToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"]),

    #("ZHtoaato4b_Incl_mA_12", ["SUSY_ZH_ZToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ZHtoaato4b_Incl_mA_15", ["SUSY_ZH_ZToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ZHtoaato4b_Incl_mA_20", ["SUSY_ZH_ZToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ZHtoaato4b_Incl_mA_25", ["SUSY_ZH_ZToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ZHtoaato4b_Incl_mA_30", ["SUSY_ZH_ZToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ZHtoaato4b_Incl_mA_35", ["SUSY_ZH_ZToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ZHtoaato4b_Incl_mA_40", ["SUSY_ZH_ZToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ZHtoaato4b_Incl_mA_45", ["SUSY_ZH_ZToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ZHtoaato4b_Incl_mA_50", ["SUSY_ZH_ZToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ZHtoaato4b_Incl_mA_55", ["SUSY_ZH_ZToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ZHtoaato4b_Incl_mA_60", ["SUSY_ZH_ZToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"]),

    #("ttHtoaato4b_Incl_mA_12", ["SUSY_TTH_TTToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ttHtoaato4b_Incl_mA_15", ["SUSY_TTH_TTToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ttHtoaato4b_Incl_mA_20", ["SUSY_TTH_TTToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ttHtoaato4b_Incl_mA_25", ["SUSY_TTH_TTToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ttHtoaato4b_Incl_mA_30", ["SUSY_TTH_TTToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ttHtoaato4b_Incl_mA_35", ["SUSY_TTH_TTToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ttHtoaato4b_Incl_mA_40", ["SUSY_TTH_TTToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ttHtoaato4b_Incl_mA_45", ["SUSY_TTH_TTToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ttHtoaato4b_Incl_mA_50", ["SUSY_TTH_TTToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ttHtoaato4b_Incl_mA_55", ["SUSY_TTH_TTToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8"]),
    #("ttHtoaato4b_Incl_mA_60", ["SUSY_TTH_TTToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"]),

    #("ggHtoaato4b_Incl_mA_All", [
    #    "SUSY_GluGluH_01J_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_GluGluH_01J_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_GluGluH_01J_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_GluGluH_01J_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_GluGluH_01J_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_GluGluH_01J_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_GluGluH_01J_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_GluGluH_01J_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_GluGluH_01J_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_GluGluH_01J_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_GluGluH_01J_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"
    #]),

    #('VBFHtoaato4b_Incl_mA_All', [
    #    "SUSY_VBFH_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_VBFH_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_VBFH_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_VBFH_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_VBFH_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_VBFH_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_VBFH_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_VBFH_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_VBFH_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_VBFH_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_VBFH_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"
    #]),

    #('WHtoaato4b_Incl_mA_All', [
    #    "SUSY_WH_WToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_WH_WToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_WH_WToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_WH_WToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_WH_WToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_WH_WToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_WH_WToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_WH_WToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_WH_WToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_WH_WToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_WH_WToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"
    #]),

    #('ZHtoaato4b_Incl_mA_All', [
    #    "SUSY_ZH_ZToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_ZH_ZToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_ZH_ZToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_ZH_ZToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_ZH_ZToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_ZH_ZToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_ZH_ZToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_ZH_ZToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_ZH_ZToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_ZH_ZToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_ZH_ZToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"
    #]),

    #('ttHtoaato4b_Incl_mA_All', [
    #    "SUSY_TTH_TTToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_TTH_TTToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_TTH_TTToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_TTH_TTToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_TTH_TTToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_TTH_TTToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_TTH_TTToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_TTH_TTToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_TTH_TTToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_TTH_TTToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8",
    #    "SUSY_TTH_TTToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"
    #]),
    
    
])

# run either "QCDIncl" or "QCDIncl_PSWeight" or both
if   QCDInclMode == 1: # "QCDIncl"
    list_tmp_ = Samples2018["QCD_Incl"]
    Samples2018.pop("QCD_Incl",          None)
    Samples2018.pop("QCD_Incl_PSWeight", None)
    Samples2018[kQCDIncl] = list_tmp_
    #del list_tmp_
    
elif QCDInclMode == 2: # "QCDIncl_PSWeight"
    list_tmp_ = Samples2018["QCD_Incl_PSWeight"]
    Samples2018.pop("QCD_Incl",          None)
    Samples2018.pop("QCD_Incl_PSWeight", None)
    Samples2018[kQCDIncl] = list_tmp_
    #del list_tmp_

