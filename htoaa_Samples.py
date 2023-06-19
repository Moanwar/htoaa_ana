#import OrderedDict as OD
from collections import OrderedDict as OD

kData         = "Data" # dict key for Datasets
kQCD_bEnrich  = "QCD_bEnrich"
kQCD_bGen     = "QCD_bGen"
kQCDIncl      = "QCD_Incl"
kZJets        = "ZJets"
kWJets        = "WJets"

QCDInclMode = 0 # 1: run "QCD_Incl", 2: run "QCD_Incl_PSWeight", 0: run both "QCD_Incl" and "QCD_Incl_PSWeight". Use 2 as "QCD_Incl_PSWeight"

Samples2018 = OD([

    (kData, [
        "JetHT_Run2018A",
        "JetHT_Run2018B",
        "JetHT_Run2018C",
        "JetHT_Run2018D"
    ]),
    
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

    ("TTJets_NLO", [
        "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8"
    ]),

    ("TTJets_Incl", [
        "TTJets_TuneCP5_13TeV-madgraphMLM-pythia8",
    ]),

    ("TTJets_HT", [
        "TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8",
        "TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8",
        "TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),
    

    ("ZJetsToQQ_HT", [
        "ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),
    
    ("WJetsToQQ_HT", [
        "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),

    ('WJetsToLNu_Incl', [
        "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
    ]),

    ('WJetsToLNu_HT', [
        #"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",
        #"W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        #"W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        #"W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        #"W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),
 
    ('W1JetsToLNu', [
        "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
    ]),   
 
    ('W2JetsToLNu', [
        "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
    ]),   
 
    ('W3JetsToLNu', [
        "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
    ]),   
 
    ('W4JetsToLNu', [
        "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
    ]),


    ("SUSY_GluGluH_01J_HToAATo4B_M-12_Incl",      ["SUSY_GluGluH_01J_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_GluGluH_01J_HToAATo4B_M-12_HPtAbv150", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8"]),
    
    ("SUSY_GluGluH_01J_HToAATo4B_M-15_Incl",      ["SUSY_GluGluH_01J_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_GluGluH_01J_HToAATo4B_M-15_HPtAbv150", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8"]),
    
    ("SUSY_GluGluH_01J_HToAATo4B_M-20_Incl",      ["SUSY_GluGluH_01J_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_GluGluH_01J_HToAATo4B_M-20_HPtAbv150", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"]),
    
    ("SUSY_GluGluH_01J_HToAATo4B_M-25_Incl",      ["SUSY_GluGluH_01J_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_GluGluH_01J_HToAATo4B_M-25_HPtAbv150", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8"]),
    
    ("SUSY_GluGluH_01J_HToAATo4B_M-30_Incl",      ["SUSY_GluGluH_01J_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_GluGluH_01J_HToAATo4B_M-30_HPtAbv150", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8"]),
    
    ("SUSY_GluGluH_01J_HToAATo4B_M-35_Incl",      ["SUSY_GluGluH_01J_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_GluGluH_01J_HToAATo4B_M-35_HPtAbv150", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8"]),
    
    ("SUSY_GluGluH_01J_HToAATo4B_M-40_Incl",      ["SUSY_GluGluH_01J_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_GluGluH_01J_HToAATo4B_M-40_HPtAbv150", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8"]),
    
    ("SUSY_GluGluH_01J_HToAATo4B_M-45_Incl",      ["SUSY_GluGluH_01J_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_GluGluH_01J_HToAATo4B_M-45_HPtAbv150", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8"]),
    
    ("SUSY_GluGluH_01J_HToAATo4B_M-50_Incl",      ["SUSY_GluGluH_01J_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_GluGluH_01J_HToAATo4B_M-50_HPtAbv150", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8"]),
    
    ("SUSY_GluGluH_01J_HToAATo4B_M-55_Incl",      ["SUSY_GluGluH_01J_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_GluGluH_01J_HToAATo4B_M-55_HPtAbv150", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8"]),
    
    ("SUSY_GluGluH_01J_HToAATo4B_M-60_Incl",      ["SUSY_GluGluH_01J_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_GluGluH_01J_HToAATo4B_M-60_HPtAbv150", ["SUSY_GluGluH_01J_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8"]),




    ("SUSY_VBFH_HToAATo4B_M-12_Incl",      ["SUSY_VBFH_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_VBFH_HToAATo4B_M-12_HPtAbv150", ["SUSY_VBFH_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_VBFH_HToAATo4B_M-15_Incl",      ["SUSY_VBFH_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_VBFH_HToAATo4B_M-15_HPtAbv150", ["SUSY_VBFH_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_VBFH_HToAATo4B_M-20_Incl",      ["SUSY_VBFH_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_VBFH_HToAATo4B_M-20_HPtAbv150", ["SUSY_VBFH_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_VBFH_HToAATo4B_M-25_Incl",      ["SUSY_VBFH_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_VBFH_HToAATo4B_M-25_HPtAbv150", ["SUSY_VBFH_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_VBFH_HToAATo4B_M-30_Incl",      ["SUSY_VBFH_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_VBFH_HToAATo4B_M-30_HPtAbv150", ["SUSY_VBFH_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_VBFH_HToAATo4B_M-35_Incl",      ["SUSY_VBFH_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_VBFH_HToAATo4B_M-35_HPtAbv150", ["SUSY_VBFH_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_VBFH_HToAATo4B_M-40_Incl",      ["SUSY_VBFH_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_VBFH_HToAATo4B_M-40_HPtAbv150", ["SUSY_VBFH_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_VBFH_HToAATo4B_M-45_Incl",      ["SUSY_VBFH_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_VBFH_HToAATo4B_M-45_HPtAbv150", ["SUSY_VBFH_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_VBFH_HToAATo4B_M-50_Incl",      ["SUSY_VBFH_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_VBFH_HToAATo4B_M-50_HPtAbv150", ["SUSY_VBFH_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_VBFH_HToAATo4B_M-55_Incl",      ["SUSY_VBFH_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_VBFH_HToAATo4B_M-55_HPtAbv150", ["SUSY_VBFH_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_VBFH_HToAATo4B_M-60_Incl",      ["SUSY_VBFH_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_VBFH_HToAATo4B_M-60_HPtAbv150", ["SUSY_VBFH_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8"]), 




    ("SUSY_WH_WToAll_HToAATo4B_M-12_Incl",      ["SUSY_WH_WToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_WH_WToAll_HToAATo4B_M-12_HPtAbv150", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_WH_WToAll_HToAATo4B_M-15_Incl",      ["SUSY_WH_WToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_WH_WToAll_HToAATo4B_M-15_HPtAbv150", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_WH_WToAll_HToAATo4B_M-20_Incl",      ["SUSY_WH_WToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_WH_WToAll_HToAATo4B_M-20_HPtAbv150", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_WH_WToAll_HToAATo4B_M-25_Incl",      ["SUSY_WH_WToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_WH_WToAll_HToAATo4B_M-25_HPtAbv150", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_WH_WToAll_HToAATo4B_M-30_Incl",      ["SUSY_WH_WToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_WH_WToAll_HToAATo4B_M-30_HPtAbv150", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_WH_WToAll_HToAATo4B_M-35_Incl",      ["SUSY_WH_WToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_WH_WToAll_HToAATo4B_M-35_HPtAbv150", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_WH_WToAll_HToAATo4B_M-40_Incl",      ["SUSY_WH_WToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_WH_WToAll_HToAATo4B_M-40_HPtAbv150", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_WH_WToAll_HToAATo4B_M-45_Incl",      ["SUSY_WH_WToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_WH_WToAll_HToAATo4B_M-45_HPtAbv150", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_WH_WToAll_HToAATo4B_M-50_Incl",      ["SUSY_WH_WToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_WH_WToAll_HToAATo4B_M-50_HPtAbv150", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_WH_WToAll_HToAATo4B_M-55_Incl",      ["SUSY_WH_WToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_WH_WToAll_HToAATo4B_M-55_HPtAbv150", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_WH_WToAll_HToAATo4B_M-60_Incl",      ["SUSY_WH_WToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_WH_WToAll_HToAATo4B_M-60_HPtAbv150", ["SUSY_WH_WToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8"]), 




    ("SUSY_ZH_ZToAll_HToAATo4B_M-12_Incl",      ["SUSY_ZH_ZToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-12_HPtAbv150", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_ZH_ZToAll_HToAATo4B_M-15_Incl",      ["SUSY_ZH_ZToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-15_HPtAbv150", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_ZH_ZToAll_HToAATo4B_M-20_Incl",      ["SUSY_ZH_ZToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-20_HPtAbv150", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_ZH_ZToAll_HToAATo4B_M-25_Incl",      ["SUSY_ZH_ZToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-25_HPtAbv150", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_ZH_ZToAll_HToAATo4B_M-30_Incl",      ["SUSY_ZH_ZToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-30_HPtAbv150", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_ZH_ZToAll_HToAATo4B_M-35_Incl",      ["SUSY_ZH_ZToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-35_HPtAbv150", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_ZH_ZToAll_HToAATo4B_M-40_Incl",      ["SUSY_ZH_ZToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-40_HPtAbv150", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_ZH_ZToAll_HToAATo4B_M-45_Incl",      ["SUSY_ZH_ZToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-45_HPtAbv150", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_ZH_ZToAll_HToAATo4B_M-50_Incl",      ["SUSY_ZH_ZToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-50_HPtAbv150", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_ZH_ZToAll_HToAATo4B_M-55_Incl",      ["SUSY_ZH_ZToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-55_HPtAbv150", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_ZH_ZToAll_HToAATo4B_M-60_Incl",      ["SUSY_ZH_ZToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-60_HPtAbv150", ["SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8"]), 


    
    ("SUSY_TTH_TTToAll_HToAATo4B_M-12_Incl",      ["SUSY_TTH_TTToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-12_HPtAbv150", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_TTH_TTToAll_HToAATo4B_M-15_Incl",      ["SUSY_TTH_TTToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-15_HPtAbv150", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_TTH_TTToAll_HToAATo4B_M-20_Incl",      ["SUSY_TTH_TTToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-20_HPtAbv150", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_TTH_TTToAll_HToAATo4B_M-25_Incl",      ["SUSY_TTH_TTToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-25_HPtAbv150", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_TTH_TTToAll_HToAATo4B_M-30_Incl",      ["SUSY_TTH_TTToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-30_HPtAbv150", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_TTH_TTToAll_HToAATo4B_M-35_Incl",      ["SUSY_TTH_TTToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-35_HPtAbv150", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_TTH_TTToAll_HToAATo4B_M-40_Incl",      ["SUSY_TTH_TTToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-40_HPtAbv150", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_TTH_TTToAll_HToAATo4B_M-45_Incl",      ["SUSY_TTH_TTToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-45_HPtAbv150", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_TTH_TTToAll_HToAATo4B_M-50_Incl",      ["SUSY_TTH_TTToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-50_HPtAbv150", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_TTH_TTToAll_HToAATo4B_M-55_Incl",      ["SUSY_TTH_TTToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-55_HPtAbv150", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8"]), 

    ("SUSY_TTH_TTToAll_HToAATo4B_M-60_Incl",      ["SUSY_TTH_TTToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8"      ]),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-60_HPtAbv150", ["SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8"]), 

    

    
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

