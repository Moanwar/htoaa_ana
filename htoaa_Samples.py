#import OrderedDict as OD
from collections import OrderedDict as OD

kData = "Data" # dict key for Datasets
kQCD_bEnrich = "QCD_bEnrich"
kQCD_bGen = "QCD_bGen"
kQCDIncl = "QCDIncl"

Samples2018 = OD([

#    (kData, [
#        "/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD",
#        "/JetHT/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD",
#        "/JetHT/Run2018C-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD",
#        "/JetHT/Run2018D-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD"
#    ]),    
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

    (kQCD_bGen,[
        "QCD_HT100to200_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT200to300_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT300to500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT700to1000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT1000to1500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT1500to2000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",
        "QCD_HT2000toInf_BGenFilter_TuneCP5_13TeV-madgraph-pythia8"
    ]),

    (kQCDIncl,[
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

#    ("QCDIncl_PSWeight", [
#        "QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#        "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
#     ]),

    ("TTJets", [
        "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8"
    ]),

    ("ZJets", [
        "ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),
    
    ("WJets", [
        "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),
    

    ("SUSY_GluGluH_01J_HToAATo4B", [
        #"SUSY_GluGluH_01J_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"
    ])
])
