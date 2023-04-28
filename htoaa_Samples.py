#import OrderedDict as OD
from collections import OrderedDict as OD

kData         = "Data" # dict key for Datasets
kQCD_bEnrich  = "QCD_bEnrich"
kQCD_bGen     = "QCD_bGen"
kQCDIncl      = "QCDIncl"
kZJets        = "ZJets"
kWJets        = "WJets"

QCDInclMode = 2 # 1: run "QCDIncl", 2: run "QCDIncl_PSWeight", 0: run both "QCDIncl" and "QCDIncl_PSWeight". Use 2 as "QCDIncl_PSWeight"

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

    ("QCDIncl", [
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

    ("QCDIncl_PSWeight", [
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

    ("TTJets", [
        "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8"
    ]),

    (kZJets, [
        "ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),
    
    (kWJets, [
        "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
        "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]),
    

    ("SUSY_GluGluH_01J_HToAATo4B", [
        #"SUSY_GluGluH_01J_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8",
        #"SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"

        # Signal sample with Higgs Pt > 150 GeV 
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8",
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8",

        # Inclusive samples
        #"SUSY_GluGluH_01J_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8",
        #"SUSY_GluGluH_01J_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8",
        #"SUSY_GluGluH_01J_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8",
        #"SUSY_GluGluH_01J_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8",
        #"SUSY_GluGluH_01J_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8",
        #"SUSY_GluGluH_01J_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8",
        #"SUSY_GluGluH_01J_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8",
        #"SUSY_GluGluH_01J_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8",
        #"SUSY_GluGluH_01J_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8",
        #"SUSY_GluGluH_01J_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8",
        #"SUSY_GluGluH_01J_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8",

    ])
])

# run either "QCDIncl" or "QCDIncl_PSWeight" or both
if   QCDInclMode == 1: # "QCDIncl"
    list_tmp_ = Samples2018["QCDIncl"]
    Samples2018.pop("QCDIncl",          None)
    Samples2018.pop("QCDIncl_PSWeight", None)
    Samples2018[kQCDIncl] = list_tmp_
    #del list_tmp_
    
elif QCDInclMode == 2: # "QCDIncl_PSWeight"
    list_tmp_ = Samples2018["QCDIncl_PSWeight"]
    Samples2018.pop("QCDIncl",          None)
    Samples2018.pop("QCDIncl_PSWeight", None)
    Samples2018[kQCDIncl] = list_tmp_
    #del list_tmp_

