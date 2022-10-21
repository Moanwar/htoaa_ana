#import OrderedDict as OD
from collections import OrderedDict as OD

kData = "Data" # dict key for Datasets

Samples2018 = OD([
    ("QCD", [
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

    ("SUSY_GluGluH_01J_HToAATo4B", [
        "SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8"
    ])
])
