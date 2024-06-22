'''
 Script to prepare Samples_2018UL.json by collecting information with 'DAQ queries'.
 To run:
        1] Set voms proxy:
            voms-proxy-init -voms cms -rfc

        2] python3 samplesList_prepare.py -era 2018
              Add -updateCrossSections to update cross-section values in the existing Samples_2018UL.json

        3] Script produces 'Samples_2018UL_v0.json' output.
           cp Samples_2018UL_v0.json Samples_2018UL.json   if you are satistied with Samples_2018UL_v0.json
'''

# List of cross-section used in this macro also with refereces can be found at 
# https://docs.google.com/spreadsheets/d/1LQDKBWGTdsT1uBumq9tz6RwRTCAVK2G0zmC6qaiqLkw/edit?usp=sharing


import os
import subprocess
import json
from collections import OrderedDict as OD
from copy import deepcopy
import argparse
import glob

from htoaa_Settings import *


printLevel = 0

sXS     = "xs"
sNameSp = "nameSp"
list_datasetAndXs_2018 = OD([
    
    ## QCD_bEnriched_HT*
    # dasgoclient --query="dataset=/QCD_bEnriched_HT*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"    
    ("/QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS: 1122000.0   }),
    ("/QCD_bEnriched_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:   79760.0   }),
    ("/QCD_bEnriched_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:   16600.0   }),
    ("/QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:    1503.0   }),    
    ("/QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:     297.4   }),
    ("/QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:      48.08  }),    
    ("/QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:       3.9510 }),
    ("/QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:       0.6957 }),
    

    ## QCD_HT*_BGenFilter
    # dasgoclient --query="dataset=/QCD_HT*_BGenFilter*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/QCD_HT100to200_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sXS: 1266000.0    }),
    ("/QCD_HT200to300_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sXS:  109900.0    }),
    ("/QCD_HT300to500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sXS:   27360.0    }),
    ("/QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sXS:    2991.0    }),
    ("/QCD_HT700to1000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",  {sXS:     731.8    }),    
    ("/QCD_HT1000to1500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:     139.3    }),
    ("/QCD_HT1500to2000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:      14.74   }),
    ("/QCD_HT2000toInf_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",  {sXS:       3.09   }),


    ## QCD_HT* PSWeights-madgraph - QCDIncl LO
    # dasgoclient --query="dataset=/QCD_HT*PSWeight*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {sXS: 187700000.0    }),
    ("/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:  23640000.0    }),
    ("/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:   1546000.0    }),
    ("/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:    321600.0    }),
    ("/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:     30310.0    }),
    ("/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",  {sXS:      6364.0    }),
    ("/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:      1117.0    }),
    ("/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:       108.4    }),
    ("/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:        22.36   }),

    
    ## QCD_HT* TuneCP5 madgraphMLM - QCDIncl LO, MatrixElement-PartonShower maching at NLO(?)
    # dasgoclient --query="dataset=/QCD_HT*TuneCP5_13TeV-madgraphMLM*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {sXS: 187700000.0   }),
    ("/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:  23500000.0   }),
    ("/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:   1552000.0   }),
    ("/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:    321100.0   }),
    ("/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:     30250.0   }),    
    ("/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:      6398.0   }),
    ("/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:      1122.0   }),
    ("/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:       109.4   }),
    ("/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:        21.74  }),

    

    ## TTbar - NLO powheg
    # dasgoclient --query="dataset=/TT*_TuneCP5_13TeV*powheg*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {sXS: 380.133 }), # 831.8 * 0.457
    ("/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {sXS: 364.328 }), # 831.8 * 0.438
    ("/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                      {sXS:  87.339 }), # 831.8 * 0.105

    ## TTbar Jets - NLO amcatnlo
    # dasgoclient --query="dataset=/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {sXS: 831.8       }),

    ## TTbar Jets - LO madgraph
    # dasgoclient --query="dataset=/TTJets_*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM" 
    ("/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                    {sXS: 831.8       }),
    ("/TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sXS:   2.4234    }),
    ("/TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sXS:   0.9818    }),
    ("/TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",      {sXS:   0.1714    }),
    ("/TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sXS:   0.001966  }),

    ## TTbar Jets To LNu - LO madgraph
    # dasgoclient --query="dataset=/TTJets_*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {sXS: 182.164 }), # 831.8 * 0.219
    ("/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 182.164 }), # 831.8 * 0.219
    ("/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",             {sXS:  87.339 }), # 831.8 * 0.105


    ## ST NLO
    # dasgoclient --query="dataset=/ST*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",     {sXS: 134.2   }),
    ("/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  80.0   }),
    ("/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",            {sXS:  39.65  }), # 79.3/2 = 39.65
    ("/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",        {sXS:  39.65  }), # 79.3/2 = 39.65
    ("/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",      {sXS:  21.63  }),  # 79.3/2*(1 - (0.6741*0.6741)) : BR(W->Hadrons)=0.6741. tW has W(from t) and W. BR(tW->NoFully hadronic)=(1 - (0.6741*0.6741))=0.546.
    ("/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:  21.63  }),  # 79.3/2*(1 - (0.6741*0.6741)) : Andrew: 21.61 
    ("/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sXS:   5.041 }),
    ("/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sXS:   4.831 }),

    ## tX
    # dasgoclient --query="dataset=/t*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:  0.0942  }),
    ("/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",          {sXS:  0.839  }), # https://cds.cern.ch/record/2227475/files/CERN-2017-002-M.pdf?version=1#page=180
    ("/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",          {sXS:  0.6008  }), # https://cds.cern.ch/record/2227475/files/CERN-2017-002-M.pdf?version=1#page=180
    

    ## ZJets
    # dasgoclient --query="dataset=/ZJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"    
    ("/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS: 1012.0   }),
    ("/ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:  114.2   }),
    ("/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:   25.34  }),
    ("/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:   12.99  }),

    # dasgoclient --query="dataset=/ZJetsToNuNu*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS: 302.96     }),  
    ("/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:  82.9239   }),  
    ("/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:  11.2574   }),  
    ("/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:   2.73349  }),  
    ("/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:   1.22321  }),  
    ("/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:   0.28526  }),  
    ("/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:   0.006366 }),  

    ## DYJetsToLL NLO
    # dasgoclient --query="dataset=/DY*Jets*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                    {sXS: 18610   }),
    ("/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                        {sXS:  6077.22   }),
    
    # DYJetsToLL LO Incl
    ("/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                     {sXS: 18610   }),
    ("/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                         {sXS:  6077.22  }),
    #("/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v1/NANOAODSIM",                    {sXS:  6077.22  }),

    # DYJetsToLL LO HT
    ("/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {sXS:   158.7445   }),
    ("/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:   159.1984   }),
    ("/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:    43.5384   }),
    ("/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:     5.9141   }),
    ("/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:     1.4377   }),
    ("/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:     0.6443   }),
    ("/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:     0.1511   }),
    ("/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:     0.003390 }),    

    
    ## WJetsToQQ
    # dasgoclient --query="dataset=/WJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS: 2549.0    }),
    ("/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:  276.5    }),
    ("/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:   59.25   }),
    ("/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:   28.75   }),


    ## WJetsToLNu
    # dasgoclient --query="dataset=/W*Jets*ToLNu*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                   {sXS: 61526.7     }), # NLO sample
    #("/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {sXS: }), # NLO sample
    #("/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {sXS: }), # NLO sample
    #("/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {sXS: }), # NLO sample    
    ("/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                    {sXS: 61526.7     }),
    ("/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",         {sXS:  1440.0     }),
    ("/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sXS:  1431.0     }),
    ("/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sXS:   382.1     }),
    ("/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sXS:    51.54    }),
    #("/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",   {sXS:    51.54    }),
    ("/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sXS:    12.49    }),
    #("/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",   {sXS:    12.49    }),
    ("/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sXS:     5.619   }),
    #("/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",  {sXS:     5.619   }),
    ("/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",      {sXS:     1.321   }),
    #("/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM", {sXS:     1.321   }),
    ("/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",       {sXS:     0.02992 }),
    ("/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {sXS: 10167.9     }),
    ("/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {sXS: 3199.45     }),
    ("/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {sXS:  941.16     }),
    ("/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                   {sXS:  439.08     }),


    ## VV, VVV
    # dasgoclient --query="dataset=/Z*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    # dasgoclient --query="dataset=/W*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                              {sXS:   16.523  }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    ("/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:   3.25 * 1.21   }), # https://indico.cern.ch/event/439995/contributions/1094416/attachments/1143460/1638648/diboson_final.pdf
    ("/ZZTo2Q2Nu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sXS:   4.07 * 1.22   }), # https://indico.cern.ch/event/439995/contributions/1094416/attachments/1143460/1638648/diboson_final.pdf    
    ("/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                    {sXS:   47.13  }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    ("/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                    {sXS:  118.7  }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV

    ("/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sXS:    0.01398  }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",     {sXS:    0.01398  }),
    ("/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sXS:    0.05565  }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",     {sXS:    0.05565  }),
    ("/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sXS:    0.1651  }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",  {sXS:    0.1651  }),
    ("/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sXS:    0.2086  }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",  {sXS:    0.2086  }),


    #("",          {sXS:    }),
    
    ## GGF HToBB
    # dasgoclient --query="dataset=/*HToBB*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    # dasgoclient --query="dataset=/*HTo*bb*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",      {sXS:  48.61 * 0.582    }), # GGF H * BR = 28.291 N3LO
    ("/GluGluHToBB_Pt-200ToInf_M-125_TuneCP5_MINLO_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:   0.2740   }),
    
    ("/VBFHToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                      {sXS:  3.766 * 0.582  }), # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
    ("/VBFHToBB_M-125_TuneCH3_13TeV-powheg-herwig/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                       {sXS:  3.766 * 0.582  }), # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
    ("/VBFHToBB_M-125_dipoleRecoilOn_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",       {sXS:  2.250  }), # HIG-21-020
    ("/VBFWH_HToBB_WToLNu_M-125_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sXS:  0.4378  }), # https://xsdb-temp.app.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=DAS%3DVBFWH_HToBB_WToLNu_M-125_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8

    ("/WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",            {sXS:  0.831 * 0.582 * 0.647    }),
    ("/WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",           {sXS:  0.831 * 0.582 * 3*0.1086 }),
    ("/WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",           {sXS:  0.527 * 0.582 * 0.647    }),
    ("/WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sXS:  0.527 * 0.582 * 3*0.1086 }),
    
    ("/ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                {sXS:  0.758 * 0.582 * 0.699     }),
    ("/ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                {sXS:  0.758 * 0.582 * 3*0.0337  }),
    ("/ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",              {sXS:  0.758 * 0.582 * 0.20      }),
    ("/ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                {sXS:  0.758 * 0.582 * 0.1512    }),

    ("/ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",              {sXS:  0.123 * 0.582 * 0.699     }),
    ("/ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",              {sXS:  0.123 * 0.582 * 3*0.0337  }),
    ("/ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",            {sXS:  0.123 * 0.582 * 0.20      }),
    ("/ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",              {sXS:  0.123 * 0.582 * 0.1512    }),

    ("/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                        {sXS:  0.5071 * 0.5824  }),
    ("/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                     {sXS:  0.5071 * (1 - 0.5824)  }),

    ("/bbHToBB_M-125_4FS_yt2_TuneCP5-13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",             {sXS:  -1  }),
    

    
    ## SUSY_GluGluH_01J_HToAATo4B_M-*   and   SUSY_GluGluH_01J_HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*GluGluH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/SUSY_GluGluH_01J_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61}),

    # SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (308) / (5379) = 5.726e-02 +- 3.168e-03	Matching efficiency = 0.5 +/- 0.0     Cross-section = 48.61 pb * 0.057	
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 }), # filter efficiency 0.05
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 }),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 }),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 }),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 }),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 }),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 }),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 }),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 }),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 }),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 }),


    ## VBF HToAATo4B_M-* and VBF HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*VBF*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/SUSY_VBFH_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 }),
    ("/SUSY_VBFH_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 }),
    ("/SUSY_VBFH_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 }),
    ("/SUSY_VBFH_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 }),
    ("/SUSY_VBFH_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 }),
    ("/SUSY_VBFH_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 }),
    ("/SUSY_VBFH_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 }),
    ("/SUSY_VBFH_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 }),
    ("/SUSY_VBFH_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 }),
    ("/SUSY_VBFH_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 }),
    ("/SUSY_VBFH_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 }),

    # SUSY_VBFH_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8		Filter efficiency (event-level)= (1730) / (10000) = 1.730e-01 +- 3.782e-03	Matching efficiency = 1.0 +/- 0.0 	Cross-section = 3.766 pb * 0.173 
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 }),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 }),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 }),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 }),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 }),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 }),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 }),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 }),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 }),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS: 3.766 * 0.173 }),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 }),


    ## WH HToAATo4B-M-* and WH HToAATo4B_Pt150_-M-*
    # dasgoclient --query="dataset=/SUSY*WH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM" 
    ("/SUSY_WH_WToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358}),

    # SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (1317) / (10000) = 1.317e-01 +- 3.382e-03	Matching efficiency = 1.0 +/- 0.0	Cross-section = 1.358 pb * 0.132 
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 }),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 }),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 }),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 }),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 }),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 }),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 }),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 }),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 }),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 }),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 }),


    ## ZH HToAATo4B-M-* and ZH HToAATo4B_Pt150_-M-*
    # dasgoclient --query="dataset=/SUSY*ZH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 }),
    
    # SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (1289) / (10000) = 1.289e-01 +- 3.351e-03	Matching efficiency = 1.0 +/- 0.0		Cross-section = 0.880 pb * 0.129 
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129  }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129  }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129  }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129  }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129  }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129  }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129  }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129  }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129  }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129  }),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129  }),


    ## ttH HToAATo4B_M-* and HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*TTH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071  }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071  }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071  }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071  }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071  }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071  }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071  }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071  }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071  }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071  }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071  }),

    # SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (2850) / (10000) = 2.850e-01 +- 4.514e-03	Matching efficiency = 1.0 +/- 0.0 Cross-section = 0.5071 pb * 0.285 
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285   }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285   }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285   }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285   }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285   }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285   }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285   }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285   }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285   }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285   }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285   }),
    
    
    


    

    ## JetHT
    # dasgoclient --query="dataset=/JetHT/*2018*UL*MiniAODv2_NanoAODv9-*/NANOAOD"
    # XS (cross-section) does not matter for data sample
    # Use JetHT/NanoAOD_GT38 instead
    #("/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD", {sXS: -1}),
    #("/JetHT/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD", {sXS: -1}),
    #("/JetHT/Run2018C-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD", {sXS: -1}),
    #("/JetHT/Run2018D-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD", {sXS: -1}),
    
    # dasgoclient --query="dataset=/JetHT/*2018*UL*MiniAODv2_NanoAODv9_GT36*/NANOAOD"
    # XS (cross-section) does not matter for data sample
    ("/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), # run 315257 to 316995
    ("/JetHT/Run2018B-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), # run 317080 to 319310
    ("/JetHT/Run2018C-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), # run 319337 to 320065
    ("/JetHT/Run2018D-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), # run 320500 to 325175

    ## SingleMuon
    # dasgoclient --query="dataset=/SingleMuon/*2018*UL*MiniAODv2_NanoAODv9_GT36*/NANOAOD"
    # XS (cross-section) does not matter for data sample    
    ("/SingleMuon/Run2018A-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}),    
    ("/SingleMuon/Run2018B-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}),
    ("/SingleMuon/Run2018C-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}),
    ("/SingleMuon/Run2018D-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}),
    
    ## EGamma
    # dasgoclient --query="dataset=/EGamma/*2018*UL*MiniAODv2_NanoAODv9_GT36*/NANOAOD"
    # XS (cross-section) does not matter for data sample    
    ("/EGamma/Run2018A-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), 
    ("/EGamma/Run2018B-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), 
    ("/EGamma/Run2018C-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), 
    ("/EGamma/Run2018D-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), 

    ## MET
    # dasgoclient --query="dataset=/MET/*2018*UL*MiniAODv2_NanoAODv9_GT36*/NANOAOD"
    # XS (cross-section) does not matter for data sample    
    ("/MET/Run2018A-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), 
    ("/MET/Run2018B-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), 
    ("/MET/Run2018C-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), 
    ("/MET/Run2018D-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {sXS: -1}), 
    

    #("", {sXS: }),    
])


sNDatasets                = "nDatasets"
sDataset                   = "dataset"
sNanoAOD_nFiles            = "nanoAOD_nFiles"
sNanoAOD                   = "nanoAOD"
sCross_section             = "cross_section"
sNEvents                   = "nEvents"
sSumEvents                 = "sumEvents"
sSkimmedNanoAOD_nFiles      = "skimmedNanoAOD_nFiles"
sSkimmedNanoAOD             = "skimmedNanoAOD"
sNFiles                    = "nFiles"
sampleDetail_dict_template = OD([
    (sCross_section,  -1.),
    (sNEvents,         0),
    (sSumEvents,       0),
    (sNDatasets,       0),
    (sDataset,        []),
    (sNanoAOD_nFiles,  0),
    (sNanoAOD,        []),
    #(sSkimmedNanoAOD_nFiles,  0),
    #(sSkimmedNanoAOD,        []),    
])

def getDatasetFiles(dataset):
    cmd1 = ['bash','-c', 'dasgoclient --query="file dataset=%s" --format=json'%(dataset)]
    if printLevel >= 10:
        print(f"cmd1: {cmd1}")
    output = subprocess.check_output(cmd1) # output in bytes
    output = output.decode("utf-8") # convert output in bytes to string
    output = json.loads(output) # convert output in 'string' to 'dict'
    nFiles = output['nresults']
    files  = []
    nEventsTotal = 0
    
    if nFiles != len(output['data']):
        print(f"nFiles != len(output['data']... something is wrong.. \t\t **** ERROR ****")
        exit(0)
        
    for iFile in range(nFiles):
        if len(output['data'][iFile]['file']) != 1:
            print(f"len(output['data'][iFile]['file']) != 1: Assumption of a single entry list 'output['data'][iFile]['file']' seems wrong... need to follow up.. \t\t **** ERROR **** ")
            exit(0)
            
        file_LFN = output['data'][iFile]['file'][0]['name']
        nEvents  = output['data'][iFile]['file'][0]['nevents']
        if printLevel >= 5:
            print(f"file_LFN: {file_LFN}, nEvents ({type(nEvents)}): {nEvents}, nEventsTotal: {nEventsTotal}  {output['data'][iFile]['file'][0]}")
            
        files.append( file_LFN )
        nEventsTotal += nEvents

    if printLevel >= 3:
        print(f"\ndataset: {dataset}, nEventsTotal: {nEventsTotal}, nFiles: {nFiles}, files: {files}")
    return nEventsTotal, nFiles, files
    

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
    parser.add_argument('-era',                 dest='era', type=str, default=Era_2018, choices=[Era_2016, Era_2017, Era_2018], required=False)
    parser.add_argument('-updateCrossSections', action='store_true', default=False, help='update cross-sections only')
    #parser.add_argument('-addSkimmedNanoAOD',   action='store_true', default=False, help='add skimmed NanoAOD files to the existing sample list')
    args=parser.parse_args()

    era                 = args.era
    updateCrossSections = args.updateCrossSections
    #addSkimmedNanoAOD   = args.addSkimmedNanoAOD
    print(f"era: {era}")
    print(f"{updateCrossSections = }")
    #print(f"{addSkimmedNanoAOD = }")


    list_datasetAndXs = None
    sFileSamplesInfo_toUse = None
    if era == Era_2018:
        list_datasetAndXs = list_datasetAndXs_2018

        
    sFileSamplesInfo_toUse = sFileSamplesInfo[era]
    sFileSamplesInfo_toUse = sFileSamplesInfo_toUse.replace('.json', '_v0.json')

    isSamplesListFromScratch = False
    samples_details = None       
    '''
    if updateCrossSections or addSkimmedNanoAOD:
        # update cross sections
        with open(sFileSamplesInfo[era]) as fSamplesInfo:
            samples_details = json.load(fSamplesInfo)
        print(f"samples_details.keys(): {samples_details.keys()}")
    else:
        samples_details = OD()
    '''
    if os.path.exists( sFileSamplesInfo[era] ):
        with open(sFileSamplesInfo[era]) as fSamplesInfo:
            samples_details = json.load(fSamplesInfo)
        print(f"samples_details.keys(): {samples_details.keys()}")
    else:
        samples_details = OD()
        isSamplesListFromScratch = True


    # Reset sNDatasets, sDataset, sSkimmedNanoAOD etc
    # Should not reset when running with updateCrossSections
    if not updateCrossSections:
        #for sampleName_ in samples_details:
        for iSample, (datasetName, datasetDetails) in enumerate(list_datasetAndXs.items()):
            datasetName_parts            = datasetName.split('/')
            sampleName                   = datasetName_parts[1]
            isMC                         = datasetName_parts[-1] == 'NANOAODSIM'
            if not isMC:
                # for data sample
                sampleName_part2 = (datasetName_parts[2]).split('-')[0] # 'Run2018A-UL2018_MiniAODv2_NanoAODv9-v2'
                sampleName = '%s_%s' % (sampleName, sampleName_part2)  # JetHT_Run2018A

            if sampleName not in samples_details: 
                samples_details[sampleName] = deepcopy(sampleDetail_dict_template)
                '''                 
                if isSamplesListFromScratch:
                    samples_details[sampleName] = deepcopy(sampleDetail_dict_template)
                else:
                    # when adding a new sample to existing samples.json, add in order as in list_datasetAndXs_2018
                    samples_details_tmp_ = samples_details.items()
                    samples_details_tmp_.insert(iSample, (sampleName, deepcopy(sampleDetail_dict_template)))
                    samples_details = samples_details_tmp_'''
            else: # reset sNDatasets, sDataset
                samples_details[sampleName][sNDatasets] = 0
                samples_details[sampleName][sDataset  ] = []
                samples_details[sampleName][sNEvents       ]         = 0
                samples_details[sampleName][sNanoAOD_nFiles]         = 0
                samples_details[sampleName][sNanoAOD       ]         = []


            samples_details[sampleName][sSkimmedNanoAOD] = OD()
            for skimName_ in sPathSkimmedNanoAODs[era]:
                samples_details[sampleName][sSkimmedNanoAOD]['%s_%s' % (skimName_, sNFiles)] = 0
                samples_details[sampleName][sSkimmedNanoAOD][skimName_                     ] = []

            # temperary fix
            #if "skimmedNanoAOD_nFiles" in samples_details[sampleName_]:
            #    samples_details[sampleName_].pop("skimmedNanoAOD_nFiles", None)


    # Now calculate..
    for datasetName, datasetDetails in list_datasetAndXs.items():
        datasetName_parts            = datasetName.split('/')
        sampleName                   = datasetName_parts[1]
        isMC                         = datasetName_parts[-1] == 'NANOAODSIM'
        #print(f"{datasetName = },  {isMC = }")

        if not isMC:
            # for data sample
            sampleName_part2 = (datasetName_parts[2]).split('-')[0] # 'Run2018A-UL2018_MiniAODv2_NanoAODv9-v2'
            sampleName = '%s_%s' % (sampleName, sampleName_part2)  # JetHT_Run2018A

        if updateCrossSections:
            if isMC:
                if sampleName in samples_details:
                    samples_details[sampleName][sCross_section] = datasetDetails[sXS]
                else:
                    print(f"Running updateCrossSections mode, but {sampleName} is not in samples_details \t **** ERROR **** \nTerminating...")
                    exit(0)
            
            continue


        # Add path to skimmed NanoAOD samples
        sDatasetExt = ''
        for sTmp1_ in (datasetName_parts[2]).split('_'): # RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2
            for sTmp2_ in sTmp1_.split('-'):
                if 'ext' in sTmp2_:
                    sDatasetExt = '_%s' % sTmp2_ # '_ext1'
        sSampleNameDir_used = '%s%s' % (datasetName_parts[1], sDatasetExt)

        sSampleTagDir_used = datasetName_parts[2]
        sSampleTagDir_used = sSampleTagDir_used.replace('_NanoAODv9', '')
        #if datasetName_parts[1] == 'SingleMuon':
        sSampleTagDir_used = sSampleTagDir_used.replace('v1', 'v*') # SingleMuon, EGamma  had DatasetTag v2/3 in MiniAOD and v1 in NanoAOD. So use wildcard charester *

        sDataType = 'MC' if isMC else 'Data'

        # Loop over multiple skim version we have
        for skimName_ in sPathSkimmedNanoAODs[era]:
            sPathSkimmedNanoAODs_toUse = sPathSkimmedNanoAODs[era][skimName_][sDataType] 
            # /eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/$SAMPLENAME/r1/PNet_*.root           
            # /eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/PNet_v1_2023_10_06/$SAMPLETAG/$SAMPLENAME/r*/PNet_*.root
            sPathSkimmedNanoAODs_toUse = sPathSkimmedNanoAODs_toUse.replace('$SAMPLENAME', sSampleNameDir_used)
            sPathSkimmedNanoAODs_toUse = sPathSkimmedNanoAODs_toUse.replace('$SAMPLETAG',  sSampleTagDir_used)

            sSkimmedNanoAODs = []
            if "*" in sPathSkimmedNanoAODs_toUse: sSkimmedNanoAODs.extend( glob.glob(sPathSkimmedNanoAODs_toUse) )
            else:                                 sSkimmedNanoAODs.append( sPathSkimmedNanoAODs_toUse )
            print(f"{sPathSkimmedNanoAODs_toUse = }, {sSkimmedNanoAODs = }")
            
            '''if sSkimmedNanoAOD_nFiles not in samples_details[sampleName]:
                samples_details[sampleName][sSkimmedNanoAOD_nFiles]  = len(sSkimmedNanoAODs)
                samples_details[sampleName][sSkimmedNanoAOD]         = sSkimmedNanoAODs
            else:
                samples_details[sampleName][sSkimmedNanoAOD_nFiles] += len(sSkimmedNanoAODs)
                samples_details[sampleName][sSkimmedNanoAOD].extend(   sSkimmedNanoAODs )'''
            samples_details[sampleName][sSkimmedNanoAOD]['%s_%s' % (skimName_, sNFiles)] += len(   sSkimmedNanoAODs )
            samples_details[sampleName][sSkimmedNanoAOD][skimName_                     ].extend(   sSkimmedNanoAODs )

        #if addSkimmedNanoAOD:
        #    continue
        

        #if sampleName not in samples_details:
        #    samples_details[sampleName] = deepcopy(sampleDetail_dict_template)

        samples_details[sampleName][sNDatasets]         += 1
        samples_details[sampleName][sDataset  ].append(    datasetName )
        
        if datasetName_parts[-1] == 'NANOAODSIM':
            # for MC sample
            samples_details[sampleName][sCross_section] = datasetDetails[sXS]
        else:
            # for data sample
            if sCross_section in samples_details[sampleName]:
                del samples_details[sampleName][sCross_section]
            if sSumEvents in samples_details[sampleName]:
                del samples_details[sampleName][sSumEvents]
        
        nEventsTotal, nFiles, files = getDatasetFiles(datasetName)
        samples_details[sampleName][sNEvents       ]         += nEventsTotal
        samples_details[sampleName][sNanoAOD_nFiles]         += nFiles
        samples_details[sampleName][sNanoAOD       ].extend(    files        )
        

    # Running the samples_prepare.py on exising Sample.json changes order of samples. 
    # Hence, to forcefully follow order of samples as in list_datasetAndXs, 
    # order samples_details keys to follow order in list_datasetAndXs
    samples_details_inOrder = OD()
    for datasetName, datasetDetails in list_datasetAndXs.items():
        datasetName_parts            = datasetName.split('/')
        sampleName                   = datasetName_parts[1]
        isMC                         = datasetName_parts[-1] == 'NANOAODSIM'
        if not isMC:
            # for data sample
            sampleName_part2 = (datasetName_parts[2]).split('-')[0] # 'Run2018A-UL2018_MiniAODv2_NanoAODv9-v2'
            sampleName = '%s_%s' % (sampleName, sampleName_part2)  # JetHT_Run2018A

        samples_details_inOrder[sampleName] = samples_details[sampleName]
    # Replace samples_details with samples_details_inOrder
    samples_details = samples_details_inOrder 


    if printLevel >= 0:
        print("\n\nsamples:: \n",json.dumps(samples_details, indent=4))

    with open(sFileSamplesInfo_toUse, "w") as fSampleInfo:
        json.dump(samples_details, fSampleInfo, indent=4)

        print(f"\n\nNew sample list wrote to {sFileSamplesInfo_toUse}")
        print(f"Copy {sFileSamplesInfo_toUse} to {sFileSamplesInfo[era]}   if you are satistied with {sFileSamplesInfo_toUse}")
        print(f"i.e.   cp {sFileSamplesInfo_toUse} {sFileSamplesInfo[era]}")
    
    

    


