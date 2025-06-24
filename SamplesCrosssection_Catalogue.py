from collections import OrderedDict as OD

sXS13TeV           = "xs13TeV"

list_XSs = OD([
    
    ## QCD_bEnriched_HT*
    # dasgoclient --query="dataset=/QCD_bEnriched_HT*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"    
    ("QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8",   {sXS13TeV: 1122000.0   , }),
    ("QCD_bEnriched_HT200to300_TuneCP5_13TeV-madgraph-pythia8",   {sXS13TeV:   79760.0   , }),
    ("QCD_bEnriched_HT300to500_TuneCP5_13TeV-madgraph-pythia8",   {sXS13TeV:   16600.0   , }),
    ("QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8",   {sXS13TeV:    1503.0   , }),    
    ("QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8",  {sXS13TeV:     297.4   , }),
    ("QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8", {sXS13TeV:      48.08  , }),    
    ("QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8", {sXS13TeV:       3.9510 , }),
    ("QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8",  {sXS13TeV:       0.6957 , }),
    

    ## QCD_HT*_BGenFilter
    # dasgoclient --query="dataset=/QCD_HT*_BGenFilter*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("QCD_HT100to200_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",   {sXS13TeV: 1266000.0    , }),
    ("QCD_HT200to300_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",   {sXS13TeV:  109900.0    , }),
    ("QCD_HT300to500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",   {sXS13TeV:   27360.0    , }),
    ("QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",   {sXS13TeV:    2991.0    , }),
    ("QCD_HT700to1000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",  {sXS13TeV:     731.8    , }),    
    ("QCD_HT1000to1500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8", {sXS13TeV:     139.3    , }),
    ("QCD_HT1500to2000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8", {sXS13TeV:      14.74   , }),
    ("QCD_HT2000toInf_BGenFilter_TuneCP5_13TeV-madgraph-pythia8",  {sXS13TeV:       3.09   , }),


    ## QCD_HT* PSWeights-madgraph - QCDIncl LO
    # dasgoclient --query="dataset=/QCD_HT*PSWeight*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8",    {sXS13TeV: 187700000.0     }),
    ("QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8",   {sXS13TeV:  23640000.0    , }),
    ("QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8",   {sXS13TeV:   1546000.0    , }),
    ("QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8",   {sXS13TeV:    321600.0    , }),
    ("QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8",   {sXS13TeV:     30310.0    , }),
    ("QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8",  {sXS13TeV:      6364.0    , }),
    ("QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8", {sXS13TeV:      1117.0    , }), ### Correction needed
    ("QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8", {sXS13TeV:       108.4    , }),
    ("QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8",  {sXS13TeV:        22.36   , }),

    
    ## QCD_HT* TuneCP5 madgraphMLM - QCDIncl LO, MatrixElement-PartonShower maching at NLO(?)
    # dasgoclient --query="dataset=/QCD_HT*TuneCP5_13TeV-madgraphMLM*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8",    {sXS13TeV: 187700000.0   }),
    ("QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8",   {sXS13TeV:  23500000.0   }),
    ("QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8",   {sXS13TeV:   1552000.0   }),
    ("QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8",   {sXS13TeV:    321100.0   }),
    ("QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8",   {sXS13TeV:     30250.0   }),    
    ("QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8",  {sXS13TeV:      6398.0   }),
    ("QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV:      1122.0   }),
    ("QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV:       109.4   }),
    ("QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8",  {sXS13TeV:        21.74  }),

    

    ## TTbar - NLO powheg
    # dasgoclient --query="dataset=/TT*_TuneCP5_13TeV*powheg*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("TTToHadronic_TuneCP5_13TeV-powheg-pythia8",                   {sXS13TeV: 380.133 , }), # 831.8 * 0.457
    ("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8",               {sXS13TeV: 364.328 , }), # 831.8 * 0.438
    ("TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",                      {sXS13TeV:  87.339 , }), # 831.8 * 0.105

    ## TTbar Jets - NLO amcatnlo
    # dasgoclient --query="dataset=/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8",                   {sXS13TeV: 831.8       }),

    ## TTbar Jets - LO madgraph
    # dasgoclient --query="dataset=/TTJets_*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM" 
    ("TTJets_TuneCP5_13TeV-madgraphMLM-pythia8",                    {sXS13TeV: 831.8       }),
    ("TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",        {sXS13TeV:   2.4234    }),
    ("TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8",       {sXS13TeV:   0.9818    }),
    ("TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8",      {sXS13TeV:   0.1714    }),
    ("TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8",       {sXS13TeV:   0.001966  }),

    ## TTbar Jets To LNu - LO madgraph
    # dasgoclient --query="dataset=/TTJets_*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8",    {sXS13TeV: 182.164 }), # 831.8 * 0.219
    ("TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV: 182.164 }), # 831.8 * 0.219
    ("TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8",             {sXS13TeV:  87.339 }), # 831.8 * 0.105


    ## ST NLO
    # dasgoclient --query="dataset=/ST*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8",     {sXS13TeV: 134.2   , }),
    ("ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8", {sXS13TeV:  80.0   , }),
    ("ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",            {sXS13TeV:  39.65  , }), # 79.3/2 = 39.65
    ("ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",        {sXS13TeV:  39.65  , }), # 79.3/2 = 39.65
    ("ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",      {sXS13TeV:  21.63  , }),  # 79.3/2*(1 - (0.6741*0.6741)) : BR(W->Hadrons)=0.6741. tW has W(from t) and W. BR(tW->NoFully hadronic)=(1 - (0.6741*0.6741))=0.546.
    ("ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",  {sXS13TeV:  21.63  , }),  # 79.3/2*(1 - (0.6741*0.6741)) : Andrew: 21.61 
    ("ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8",        {sXS13TeV:   5.041 , }),
    ("ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8",          {sXS13TeV:   4.831 , }),

    ## tX            
    # dasgoclient --query="dataset=/t*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8",   {sXS13TeV:  0.0942 , }),
    ("ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8",          {sXS13TeV:  0.839  , }), # https://cds.cern.ch/record/2227475/files/CERN-2017-002-M.pdf?version=1#page=180
    ("ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8",          {sXS13TeV:  0.6008 , }), # https://cds.cern.ch/record/2227475/files/CERN-2017-002-M.pdf?version=1#page=180
    

    ## ZJets
    # dasgoclient --query="dataset=/ZJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"    
    ("ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV: 1012.0   , }),
    ("ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV:  114.2   , }),
    ("ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV:   25.34  , }), ### Correction needed
    ("ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV:   12.99  , }), ### Correction needed

    ## DYJetsToLL NLO
    # dasgoclient --query="dataset=/DY*Jets*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8",  {sXS13TeV: 18610    }),
    ("DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8",      {sXS13TeV:  6077.22 }),
    
    # DYJetsToLL LO Incl
    ("DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8",                     {sXS13TeV: 18610    , }),
    ("DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",                         {sXS13TeV:  6077.22 , }),
    #("DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",                    {sXS13TeV:  6077.22  }),

    # DYJetsToLL LO HT
    ("DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",    {sXS13TeV:   158.7445   , }),
    ("DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",   {sXS13TeV:   159.1984   , }),
    ("DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",   {sXS13TeV:    43.5384   , }),
    ("DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",   {sXS13TeV:     5.9141   , }),
    ("DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",   {sXS13TeV:     1.4377   , }),
    ("DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",  {sXS13TeV:     0.6443   , }),
    ("DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8", {sXS13TeV:     0.1511   , }),
    ("DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",  {sXS13TeV:     0.003390 , }),    

    # dasgoclient --query="dataset=/ZJetsToNuNu*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",   {sXS13TeV: 302.96     , }),  
    ("ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",   {sXS13TeV:  82.9239   , }),  
    ("ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",   {sXS13TeV:  11.2574   , }),  
    ("ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",   {sXS13TeV:   2.73349  , }),  
    ("ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",  {sXS13TeV:   1.22321  , }),  
    ("ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV:   0.28526  , }),  
    ("ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",  {sXS13TeV:   0.006366 , }),  
    
    ## WJetsToQQ
    # dasgoclient --query="dataset=/WJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV: 2549.0    , }),
    ("WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV:  276.5    , }),
    ("WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV:   59.25   }),
    ("WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8", {sXS13TeV:   28.75   }),


    ## WJetsToLNu
    # dasgoclient --query="dataset=/W*Jets*ToLNu*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8",                   {sXS13TeV: 61526.7     , }), # NLO sample
    #("WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8",               {sXS13TeV: }), # NLO sample
    #("WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8",               {sXS13TeV: }), # NLO sample
    #("WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8",               {sXS13TeV: }), # NLO sample    
    ("WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",                    {sXS13TeV: 61526.7      }),
    ("WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8",         {sXS13TeV:  1440.0     ,  }),
    ("WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",        {sXS13TeV:  1431.0     ,  }),
    ("WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",        {sXS13TeV:   382.1     ,  }),
    ("WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",        {sXS13TeV:    51.54    ,  }),
    ("WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",        {sXS13TeV:    12.49    ,  }),
    ("WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",       {sXS13TeV:     5.619   ,  }), 
    ("WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",      {sXS13TeV:     1.321   ,  }),
    ("WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",       {sXS13TeV:     0.02992 ,  }),
    #("W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",                   {sXS13TeV: 10167.9     }),
    #("W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",                   {sXS13TeV: 3199.45     }),
    #("W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",                   {sXS13TeV:  941.16     }),
    #("W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",                   {sXS13TeV:  439.08     }),


    ## VV, VVV
    # dasgoclient --query="dataset=/Z*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    # dasgoclient --query="dataset=/W*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("ZZ_TuneCP5_13TeV-pythia8",                              {sXS13TeV:   16.523  ,     }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    ("ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8", {sXS13TeV:   3.25 * 1.21,  }), # https://indico.cern.ch/event/439995/contributions/1094416/attachments/1143460/1638648/diboson_final.pdf
    ("ZZTo2Q2Nu_TuneCP5_13TeV-amcatnloFXFX-pythia8",          {sXS13TeV:   4.07 * 1.22   }), # https://indico.cern.ch/event/439995/contributions/1094416/attachments/1143460/1638648/diboson_final.pdf    
    ("WZ_TuneCP5_13TeV-pythia8",                              {sXS13TeV:   47.13  ,      }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    ("WW_TuneCP5_13TeV-pythia8",                              {sXS13TeV:  118.7  ,       }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV

    ("ZZZ_TuneCP5_13TeV-amcatnlo-pythia8",          {sXS13TeV:    0.01398  ,  }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    ("WZZ_TuneCP5_13TeV-amcatnlo-pythia8",          {sXS13TeV:    0.05565  ,  }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    ("WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8",       {sXS13TeV:    0.1651  ,  }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    ("WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8",       {sXS13TeV:    0.2086  ,  }), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    

    #("",          {sXS13TeV:    }), 
    
    ## GGF HToBB
    # dasgoclient --query="dataset=/*HToBB*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    # dasgoclient --query="dataset=/*HTo*bb*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8",      {sXS13TeV:  48.61 * 0.582    }), # GGF H * BR = 28.291 N3LO
    ("GluGluHToBB_Pt-200ToInf_M-125_TuneCP5_MINLO_13TeV-powheg-pythia8", {sXS13TeV:   0.2740   ,  }),
    
    ("VBFHToBB_M-125_TuneCP5_13TeV-powheg-pythia8",                      {sXS13TeV:  3.766 * 0.582  }), # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
    ("VBFHToBB_M-125_TuneCH3_13TeV-powheg-herwig",                       {sXS13TeV:  3.766 * 0.582  }), # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
    ("VBFHToBB_M-125_dipoleRecoilOn_TuneCP5_13TeV-powheg-pythia8",       {sXS13TeV:  2.250  }), # HIG-21-020
    ("VBFWH_HToBB_WToLNu_M-125_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8",       {sXS13TeV:  0.4378  }), # https://xsdb-temp.app.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=DAS%3DVBFWH_HToBB_WToLNu_M-125_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8
    ("VBFHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8",                   {sXS13TeV:  1  }),

    ("WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8",            {sXS13TeV:  0.831 * 0.582 * 0.647    }),
    ("WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8",           {sXS13TeV:  0.831 * 0.582 * 3*0.1086 ,  }),
    ("WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8",           {sXS13TeV:  0.527 * 0.582 * 0.647    }),
    ("WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8",          {sXS13TeV:  0.527 * 0.582 * 3*0.1086 ,  }),
    ("WHToMuMuG_M125_Dalitz_012j_TuneCP5_13TeV_amcatnloFXFX_pythia8",    {sXS13TeV:  1    }),
    ("WplusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8",                 {sXS13TeV:  1    }),
    ("WminusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8",                {sXS13TeV:  1    }),    

    ("ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8",                {sXS13TeV:  0.758 * 0.582 * 0.699     }),
    ("ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8",                {sXS13TeV:  0.758 * 0.582 * 3*0.0337 ,  }),
    ("ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8",              {sXS13TeV:  0.758 * 0.582 * 0.20      }),
    ("ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8",                {sXS13TeV:  0.758 * 0.582 * 0.1512   , }),
    ("ZHToMuMuG_M125_Dalitz_012j_TuneCP5_13TeV_amcatnloFXFX_pythia8",    {sXS13TeV:  1    }),
    ("ZHToTauTau_M125_CP5_13TeV-powheg-pythia8",                         {sXS13TeV:  1    }),    

    ("ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8",              {sXS13TeV:  0.123 * 0.582 * 0.699     }),
    ("ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8",              {sXS13TeV:  0.123 * 0.582 * 3*0.0337  }),
    ("ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8",            {sXS13TeV:  0.123 * 0.582 * 0.20      }),
    ("ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8",              {sXS13TeV:  0.123 * 0.582 * 0.1512    }),

    ("ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8",                        {sXS13TeV:  0.5071 * 0.5824  , }),
    ("ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8",                     {sXS13TeV:  0.5071 * (1 - 0.5824)  }),
    ("ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8",                    {sXS13TeV:  1  }),

    ("bbHToBB_M-125_4FS_yt2_TuneCP5-13TeV-amcatnlo-pythia8",             {sXS13TeV:  1  ,               }),
    

    # , 
    ## SUSY_GluGluH_01J_HToAATo4B_M-*   and   SUSY_GluGluH_01J_HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*GluGluH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("SUSY_GluGluH_01J_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61}),
    ("SUSY_GluGluH_01J_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61}),
    ("SUSY_GluGluH_01J_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61}),
    ("SUSY_GluGluH_01J_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61}),
    ("SUSY_GluGluH_01J_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61}),
    ("SUSY_GluGluH_01J_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61}),
    ("SUSY_GluGluH_01J_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61}),
    ("SUSY_GluGluH_01J_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61}),
    ("SUSY_GluGluH_01J_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61}),
    ("SUSY_GluGluH_01J_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61}),
    ("SUSY_GluGluH_01J_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61}),

    # SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (308) / (5379) = 5.726e-02 +- 3.168e-03	Matching efficiency = 0.5 +/- 0.0     Cross-section = 48.61 pb * 0.057	
    ("SUSY_GluGluH_01J_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61 * 0.057 ,  }), # filter efficiency 0.05
    ("SUSY_GluGluH_01J_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61 * 0.057 ,  }),
    ("SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61 * 0.057 ,  }),
    ("SUSY_GluGluH_01J_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61 * 0.057 ,  }),
    ("SUSY_GluGluH_01J_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61 * 0.057 ,  }),
    ("SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61 * 0.057 ,  }),
    ("SUSY_GluGluH_01J_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61 * 0.057 ,  }),
    ("SUSY_GluGluH_01J_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61 * 0.057 ,  }),
    ("SUSY_GluGluH_01J_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61 * 0.057 ,  }),
    ("SUSY_GluGluH_01J_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61 * 0.057 ,  }),
    ("SUSY_GluGluH_01J_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 48.61 * 0.057 ,  }),


    ## VBF HToAATo4B_M-* and VBF HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*VBF*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("SUSY_VBFH_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 }),
    ("SUSY_VBFH_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 }),
    ("SUSY_VBFH_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 }),
    ("SUSY_VBFH_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 }),
    ("SUSY_VBFH_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 }),
    ("SUSY_VBFH_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 }),
    ("SUSY_VBFH_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 }),
    ("SUSY_VBFH_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 }),
    ("SUSY_VBFH_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 }),
    ("SUSY_VBFH_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 }),
    ("SUSY_VBFH_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 }),

    # SUSY_VBFH_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8		Filter efficiency (event-level)= (1730) / (10000) = 1.730e-01 +- 3.782e-03	Matching efficiency = 1.0 +/- 0.0 	Cross-section = 3.766 pb * 0.173 
    ("SUSY_VBFH_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 * 0.173 ,  }),
    ("SUSY_VBFH_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 * 0.173 ,  }),
    ("SUSY_VBFH_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 * 0.173 ,  }),
    ("SUSY_VBFH_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 * 0.173 ,  }),
    ("SUSY_VBFH_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 * 0.173 ,  }),
    ("SUSY_VBFH_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 * 0.173 ,  }),
    ("SUSY_VBFH_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 * 0.173 ,  }),
    ("SUSY_VBFH_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 * 0.173 ,  }),
    ("SUSY_VBFH_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 * 0.173 ,  }),
    ("SUSY_VBFH_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 * 0.173 ,  }),
    ("SUSY_VBFH_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 3.766 * 0.173 ,  }),


    ## WH HToAATo4B-M-* and WH HToAATo4B_Pt150_-M-*
    # dasgoclient --query="dataset=/SUSY*WH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM" 
    ("SUSY_WH_WToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358}),
    ("SUSY_WH_WToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358}),
    ("SUSY_WH_WToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358}),
    ("SUSY_WH_WToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358}),
    ("SUSY_WH_WToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358}),
    ("SUSY_WH_WToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358}),
    ("SUSY_WH_WToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358}),
    ("SUSY_WH_WToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358}),
    ("SUSY_WH_WToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358}),
    ("SUSY_WH_WToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358}),
    ("SUSY_WH_WToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358}),

    # SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (1317) / (10000) = 1.317e-01 +- 3.382e-03	Matching efficiency = 1.0 +/- 0.0	Cross-section = 1.358 pb * 0.132 
    ("SUSY_WH_WToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358 * 0.132 ,  }),
    ("SUSY_WH_WToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358 * 0.132 ,  }),
    ("SUSY_WH_WToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358 * 0.132 ,  }),
    ("SUSY_WH_WToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358 * 0.132 ,  }),
    ("SUSY_WH_WToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358 * 0.132 ,  }),
    ("SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358 * 0.132 ,  }),
    ("SUSY_WH_WToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358 * 0.132 ,  }),
    ("SUSY_WH_WToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358 * 0.132 ,  }),
    ("SUSY_WH_WToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358 * 0.132 ,  }),
    ("SUSY_WH_WToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358 * 0.132 ,  }),
    ("SUSY_WH_WToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  1.358 * 0.132 ,  }),


    ## ZH HToAATo4B-M-* and ZH HToAATo4B_Pt150_-M-*
    # dasgoclient --query="dataset=/SUSY*ZH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("SUSY_ZH_ZToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 }),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 }),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 }),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 }),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 }),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 }),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 }),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 }),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 }),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 }),
    ("SUSY_ZH_ZToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 }),
    
    # SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (1289) / (10000) = 1.289e-01 +- 3.351e-03	Matching efficiency = 1.0 +/- 0.0		Cross-section = 0.880 pb * 0.129 
    ("SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 * 0.129 ,  }),
    ("SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 * 0.129 ,  }),
    ("SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 * 0.129 ,  }),
    ("SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 * 0.129 ,  }),
    ("SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 * 0.129 ,  }),
    ("SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 * 0.129 ,  }),
    ("SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 * 0.129 ,  }),
    ("SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 * 0.129 ,  }),
    ("SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 * 0.129 ,  }),
    ("SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 * 0.129 ,  }),
    ("SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV:  0.880 * 0.129 ,  }),


    ## ttH HToAATo4B_M-* and HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*TTH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("SUSY_TTH_TTToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071  }),

    # SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (2850) / (10000) = 2.850e-01 +- 4.514e-03	Matching efficiency = 1.0 +/- 0.0 Cross-section = 0.5071 pb * 0.285 
    ("SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071 * 0.285 ,  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071 * 0.285 ,  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071 * 0.285 ,  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071 * 0.285 ,  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071 * 0.285 ,  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071 * 0.285 ,  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071 * 0.285 ,  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071 * 0.285 ,  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071 * 0.285 ,  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071 * 0.285 ,  }),
    ("SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8", {sXS13TeV: 0.5071 * 0.285 ,  }),
    

    ## H->aa->4Tau
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYGluGluHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    
    ("SUSYVBFHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVBFHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    
    ("SUSYVHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYVHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),

    ("SUSYttHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),
    ("SUSYttHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8", {sXS13TeV: 1   }),



    #("", {sXS13TeV: }),    
])