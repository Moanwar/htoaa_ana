from collections import OrderedDict as OD

sXS           = "xs"
sNameSp       = "nameSp"
sNanoSkimv2   = "skim_v2"
sNEvents      = "nEvents"
sSumEvents    = "sumEvents"
sNEvtSkimv2   = "%s_%s" % (sNanoSkimv2, sNEvents)
sSumEvtSkimv2 = "%s_%s" % (sNanoSkimv2, sSumEvents)

list_datasetAndXs_2018 = OD([
    
    ## QCD_bEnriched_HT*
    # dasgoclient --query="dataset=/QCD_bEnriched_HT*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"    
    ("/QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS: 1122000.0   , sNEvtSkimv2: 36118282, sSumEvtSkimv2: 36118282}),
    ("/QCD_bEnriched_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:   79760.0   , sNEvtSkimv2: 18462183, sSumEvtSkimv2: 18462183}),
    ("/QCD_bEnriched_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:   16600.0   , sNEvtSkimv2: 11197722, sSumEvtSkimv2: 11197722}),
    ("/QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:    1503.0   , sNEvtSkimv2:  9246898, sSumEvtSkimv2:  9246898}),    
    ("/QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:     297.4   , sNEvtSkimv2:  1844165, sSumEvtSkimv2:  1844165}),
    ("/QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:      48.08  , sNEvtSkimv2:  1330829, sSumEvtSkimv2:  1330829}),    
    ("/QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:       3.9510 , sNEvtSkimv2: 1431254, sSumEvtSkimv2:  1431254}),
    ("/QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:       0.6957 , sNEvtSkimv2: 1429280, sSumEvtSkimv2:  1429280}),
    

    ## QCD_HT*_BGenFilter
    # dasgoclient --query="dataset=/QCD_HT*_BGenFilter*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/QCD_HT100to200_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sXS: 1266000.0    , sNEvtSkimv2: 35884072, sSumEvtSkimv2: 35884072}),
    ("/QCD_HT200to300_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sXS:  109900.0    , sNEvtSkimv2: 14830551, sSumEvtSkimv2: 14830551}),
    ("/QCD_HT300to500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sXS:   27360.0    , sNEvtSkimv2: 14144826, sSumEvtSkimv2: 14144826}),
    ("/QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sXS:    2991.0    , sNEvtSkimv2:  8015859, sSumEvtSkimv2:  8015859}),
    ("/QCD_HT700to1000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",  {sXS:     731.8    , sNEvtSkimv2:  4642245, sSumEvtSkimv2:  4642245}),    
    ("/QCD_HT1000to1500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:     139.3    , sNEvtSkimv2:  1537452, sSumEvtSkimv2:  1537452}),
    ("/QCD_HT1500to2000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:      14.74   , sNEvtSkimv2:  1263157, sSumEvtSkimv2:  1263157}),
    ("/QCD_HT2000toInf_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",  {sXS:       3.09   , sNEvtSkimv2:  1300672, sSumEvtSkimv2:  1300672}),


    ## QCD_HT* PSWeights-madgraph - QCDIncl LO
    # dasgoclient --query="dataset=/QCD_HT*PSWeight*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {sXS: 187700000.0     }),
    ("/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:  23640000.0    , sNEvtSkimv2: 84461486, sSumEvtSkimv2: 84461486}),
    ("/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:   1546000.0    , sNEvtSkimv2: 57336623, sSumEvtSkimv2: 57336623}),
    ("/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:    321600.0    , sNEvtSkimv2: 61705174, sSumEvtSkimv2: 61705174}),
    ("/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:     30310.0    , sNEvtSkimv2: 49184771, sSumEvtSkimv2: 49184771}),
    ("/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",  {sXS:      6364.0    , sNEvtSkimv2: 48506751, sSumEvtSkimv2: 48506751}),
    ("/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:      1117.0    , sNEvtSkimv2: 13026372, sSumEvtSkimv2: 13026372}), ### Correction needed
    ("/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:       108.4    , sNEvtSkimv2: 10871473, sSumEvtSkimv2: 10871473}),
    ("/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:        22.36   , sNEvtSkimv2:  5374711, sSumEvtSkimv2:  5374711}),

    
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
    ("/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {sXS: 380.133 , sNEvtSkimv2: 343248000, sSumEvtSkimv2: 338171646}), # 831.8 * 0.457
    ("/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {sXS: 364.328 , sNEvtSkimv2: 478982000, sSumEvtSkimv2: 471724160}), # 831.8 * 0.438
    ("/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                      {sXS:  87.339 , sNEvtSkimv2: 146010000, sSumEvtSkimv2: 143887383}), # 831.8 * 0.105

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
    ("/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",     {sXS: 134.2   , sNEvtSkimv2: 198569000, sSumEvtSkimv2: 197426478}),
    ("/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  80.0   , sNEvtSkimv2:  99264000, sSumEvtSkimv2:  98724306}),
    ("/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",            {sXS:  39.65  , sNEvtSkimv2:   7956000, sSumEvtSkimv2:   7953729}), # 79.3/2 = 39.65
    ("/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",        {sXS:  39.65  , sNEvtSkimv2:   7749000, sSumEvtSkimv2:   7747316}), # 79.3/2 = 39.65
    ("/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",      {sXS:  21.63  , sNEvtSkimv2:  11270430, sSumEvtSkimv2:  11269836}),  # 79.3/2*(1 - (0.6741*0.6741)) : BR(W->Hadrons)=0.6741. tW has W(from t) and W. BR(tW->NoFully hadronic)=(1 - (0.6741*0.6741))=0.546.
    ("/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:  21.63  , sNEvtSkimv2:  11015956, sSumEvtSkimv2:  11012326}),  # 79.3/2*(1 - (0.6741*0.6741)) : Andrew: 21.61 
    ("/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sXS:   5.041 , sNEvtSkimv2:  16391000, sSumEvtSkimv2:  10793026}),
    ("/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sXS:   4.831 , sNEvtSkimv2:  19365999, sSumEvtSkimv2:  12607741}),

    ## tX            
    # dasgoclient --query="dataset=/t*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:  0.0942 , sNEvtSkimv2:  11916000, sSumEvtSkimv2:  3229640}),
    ("/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",          {sXS:  0.839  , sNEvtSkimv2:  32793815, sSumEvtSkimv2: 32793815}), # https://cds.cern.ch/record/2227475/files/CERN-2017-002-M.pdf?version=1#page=180
    ("/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",          {sXS:  0.6008 , sNEvtSkimv2:  27686862, sSumEvtSkimv2: 27686862}), # https://cds.cern.ch/record/2227475/files/CERN-2017-002-M.pdf?version=1#page=180
    

    ## ZJets
    # dasgoclient --query="dataset=/ZJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"    
    ("/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS: 1012.0   , sNEvtSkimv2: 14808858, sSumEvtSkimv2: 14808858}),
    ("/ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:  114.2   , sNEvtSkimv2: 13930474, sSumEvtSkimv2: 13930474}),
    ("/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:   25.34  , sNEvtSkimv2: 12029507, sSumEvtSkimv2: 11668622}), ### Correction needed
    ("/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:   12.99  , sNEvtSkimv2:  9681521, sSumEvtSkimv2:  9003815}), ### Correction needed

    ## DYJetsToLL NLO
    # dasgoclient --query="dataset=/DY*Jets*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS: 18610    }),
    ("/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",      {sXS:  6077.22 }),
    
    # DYJetsToLL LO Incl
    ("/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                     {sXS: 18610    , sNEvtSkimv2: 99288125, sSumEvtSkimv2: 99288125}),
    ("/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                         {sXS:  6077.22 , sNEvtSkimv2: 96233328, sSumEvtSkimv2: 96233328}),
    #("/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v1/NANOAODSIM",                    {sXS:  6077.22  }),

    # DYJetsToLL LO HT
    ("/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {sXS:   158.7445   , sNEvtSkimv2: 17004433, sSumEvtSkimv2: 17004433}),
    ("/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:   159.1984   , sNEvtSkimv2: 26202328, sSumEvtSkimv2: 26202328}),
    ("/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:    43.5384   , sNEvtSkimv2: 18455718, sSumEvtSkimv2: 18455718}),
    ("/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:     5.9141   , sNEvtSkimv2:  8908406, sSumEvtSkimv2:  8908406}),
    ("/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:     1.4377   , sNEvtSkimv2:  7035971, sSumEvtSkimv2:  7035971}),
    ("/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:     0.6443   , sNEvtSkimv2:  6678036, sSumEvtSkimv2:  6678036}),
    ("/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:     0.1511   , sNEvtSkimv2:  6166852, sSumEvtSkimv2:  6166852}),
    ("/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:     0.003390 , sNEvtSkimv2:  1978203, sSumEvtSkimv2:  1978203}),    

    # dasgoclient --query="dataset=/ZJetsToNuNu*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS: 302.96     , sNEvtSkimv2: 28876062, sSumEvtSkimv2: 28876062}),  
    ("/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:  82.9239   , sNEvtSkimv2: 22749608, sSumEvtSkimv2: 22749608}),  
    ("/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:  11.2574   , sNEvtSkimv2: 19810491, sSumEvtSkimv2: 19810491}),  
    ("/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sXS:   2.73349  , sNEvtSkimv2: 5067605, sSumEvtSkimv2: 5067605}),  
    ("/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:   1.22321  , sNEvtSkimv2: 1948816, sSumEvtSkimv2: 1948816}),  
    ("/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:   0.28526  , sNEvtSkimv2: 381695, sSumEvtSkimv2: 381695}),  
    ("/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sXS:   0.006366 , sNEvtSkimv2: 268224, sSumEvtSkimv2: 268224}),  
    
    ## WJetsToQQ
    # dasgoclient --query="dataset=/WJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS: 2549.0    , sNEvtSkimv2: 14494966, sSumEvtSkimv2: 14494966}),
    ("/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:  276.5    , sNEvtSkimv2:  9335298, sSumEvtSkimv2:  9335298}),
    ("/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:   59.25   }),
    ("/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS:   28.75   }),


    ## WJetsToLNu
    # dasgoclient --query="dataset=/W*Jets*ToLNu*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                   {sXS: 61526.7     , sNEvtSkimv2: 29117828, sSumEvtSkimv2: 12361080}), # NLO sample
    #("/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {sXS: }), # NLO sample
    #("/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {sXS: }), # NLO sample
    #("/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {sXS: }), # NLO sample    
    #("/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                    {sXS: 61526.7      }),
    ("/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",         {sXS:  1440.0     , sNEvtSkimv2:  66569448, sSumEvtSkimv2: 66569448}),
    ("/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sXS:  1431.0     , sNEvtSkimv2:  51541593, sSumEvtSkimv2: 51541593}),
    ("/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sXS:   382.1     , sNEvtSkimv2:  58331446, sSumEvtSkimv2: 58331446}),
    ("/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sXS:    51.54    , sNEvtSkimv2:   7444030, sSumEvtSkimv2:  7444030}),
    #("/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",   {sXS:    51.54    }),
    ("/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sXS:    12.49    , sNEvtSkimv2:   7718765, sSumEvtSkimv2:  7718765}),
    #("/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",   {sXS:    12.49    }),
    ("/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sXS:     5.619   , sNEvtSkimv2:   7306187, sSumEvtSkimv2:  7296689}), ### COrrection needed
    #("/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",  {sXS:     5.619   }),
    ("/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",      {sXS:     1.321   , sNEvtSkimv2:   6481518, sSumEvtSkimv2:  6481518}),
    #("/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM", {sXS:     1.321   }),
    ("/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",       {sXS:     0.02992 , sNEvtSkimv2:   2097648, sSumEvtSkimv2:  2097648}),
    #("/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {sXS: 10167.9     }),
    #("/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {sXS: 3199.45     }),
    #("/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {sXS:  941.16     }),
    #("/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                   {sXS:  439.08     }),


    ## VV, VVV
    # dasgoclient --query="dataset=/Z*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    # dasgoclient --query="dataset=/W*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                              {sXS:   16.523  ,    sNEvtSkimv2:   3526000, sSumEvtSkimv2:  3526000}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    ("/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:   3.25 * 1.21, sNEvtSkimv2:  29357938, sSumEvtSkimv2: 15840251}), # https://indico.cern.ch/event/439995/contributions/1094416/attachments/1143460/1638648/diboson_final.pdf
    ("/ZZTo2Q2Nu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sXS:   4.07 * 1.22   }), # https://indico.cern.ch/event/439995/contributions/1094416/attachments/1143460/1638648/diboson_final.pdf    
    ("/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                              {sXS:   47.13  ,     sNEvtSkimv2:   7940000, sSumEvtSkimv2:  7940000}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    ("/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                              {sXS:  118.7  ,      sNEvtSkimv2:  15679000, sSumEvtSkimv2: 15679000}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV

    ("/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sXS:    0.01398  , sNEvtSkimv2:  240000, sSumEvtSkimv2: 210560}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",     {sXS:    0.01398  }),
    ("/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sXS:    0.05565  , sNEvtSkimv2:  248000, sSumEvtSkimv2: 229607}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",     {sXS:    0.05565  }),
    ("/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sXS:    0.1651  , sNEvtSkimv2:  300000, sSumEvtSkimv2: 275080}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",  {sXS:    0.1651  }),
    ("/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sXS:    0.2086  , sNEvtSkimv2:  250000, sSumEvtSkimv2: 237222}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",  {sXS:    0.2086  }),


    #("",          {sXS:    }), 
    
    ## GGF HToBB
    # dasgoclient --query="dataset=/*HToBB*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    # dasgoclient --query="dataset=/*HTo*bb*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",      {sXS:  48.61 * 0.582    }), # GGF H * BR = 28.291 N3LO
    ("/GluGluHToBB_Pt-200ToInf_M-125_TuneCP5_MINLO_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:   0.2740   , sNEvtSkimv2:  497000, sSumEvtSkimv2: 462530}),
    
    ("/VBFHToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                      {sXS:  3.766 * 0.582  }), # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
    ("/VBFHToBB_M-125_TuneCH3_13TeV-powheg-herwig/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                       {sXS:  3.766 * 0.582  }), # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
    ("/VBFHToBB_M-125_dipoleRecoilOn_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",       {sXS:  2.250  }), # HIG-21-020
    ("/VBFWH_HToBB_WToLNu_M-125_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sXS:  0.4378  }), # https://xsdb-temp.app.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=DAS%3DVBFWH_HToBB_WToLNu_M-125_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8
    ("/VBFHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {sXS:  1  }),

    ("/WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",            {sXS:  0.831 * 0.582 * 0.647    }),
    ("/WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",           {sXS:  0.831 * 0.582 * 3*0.1086 , sNEvtSkimv2:  4853959, sSumEvtSkimv2: 4562550}),
    ("/WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",           {sXS:  0.527 * 0.582 * 0.647    }),
    ("/WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sXS:  0.527 * 0.582 * 3*0.1086 , sNEvtSkimv2:  4870968, sSumEvtSkimv2: 4870968}),
    ("/WHToMuMuG_M125_Dalitz_012j_TuneCP5_13TeV_amcatnloFXFX_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {sXS:  1    }),
    ("/WplusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                 {sXS:  1    }),
    ("/WminusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                {sXS:  1    }),    

    ("/ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                {sXS:  0.758 * 0.582 * 0.699     }),
    ("/ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                {sXS:  0.758 * 0.582 * 3*0.0337 , sNEvtSkimv2: 4885835, sSumEvtSkimv2: 4592693 }),
    ("/ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",              {sXS:  0.758 * 0.582 * 0.20      }),
    ("/ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                {sXS:  0.758 * 0.582 * 0.1512   , sNEvtSkimv2: 4776217, sSumEvtSkimv2: 4494803}),
    ("/ZHToMuMuG_M125_Dalitz_012j_TuneCP5_13TeV_amcatnloFXFX_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {sXS:  1    }),
    ("/ZHToTauTau_M125_CP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                         {sXS:  1    }),    

    ("/ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",              {sXS:  0.123 * 0.582 * 0.699     }),
    ("/ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",              {sXS:  0.123 * 0.582 * 3*0.0337  }),
    ("/ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",            {sXS:  0.123 * 0.582 * 0.20      }),
    ("/ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",              {sXS:  0.123 * 0.582 * 0.1512    }),

    ("/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                        {sXS:  0.5071 * 0.5824  , sNEvtSkimv2: 9599000, sSumEvtSkimv2: 9399658}),
    ("/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                     {sXS:  0.5071 * (1 - 0.5824)  }),
    ("/ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                    {sXS:  1  }),

    ("/bbHToBB_M-125_4FS_yt2_TuneCP5-13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",             {sXS:  1  ,               sNEvtSkimv2: 9558000, sSumEvtSkimv2: 3200428}),
    

    # , sNEvtSkimv2: , sSumEvtSkimv2: 
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
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 , sNEvtSkimv2:  468175, sSumEvtSkimv2: 468175}), # filter efficiency 0.05
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 , sNEvtSkimv2:  468518, sSumEvtSkimv2: 468518}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 , sNEvtSkimv2:  452529, sSumEvtSkimv2: 452529}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 , sNEvtSkimv2:  503270, sSumEvtSkimv2: 503270}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 , sNEvtSkimv2:  464480, sSumEvtSkimv2: 464480}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 , sNEvtSkimv2:  463836, sSumEvtSkimv2: 463836}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 , sNEvtSkimv2:  480773, sSumEvtSkimv2: 480773}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 , sNEvtSkimv2:  512854, sSumEvtSkimv2: 512854}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 , sNEvtSkimv2:  460142, sSumEvtSkimv2: 460142}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 , sNEvtSkimv2:  477681, sSumEvtSkimv2: 477681}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 48.61 * 0.057 , sNEvtSkimv2:  460168, sSumEvtSkimv2: 460168}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-11.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 344458, sSumEvtSkimv2: 344458}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-11.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 320198, sSumEvtSkimv2: 320198}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-12.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 319402, sSumEvtSkimv2: 319402}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-13.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 320347, sSumEvtSkimv2: 320347}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-13.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 333902, sSumEvtSkimv2: 333902}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-14.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 322135, sSumEvtSkimv2: 322135}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-16.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 310604, sSumEvtSkimv2: 310604}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-17.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 316949, sSumEvtSkimv2: 316949}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 321407, sSumEvtSkimv2: 321407}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-21.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 313342, sSumEvtSkimv2: 313342}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-23.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 350144, sSumEvtSkimv2: 350144}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-27.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 312697, sSumEvtSkimv2: 312697}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-32.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 335356, sSumEvtSkimv2: 335356}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-37.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 333553, sSumEvtSkimv2: 333553}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-42.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 304082, sSumEvtSkimv2: 304082}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-47.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 339873, sSumEvtSkimv2: 339873}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-52.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 331319, sSumEvtSkimv2: 331319}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-57.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 314106, sSumEvtSkimv2: 314106}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-62.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 317541, sSumEvtSkimv2: 317541}),


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
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 , sNEvtSkimv2:  194401, sSumEvtSkimv2: 194401}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 , sNEvtSkimv2:  195823, sSumEvtSkimv2: 195823}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 , sNEvtSkimv2:  197952, sSumEvtSkimv2: 197952}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 , sNEvtSkimv2:  199104, sSumEvtSkimv2: 199104}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 , sNEvtSkimv2:  183328, sSumEvtSkimv2: 183328}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 , sNEvtSkimv2:  192347, sSumEvtSkimv2: 192347}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 , sNEvtSkimv2:  186254, sSumEvtSkimv2: 186254}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 , sNEvtSkimv2:  191943, sSumEvtSkimv2: 191943}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 , sNEvtSkimv2:  207777, sSumEvtSkimv2: 207777}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS: 3.766 * 0.173 , sNEvtSkimv2:  202133, sSumEvtSkimv2: 202133}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 3.766 * 0.173 , sNEvtSkimv2:  203354, sSumEvtSkimv2: 203354}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-11.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 284704, sSumEvtSkimv2: 284704}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-11.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 224593, sSumEvtSkimv2: 224593}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-12.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 286966, sSumEvtSkimv2: 286966}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-13.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 284253, sSumEvtSkimv2: 284253}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-13.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 272918, sSumEvtSkimv2: 272918}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-14.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 283886, sSumEvtSkimv2: 283886}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-16.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 284520, sSumEvtSkimv2: 284520}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-17.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 293736, sSumEvtSkimv2: 293736}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 279919, sSumEvtSkimv2: 279919}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-21.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 289270, sSumEvtSkimv2: 289270}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-23.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 278222, sSumEvtSkimv2: 278222}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-27.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 262841, sSumEvtSkimv2: 262841}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-32.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 225766, sSumEvtSkimv2: 225766}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-37.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 225349, sSumEvtSkimv2: 225349}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-42.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 226984, sSumEvtSkimv2: 226984}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-47.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 227532, sSumEvtSkimv2: 227532}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-52.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 241612, sSumEvtSkimv2: 241612}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-57.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 228574, sSumEvtSkimv2: 228574}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-62.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 240552, sSumEvtSkimv2: 240552}),


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
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 , sNEvtSkimv2:  96147, sSumEvtSkimv2: 96147}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 , sNEvtSkimv2:  103447, sSumEvtSkimv2: 103447}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 , sNEvtSkimv2:  86681, sSumEvtSkimv2: 86681}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 , sNEvtSkimv2:  100712, sSumEvtSkimv2: 100712}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 , sNEvtSkimv2:  104740, sSumEvtSkimv2: 104740}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 , sNEvtSkimv2:  100347, sSumEvtSkimv2: 100347}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 , sNEvtSkimv2:  101174, sSumEvtSkimv2: 101174}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 , sNEvtSkimv2:  101364, sSumEvtSkimv2: 101364}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 , sNEvtSkimv2:  100499, sSumEvtSkimv2: 100499}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 , sNEvtSkimv2:  104342, sSumEvtSkimv2: 104342}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  1.358 * 0.132 , sNEvtSkimv2:  104636, sSumEvtSkimv2: 104636}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-11.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 142539, sSumEvtSkimv2: 142539}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-11.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 133026, sSumEvtSkimv2: 133026}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-12.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 130865, sSumEvtSkimv2: 130865}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-13.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 124564, sSumEvtSkimv2: 124564}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-13.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 121600, sSumEvtSkimv2: 121600}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-14.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 144090, sSumEvtSkimv2: 144090}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-16.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 112581, sSumEvtSkimv2: 112581}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-17.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 149426, sSumEvtSkimv2: 149426}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 123524, sSumEvtSkimv2: 123524}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-21.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 145756, sSumEvtSkimv2: 145756}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-23.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 139194, sSumEvtSkimv2: 139194}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-27.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 148004, sSumEvtSkimv2: 148004}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-32.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 112694, sSumEvtSkimv2: 112694}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-37.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 148351, sSumEvtSkimv2: 148351}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-42.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 115271, sSumEvtSkimv2: 115271}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-47.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 116578, sSumEvtSkimv2: 116578}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-52.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 123945, sSumEvtSkimv2: 123945}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-57.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 118857, sSumEvtSkimv2: 118857}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-62.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 130501, sSumEvtSkimv2: 130501}),


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
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129 , sNEvtSkimv2:  101633, sSumEvtSkimv2:  101633}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129 , sNEvtSkimv2:  106743, sSumEvtSkimv2:  106743}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129 , sNEvtSkimv2:  97716, sSumEvtSkimv2:  97716}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129 , sNEvtSkimv2:  100766, sSumEvtSkimv2:  100766}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129 , sNEvtSkimv2:  101445, sSumEvtSkimv2:  101445}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129 , sNEvtSkimv2:  102774, sSumEvtSkimv2:  102774}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129 , sNEvtSkimv2:  100835, sSumEvtSkimv2:  100835}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129 , sNEvtSkimv2:  102724, sSumEvtSkimv2:  102724}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129 , sNEvtSkimv2:  99612, sSumEvtSkimv2:  99612}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129 , sNEvtSkimv2:  95365, sSumEvtSkimv2:  95365}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS:  0.880 * 0.129 , sNEvtSkimv2:  100569, sSumEvtSkimv2:  100569}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-11.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 126474, sSumEvtSkimv2: 126474}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-11.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 126889, sSumEvtSkimv2: 126889}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-12.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 122989, sSumEvtSkimv2: 122989}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-13.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 126016, sSumEvtSkimv2: 126016}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-13.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 128145, sSumEvtSkimv2: 128145}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-14.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 126910, sSumEvtSkimv2: 126910}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-16.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 122854, sSumEvtSkimv2: 122854}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-17.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 124177, sSumEvtSkimv2: 124177}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 126788, sSumEvtSkimv2: 126788}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-21.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 130247, sSumEvtSkimv2: 130247}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-23.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 133181, sSumEvtSkimv2: 133181}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-27.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 132964, sSumEvtSkimv2: 132964}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-32.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 133549, sSumEvtSkimv2: 133549}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-37.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 131627, sSumEvtSkimv2: 131627}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-42.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 129405, sSumEvtSkimv2: 129405}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-47.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 132768, sSumEvtSkimv2: 132768}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-52.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 132143, sSumEvtSkimv2: 132143}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-57.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 125970, sSumEvtSkimv2: 125970}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-62.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 129540, sSumEvtSkimv2: 129540}),


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
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285 , sNEvtSkimv2:  199812, sSumEvtSkimv2:  199812 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285 , sNEvtSkimv2:  188898, sSumEvtSkimv2:  188898}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285 , sNEvtSkimv2:  188884, sSumEvtSkimv2:  188884 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285 , sNEvtSkimv2:  200687, sSumEvtSkimv2:  200687 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285 , sNEvtSkimv2:  196711, sSumEvtSkimv2:  196711 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285 , sNEvtSkimv2:  197161, sSumEvtSkimv2:  197161 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285 , sNEvtSkimv2:  199334, sSumEvtSkimv2:  199334 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285 , sNEvtSkimv2:  191927, sSumEvtSkimv2:  191927 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285 , sNEvtSkimv2:  194429, sSumEvtSkimv2:  194429 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285 , sNEvtSkimv2:  198645, sSumEvtSkimv2:  198645 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 0.5071 * 0.285 , sNEvtSkimv2:  203990, sSumEvtSkimv2:  203990 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-11.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 237248, sSumEvtSkimv2: 237248}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-11.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 245178, sSumEvtSkimv2: 245178}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-12.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 244365, sSumEvtSkimv2: 244365}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-13.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 242777, sSumEvtSkimv2: 242777}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-13.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 237450, sSumEvtSkimv2: 237450}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-14.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 237302, sSumEvtSkimv2: 237302}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-16.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 235117, sSumEvtSkimv2: 235117}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-17.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 233498, sSumEvtSkimv2: 233498}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 242612, sSumEvtSkimv2: 242612}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-21.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 235910, sSumEvtSkimv2: 235910}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-23.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 233671, sSumEvtSkimv2: 233671}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-27.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 231705, sSumEvtSkimv2: 231705}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-32.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 231812, sSumEvtSkimv2: 231812}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-37.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 229240, sSumEvtSkimv2: 229240}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-42.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 238661, sSumEvtSkimv2: 238661}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-47.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 239480, sSumEvtSkimv2: 239480}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-52.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 239845, sSumEvtSkimv2: 239845}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-57.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 237874, sSumEvtSkimv2: 237874}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-62.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 235987, sSumEvtSkimv2: 235987}),


    ## H->aa->4Tau
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    
    ("/SUSYVHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYVHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),

    ("/SUSYttHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),
    ("/SUSYttHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1   }),

    

    
    


    

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