from collections import OrderedDict as OD

#sXS           = "xs"
sNameSp       = "nameSp"
sNanoSkimv2   = "skim_v2"
sNEvents      = "nEvents"
sSumEvents    = "sumEvents"
sNEvtSkimv2   = "%s_%s" % (sNanoSkimv2, sNEvents)
sSumEvtSkimv2 = "%s_%s" % (sNanoSkimv2, sSumEvents)

list_datasets_2018 = OD([
    
    ## QCD_bEnriched_HT*
    # dasgoclient --query="dataset=/QCD_bEnriched_HT*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"    
    ("/QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 36118282, sSumEvtSkimv2: 36118282}),
    ("/QCD_bEnriched_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 18462183, sSumEvtSkimv2: 18462183}),
    ("/QCD_bEnriched_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 11197722, sSumEvtSkimv2: 11197722}),
    ("/QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2:  9246898, sSumEvtSkimv2:  9246898}),    
    ("/QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sNEvtSkimv2:  1844165, sSumEvtSkimv2:  1844165}),
    ("/QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  1330829, sSumEvtSkimv2:  1330829}),    
    ("/QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 1431254, sSumEvtSkimv2:  1431254}),
    ("/QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sNEvtSkimv2: 1429280, sSumEvtSkimv2:  1429280}),
    

    ## QCD_HT*_BGenFilter
    # dasgoclient --query="dataset=/QCD_HT*_BGenFilter*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/QCD_HT100to200_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sNEvtSkimv2: 35884072, sSumEvtSkimv2: 35884072}),
    ("/QCD_HT200to300_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sNEvtSkimv2: 14830551, sSumEvtSkimv2: 14830551}),
    ("/QCD_HT300to500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sNEvtSkimv2: 14144826, sSumEvtSkimv2: 14144826}),
    ("/QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",   {sNEvtSkimv2:  8015859, sSumEvtSkimv2:  8015859}),
    ("/QCD_HT700to1000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",  {sNEvtSkimv2:  4642245, sSumEvtSkimv2:  4642245}),    
    ("/QCD_HT1000to1500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sNEvtSkimv2:  1537452, sSumEvtSkimv2:  1537452}),
    ("/QCD_HT1500to2000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sNEvtSkimv2:  1263157, sSumEvtSkimv2:  1263157}),
    ("/QCD_HT2000toInf_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",  {sNEvtSkimv2:  1300672, sSumEvtSkimv2:  1300672}),


    ## QCD_HT* PSWeights-madgraph - QCDIncl LO
    # dasgoclient --query="dataset=/QCD_HT*PSWeight*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {}),
    ("/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 84461486, sSumEvtSkimv2: 84461486}),
    ("/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 57336623, sSumEvtSkimv2: 57336623}),
    ("/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 61705174, sSumEvtSkimv2: 61705174}),
    ("/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 49184771, sSumEvtSkimv2: 49184771}),
    ("/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",  {sNEvtSkimv2: 48506751, sSumEvtSkimv2: 48506751}),
    ("/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 13026372, sSumEvtSkimv2: 13026372}), ### Correction needed
    ("/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 10871473, sSumEvtSkimv2: 10871473}),
    ("/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sNEvtSkimv2:  5374711, sSumEvtSkimv2:  5374711}),

    
    ## QCD_HT* TuneCP5 madgraphMLM - QCDIncl LO, MatrixElement-PartonShower maching at NLO(?)
    # dasgoclient --query="dataset=/QCD_HT*TuneCP5_13TeV-madgraphMLM*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {}),
    ("/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {}),
    ("/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {}),
    ("/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {}),
    ("/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {}),    
    ("/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {}),
    ("/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {}),

    

    ## TTbar - NLO powheg
    # dasgoclient --query="dataset=/TT*_TuneCP5_13TeV*powheg*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {sNEvtSkimv2: 343248000, sSumEvtSkimv2: 338171646}), # 831.8 * 0.457
    ("/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {sNEvtSkimv2: 478982000, sSumEvtSkimv2: 471724160}), # 831.8 * 0.438
    ("/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                      {sNEvtSkimv2: 146010000, sSumEvtSkimv2: 143887383}), # 831.8 * 0.105

    ## TTbar Jets - NLO amcatnlo
    # dasgoclient --query="dataset=/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {}),

    ## TTbar Jets - LO madgraph
    # dasgoclient --query="dataset=/TTJets_*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM" 
    ("/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                    {}),
    ("/TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {}),
    ("/TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {}),
    ("/TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",      {}),
    ("/TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {}),

    ## TTbar Jets To LNu - LO madgraph
    # dasgoclient --query="dataset=/TTJets_*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {}), # 831.8 * 0.219
    ("/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}), # 831.8 * 0.219
    ("/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",             {}), # 831.8 * 0.105


    ## ST NLO
    # dasgoclient --query="dataset=/ST*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",     {sNEvtSkimv2: 198569000, sSumEvtSkimv2: 197426478}),
    ("/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  99264000, sSumEvtSkimv2:  98724306}),
    ("/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",            {sNEvtSkimv2:   7956000, sSumEvtSkimv2:   7953729}), # 79.3/2 = 39.65
    ("/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",        {sNEvtSkimv2:   7749000, sSumEvtSkimv2:   7747316}), # 79.3/2 = 39.65
    ("/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",      {sNEvtSkimv2:  11270430, sSumEvtSkimv2:  11269836}),  # 79.3/2*(1 - (0.6741*0.6741)) : BR(W->Hadrons)=0.6741. tW has W(from t) and W. BR(tW->NoFully hadronic)=(1 - (0.6741*0.6741))=0.546.
    ("/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sNEvtSkimv2:  11015956, sSumEvtSkimv2:  11012326}),  # 79.3/2*(1 - (0.6741*0.6741)) : Andrew: 21.61 
    ("/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sNEvtSkimv2:  16391000, sSumEvtSkimv2:  10793026}),
    ("/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sNEvtSkimv2:  19365999, sSumEvtSkimv2:  12607741}),

    ## tX            
    # dasgoclient --query="dataset=/t*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2:  11916000, sSumEvtSkimv2:  3229640}),
    ("/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",          {sNEvtSkimv2:  32793815, sSumEvtSkimv2: 32793815}), # https://cds.cern.ch/record/2227475/files/CERN-2017-002-M.pdf?version=1#page=180
    ("/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",          {sNEvtSkimv2:  27686862, sSumEvtSkimv2: 27686862}), # https://cds.cern.ch/record/2227475/files/CERN-2017-002-M.pdf?version=1#page=180
    

    ## ZJets
    # dasgoclient --query="dataset=/ZJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"    
    ("/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sNEvtSkimv2: 14808858, sSumEvtSkimv2: 14808858}),
    ("/ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sNEvtSkimv2: 13930474, sSumEvtSkimv2: 13930474}),
    ("/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sNEvtSkimv2: 12029507, sSumEvtSkimv2: 11668622}), ### Correction needed
    ("/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sNEvtSkimv2:  9681521, sSumEvtSkimv2:  9003815}), ### Correction needed

    ## DYJetsToLL NLO
    # dasgoclient --query="dataset=/DY*Jets*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {}),
    ("/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",      {}),
    
    # DYJetsToLL LO Incl
    ("/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                     {sNEvtSkimv2: 99288125, sSumEvtSkimv2: 99288125}),
    ("/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                         {sNEvtSkimv2: 96233328, sSumEvtSkimv2: 96233328}),
    #("/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v1/NANOAODSIM",                    {}),

    # DYJetsToLL LO HT
    ("/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {sNEvtSkimv2: 17004433, sSumEvtSkimv2: 17004433}),
    ("/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 26202328, sSumEvtSkimv2: 26202328}),
    ("/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 18455718, sSumEvtSkimv2: 18455718}),
    ("/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2:  8908406, sSumEvtSkimv2:  8908406}),
    ("/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2:  7035971, sSumEvtSkimv2:  7035971}),
    ("/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sNEvtSkimv2:  6678036, sSumEvtSkimv2:  6678036}),
    ("/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  6166852, sSumEvtSkimv2:  6166852}),
    ("/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sNEvtSkimv2:  1978203, sSumEvtSkimv2:  1978203}),    

    # dasgoclient --query="dataset=/ZJetsToNuNu*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 28876062, sSumEvtSkimv2: 28876062}),  
    ("/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 22749608, sSumEvtSkimv2: 22749608}),  
    ("/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 19810491, sSumEvtSkimv2: 19810491}),  
    ("/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",   {sNEvtSkimv2: 5067605, sSumEvtSkimv2: 5067605}),  
    ("/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sNEvtSkimv2: 1948816, sSumEvtSkimv2: 1948816}),  
    ("/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2: 381695, sSumEvtSkimv2: 381695}),  
    ("/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",  {sNEvtSkimv2: 268224, sSumEvtSkimv2: 268224}),  
    
    ## WJetsToQQ
    # dasgoclient --query="dataset=/WJetsToQQ*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sNEvtSkimv2: 14494966, sSumEvtSkimv2: 14494966}),
    ("/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sNEvtSkimv2:  9335298, sSumEvtSkimv2:  9335298}),
    ("/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {}),
    ("/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {}),


    ## WJetsToLNu
    # dasgoclient --query="dataset=/W*Jets*ToLNu*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                   {sNEvtSkimv2: 29117828, sSumEvtSkimv2: 12361080}), # NLO sample
    #("/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {}), # NLO sample
    #("/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {}), # NLO sample
    #("/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",               {}), # NLO sample    
    #("/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                    {}),
    ("/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",         {sNEvtSkimv2:  66569448, sSumEvtSkimv2: 66569448}),
    ("/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sNEvtSkimv2:  51541593, sSumEvtSkimv2: 51541593}),
    ("/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sNEvtSkimv2:  58331446, sSumEvtSkimv2: 58331446}),
    ("/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sNEvtSkimv2:   7444030, sSumEvtSkimv2:  7444030}),
    #("/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",   {}),
    ("/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",        {sNEvtSkimv2:   7718765, sSumEvtSkimv2:  7718765}),
    #("/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",   {}),
    ("/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sNEvtSkimv2:   7306187, sSumEvtSkimv2:  7296689}), ### COrrection needed
    #("/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",  {}),
    ("/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",      {sNEvtSkimv2:   6481518, sSumEvtSkimv2:  6481518}),
    #("/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM", {}),
    ("/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",       {sNEvtSkimv2:   2097648, sSumEvtSkimv2:  2097648}),
    #("/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {}),
    #("/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {}),
    #("/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {}),
    #("/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                   {}),


    ## VV, VVV
    # dasgoclient --query="dataset=/Z*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    # dasgoclient --query="dataset=/W*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                              {sNEvtSkimv2:   3526000, sSumEvtSkimv2:  3526000}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    ("/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  29357938, sSumEvtSkimv2: 15840251}), # https://indico.cern.ch/event/439995/contributions/1094416/attachments/1143460/1638648/diboson_final.pdf
    ("/ZZTo2Q2Nu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {}), # https://indico.cern.ch/event/439995/contributions/1094416/attachments/1143460/1638648/diboson_final.pdf    
    ("/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                              {sNEvtSkimv2:   7940000, sSumEvtSkimv2:  7940000}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    ("/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                              {sNEvtSkimv2:  15679000, sSumEvtSkimv2: 15679000}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV

    ("/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sNEvtSkimv2:  240000, sSumEvtSkimv2: 210560}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",     {}),
    ("/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sNEvtSkimv2:  248000, sSumEvtSkimv2: 229607}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",     {}),
    ("/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sNEvtSkimv2:  300000, sSumEvtSkimv2: 275080}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",  {}),
    ("/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {sNEvtSkimv2:  250000, sSumEvtSkimv2: 237222}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    #("/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2/NANOAODSIM",  {}),


    #("",          {   }), 
    
    ## GGF HToBB
    # dasgoclient --query="dataset=/*HToBB*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    # dasgoclient --query="dataset=/*HTo*bb*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",      {}), # GGF H * BR = 28.291 N3LO
    ("/GluGluHToBB_Pt-200ToInf_M-125_TuneCP5_MINLO_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  497000, sSumEvtSkimv2: 462530}),
    
    ("/VBFHToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                      {}), # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
    ("/VBFHToBB_M-125_TuneCH3_13TeV-powheg-herwig/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                       {}), # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
    ("/VBFHToBB_M-125_dipoleRecoilOn_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",       {}), # HIG-21-020
    ("/VBFWH_HToBB_WToLNu_M-125_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",       {}), # https://xsdb-temp.app.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=DAS%3DVBFWH_HToBB_WToLNu_M-125_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8
    ("/VBFHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                   {}),

    ("/WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",            {}),
    ("/WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",           {sNEvtSkimv2:  4853959, sSumEvtSkimv2: 4562550}),
    ("/WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",           {}),
    ("/WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",          {sNEvtSkimv2:  4870968, sSumEvtSkimv2: 4870968}),
    ("/WHToMuMuG_M125_Dalitz_012j_TuneCP5_13TeV_amcatnloFXFX_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {}),
    ("/WplusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                 {}),
    ("/WminusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                {}),    

    ("/ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                {}),
    ("/ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                {sNEvtSkimv2: 4885835, sSumEvtSkimv2: 4592693 }),
    ("/ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",              {}),
    ("/ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                {sNEvtSkimv2: 4776217, sSumEvtSkimv2: 4494803}),
    ("/ZHToMuMuG_M125_Dalitz_012j_TuneCP5_13TeV_amcatnloFXFX_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",    {}),
    ("/ZHToTauTau_M125_CP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                         {}),    

    ("/ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",              {}),
    ("/ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",              {}),
    ("/ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",            {}),
    ("/ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",              {}),

    ("/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                        {sNEvtSkimv2: 9599000, sSumEvtSkimv2: 9399658}),
    ("/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM",                     {}),
    ("/ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",                    {}),

    ("/bbHToBB_M-125_4FS_yt2_TuneCP5-13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",             {sNEvtSkimv2: 9558000, sSumEvtSkimv2: 3200428}),
    

    # , sNEvtSkimv2: , sSumEvtSkimv2: 
    ## SUSY_GluGluH_01J_HToAATo4B_M-*   and   SUSY_GluGluH_01J_HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*GluGluH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/SUSY_GluGluH_01J_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),

    # SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (308) / (5379) = 5.726e-02 +- 3.168e-03	Matching efficiency = 0.5 +/- 0.0     Cross-section = 48.61 pb * 0.057	
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  468175, sSumEvtSkimv2: 468175}), # filter efficiency 0.05
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  468518, sSumEvtSkimv2: 468518}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  452529, sSumEvtSkimv2: 452529}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  503270, sSumEvtSkimv2: 503270}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  464480, sSumEvtSkimv2: 464480}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  463836, sSumEvtSkimv2: 463836}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  480773, sSumEvtSkimv2: 480773}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  512854, sSumEvtSkimv2: 512854}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  460142, sSumEvtSkimv2: 460142}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  477681, sSumEvtSkimv2: 477681}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  460168, sSumEvtSkimv2: 460168}),


    ## VBF HToAATo4B_M-* and VBF HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*VBF*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/SUSY_VBFH_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),

    # SUSY_VBFH_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8		Filter efficiency (event-level)= (1730) / (10000) = 1.730e-01 +- 3.782e-03	Matching efficiency = 1.0 +/- 0.0 	Cross-section = 3.766 pb * 0.173 
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  194401, sSumEvtSkimv2: 194401}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  195823, sSumEvtSkimv2: 195823}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  197952, sSumEvtSkimv2: 197952}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  199104, sSumEvtSkimv2: 199104}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  183328, sSumEvtSkimv2: 183328}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  192347, sSumEvtSkimv2: 192347}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  186254, sSumEvtSkimv2: 186254}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  191943, sSumEvtSkimv2: 191943}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  207777, sSumEvtSkimv2: 207777}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sNEvtSkimv2:  202133, sSumEvtSkimv2: 202133}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  203354, sSumEvtSkimv2: 203354}),


    ## WH HToAATo4B-M-* and WH HToAATo4B_Pt150_-M-*
    # dasgoclient --query="dataset=/SUSY*WH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM" 
    ("/SUSY_WH_WToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),

    # SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (1317) / (10000) = 1.317e-01 +- 3.382e-03	Matching efficiency = 1.0 +/- 0.0	Cross-section = 1.358 pb * 0.132 
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  96147, sSumEvtSkimv2: 96147}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  103447, sSumEvtSkimv2: 103447}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  86681, sSumEvtSkimv2: 86681}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  100712, sSumEvtSkimv2: 100712}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  104740, sSumEvtSkimv2: 104740}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  100347, sSumEvtSkimv2: 100347}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  101174, sSumEvtSkimv2: 101174}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  101364, sSumEvtSkimv2: 101364}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  100499, sSumEvtSkimv2: 100499}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  104342, sSumEvtSkimv2: 104342}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  104636, sSumEvtSkimv2: 104636}),


    ## ZH HToAATo4B-M-* and ZH HToAATo4B_Pt150_-M-*
    # dasgoclient --query="dataset=/SUSY*ZH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    
    # SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (1289) / (10000) = 1.289e-01 +- 3.351e-03	Matching efficiency = 1.0 +/- 0.0		Cross-section = 0.880 pb * 0.129 
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  101633, sSumEvtSkimv2:  101633}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  106743, sSumEvtSkimv2:  106743}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  97716, sSumEvtSkimv2:  97716}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  100766, sSumEvtSkimv2:  100766}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  101445, sSumEvtSkimv2:  101445}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  102774, sSumEvtSkimv2:  102774}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  100835, sSumEvtSkimv2:  100835}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  102724, sSumEvtSkimv2:  102724}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  99612, sSumEvtSkimv2:  99612}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  95365, sSumEvtSkimv2:  95365}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  100569, sSumEvtSkimv2:  100569}),


    ## ttH HToAATo4B_M-* and HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*TTH*HToAATo4B*M*/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v*/NANOAODSIM"
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),

    # SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (2850) / (10000) = 2.850e-01 +- 4.514e-03	Matching efficiency = 1.0 +/- 0.0 Cross-section = 0.5071 pb * 0.285 
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  199812, sSumEvtSkimv2:  199812 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  188898, sSumEvtSkimv2:  188898}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  188884, sSumEvtSkimv2:  188884 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  200687, sSumEvtSkimv2:  200687 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  196711, sSumEvtSkimv2:  196711 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  197161, sSumEvtSkimv2:  197161 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  199334, sSumEvtSkimv2:  199334 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  191927, sSumEvtSkimv2:  191927 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  194429, sSumEvtSkimv2:  194429 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  198645, sSumEvtSkimv2:  198645 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sNEvtSkimv2:  203990, sSumEvtSkimv2:  203990 }),
    

    ## H->aa->4Tau
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    
    ("/SUSYVHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),

    ("/SUSYttHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {}),

    

    
    


    

    ## JetHT
    # dasgoclient --query="dataset=/JetHT/*2018*UL*MiniAODv2_NanoAODv9-*/NANOAOD"
    # XS (cross-section) does not matter for data sample
    # Use JetHT/NanoAOD_GT38 instead
    #("/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD", {}),
    #("/JetHT/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD", {}),
    #("/JetHT/Run2018C-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD", {}),
    #("/JetHT/Run2018D-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD", {}),
    
    # dasgoclient --query="dataset=/JetHT/*2018*UL*MiniAODv2_NanoAODv9_GT36*/NANOAOD"
    # XS (cross-section) does not matter for data sample
    ("/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), # run 315257 to 316995
    ("/JetHT/Run2018B-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), # run 317080 to 319310
    ("/JetHT/Run2018C-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), # run 319337 to 320065
    ("/JetHT/Run2018D-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), # run 320500 to 325175

    ## SingleMuon
    # dasgoclient --query="dataset=/SingleMuon/*2018*UL*MiniAODv2_NanoAODv9_GT36*/NANOAOD"
    # XS (cross-section) does not matter for data sample    
    ("/SingleMuon/Run2018A-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}),    
    ("/SingleMuon/Run2018B-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}),
    ("/SingleMuon/Run2018C-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}),
    ("/SingleMuon/Run2018D-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}),
    
    ## EGamma
    # dasgoclient --query="dataset=/EGamma/*2018*UL*MiniAODv2_NanoAODv9_GT36*/NANOAOD"
    # XS (cross-section) does not matter for data sample    
    ("/EGamma/Run2018A-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), 
    ("/EGamma/Run2018B-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), 
    ("/EGamma/Run2018C-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), 
    ("/EGamma/Run2018D-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), 

    ## MET
    # dasgoclient --query="dataset=/MET/*2018*UL*MiniAODv2_NanoAODv9_GT36*/NANOAOD"
    # XS (cross-section) does not matter for data sample    
    ("/MET/Run2018A-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), 
    ("/MET/Run2018B-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), 
    ("/MET/Run2018C-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), 
    ("/MET/Run2018D-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}), 
    

    #("", {}),    
])