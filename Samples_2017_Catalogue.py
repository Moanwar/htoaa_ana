from collections import OrderedDict as OD

#sXS           = "xs"
sNameSp       = "nameSp"
sNanoSkimv2   = "skim_v2"
sNEvents      = "nEvents"
sSumEvents    = "sumEvents"
sNEvtSkimv2   = "%s_%s" % (sNanoSkimv2, sNEvents)
sSumEvtSkimv2 = "%s_%s" % (sNanoSkimv2, sSumEvents)

list_datasets_2017 = OD([
    
    ## QCD_bEnriched_HT*
    # dasgoclient --query="dataset=/QCD_bEnriched_HT*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"    
    ("/QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 37427427, sSumEvtSkimv2: 37427427}),
    ("/QCD_bEnriched_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 19844424, sSumEvtSkimv2: 19844424}),
    ("/QCD_bEnriched_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 11312350, sSumEvtSkimv2: 11312350}),
    ("/QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2:  10203561, sSumEvtSkimv2:  10203561}),    
    ("/QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",  {sNEvtSkimv2:  1881618, sSumEvtSkimv2:  1881618}),
    ("/QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  1385631, sSumEvtSkimv2:  1385631}),    
    ("/QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 1458069, sSumEvtSkimv2:  1458069}),
    ("/QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",  {sNEvtSkimv2: 1408971, sSumEvtSkimv2:  1408971}),

    ## QCD_HT*_BGenFilter
    # dasgoclient --query="dataset=/QCD_HT*_BGenFilter*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/QCD_HT100to200_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",   {sNEvtSkimv2: 30667858, sSumEvtSkimv2: 30667858}),
    ("/QCD_HT200to300_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",   {sNEvtSkimv2: 15085476, sSumEvtSkimv2: 15085476}),
    ("/QCD_HT300to500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",   {sNEvtSkimv2: 14259151, sSumEvtSkimv2: 14259151}),
    ("/QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",   {sNEvtSkimv2:  5969789, sSumEvtSkimv2:  5969789}),
    ("/QCD_HT700to1000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",  {sNEvtSkimv2:  4892427, sSumEvtSkimv2:  4892427}),    
    ("/QCD_HT1000to1500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {sNEvtSkimv2:  1642413, sSumEvtSkimv2:  1642413}),
    ("/QCD_HT1500to2000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {sNEvtSkimv2:  1263823, sSumEvtSkimv2:  1263823}),
    ("/QCD_HT2000toInf_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",  {sNEvtSkimv2:  1293056, sSumEvtSkimv2:  1293056}),

    ## QCD_HT* PSWeights-madgraph - QCDIncl LO
    # dasgoclient --query="dataset=/QCD_HT*PSWeight*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",    {}),
    ("/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 54760426, sSumEvtSkimv2: 54760426}),
    ("/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 42714435, sSumEvtSkimv2: 42714435}),
    ("/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 43589739, sSumEvtSkimv2: 43589739}),
    ("/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 36194860, sSumEvtSkimv2: 36194860}),
    ("/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",  {sNEvtSkimv2: 34051754, sSumEvtSkimv2: 34051754}),
    ("/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 10256089, sSumEvtSkimv2: 10256089}), ### Correction needed
    ("/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  7701876, sSumEvtSkimv2:  7701876}),
    ("/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",  {sNEvtSkimv2:  4112573, sSumEvtSkimv2:  4112573}),
    

    

    ## TTbar - NLO powheg
    # dasgoclient --query="dataset=/TT*_TuneCP5_13TeV*powheg*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                   {sNEvtSkimv2: 235719999, sSumEvtSkimv2: 231910835}), # 831.8 * 0.457
    ("/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",               {sNEvtSkimv2: 355332000, sSumEvtSkimv2: 352462632}), # 831.8 * 0.438
    ("/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                      {sNEvtSkimv2: 106724000, sSumEvtSkimv2: 105859990}), # 831.8 * 0.105

    ## ST NLO
    # dasgoclient --query="dataset=/ST*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",     {sNEvtSkimv2: 142078000, sSumEvtSkimv2: 141260698}),
    ("/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {sNEvtSkimv2:  70203000, sSumEvtSkimv2:  69821018}),
    ("/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",            {sNEvtSkimv2:   5649000, sSumEvtSkimv2:   5648712}), # 79.3/2 = 39.65
    ("/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",        {sNEvtSkimv2:   5674000, sSumEvtSkimv2:   5673700}), # 79.3/2 = 39.65
    ("/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",      {sNEvtSkimv2:   8507203, sSumEvtSkimv2:   8506765}),  # 79.3/2*(1 - (0.6741*0.6741)) : BR(W->Hadrons)=0.6741. tW has W(from t) and W. BR(tW->NoFully hadronic)=(1 - (0.6741*0.6741))=0.546.
    ("/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",  {sNEvtSkimv2:   8433998, sSumEvtSkimv2:   8433562}),  # 79.3/2*(1 - (0.6741*0.6741)) : Andrew: 21.61 
    ("/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",        {sNEvtSkimv2:  11696999, sSumEvtSkimv2:   7615207}),
    ("/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",          {sNEvtSkimv2:  13882000, sSumEvtSkimv2:   9037288}),

    ## tX            
    # dasgoclient --query="dataset=/t*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2:  9530000, sSumEvtSkimv2:  2578040}),
    ("/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",          {sNEvtSkimv2:  32503711, sSumEvtSkimv2: 32503711}), # https://cds.cern.ch/record/2227475/files/CERN-2017-002-M.pdf?version=1#page=180
    ("/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",          {sNEvtSkimv2:  27662138, sSumEvtSkimv2: 27662138}), # https://cds.cern.ch/record/2227475/files/CERN-2017-002-M.pdf?version=1#page=180
    
    ## ZJets
    # dasgoclient --query="dataset=/ZJetsToQQ*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"    
    ("/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {sNEvtSkimv2: 15256299, sSumEvtSkimv2: 15256299}),
    ("/ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {sNEvtSkimv2: 14884962, sSumEvtSkimv2: 14884962}),
    ("/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {sNEvtSkimv2: 11702567, sSumEvtSkimv2: 11702567}), ### Correction needed
    ("/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {sNEvtSkimv2:  9384525, sSumEvtSkimv2:  9384525}), ### Correction needed
    
    # DYJetsToLL LO Incl
    # dasgoclient --query="dataset=/DY*Jets*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                     {sNEvtSkimv2:  68480179, sSumEvtSkimv2:  68480179}),
    ("/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                         {sNEvtSkimv2: 103344974, sSumEvtSkimv2: 103344974}),

    # DYJetsToLL LO HT
    ("/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",    {sNEvtSkimv2: 12205958, sSumEvtSkimv2: 12205958}),
    ("/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 18955253, sSumEvtSkimv2: 18955253}),
    ("/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 12513057, sSumEvtSkimv2: 12513057}),
    ("/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2:  5543804, sSumEvtSkimv2:  5543804}),
    ("/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2:  5278417, sSumEvtSkimv2:  5278417}),
    ("/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",  {sNEvtSkimv2:  4506887, sSumEvtSkimv2:  4506887}),
    ("/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  4802716, sSumEvtSkimv2:  4802716}),
    ("/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",  {sNEvtSkimv2:  1480047, sSumEvtSkimv2:  1480047}),    

    # dasgoclient --query="dataset=/ZJetsToNuNu*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 18983897, sSumEvtSkimv2: 18983897}),  
    ("/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 17349597, sSumEvtSkimv2: 17349597}),  
    ("/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2: 13963690, sSumEvtSkimv2: 13963690}),  
    ("/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",   {sNEvtSkimv2:  4418971, sSumEvtSkimv2:  4418971}),  
    ("/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",  {sNEvtSkimv2:  1513585, sSumEvtSkimv2:  1513585}),  
    ("/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:   267125, sSumEvtSkimv2:   267125}),  
    ("/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",  {sNEvtSkimv2:   176201, sSumEvtSkimv2:   176201}),  

    ## WJetsToQQ
    # dasgoclient --query="dataset=/WJetsToQQ*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {sNEvtSkimv2: 15968057, sSumEvtSkimv2: 15968057}),
    ("/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {sNEvtSkimv2:  9927793, sSumEvtSkimv2:  9927793}),
    ("/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {sNEvtSkimv2: 14667933, sSumEvtSkimv2: 14667933}),
    ("/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {sNEvtSkimv2: 14722417, sSumEvtSkimv2: 14722417}),

    ## WJetsToLNu
    # dasgoclient --query="dataset=/W*Jets*ToLNu*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",                   {sNEvtSkimv2:  26838734, sSumEvtSkimv2: 18528826}), # NLO sample
    ("/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                    {sNEvtSkimv2:  78981243, sSumEvtSkimv2: 78981243}),
    ("/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",         {sNEvtSkimv2:  44736228, sSumEvtSkimv2: 44736228}),
    ("/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",        {sNEvtSkimv2:  47424468, sSumEvtSkimv2: 47424468}),
    ("/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",        {sNEvtSkimv2:  42602407, sSumEvtSkimv2: 42602407}),
    ("/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",        {sNEvtSkimv2:   5468473, sSumEvtSkimv2:  5468473}),
    ("/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",        {sNEvtSkimv2:   5545298, sSumEvtSkimv2:  5545298}),
    ("/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v3/NANOAODSIM",       {sNEvtSkimv2:   5088483, sSumEvtSkimv2:  5088483}), ### COrrection needed
    ("/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",      {sNEvtSkimv2:   4955636, sSumEvtSkimv2:  4955636}),
    ("/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",       {sNEvtSkimv2:   1185699, sSumEvtSkimv2:  1185699}),

    ## VV, VVV
    # dasgoclient --query="dataset=/Z*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    # dasgoclient --query="dataset=/W*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                              {sNEvtSkimv2:   2706000, sSumEvtSkimv2:  2706000}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    ("/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  29521496, sSumEvtSkimv2: 19134900}), # https://indico.cern.ch/event/439995/contributions/1094416/attachments/1143460/1638648/diboson_final.pdf
    #("/ZZTo2Q2Nu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",          {}), # https://indico.cern.ch/event/439995/contributions/1094416/attachments/1143460/1638648/diboson_final.pdf    
    ("/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                              {sNEvtSkimv2:   7889000, sSumEvtSkimv2:  7889000}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    ("/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                              {sNEvtSkimv2:  15634000, sSumEvtSkimv2: 15634000}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV

    ("/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",          {sNEvtSkimv2: 178000, sSumEvtSkimv2: 156384}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    ("/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",          {sNEvtSkimv2: 298000, sSumEvtSkimv2: 271070}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    ("/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",       {sNEvtSkimv2: 178000, sSumEvtSkimv2: 162156}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    ("/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",       {sNEvtSkimv2: 171000, sSumEvtSkimv2: 156692}), # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Triboson
    

    #("",          {   }), sNEvtSkimv2: , sSumEvtSkimv2: 
    
    ## GGF HToBB
    # dasgoclient --query="dataset=/*HToBB*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    # dasgoclient --query="dataset=/*HTo*bb*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",      {}), # GGF H * BR = 28.291 N3LO
    ("/GluGluHToBB_Pt-200ToInf_M-125_TuneCP5_MINLO_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  497000, sSumEvtSkimv2: 496728}),
    
    ("/VBFHToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                      {sNEvtSkimv2: 7352000, sSumEvtSkimv2: 7339234}), # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
    ("/VBFHToBB_M-125_TuneCH3_13TeV-powheg-herwig/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                       {}), # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNHLHE2019
    ("/VBFHToBB_M-125_dipoleRecoilOn_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",       {sNEvtSkimv2: 7692000, sSumEvtSkimv2: 7679477}), # HIG-21-020
    ("/VBFWH_HToBB_WToLNu_M-125_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",       {}), # https://xsdb-temp.app.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=DAS%3DVBFWH_HToBB_WToLNu_M-125_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8
    ("/VBFHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                   {}),

    ("/WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",            {sNEvtSkimv2: 7940019, sSumEvtSkimv2: 7683452}),
    ("/WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",           {sNEvtSkimv2:  4862115, sSumEvtSkimv2: 4862115}),
    ("/WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",           {sNEvtSkimv2:  7951919, sSumEvtSkimv2: 7684812}),
    ("/WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",          {sNEvtSkimv2:  4824508, sSumEvtSkimv2: 4824508}),
    ("/WHToMuMuG_M125_Dalitz_012j_TuneCP5_13TeV_amcatnloFXFX_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",    {}),
    ("/WplusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                 {}),
    ("/WminusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                {}),    

    ("/ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                {sNEvtSkimv2: 9991620, sSumEvtSkimv2: 9659256}),
    ("/ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                {sNEvtSkimv2: 4757708, sSumEvtSkimv2: 4472444 }),
    ("/ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",              {sNEvtSkimv2: 4869506, sSumEvtSkimv2: 4700684}),
    ("/ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",                {sNEvtSkimv2: 4734791, sSumEvtSkimv2: 4456443}),
    ("/ZHToMuMuG_M125_Dalitz_012j_TuneCP5_13TeV_amcatnloFXFX_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",    {}),
    ("/ZHToTauTau_M125_CP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",                         {}),    

    ("/ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",              {}),
    ("/ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",              {}),
    ("/ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",            {}),
    ("/ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",              {}),

    ("/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",                        {sNEvtSkimv2: 7825000, sSumEvtSkimv2: 7661778}),
    ("/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM",                     {}),
    ("/ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",                    {}),

    #("/bbHToBB_M-125_4FS_yt2_TuneCP5-13TeV-amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM",             {sNEvtSkimv2: 0, sSumEvtSkimv2: 0}),
    

    # , sNEvtSkimv2: , sSumEvtSkimv2: 
    ## SUSY_GluGluH_01J_HToAATo4B_M-*   and   SUSY_GluGluH_01J_HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*GluGluH*HToAATo4B*M*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/SUSY_GluGluH_01J_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_GluGluH_01J_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),

    # SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (308) / (5379) = 5.726e-02 +- 3.168e-03	Matching efficiency = 0.5 +/- 0.0     Cross-section = 48.61 pb * 0.057	
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  238687, sSumEvtSkimv2: 238687}), # filter efficiency 0.05
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  254262, sSumEvtSkimv2: 254262}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  240973, sSumEvtSkimv2: 240973}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  251887, sSumEvtSkimv2: 251887}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  248978, sSumEvtSkimv2: 248978}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  252264, sSumEvtSkimv2: 252264}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  251418, sSumEvtSkimv2: 251418}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  248045, sSumEvtSkimv2: 248045}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  249671, sSumEvtSkimv2: 249671}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  228347, sSumEvtSkimv2: 228347}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  244590, sSumEvtSkimv2: 244590}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-11.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 156385, sSumEvtSkimv2: 156385}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-11.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 163815, sSumEvtSkimv2: 163815}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-12.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 159084, sSumEvtSkimv2: 159084}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-13.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 177337, sSumEvtSkimv2: 177337}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-13.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 175432, sSumEvtSkimv2: 175432}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-14.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 156058, sSumEvtSkimv2: 156058}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-16.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 151409, sSumEvtSkimv2: 151409}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-17.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 177947, sSumEvtSkimv2: 177947}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 181600, sSumEvtSkimv2: 181600}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-21.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 184671, sSumEvtSkimv2: 184671}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-23.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 172630, sSumEvtSkimv2: 172630}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-27.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 169873, sSumEvtSkimv2: 169873}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-32.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 167890, sSumEvtSkimv2: 167890}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-37.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 168514, sSumEvtSkimv2: 168514}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-42.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 154175, sSumEvtSkimv2: 154175}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-47.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 162892, sSumEvtSkimv2: 162892}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-52.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 161334, sSumEvtSkimv2: 161334}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-57.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 170905, sSumEvtSkimv2: 170905}),
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-62.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 171029, sSumEvtSkimv2: 171029}),


    ## VBF HToAATo4B_M-* and VBF HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*VBF*HToAATo4B*M*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/SUSY_VBFH_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_VBFH_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),

    # SUSY_VBFH_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8		Filter efficiency (event-level)= (1730) / (10000) = 1.730e-01 +- 3.782e-03	Matching efficiency = 1.0 +/- 0.0 	Cross-section = 3.766 pb * 0.173 
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  101503, sSumEvtSkimv2: 101503}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  102251, sSumEvtSkimv2: 102251}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:   96201, sSumEvtSkimv2:  96201}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:   97487, sSumEvtSkimv2:  97487}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:   97312, sSumEvtSkimv2:  97312}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:   96089, sSumEvtSkimv2:  96089}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:   97959, sSumEvtSkimv2:  97959}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  104599, sSumEvtSkimv2: 104599}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  105876, sSumEvtSkimv2: 105876}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  101319, sSumEvtSkimv2: 101319}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  103443, sSumEvtSkimv2: 103443}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-11.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 121143, sSumEvtSkimv2: 121143}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-11.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 127149, sSumEvtSkimv2: 127149}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-12.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 134160, sSumEvtSkimv2: 134160}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-13.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 114617, sSumEvtSkimv2: 114617}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-13.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 128255, sSumEvtSkimv2: 128255}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-14.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 126285, sSumEvtSkimv2: 126285}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-16.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 136814, sSumEvtSkimv2: 136814}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-17.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 129836, sSumEvtSkimv2: 129836}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 117232, sSumEvtSkimv2: 117232}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-21.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 126103, sSumEvtSkimv2: 126103}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-23.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 127154, sSumEvtSkimv2: 127154}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-27.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 122972, sSumEvtSkimv2: 122972}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-32.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 127000, sSumEvtSkimv2: 127000}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-37.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 134478, sSumEvtSkimv2: 134478}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-42.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 130092, sSumEvtSkimv2: 130092}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-47.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 133929, sSumEvtSkimv2: 133929}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-52.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 132667, sSumEvtSkimv2: 132667}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-57.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 124337, sSumEvtSkimv2: 124337}),
    ("/SUSY_VBFH_HToAATo4B_Pt150_M-62.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 129561, sSumEvtSkimv2: 129561}),


    ## WH HToAATo4B-M-* and WH HToAATo4B_Pt150_-M-*
    # dasgoclient --query="dataset=/SUSY*WH*HToAATo4B*M*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM" 
    ("/SUSY_WH_WToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_WH_WToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),

    # SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (1317) / (10000) = 1.317e-01 +- 3.382e-03	Matching efficiency = 1.0 +/- 0.0	Cross-section = 1.358 pb * 0.132 
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  51724, sSumEvtSkimv2: 51724}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  51935, sSumEvtSkimv2: 51935}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  50283, sSumEvtSkimv2: 50283}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  49781, sSumEvtSkimv2: 49781}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  49700, sSumEvtSkimv2: 49700}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  50800, sSumEvtSkimv2: 50800}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  49376, sSumEvtSkimv2: 49376}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  52608, sSumEvtSkimv2: 52608}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  50209, sSumEvtSkimv2: 50209}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  52558, sSumEvtSkimv2: 52558}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  50331, sSumEvtSkimv2: 50331}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-11.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 73439, sSumEvtSkimv2: 73439}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-11.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 73264, sSumEvtSkimv2: 73264}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-12.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 67679, sSumEvtSkimv2: 67679}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-13.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 67789, sSumEvtSkimv2: 67789}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-13.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 70536, sSumEvtSkimv2: 70536}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-14.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 71647, sSumEvtSkimv2: 71647}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-16.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 67214, sSumEvtSkimv2: 67214}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-17.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 67714, sSumEvtSkimv2: 67714}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 68153, sSumEvtSkimv2: 68153}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-21.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 65740, sSumEvtSkimv2: 65740}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-23.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 73737, sSumEvtSkimv2: 73737}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-27.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 71324, sSumEvtSkimv2: 71324}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-32.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 73450, sSumEvtSkimv2: 73450}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-37.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 67496, sSumEvtSkimv2: 67496}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-42.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 79965, sSumEvtSkimv2: 79965}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-47.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 65112, sSumEvtSkimv2: 65112}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-52.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 59627, sSumEvtSkimv2: 59627}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-57.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 77635, sSumEvtSkimv2: 77635}),
    ("/SUSY_WH_WToAll_HToAATo4B_Pt150_M-62.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 79193, sSumEvtSkimv2: 79193}),


    ## ZH HToAATo4B-M-* and ZH HToAATo4B_Pt150_-M-*
    # dasgoclient --query="dataset=/SUSY*ZH*HToAATo4B*M*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    
    # SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (1289) / (10000) = 1.289e-01 +- 3.351e-03	Matching efficiency = 1.0 +/- 0.0		Cross-section = 0.880 pb * 0.129 
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  51149, sSumEvtSkimv2:  51149}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  47582, sSumEvtSkimv2:  47582}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  49243, sSumEvtSkimv2:  49243}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  52568, sSumEvtSkimv2:  52568}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  49151, sSumEvtSkimv2:  49151}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  48670, sSumEvtSkimv2:  48670}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  50303, sSumEvtSkimv2:  50303}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  51097, sSumEvtSkimv2:  51097}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  50786, sSumEvtSkimv2:  50786}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  48896, sSumEvtSkimv2:  48896}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  50888, sSumEvtSkimv2:  50888}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-11.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 74086, sSumEvtSkimv2: 74086}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-11.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 73223, sSumEvtSkimv2: 73223}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-12.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 65146, sSumEvtSkimv2: 65146}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-13.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 75089, sSumEvtSkimv2: 75089}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-13.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 68796, sSumEvtSkimv2: 68796}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-14.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 70436, sSumEvtSkimv2: 70436}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-16.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 67714, sSumEvtSkimv2: 67714}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-17.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 60862, sSumEvtSkimv2: 60862}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 66558, sSumEvtSkimv2: 66558}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-21.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 70018, sSumEvtSkimv2: 70018}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-23.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 60642, sSumEvtSkimv2: 60642}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-27.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 60472, sSumEvtSkimv2: 60472}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-32.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 60881, sSumEvtSkimv2: 60881}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-37.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 65110, sSumEvtSkimv2: 65110}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-42.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 68121, sSumEvtSkimv2: 68121}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-47.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 64816, sSumEvtSkimv2: 64816}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-52.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 61659, sSumEvtSkimv2: 61659}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-57.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 60123, sSumEvtSkimv2: 60123}),
    ("/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-62.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 66736, sSumEvtSkimv2: 66736}),


    ## ttH HToAATo4B_M-* and HToAATo4B_Pt150_M-*
    # dasgoclient --query="dataset=/SUSY*TTH*HToAATo4B*M*/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v*/NANOAODSIM"
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {}),

    # SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8	Filter efficiency (event-level)= (2850) / (10000) = 2.850e-01 +- 4.514e-03	Matching efficiency = 1.0 +/- 0.0 Cross-section = 0.5071 pb * 0.285 
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  100110, sSumEvtSkimv2:  100110 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:   95752, sSumEvtSkimv2:   95752}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  100480, sSumEvtSkimv2:  100480 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-25_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  105083, sSumEvtSkimv2:  105083 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  102876, sSumEvtSkimv2:  102876 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-35_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:   99671, sSumEvtSkimv2:   99671 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:   99735, sSumEvtSkimv2:   99735 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  100078, sSumEvtSkimv2:  100078 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-50_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:   99284, sSumEvtSkimv2:   99284 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-55_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:  100285, sSumEvtSkimv2:  100285 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-60_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2:   97782, sSumEvtSkimv2:   97782 }),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-11.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 142055, sSumEvtSkimv2: 142055}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-11.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 142562, sSumEvtSkimv2: 142562}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-12.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 143167, sSumEvtSkimv2: 143167}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-13.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 132053, sSumEvtSkimv2: 132053}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-13.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 119337, sSumEvtSkimv2: 119337}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-14.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 114879, sSumEvtSkimv2: 114879}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-16.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 145289, sSumEvtSkimv2: 145289}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-17.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 144755, sSumEvtSkimv2: 144755}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 141231, sSumEvtSkimv2: 141231}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-21.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 147211, sSumEvtSkimv2: 147211}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-23.0_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 125587, sSumEvtSkimv2: 125587}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-27.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 137587, sSumEvtSkimv2: 137587}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-32.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 158245, sSumEvtSkimv2: 158245}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-37.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 156246, sSumEvtSkimv2: 156246}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-42.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 151162, sSumEvtSkimv2: 151162}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-47.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 141611, sSumEvtSkimv2: 141611}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-52.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 122382, sSumEvtSkimv2: 122382}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-57.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 117039, sSumEvtSkimv2: 117039}),
    ("/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-62.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM", {sNEvtSkimv2: 135902, sSumEvtSkimv2: 135902}),


    ## H->aa->4Tau
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYGluGluHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVBFHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    
    ("/SUSYVHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYVHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),

    ("/SUSYttHToAA_AToTauTau_M-125_M-4_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-5_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-6_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-7_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-8_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-9_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-10_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-11_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-12_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-13_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-14_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-15_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-16_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-17_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-18_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-19_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-20_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),
    ("/SUSYttHToAA_AToTauTau_M-125_M-21_TuneCP5_13TeV_PSWeights_pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM", {}),

    

    
    


    

    ## JetHT
    # dasgoclient --query="dataset=/JetHT/*2018*UL*MiniAODv2_NanoAODv9-*/NANOAOD"
    # XS (cross-section) does not matter for data sample
    # Use JetHT/NanoAOD_GT38 instead
    #("/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD", {}),
    #("/JetHT/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD", {}),
    #("/JetHT/Run2018C-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD", {}),
    #("/JetHT/Run2018D-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD", {}),
    
    # dasgoclient --query="dataset=/JetHT/*2017*UL*MiniAODv2_NanoAODv9*/NANOAOD"
    # XS (cross-section) does not matter for data sample
    ("/JetHT/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), # run 297047 to 299329
    ("/JetHT/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), # run 299368 to 302029
    ("/JetHT/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), # run 302031 to 302663
    ("/JetHT/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), # run 303824 to 304797
    ("/JetHT/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), # run 305040 to 306460
    
    ## SingleMuon
    # dasgoclient --query="dataset=/SingleMuon/*2017*UL*MiniAODv2_NanoAODv9_GT36*/NANOAOD"
    # XS (cross-section) does not matter for data sample    
    ("/SingleMuon/Run2017B-UL2017_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}),    
    ("/SingleMuon/Run2017C-UL2017_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}),
    ("/SingleMuon/Run2017D-UL2017_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}),
    ("/SingleMuon/Run2017E-UL2017_MiniAODv2_NanoAODv9_GT36-v2/NANOAOD", {}),
    ("/SingleMuon/Run2017F-UL2017_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD", {}),
    
    ## EGamma
    # dasgoclient --query="dataset=/SingleElectron/*2017*UL*MiniAODv2_NanoAODv9*/NANOAOD"
    # XS (cross-section) does not matter for data sample    
    ("/SingleElectron/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/SingleElectron/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/SingleElectron/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/SingleElectron/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/SingleElectron/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}),

    ## MET
    # dasgoclient --query="dataset=/MET/*2017*UL*MiniAODv2_NanoAODv9*/NANOAOD"
    # XS (cross-section) does not matter for data sample    
    ("/MET/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/MET/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/MET/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/MET/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/MET/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}),
    
    ## BTagCSV
    # dasgoclient --query="dataset=/BTagCSV/Run2017*UL2017_MiniAODv2_NanoAODv9*/NANOAOD"
    ("/BTagCSV/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/BTagCSV/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/BTagCSV/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/BTagCSV/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}), 
    ("/BTagCSV/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD", {}),

    #("", {}),    
])