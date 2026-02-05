from glob import glob
import ROOT as R

def getTreeEntries(sFile, sTreeName="Events"):
    tChain = R.TChain(sTreeName)
    tChain.Add(sFile)
    #print(f"{sFile = }, {tChain.GetEntries() = }")
    return tChain.GetEntries()

def countTotalEntriesFromTrees(sFile, sTreeName="Events"):
    files_list = glob(sFile)
    nEventsTotal = 0
    for sFile in files_list:
        #nEvents_i = getTreeEntries(sFile)
        #nEventsTotal += nEvents_i
        nEventsTotal += getTreeEntries(sFile)
    return nEventsTotal



if __name__ == "__main__":
    
    # TChain does not work with wild-card * for directories and filenames together

    ''''
    sFName = "/eos/cms/store/group/phys_susy/HToaaTo4b/MiniAOD/2018/MC/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-62.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18/*/MiniAODv2_*.root"
    files_list = glob(sFName)
    nEventsTotal = 0
    for sFile in files_list:
        nEvents_i = getTreeEntries(sFile)
        nEventsTotal += nEvents_i
        #print(f"\t {sFile}: {nEvents_i},    {nEventsTotal}")

    print(f"{len(files_list) = }, {nEventsTotal = }")

    print(f"\n\n TChain with *::\n{sFName}: {getTreeEntries(sFName) = }")
    '''

    ## Higgs production mode
    prodmodes=[
        "SUSY_GluGluH_01J_HToAATo4B",
        "SUSY_VBFH_HToAATo4B",
        "SUSY_WH_WToAll_HToAATo4B",
        "SUSY_ZH_ZToAll_HToAATo4B",
        "SUSY_TTH_TTToAll_HToAATo4B",
    ]
    #prodmodes=["SUSY_TTH_TTToAll_HToAATo4B"]

    mApoints=["11.0", "11.5", "12.5", "13.0", "13.5", "14.0", "16.0", "17.0", "18.5", "21.5", "23.0", "27.5", "32.5", "37.5", "42.5", "47.5", "52.5", "57.5", "62.5"]
    #mApoints=["11.0"]

    Eras=[
        "RunIISummer20UL18",
        "RunIISummer20UL17",
        "RunIISummer20UL16",
        "RunIISummer20UL16APV",
    ]
    #Eras=["RunIISummer20UL16",]
    #Eras=["RunIISummer20UL18",]

    countEntriesNanoAOD = False # Default: True

    #sFNEvts="nEvents_SignalIntermediateMassPoints.txt"
    #print(f"\n rm ${sFNEvts} : \n")
    #rm ${sFParams}

    nEvents_dict = {}
    for Era in Eras:
        EraYear='2018'
        if   "16" in Era and "APV" in Era:  EraYear = "2016APV"
        elif "16" in Era:                   EraYear = "2016"
        elif "17" in Era:                   EraYear = "2017"
        elif "18" in Era:                   EraYear = "2018"
        print(f"\n\n{Era = }, {EraYear = }")
        nEvents_dict[EraYear] = {}

        for prodmode in prodmodes:
            print(f"\n\n{prodmode = }")
            nEvents_dict[EraYear][prodmode] = {}

            for mA in mApoints:

                # /eos/cms/store/group/phys_susy/HToaaTo4b/MiniAOD/2018/MC/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-11.5_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18/0018/MiniAODv2_1898_nEvents500.root
                sMiniAOD = f"/eos/cms/store/group/phys_susy/HToaaTo4b/MiniAOD/{EraYear}/MC/{prodmode}_Pt150_M-{mA}_TuneCP5_13TeV_madgraph_pythia8/{Era}/*/MiniAODv2_*.root"

                if EraYear == "2018":
                    crabDir_timeStamp="20251107_000000"
                else:
                    crabDir_timeStamp="20260121_000000"
                # /eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v2_2024_11_22/SUSY_VBFH_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/r1/20251107_000000/0000/PNet_v1_14_Skim.root
                sNanoAOD_0 = f"/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/{EraYear}/MC/PNet_v2_2024_11_22/{prodmode}_Pt150_M-{mA}_TuneCP5_13TeV_madgraph_pythia8/r1/{crabDir_timeStamp}/*/PNet_*_Skim.root"
                # /eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v2_2024_11_22/SUSY_VBFH_HToAATo4B_Pt150_M-18.5_TuneCP5_13TeV_madgraph_pythia8/r1/PNet_v1_Skim_0.root
                sNanoAOD_hadded = f"/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/{EraYear}/MC/PNet_v2_2024_11_22/{prodmode}_Pt150_M-{mA}_TuneCP5_13TeV_madgraph_pythia8/r1/PNet_v1_Skim_*.root"


                nEvents_MiniAOD   = countTotalEntriesFromTrees(sMiniAOD)
                nEvents_NanoAOD_0 = -1
                nEvents_NanoAOD   = -1
                if countEntriesNanoAOD:
                    nEvents_NanoAOD_0 = countTotalEntriesFromTrees(sNanoAOD_0)
                    nEvents_NanoAOD   = countTotalEntriesFromTrees(sNanoAOD_hadded)
                print(f"\t {mA = }: \t {nEvents_MiniAOD} \t {nEvents_NanoAOD_0} \t {nEvents_NanoAOD}", flush=True)

                nEvents_dict[EraYear][prodmode][mA] = {'MiniAOD': nEvents_MiniAOD, 'NanoAOD': nEvents_NanoAOD, 'NanoAOD_0': nEvents_NanoAOD}

    print(f"\n\n\n\nPrint Sample: \t nEvents_MiniAOD \t nEvents_NanoAOD \t nEvents_NanoAOD_0")
    for Era in Eras:
        EraYear='2018'
        if   "16" in Era and "APV" in Era:  EraYear = "2016APV"
        elif "16" in Era:                   EraYear = "2016"
        elif "17" in Era:                   EraYear = "2017"
        elif "18" in Era:                   EraYear = "2018"
        
        print(f"\n{EraYear:}")
        for prodmode in prodmodes:                        
            for mA in mApoints:
                print(f"{prodmode}_Pt150_M-{mA} \t {nEvents_dict[EraYear][prodmode][mA]['MiniAOD']} \t {nEvents_dict[EraYear][prodmode][mA]['NanoAOD']} \t {nEvents_dict[EraYear][prodmode][mA]['NanoAOD_0']}", flush=True)



    print(f"\n\n\n\nPrint for sample dictionary")
    for Era in Eras:
        EraYear='2018'
        if   "16" in Era and "APV" in Era:  EraYear = "2016APV"
        elif "16" in Era:                   EraYear = "2016"
        elif "17" in Era:                   EraYear = "2017"
        elif "18" in Era:                   EraYear = "2018"

        if   "16" in Era and "APV" in Era:  DatasetNamePart2 = "RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1"
        elif "16" in Era:                   DatasetNamePart2 = "RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1"
        elif "17" in Era:                   DatasetNamePart2 = "RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1"
        elif "18" in Era:                   DatasetNamePart2 = "RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1"

        print(f"\n{EraYear:}")
        for prodmode in prodmodes:                        
            for mA in mApoints:
                print('    ("/%s_Pt150_M-%s_TuneCP5_13TeV_madgraph_pythia8/%s/NANOAODSIM", {sNEvtSkimv2: %d, sSumEvtSkimv2: %d}),'%(prodmode,mA,DatasetNamePart2, nEvents_dict[EraYear][prodmode][mA]['MiniAOD'],nEvents_dict[EraYear][prodmode][mA]['MiniAOD']));
            print(' ')

                
        




