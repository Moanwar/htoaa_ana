import os
import sys

import uproot

if __name__ == '__main__':
    xrootd_redirectorName = "root://xrootd-cms.infn.it//"
    sInputFiles = [
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/BC10B8D4-88BC-CB47-9746-7E93AE25CD7B.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/51E2F535-922D-BC46-92E2-C01060B5D992.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/5E60B5D0-7718-7343-B1EE-19704A647E4E.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/CBD62D40-7933-4E4F-9D27-59F40372E3FD.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/80000/9D1C2E44-61D4-3D4C-9A3A-10D66C8532C8.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/0B8D33AD-A7A4-2C44-9359-EC7827A73B61.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/1C334E46-52FB-C34B-BC69-726312AECF02.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/91531500-CC2D-9648-87B5-12546D14890C.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/4670D738-2AD4-4743-86CF-8633A7DB8129.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/B2B49191-D807-CC41-B8DA-5AE38DF1D461.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/B87227FB-E788-2547-AB87-F88F954A0E25.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/749F8769-592F-E74E-9404-4DB95863FEFE.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/03E733AE-EEAF-674A-BDF5-5EA941F7D7E9.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/C37B2077-8142-3244-9AA3-4C4E4932569E.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/8E1F11D3-1CD8-3D44-AFA8-49E6C38BB54B.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/206062EF-523E-EE46-9746-869C5852EB30.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/00EA8466-C151-3E4C-ACF9-A908DB67D610.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2560000/0FC587FA-1940-E847-9B29-72D9F3CE7EB0.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/80000/7855F5E5-8C0B-5A41-8A07-C1CC40277C58.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/80000/52847EE8-4FC6-3F42-8F24-976B09110256.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2430000/C281F483-01E3-3842-9298-74E78BB7279B.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/F25278F5-7C34-CD4E-ACD4-523F3E964BD1.root",
	    "/store/mc/RunIISummer20UL18NanoAODv9/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2810000/7F6E2AB9-A31A-8C46-8772-3E1B3AC0700F.root"
        ]

    
    for iFile in range(len(sInputFiles)):
        if sInputFiles[iFile].startswith("/store/"): # LFN: Logical File Name
            sInputFiles[iFile] = xrootd_redirectorName + sInputFiles[iFile]
    print(f"sInputFiles ({len(sInputFiles)}): {sInputFiles}"); sys.stdout.flush()

    nEventsTot = 0
    for sInputFile in sInputFiles:
        fInputFile = None
        try:
            fInputFile = uproot.open(sInputFile)
        except:
            print(f"File {sInputFile} could not open: {fInputFile}")
        else:
        #with uproot.open(sInputFile) as fInputFile:
            print(f"File {sInputFile} is open: {fInputFile}")
            nEventsTot += fInputFile['Events'].num_entries
            print(f"fInputFile['Events']: {fInputFile['Events']}, nEvents: {fInputFile['Events'].num_entries}, nEventsTot: {nEventsTot}")

    print(f"nEventsTot: {nEventsTot}")
