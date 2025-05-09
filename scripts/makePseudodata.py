import ROOT
import os
import re
import shutil


def merge_histograms(path, signal):
    sTagExtra = '' # '_Xto4bv2': for Mohamed,  '': for Siddhesh

    sOpDataOriginalFile = path + f"{signal+sTagExtra}_Data_2018.root"
    sOpDataBackupFile   = path + f"{signal+sTagExtra}_Data_2018_backup.root"
    sOpPseudodataFile   = path + f"{signal+sTagExtra}_Data_2018.root"

    print(f"\n\n{path = }")
    #print(f"Initially {os.listdir(path) = }", flush=True)

    # Keep back up of data
    #shutil.copy(sOpDataOriginalFile, sOpDataBackupFile)
    shutil.move(sOpDataOriginalFile, sOpDataBackupFile)
    #print(f"\nStage 1   {os.listdir(path) = }", flush=True)

    # Open i/p files
    fData       = ROOT.TFile(sOpDataBackupFile, "READ")
    fPseudodata = ROOT.TFile(sOpPseudodataFile, "RECREATE")
    fBkgs_dict  = { f:ROOT.TFile(path+f, "READ")  for f in os.listdir(path) if f.endswith(".root") and "Htoaato4b" not in f and "Data" not in f}
    #print(f"\nStage 2   {os.listdir(path) = }", flush=True)
    print(f"\n{fBkgs_dict.keys() = }")

    # Loop over histogram-name combinations
    for MASSH in MASSH_list:
        for WP in WPs:
            for Region in Regions:
                Systematics = 'Nom'
                print(f"\n\n{MASSH = }, {WP = }, {Region = }, {Systematics = }")

                hData = None
                hPseudodata = None
                sPseudodata = None
                
                # Read hData
                for histNameInFile in fData.GetListOfKeys():
                    histNameInFile = histNameInFile.GetName()
                    # select histogram name under consideration
                    if MASSH in histNameInFile and WP in histNameInFile and Region in histNameInFile and Systematics in histNameInFile:
                        sPseudodata = histNameInFile
                        hData = fData.Get(histNameInFile)
                        hData.SetName('%s_backup' % (hData.GetName()))
                        hData.SetTitle('%s_backup' % (hData.GetTitle()))                        
                print(f"{sPseudodata = }, hData: {hData.Integral()}")

                # hadd all backgrounds
                for sBkg, fBkg in fBkgs_dict.items():
                    # Read hBkg
                    for histNameInFile in fBkg.GetListOfKeys():
                        histNameInFile = histNameInFile.GetName()
                        # select histogram name under consideration
                        if MASSH in histNameInFile and WP in histNameInFile and Region in histNameInFile and Systematics in histNameInFile:
                            hBkg_i = fBkg.Get(histNameInFile)
                            if hPseudodata is None:
                                hPseudodata = hBkg_i.Clone(sPseudodata)
                                hPseudodata.SetDirectory(0)
                            else:
                                hPseudodata.Add(hBkg_i)
                            print(f"\t {sBkg} \t:{histNameInFile}, \t hBkg_i: {hBkg_i.Integral()}, \t hPseudodata: {hPseudodata.Integral()}")
                print(f"{sPseudodata = }, hData: {hData.Integral()}, \t hPseudodata: {hPseudodata.Integral()}")

                # Scale pseudodata to match number of data events
                scale_factor = hData.Integral() / hPseudodata.Integral() if hPseudodata.Integral() > 0 else 1
                hPseudodata.Scale(scale_factor)
                print(f"{sPseudodata = }, hData: {hData.Integral()}, \t hPseudodata: {hPseudodata.Integral()} finally <<<<<<<<<<<<<<<<<<<<")

                fPseudodata.cd()
                hPseudodata.Write()

    fData.Close()
    fPseudodata.Close()
    fBkgs_dict.clear()




                        
                        
MASSH_list   = ['pnet_vs_massA34a']    ## Higgs mass regression (mass, msoft, pnet)
WPs          = ['WP40', 'WP60']
Regions      = ['Pass', 'Fail']
Categories   = ['gg0lHi', 'gg0lLo', 'gg0lIncl']

for signal in Categories:
    #pathOriginal   = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_gg0l_FullSyst/2018/2DAlphabet_inputFiles/%s/' % (signal) #'/eos/user/m/moanwar/htoaa/analysis/VBF_channel/2DAlphabetfiles_VBF_BKgIncl/VBFHi_Xto4bv2/'
    #pathPseudoData = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_gg0l_FullSyst/2018/2DAlphabet_inputFiles_pseudodata_Correct_1/%s/' % (signal)
    pathOriginal   = '/eos/cms/store/user/ssawant/htoaa/analysis/20250502_gg0l_FullSyst/2018/2DAlphabet_inputFiles/%s/' % (signal) #'/eos/user/m/moanwar/htoaa/analysis/VBF_channel/2DAlphabetfiles_VBF_BKgIncl/VBFHi_Xto4bv2/'
    pathPseudoData = '/eos/cms/store/user/ssawant/htoaa/analysis/20250502_gg0l_FullSyst/2018/2DAlphabet_inputFiles_pseudodata/%s/' % (signal)
    

    # Copy the original directory and work with the copy
    if os.path.isdir(pathPseudoData): shutil.rmtree(pathPseudoData)
    shutil.copytree(pathOriginal, pathPseudoData, dirs_exist_ok=True)
    merge_histograms(pathPseudoData, signal)
