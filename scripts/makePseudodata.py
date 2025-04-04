import ROOT
import os
import re
import shutil

MASSH   = 'pnet_vs_massA34a'    ## Higgs mass regression (mass, msoft, pnet)


def merge_histograms(path, signal):
    sTagExtra = '' # '_Xto4bv2': for Mohamed,  '': for Siddhesh

    sOpDataOriginalFile = path + f"{signal+sTagExtra}_Data_2018.root"
    sOpDataBackupFile = path + f"{signal+sTagExtra}_Data_2018_backup.root"
    sOpPseudodataFile = path + f"{signal+sTagExtra}_Data_2018.root"
    
    # Keep back up of data
    shutil.copy(sOpDataOriginalFile, sOpDataBackupFile)
    print(f"{os.listdir(path) = }", flush=True)
    

    # Output file for merged histograms
    #output_file = ROOT.TFile(path + f"{signal}_Xto4bv2_Data_2018.root", "RECREATE")
    output_file = ROOT.TFile(sOpPseudodataFile, "RECREATE")

    # Get all ROOT files except those containing "Htoaato4b" or "Data"
    root_files = [path+f for f in os.listdir(path) if f.endswith(".root") and "Htoaato4b" not in f and "Data" not in f]

    # Initialize histograms for merging
    h_pass = None
    h_fail = None

    # Loop over ROOT files
    for file_name in root_files:
        input_file = ROOT.TFile.Open(file_name, "READ")
        if not input_file or input_file.IsZombie():
            print(f"Error opening file: {file_name}")
            continue

        # Loop over histograms in the file
        for key in input_file.GetListOfKeys():
            hist_name = key.GetName()
            print("Processing histogram:", hist_name)

            hist = input_file.Get(hist_name)
            if not hist:
                continue  # Skip if histogram does not exist

            # Merge "Pass" histograms
            if 'Pass' in hist_name:
                if h_pass is None:
                    h_pass = hist.Clone("histogram_Pass")
                    h_pass.SetDirectory(0)
                else:
                    h_pass.Add(hist)

            # Merge "Fail" histograms
            elif 'Fail' in hist_name:
                if h_fail is None:
                    h_fail = hist.Clone("histogram_Fail")
                    h_fail.SetDirectory(0)
                else:
                    h_fail.Add(hist)

        input_file.Close()

    # Save merged histograms
    output_file.cd()
    if h_pass:
        h_pass.Write()
    if h_fail:
        h_fail.Write()
    output_file.Close()

    print(f"Merged histograms saved in {sOpPseudodataFile}")

    Bkg_file = ROOT.TFile.Open(sOpPseudodataFile, "UPDATE")
    if not Bkg_file or Bkg_file.IsZombie():
        print("Error opening background file")
        return

    h_pass = Bkg_file.Get("histogram_Pass")
    h_fail = Bkg_file.Get("histogram_Fail")

    data_file = ROOT.TFile.Open(sOpDataBackupFile, "READ")
    if not data_file or data_file.IsZombie():
        print("Error opening data.root")
        return

    h_data_pass = data_file.Get(f"{signal+sTagExtra}_Data_2018_{MASSH}_WP40_Pass_Nom")
    h_data_fail = data_file.Get(f"{signal+sTagExtra}_Data_2018_{MASSH}_WP40_Fail_Nom")

    # Compute and apply scale factors
    if h_pass and h_data_pass:
        integral_pass = h_pass.Integral() if h_pass else 0
        integral_data_pass = h_data_pass.Integral() if h_data_pass else 0
        if integral_pass > 0 and integral_data_pass > 0:
            scale_pass = integral_data_pass / integral_pass
            h_pass.Scale(scale_pass)
            print(f"Scaling Pass histograms by {scale_pass}")

    if h_fail and h_data_fail:
        integral_fail = h_fail.Integral() if h_fail else 0
        integral_data_fail = h_data_fail.Integral() if h_data_fail else 0
        if integral_fail > 0 and integral_data_fail > 0:
            scale_fail = integral_data_fail / integral_fail
            h_fail.Scale(scale_fail)
            print(f"Scaling Fail histograms by {scale_fail}")

    # Save scaled histograms
    Bkg_file.cd()
    if h_pass:
        h_pass.Write(f"{signal+sTagExtra}_Data_2018_{MASSH}_WP40_Pass_Nom")
    if h_fail:
        h_fail.Write(f"{signal+sTagExtra}_Data_2018_{MASSH}_WP40_Fail_Nom")
    Bkg_file.Close()
    data_file.Close()

    print(f"Scaled histograms saved in {signal+sTagExtra}_Data_2018.root")

# Run the function
pathOriginal   = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_gg0l_FullSyst/2018/2DAlphabet_inputFiles_pseudodata/gg0lIncl/' #'/eos/user/m/moanwar/htoaa/analysis/VBF_channel/2DAlphabetfiles_VBF_BKgIncl/VBFHi_Xto4bv2/'
pathPseudoData = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_gg0l_FullSyst/2018/2DAlphabet_inputFiles_pseudodata_0/gg0lIncl/'
signal = 'gg0lIncl' #'VBFjjHi'

# Copy the original directory and work with the copy
shutil.copytree(pathOriginal, pathPseudoData, dirs_exist_ok=True)
merge_histograms(pathPseudoData, signal)
