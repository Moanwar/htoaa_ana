import os
import sys
from datetime import datetime
import subprocess
import shlex
import logging
import json
import numpy as np
import math
import awkward as ak
import uproot as uproot
from coffea import hist as coffea_hist
import hist as hist
#import ROOT as R
from parse import *
import logging
import correctionlib
from coffea import util
from coffea.jetmet_tools import CorrectedJetsFactory, JECStack
from coffea.lookup_tools import extractor

from htoaa_Settings import * 
from htoaa_Samples import (
    kData, kQCD_bEnrich, kQCD_bGen, kQCDIncl, kZJets, kWJets
)
#from numba import jit

def calculate_lumiScale(luminosity, crossSection, sumEvents):
    lumiScale = 1
    # as crosssection is in pb and luminosity in fb
    pb_to_fb_conversionFactor = 1000
    
    if sumEvents != 0: lumiScale = luminosity * crossSection * pb_to_fb_conversionFactor / sumEvents
    return lumiScale


def getSampleHTRange(sample_datasetNameFull):
    sample_HT_Min = sample_HT_Max = None
    # for e.g. sample_dataset: "QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8"
    for sample_dataset_parts in sample_datasetNameFull.split('_'):
        if 'HT' in sample_dataset_parts and 'to' in sample_dataset_parts:
            # HT100to200
            sample_HT_Min = int(sample_dataset_parts.split('HT')[1].split('to')[0])
            sample_HT_Max = sample_dataset_parts.split('HT')[1].split('to')[1]
            if sample_HT_Max == 'Inf':
                sample_HT_Max = kLHE_HT_Max
            else:
                sample_HT_Max = int(sample_HT_Max)            
#            try:
#                sample_HT_Max = int(sample_HT_Max) # to take care of sample_HT_Max = 'Inf'
#            except:
#                sample_HT_Max = -1
            break
    return sample_HT_Min, sample_HT_Max
    
'''
def update_crosssection(sample_category, sample_dataset, sample_crossSection):
    if sample_category not in [kQCD_bGen]: return sample_crossSection

    # HTSamplesStitch SF -------------------------------------------------------------------
    sample_HT_Min, sample_HT_Max = getSampleHTRange(sample_dataset)
    sample_HT_toUse = sample_HT_Min # int(sample_HT_Min)
    
    sIpFile_HTSamplesStitchSF        = Corrections['HTSamplesStitch']['inputFile']
    sHistogramName_HTSamplesStitchSF = Corrections['HTSamplesStitch']['histogramName']
    sHistogramName_HTSamplesStitchSF = sHistogramName_HTSamplesStitchSF.replace('$SAMPLECATEGORY', sample_category)

    HTSamplesStitchSF = None
    ipFile_HTSamplesStitchSF = R.TFile(sIpFile_HTSamplesStitchSF)
    if not ipFile_HTSamplesStitchSF.IsOpen():
        logging.error   ("update_crosssection(): Colud not open inputfile %s ." % (sIpFile_HTSamplesStitchSF), exc_info=True)        
        exit(0)

    hHTSamplesStitchSF = None
    hHTSamplesStitchSF = ipFile_HTSamplesStitchSF.Get(sHistogramName_HTSamplesStitchSF)
    if not hHTSamplesStitchSF:
        logging.error   ("update_crosssection(): Histogram %s could not read from inputfile %s ." % (sHistogramName_HTSamplesStitchSF, sIpFile_HTSamplesStitchSF), exc_info=True)
        exit(0)
        
    try:
        HTSamplesStitchSF = hHTSamplesStitchSF.GetBinContent( hHTSamplesStitchSF.FindBin(sample_HT_toUse) )
    except:
        logging.error   ("update_crosssection(): Could not read SF @HT %g from histogram %s could not open." % (sample_HT_toUse, sHistogramName_HTSamplesStitchSF), exc_info=True)
        exit(0)

    sample_crossSection_corr = sample_crossSection * HTSamplesStitchSF
    print(f"update_crosssection():: sample_category: {sample_category}, sample_dataset: {sample_dataset}, sample_crossSection (original): {sample_crossSection}, HTSamplesStitchSF(@HT {sample_HT_toUse}): {HTSamplesStitchSF}, sample_crossSection_corr: {sample_crossSection_corr}")
    ipFile_HTSamplesStitchSF.Close()
    # ----------------------------------------------------------------------------------------
    
    return sample_crossSection_corr
'''    

def getTH1BinContent(histo, xValue):
    # histo: uprootTH1D.to_hist()
    xBin = None
    for iBin in range(len(histo.axes[0])):
        if xValue >= histo.axes[0][iBin][0] and xValue < histo.axes[0][iBin][1]:
            xBin = iBin
            break

    #print(f"htoaa_CommonTools.py::getTH1BinContent():: {histo = },   {xBin = },    {histo[xBin] = },   {histo[xBin].value = }")
    return histo[xBin]

    
def getNanoAODFile(
        fileName, 
        useLocalFileIfExists = True, 
        downloadFile = True, 
        fileNameLocal = './inputFiles/fLocal.root', 
        nTriesToDownload = 3, 
        server = 'lxplus'
        ):
    # MC:
    # DAS file: "/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/2BBE7B3F-C5A7-0D48-A384-FAD06B127FD8.root"
    # eos file: "/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9/2BBE7B3F-C5A7-0D48-A384-FAD06B127FD8.root"
    # Data:
    # DAS file: "/store/data/Run2018A/JetHT/NANOAOD/UL2018_MiniAODv2_NanoAODv9_GT36-v1/2820000/97F68EC0-0E12-C04C-A5D6-2B7A7C6688F8.root"
    # eos file: "/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9_GT36-v1/97F68EC0-0E12-C04C-A5D6-2B7A7C6688F8.root"

    if downloadFile  and  os.path.exists(fileNameLocal):
        # local copy of the i/p file exists
        print(f"{fileNameLocal = } exists")
        return fileNameLocal, True

    fileName_toUse = fileName
    cp_command = 'eos cp' if server in ['lxplus'] else 'xrdcp'
    print(f"htoaa_CommonTools::getNanoAODFile() here1 {datetime.now() = } {fileName = }"); sys.stdout.flush()


    if useLocalFileIfExists and server in ['lxplus']: # 'eos cp' works on lxplus
        fileName_EOS = None

        # check if Central NanoAOD exists in /eos area
        if fileName.startswith("/store/"):
            fileNameTemplate_DAS = "/store/{IsMC}/{SampleProductionCampaign}/{SampleName}/{DatasetTier}/{GT}/{SampleDir}/{SampleFileName}"
            fileNameTemplate_EOS = "/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9/2BBE7B3F-C5A7-0D48-A384-FAD06B127FD8.root"

            r_ = parse(fileNameTemplate_DAS, fileName)
            IsMC                       = 'MC' if r_['IsMC'].lower() == 'mc' else 'data'
            SampleProductionCampaign   = r_['SampleProductionCampaign'] if r_['IsMC'].lower() == 'mc' else '%s-%s' % (r_['SampleProductionCampaign'], r_['GT'])
            SampleName                 = r_['SampleName']
            DatasetTier                = r_['DatasetTier']
            GT                         = r_['GT']
            SampleDir                  = r_['SampleDir']
            SampleFileName             = r_['SampleFileName'] 
            Era                        = None
            if   'UL16' in SampleProductionCampaign or 'UL2016' in SampleProductionCampaign:
                Era = '2016'
            elif 'UL17' in SampleProductionCampaign or 'UL2017' in SampleProductionCampaign:
                Era = '2017'
            elif 'UL18' in SampleProductionCampaign or 'UL2018' in SampleProductionCampaign:
                Era = '2018'       
            fileName_EOS = f"/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/{Era}/{IsMC}/{SampleName}/{SampleProductionCampaign}/{SampleFileName}"

        
        if fileName.startswith("/eos/"): # File is stored on /eos/ area on lxplus
            fileName_EOS = fileName

        print(f"Checking for eos file: {fileName_EOS = }") 
        print(f"htoaa_CommonTools::getNanoAODFile() here2 {datetime.now() = }"); sys.stdout.flush()
        print(f"{fileName_EOS = }: {os.path.exists(fileName_EOS) = } ")

        if (not downloadFile) and (os.path.exists(fileName_EOS)):
            print(f"htoaa_CommonTools::getNanoAODFile() here2.1 {datetime.now() = }: Reading directly {fileName_EOS = }"); sys.stdout.flush()
            return fileName_EOS, True

        # os.path.exists() for files on /eos are always return False. So try 'eos cp' to check if the file on eos exists or not
        if  xrdcpFile(fileName_EOS, fileNameLocal, nTry = 3, cp_command = 'eos cp'):
            print(f"Forced xrdcp for {fileName_EOS = } successful.")
            print(f"htoaa_CommonTools::getNanoAODFile() here3 {datetime.now() = }"); sys.stdout.flush()
        print(f"{fileNameLocal = }: {os.path.exists(fileNameLocal) = } ")            
        print(f"List directory {os.path.dirname(fileNameLocal) = }:  {os.listdir(os.path.dirname(fileNameLocal)) = }")
        if os.path.exists(fileNameLocal):
            return fileNameLocal, True

    print(f"htoaa_CommonTools::getNanoAODFile() here4 {datetime.now() = }"); sys.stdout.flush()
    if fileName_toUse.startswith("/store/"):
        # Copy of the NanoAOD file is not on /eos
        isReadingSuccessful = False
        for redirector in xrootd_redirectorNames:
            fileName_toUse_i = redirector + fileName_toUse
            print(f"getNanoAODFile():: Checking {fileName_toUse_i}"); sys.stdout.flush()
            print(f"htoaa_CommonTools::getNanoAODFile() here5 {datetime.now() = }"); sys.stdout.flush()

            if downloadFile:
                if  xrdcpFile(fileName_toUse_i, fileNameLocal, nTry = 3):
                    print(f"{fileNameLocal = } xrdcp successfully")
                    fileName_toUse = fileNameLocal
                    print(f"htoaa_CommonTools::getNanoAODFile() here6 {datetime.now() = }"); sys.stdout.flush()
                    isReadingSuccessful = True
                    break
                
            else:
                file1 = None
                try:
                    file1 = uproot.open(fileName_toUse_i)
                except:
                    print(f"getNanoAODFile():: File open {fileName_toUse_i} failed"); sys.stdout.flush()
                else:
                    #print(f"\n{redirector + fileName}: file1.keys(): {file1.keys()}")
                    #print(f"\n{redirector + fileName}: file1.keys(): {file1['Events'].numentries}"); sys.stdout.flush()

                    nEntries = file1['Events'].numentries
                    file1.close()
                    #if file1['Events'].num_entries > 0:
                    if nEntries > 0:
                        print(f"{fileName_toUse_i}: {nEntries}"); sys.stdout.flush()
                        fileName_toUse = fileName_toUse_i
                        isReadingSuccessful = True
                        break
            print(f"htoaa_CommonTools::getNanoAODFile() here7 {datetime.now() = }"); sys.stdout.flush()
        print(f"htoaa_CommonTools::getNanoAODFile() here8 {datetime.now() = }"); sys.stdout.flush()
        return fileName_toUse, True




def setXRootDRedirector(fileName, useLocalFileIfExists = True):
    # DAS file: "/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/2BBE7B3F-C5A7-0D48-A384-FAD06B127FD8.root"
    # eos file: "/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9/2BBE7B3F-C5A7-0D48-A384-FAD06B127FD8.root"

    # DAS file: "/store/data/Run2018A/JetHT/NANOAOD/UL2018_MiniAODv2_NanoAODv9_GT36-v1/2820000/97F68EC0-0E12-C04C-A5D6-2B7A7C6688F8.root"
    # eos file: "/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9_GT36-v1/97F68EC0-0E12-C04C-A5D6-2B7A7C6688F8.root"
    
    if not fileName.startswith("/store/"):
        return fileName

    if useLocalFileIfExists:
        fileNameTemplate_DAS = "/store/{IsMC}/{SampleProductionCampaign}/{SampleName}/{DatasetTier}/{GT}/{SampleDir}/{SampleFileName}"
        fileNameTemplate_EOS = "/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/QCD_HT500to700_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9/2BBE7B3F-C5A7-0D48-A384-FAD06B127FD8.root"

        r_ = parse(fileNameTemplate_DAS, fileName)
        IsMC = 'MC' if 'mc' in r_['IsMC'].lower() else 'data'
        SampleProductionCampaign = r_['SampleProductionCampaign']
        SampleName = r_['SampleName']
        DatasetTier = r_['DatasetTier']
        GT = r_['GT']
        SampleDir = r_['SampleDir']
        SampleFileName = r_['SampleFileName']
        Era = None
        if 'UL16' in SampleProductionCampaign:
            Era = '2016'
        elif 'UL17' in SampleProductionCampaign:
            Era = '2017'
        elif 'UL18' in SampleProductionCampaign:
            Era = '2018'
           
        fileName_EOS = f"/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/{Era}/{IsMC}/{SampleName}/{SampleProductionCampaign}/{SampleFileName}"
    
    redirector_toUse = None
    for redirector in xrootd_redirectorNames:
        print(f"setXRootDRedirector():: Checking {redirector + fileName}"); sys.stdout.flush()
        '''
        with uproot.open(redirector + fileName) as file1:
            #print(f"\n{redirector + fileName}: file1.keys(): {file1.keys()}")
            #print(f"\n{redirector + fileName}: file1.keys(): {file1['Events'].numentries}"); sys.stdout.flush()

            #if file1['Events'].num_entries > 0:
            if file1['Events'].numentries > 0:
                print(f"{redirector + fileName}: file1.keys(): {file1['Events'].numentries}"); sys.stdout.flush()
                redirector_toUse = redirector
                break
        '''
        file1 = None
        try:
            file1 = uproot.open(redirector + fileName)
        except:
            print(f"setXRootDRedirector():: File open {redirector + fileName} failed"); sys.stdout.flush()
        else:
            #print(f"\n{redirector + fileName}: file1.keys(): {file1.keys()}")
            #print(f"\n{redirector + fileName}: file1.keys(): {file1['Events'].numentries}"); sys.stdout.flush()

            nEntries = file1['Events'].numentries
            file1.close()
            #if file1['Events'].num_entries > 0:
            if nEntries > 0:
                print(f"{redirector + fileName}: {nEntries}"); sys.stdout.flush()
                redirector_toUse = redirector
                break
            
    #print(f"redirector_toUse: {redirector_toUse}")
    
    return redirector_toUse + fileName

def xrdcpFile(sFileName, sFileNameLocal, nTry = 3, cp_command = 'xrdcp'):
    command_ = "time %s %s %s" % (cp_command, sFileName, sFileNameLocal)
    command_list_ = command_.split(" ")
    dirName_ = os.path.dirname(sFileNameLocal)
    os.makedirs(dirName_, exist_ok=True)
    print(f"{dirName_ = }: {os.path.exists(dirName_)} ")
    print(f"{command_ = }")    
    for iTry in range(nTry):
        print(f"htoaa_CommonTools::xrdcpFile() here1 {iTry = } {datetime.now() = }"); sys.stdout.flush()
        process = subprocess.Popen(command_list_,
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True
                                   )
        stdout, stderr = process.communicate()
        print(f"  {iTry = } {stdout = }, {stderr = }");  sys.stdout.flush()
        stderr_lower = stderr.lower()
        fileSize = 0
        if os.path.exists(sFileNameLocal):
            try:
                fileSize = os.path.getsize(sFileNameLocal) / (1024 * 1024) # file size in MB
                print(f"sFileNameLocal: {sFileNameLocal} ({fileSize} MB) ")
            except  FileNotFoundError:
                print(f"sFileNameLocal: {sFileNameLocal} file not found.")
            except OSError: 
                print(f"sFileNameLocal: {sFileNameLocal} OS error occurred.")
        #if 'FATAL' not in stderr and 'ERROR' not in stderr : # download was successful
#        if ('fatal' not in stderr_lower and 'error' not in stderr_lower) or \
#            (os.path.exists(sFileNameLocal) and fileSize > NanoAODFileSize_Min) : # download was successful
        if ('fatal' not in stderr_lower and 'error' not in stderr_lower):
            print(f"xrdcpFile():: successful")
            return True

    return False
    


def selectRunLuminosityBlock_ApprochEventBase(golden_json_path, runNumber_list, luminosityBlock_list):
    print(f" selectRunLuminosityBlock "); sys.stdout.flush();
    with open(golden_json_path) as fDataGoldenJSON:
        dataLSSelGoldenJSON = json.load(fDataGoldenJSON)    
        #print(f"0 : {dataLSSelGoldenJSON = }")
        dataLSSelGoldenJSON = {int(k): v for k, v in dataLSSelGoldenJSON.items()}

    #mask_run_ls = []
    mask_run_ls = np.full(len(runNumber_list), False, dtype=bool)
    #print(f"zip(runNumber_list, luminosityBlock_list) ({type(zip(runNumber_list, luminosityBlock_list))}) ({len(zip(runNumber_list, luminosityBlock_list))}): {zip(runNumber_list, luminosityBlock_list)} ")
    for idx_, (r,ls) in enumerate(zip(runNumber_list, luminosityBlock_list)):
        #mask_ = False
        if r in dataLSSelGoldenJSON:
            for ls_range in dataLSSelGoldenJSON[r]:
                if ls >= ls_range[0] and ls <= ls_range[1]:
                    #mask_ = True
                    mask_run_ls[idx_] = True
                    break
        #mask_run_ls.append(mask_)

    #print(f"mask_run_ls ({type(mask_run_ls)}) ({len(mask_run_ls)}){' '*6}: {mask_run_ls}")
    return mask_run_ls

def selectRunLuminosityBlock(dataLSSelGoldenJSON, runNumber_list, luminosityBlock_list):
    mask_run_ls = np.full(len(runNumber_list), False, dtype=bool)
    #print(f"zip(runNumber_list, luminosityBlock_list) ({type(zip(runNumber_list, luminosityBlock_list))}) ({len(zip(runNumber_list, luminosityBlock_list))}): {zip(runNumber_list, luminosityBlock_list)} ")
    for idx_, (r,ls) in enumerate(zip(runNumber_list, luminosityBlock_list)):
        #mask_ = False
        if r in dataLSSelGoldenJSON:
            for ls_range in dataLSSelGoldenJSON[r]:
                if ls >= ls_range[0] and ls <= ls_range[1]:
                    #mask_ = True
                    mask_run_ls[idx_] = True
                    break
    return mask_run_ls


def selectMETFilters(flags_list, era, isMC):
    sFLagDataOrMC = "MC" if isMC else "Data"

    mask_METFilters = np.full(len(flags_list), True, dtype=bool)

    if "goodVertices" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.goodVertices
    
    if "globalSuperTightHalo2016Filter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.globalSuperTightHalo2016Filter
    
    if "HBHENoiseFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.HBHENoiseFilter
    
    if "HBHENoiseIsoFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.HBHENoiseIsoFilter
    
    if "EcalDeadCellTriggerPrimitiveFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.EcalDeadCellTriggerPrimitiveFilter
    
    if "BadPFMuonFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.BadPFMuonFilter
    
    if "BadPFMuonDzFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.BadPFMuonDzFilter
    
    if "hfNoisyHitsFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.hfNoisyHitsFilter
    
    if "eeBadScFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.eeBadScFilter
    
    if "ecalBadCalibFilter" in MET_Filters[era][sFLagDataOrMC]:
        mask_METFilters = mask_METFilters & flags_list.ecalBadCalibFilter

    return mask_METFilters
    

def selectAK4Jets(Jets, era, pT_Thsh=0):
    # Not sure what to refer?
    #   1) https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVUL
    #   2) https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Jets

    '''
    # 1) https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVUL
    if '2018' in era or '2017' in era:
        # AK4CHS jets
        maskJetsSelected_HB = (
            (abs(Jets.eta)      <= 2.6) & 
            (Jets.neHEF         < 0.90) & # Neutral Hadron Fraction: Jet_neHEF	Float_t	neutral Hadron Energy Fraction
            (Jets.neEmEF        < 0.90) & # Neutral EM Fraction: Jet_neEmEF	Float_t	neutral Electromagnetic Energy Fraction
            (Jets.nConstituents > 1)    & # Number of Constituents: Jet_nConstituents	UChar_t	Number of particles in the jet
            (Jets.muEF          < 0.80) & # Muon Fraction: Jet_muEF	Float_t	muon Energy Fraction
            (Jets.chHEF         > 0)    & # Charged Hadron Fraction: Jet_chHEF	Float_t	charged Hadron Energy Fraction
            #(                     )    & # Charged Multiplicity > 0 ??
            (Jets.chEmEF        < 0.80)   # Charged EM Fraction: Jet_chEmEF	Float_t	charged Electromagnetic Energy Fraction
        )
        maskJetsSelected_HE1 = (
            (abs(Jets.eta)      > 2.6)  & (abs(Jets.eta)      <= 2.7) & 
            (Jets.neHEF         < 0.90) & # Neutral Hadron Fraction: Jet_neHEF	Float_t	neutral Hadron Energy Fraction
            (Jets.neEmEF        < 0.99) & # Neutral EM Fraction: Jet_neEmEF	Float_t	neutral Electromagnetic Energy Fraction
            (Jets.muEF          < 0.80) & # Muon Fraction: Jet_muEF	Float_t	muon Energy Fraction
            #() & # Charged Multiplicity > 0 ??
            (Jets.chEmEF        < 0.80)   # Charged EM Fraction: Jet_chEmEF	Float_t	charged Electromagnetic Energy Fraction
        )
        maskJetsSelected_HE2 = (
            (abs(Jets.eta)      > 2.7)  & (abs(Jets.eta)      <= 3.0) & 
            (Jets.neEmEF        > 0.01) & (Jets.neEmEF        < 0.99)   # Neutral EM Fraction: Jet_neEmEF	Float_t	neutral Electromagnetic Energy Fraction
            #() & # Number of Neutral Particles: 
        )
        maskJetsSelected_HF = (
            (abs(Jets.eta)      > 3.0)  & (abs(Jets.eta)      <= 5.0) & 
            (Jets.neHEF         < 0.20) & # Neutral Hadron Fraction: Jet_neHEF	Float_t	neutral Hadron Energy Fraction
            (Jets.neEmEF        < 0.90)   # Neutral EM Fraction: Jet_neEmEF	Float_t	neutral Electromagnetic Energy Fraction
            #() & # Number of Neutral Particles: 
        )
        maskJetsSelected = ( maskJetsSelected_HB | maskJetsSelected_HE1 | maskJetsSelected_HE2 | maskJetsSelected_HF )
    '''

    # 2) https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Jets
    maskJetsSelected = (
        (Jets.jetId >= 6) & 
        ( (Jets.pt > 50) | (Jets.puId >= 4 ) )
    )

    return Jets[maskJetsSelected & (Jets.pt > pT_Thsh)]
        

def selectMuons(eventsObj, pT_Thsh=10, MVAId=3, MiniIsoId=3, MVATTHThsh=0.5):
    # MuonMVAId    : (1=MvaLoose, 2=MvaMedium, 3=MvaTight, 4=MvaVTight, 5=MvaVVTight)
    # MuonMiniIsoId: (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)
    maskSelMuons = (
        (eventsObj.pt > pT_Thsh) &
        (abs(eventsObj.eta) < 2.4) & 
        (eventsObj.mvaId >= MVAId) &
        (eventsObj.miniIsoId >= MiniIsoId) & 
        (eventsObj.mvaTTH > MVATTHThsh)
    )
    return eventsObj[maskSelMuons]
    

def selectElectrons(eventsObj, pT_Thsh=10, MVAId='mvaFall17V2Iso_WP80', MVATTHThsh=0.3):
    # ElectronMVAId: 'mvaFall17V2Iso_WP80', 'mvaFall17V2Iso_WP90' 'mvaFall17V2Iso_WPL'
    maskSelElectrons = (
        (eventsObj.pt > pT_Thsh) &
        (abs(eventsObj.eta) < 2.3) & 
        (eventsObj[MVAId] > 0) &
        (eventsObj.mvaTTH > MVATTHThsh)
    )
    return eventsObj[maskSelElectrons]


def getLumiScaleForPhSpOverlapRewgtMode(
        hLumiScale ,
        sample_category ,
        sample_HT_value ,
        mask_PhSp_dict ):
    lumiScale = None
    
    if 'QCD' in sample_category:
        xBin = None
        for iBin in range(len(hLumiScale.axes[0])):
            if sample_HT_value >= hLumiScale.axes[0][iBin][0] and sample_HT_value < hLumiScale.axes[0][iBin][1]:
                xBin = iBin
                break

        nEvents = len(list(mask_PhSp_dict.values())[0])
        lumiScale = np.ones(nEvents)
        for PhSpName, mask_PhSp in mask_PhSp_dict.items():
            lumiScale_value = hLumiScale[xBin, PhSpName]
            lumiScale = np.where(mask_PhSp, np.ones(nEvents) * lumiScale_value, lumiScale)
            #print(f'htoaa_CommonTools::getLumiScaleForPhSpOverlapRewgtMode():: {sample_category = }, {sample_HT_value = }, {xBin = }, {PhSpName = }, {lumiScale_value = }, ')
                            
    else:
        logging.error(f'htoaa_CommonTools::getLumiScaleForPhSpOverlapRewgtMode():: {sample_category = } not implemented')
        exit(0)

    return lumiScale



def getTopPtRewgt(eventsGenPart, isPythiaTuneCP5):
    '''
    Top pT reweights calculated using NNLO calculation and POWHEG TuneCP5 sample. Additional corrections for POWHEG TuneCUETP samples.
    Top pT reweights = sqrt( SF(top pT, TuneCP5) * SF(top pT, TuneCUETP) * SF(anti-top pT, TuneCP5) * SF(anti-top pT, TuneCUETP) )
    Ref: https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf#page=12
    '''
    #printVariable('eventsGenPart.pt', eventsGenPart.pt)
    #print(f"{isPythiaTuneCP5 = }")

    list_pT_top_antiTop = []
    list_pT_top_antiTop.append( ak.firsts( eventsGenPart[(eventsGenPart.pdgId ==       PDGID_TopQuark )].pt ) ) # Convert events.pt: [[1.1], [1.09], ...] -> ak.firsts() -> [1.1, 1.09, ..]
    list_pT_top_antiTop.append( ak.firsts( eventsGenPart[(eventsGenPart.pdgId == (-1 * PDGID_TopQuark))].pt ) )

    fitRangeMin         = Corrections["TopPtRewgt"]["TuneCP5"]["FitRange"][0]
    fitRangeMax         = Corrections["TopPtRewgt"]["TuneCP5"]["FitRange"][1]   
    sFitFunctionFormat  = Corrections["TopPtRewgt"]["TuneCP5"]["FitFunctionFormat"]
    sFitFunction        = Corrections["TopPtRewgt"]["TuneCP5"]["FitFunction"] 
    if sFitFunctionFormat == "exp( {a} + ({b} * x) + ({c} * x * x) + ({d}/(x + {e})) )":
        fitResult_        = parse(sFitFunctionFormat, sFitFunction) # https://pypi.org/project/parse/
        pTuneCP5_a = float(fitResult_['a'])
        pTuneCP5_b = float(fitResult_['b'])
        pTuneCP5_c = float(fitResult_['c'])
        pTuneCP5_d = float(fitResult_['d'])
        pTuneCP5_e = float(fitResult_['e'])
        #print(f"sFitFunction: {sFitFunction}: a {pTuneCP5_a}, b {pTuneCP5_b}, c {pTuneCP5_c}, d {pTuneCP5_d}, e {pTuneCP5_e} ")
    else:
        logging.critical(f'htoaa_CommonTools::getTopPtRewgt():: The current TopPtReweigting function is not implemented.. Fix it..')
        exit(0)

    if not isPythiaTuneCP5:
        if not math.isclose(Corrections["TopPtRewgt"]["TuneCUETP"]["FitRange"][0], fitRangeMin, rel_tol=1e-5):
            logging.critical(f'htoaa_CommonTools::getTopPtRewgt():: {Corrections["TopPtRewgt"]["TuneCP5"]["FitRange"][0] = } is not equal to {Corrections["TopPtRewgt"]["TuneCUETP"]["FitRange"][0]}. Code is not able to handle it... Fix it..')
            exit(0)
        if not math.isclose(Corrections["TopPtRewgt"]["TuneCUETP"]["FitRange"][1], fitRangeMax, rel_tol=1e-5):
            logging.critical(f'htoaa_CommonTools::getTopPtRewgt():: {Corrections["TopPtRewgt"]["TuneCP5"]["FitRange"][1] = } is not equal to {Corrections["TopPtRewgt"]["TuneCUETP"]["FitRange"][1]}. Code is not able to handle it... Fix it..')
            exit(0)    
        sFitFunctionFormat    = Corrections["TopPtRewgt"]["TuneCUETP"]["FitFunctionFormat"]
        sFitFunction          = Corrections["TopPtRewgt"]["TuneCUETP"]["FitFunction"]
        if sFitFunctionFormat == "{a} + ({b} * TanH({c} + ({d} * x) )":
            fitResult_        = parse(sFitFunctionFormat, sFitFunction) # https://pypi.org/project/parse/
            pTuneCUETP_a = float(fitResult_['a'])
            pTuneCUETP_b = float(fitResult_['b'])
            pTuneCUETP_c = float(fitResult_['c'])
            pTuneCUETP_d = float(fitResult_['d'])
            #print(f"sFitFunction: {sFitFunction}: a {pTuneCUETP_a}, b {pTuneCUETP_b}, c {pTuneCUETP_c}, d {pTuneCUETP_d} ")
        else:
            logging.critical(f'htoaa_CommonTools::getTopPtRewgt():: The current TopPtReweigting function is not implemented.. Fix it..')
            exit(0)



    #wgt_TopPtRewgt =  np.ones(len(eventsGenPart))
    wgt_TopPtRewgt =  np.full(len(eventsGenPart), 1.0)
    pTMax_list = np.full(len(eventsGenPart), fitRangeMax )
    #printVariable('wgt_TopPtRewgt ', wgt_TopPtRewgt)
    for pT_list in list_pT_top_antiTop: # loop over top pT and antitop pT
        #printVariable('pT_list ', pT_list)
        #printVariable('ak.firsts(pT_list) ', ak.firsts(pT_list))

        if Corrections["TopPtRewgt"]["TuneCP5"]["FitFunctionFormat"] == "exp( {a} + ({b} * x) + ({c} * x * x) + ({d}/(x + {e})) )":
            wgt_withinFitRange  = np.exp( pTuneCP5_a + (pTuneCP5_b * pT_list   ) + (pTuneCP5_c * pT_list    * pT_list   ) + (pTuneCP5_d / (pT_list    + pTuneCP5_e)) )
            wgt_outsideFitRange = np.exp( pTuneCP5_a + (pTuneCP5_b * pTMax_list) + (pTuneCP5_c * pTMax_list * pTMax_list) + (pTuneCP5_d / (pTMax_list + pTuneCP5_e)) )
            
            wgt_TopPtRewgt = wgt_TopPtRewgt * ak.where(
                np.less_equal(pT_list, fitRangeMax),
                wgt_withinFitRange,
                wgt_outsideFitRange
            ) 
        #printVariable('wgt_TopPtRewgt ', wgt_TopPtRewgt)

        if not isPythiaTuneCP5:
            if Corrections["TopPtRewgt"]["TuneCUETP"]["FitFunctionFormat"] == "{a} + ({b} * TanH({c} + ({d} * x) )":
                wgt_withinFitRange  = pTuneCUETP_a + (pTuneCUETP_b * np.tanh( pTuneCUETP_c + (pTuneCUETP_d * pT_list   ) ))
                wgt_outsideFitRange = pTuneCUETP_a + (pTuneCUETP_b * np.tanh( pTuneCUETP_c + (pTuneCUETP_d * pTMax_list) ))

                wgt_TopPtRewgt = wgt_TopPtRewgt * ak.where(
                    np.less_equal(pT_list, fitRangeMax),
                    wgt_withinFitRange,
                    wgt_outsideFitRange
                )
        #printVariable('wgt_TopPtRewgt ', wgt_TopPtRewgt)

    wgt_TopPtRewgt = np.sqrt(wgt_TopPtRewgt)
    #printVariable('wgt_TopPtRewgt ', wgt_TopPtRewgt)
    return wgt_TopPtRewgt


def getPURewgts(PU_list, hPURewgt):
    # hPURewgt: Hist() object

    wgt_PU = np.ones(len(PU_list))
    for iBin in range(len(hPURewgt.values())):
        xBin_edge_low = hPURewgt.axes[0][iBin][0]
        xBin_edge_up  = hPURewgt.axes[0][iBin][1]
        xBin_value    = hPURewgt.values()[iBin]

        wgt_PU = ak.where(
            np.logical_and( np.greater_equal(PU_list, xBin_edge_low), np.less(PU_list, xBin_edge_up) ),
            np.full(len(PU_list), xBin_value),
            wgt_PU
        )

    #print(f"PU_list ({len(PU_list)}): {PU_list}")
    #print(f"wgt_PU ({len(wgt_PU)}): {wgt_PU}")
    return wgt_PU


def getPURewgts_variation(events, year):
    ## json files from: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/LUM
    fname = "data/correction/mc/PURewgt/{0}_UL/puWeights.json.gz".format(year)
    hname = {
        "2016APV": "Collisions16_UltraLegacy_goldenJSON",
        "2016"   : "Collisions16_UltraLegacy_goldenJSON",
        "2017"   : "Collisions17_UltraLegacy_goldenJSON",
        "2018"   : "Collisions18_UltraLegacy_goldenJSON"
    }
    evaluator = correctionlib.CorrectionSet.from_file(fname)

    puUp = evaluator[hname[str(year)]].evaluate(np.array(events.Pileup.nTrueInt), "up")
    puDown = evaluator[hname[str(year)]].evaluate(np.array(events.Pileup.nTrueInt), "down")
    puNom = evaluator[hname[str(year)]].evaluate(np.array(events.Pileup.nTrueInt), "nominal")

    return [puNom, puUp, puDown]


def getHiggsPtRewgtForGGToHToAATo4B(GenHiggsPt_list): # GenHiggsPt_list
    # Used in Brook's analysis
    #wgt_HiggsPt = (3.9 - (0.4 * np.log2(pT)))
    #wgt_HiggsPt = np.maximum(wgt_HiggsPt, np.full(len(pT), 0.1) )

    # https://indico.cern.ch/event/1348321/#19-siddhesh-sawant
    # min(max(1.45849 + -0.00400668*x + 4.02577e-06*pow(x, 2) + -1.38804e-09*pow(x, 3), 0.09), 1.02)
    wgt_HiggsPt = 1.45849 - 0.00400668*GenHiggsPt_list + 4.02577e-06*GenHiggsPt_list**2 - 1.38804e-09*GenHiggsPt_list**3 
    wgt_HiggsPt = np.maximum(wgt_HiggsPt, np.full(len(GenHiggsPt_list), 0.09) )
    wgt_HiggsPt = np.minimum(wgt_HiggsPt, np.full(len(GenHiggsPt_list), 1.02) )
    return wgt_HiggsPt


def getHTReweight(HT_list, sFitFunctionFormat, sFitFunction, sFitFunctionRange):
    wgt_HT = None
    
    # 'Corrections' variable defined in htoaa_Settings
    if sFitFunctionFormat == "{p0} + ({p1} * (x - {HTBinMin}))":
        fitResult_        = parse(sFitFunctionFormat, sFitFunction) # https://pypi.org/project/parse/
        fitParam_p0       = float(fitResult_['p0'])
        fitParam_p1       = float(fitResult_['p1'])
        fitConst_HTBinMin = float(fitResult_['HTBinMin'])

        # for e.g. HT500to700
        fitResult_  = parse("HT{FitRangeMin}to{FitRangeMax}", sFitFunctionRange)
        fitRangeMin = float(fitResult_['FitRangeMin'])
        fitRangeMax = float(fitResult_['FitRangeMax'])

        wgt_HT_cal = fitParam_p0 + (fitParam_p1 * (HT_list - fitConst_HTBinMin))
        wgt_HT = np.where(
            np.logical_and( np.greater_equal(HT_list, fitRangeMin), np.less(HT_list, fitRangeMax) ),
            wgt_HT_cal,
            np.ones(len(HT_list))
        )

    else:
        print(f'htoaa_CommonTools.py::getHTReweight():: {Corrections["HTRewgt"]["QCD_bGen"]["FitFunctionFormat"] = } is not implemented \t\t **** ERROR **** \n')
        exit(0)

    return wgt_HT
    


def get_PSWeight(events, dataset):
    """
    Parton Shower Weights (FSR and ISR)
    "Default" variation: https://twiki.cern.ch/twiki/bin/view/CMS/HowToPDF#Which_set_of_weights_to_use
    i.e. scaling ISR up and down

    PS weights (w_var / w_nominal);   [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2; [2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5
    """
    nweights = len(events)
    nom = np.ones(nweights)

    up_isr   = np.ones(nweights)
    down_isr = np.ones(nweights)
    up_fsr   = np.ones(nweights)
    down_fsr = np.ones(nweights)

    if hasattr(events, 'PSWeight') and "HToAATo4B" in dataset:
        if len(events.PSWeight[0]) == 4:
            up_isr   = events.PSWeight[:, 0]  # ISR=2, FSR=1
            down_isr = events.PSWeight[:, 2]  # ISR=0.5, FSR=1

            up_fsr   = events.PSWeight[:, 1]  # ISR=1, FSR=2
            down_fsr = events.PSWeight[:, 3]  # ISR=1, FSR=0.5

        elif len(events.PSWeight[0]) > 1:
            print("PS weight vector has length ", len(events.PSWeight[0]))

    return [nom, up_isr, down_isr, up_fsr, down_fsr]


def add_pdf_as_weight(events, pdf_weights,dataset):


    nom = np.ones(len(events))
    up_pdfas   = up_aS   = up_pdf    = np.ones(len(events))
    down_pdfas = down_aS = down_pdf  = np.ones(len(events))

    docstring = pdf_weights.__doc__

    # NNPDF31_nnlo_as_0118_nf_4_mc_hessian
    # https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_as_0118_nf_4_mc_hessian/NNPDF31_nnlo_as_0118_nf_4_mc_hessian.info
    # if True: #LHA IDs "325500 - 325600" in docstring:
    # Hessian PDF weights
    # Eq. 21 of https://arxiv.org/pdf/1510.03865v1.pdf                                   
    #print (" no. of PDF column  ",len(pdf_weights[0]))
    if "HToAATo4B" in dataset:
        arg = pdf_weights[:,1:]-np.ones((len(events),100)) #np.ones((len(events),100))
        summed = ak.sum(np.square(arg),axis=1)
        #pdf_unc = np.sqrt( (1./99.) * summed )
        pdf_unc = np.sqrt( summed )
        up_pdf   = nom + pdf_unc
        down_pdf = nom - pdf_unc

    #anther pdf unc definition 
    #pdfUnc = ak.std(events.LHEPdfWeight,axis=1)/ak.mean(events.LHEPdfWeight,axis=1) 
    #pdfUnc = ak.fill_none(pdfUnc, 0.00)
    #up_pdf = nom + pdfUnc
    #down_pdf = nom - pdfUnc

    # alpha_S weights
    # Eq. 27 of same ref
    #as_unc = 0.5*(pdf_weights[:,102] - pdf_weights[:,101])
    #up_pdf   = nom + as_unc
    #down_pdf = nom - as_unc
    
    
    # PDF + alpha_S weights
    # Eq. 28 of same ref
    #pdfas_unc = np.sqrt( np.square(pdf_unc) + np.square(as_unc) )
    #weights.add('PDFaS_weight', nom, pdfas_unc + nom) 
    #up_pdfas   = nom + pdfas_unc
    #down_pdfas = nom - pdfas_unc

    return [nom, up_pdf, down_pdf]#, up_aS, down_aS, up_pdfas, down_pdfas, up_pdfas, down_pdfas]


def get_QCDScaleWeight(events, dataset):
    nEvents = len(events)
    nom  = renorm_up = renorm_down = factr_up = factr_down = np.ones(nEvents)

    if hasattr(events, 'LHEScaleWeight') and "HToAATo4B" in dataset:
        if len(events.LHEScaleWeight[0]) == 9:
            # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv9/2018UL/doc_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1.html#LHEPdfWeight
            # LHEScaleWeight	Float_t	LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0
            # [1] is renscfact=0.5d0 facscfact=1d0.      [7] is renscfact=2d0 facscfact=1d0
            # [3] is renscfact=1d0 facscfact=0.5d0.      [5] is renscfact=1d0 facscfact=2d0 ;
            # renorm_up = 1, down = 7.   fact_up = 3, down = 5
            renorm_up   = events.LHEScaleWeight[:, 1]
            renorm_down = events.LHEScaleWeight[:, 7]
            factr_up    = events.LHEScaleWeight[:, 3]
            factr_down  = events.LHEScaleWeight[:, 5]
        
        elif len(events.nLHEScaleWeight[0]) > 1:
            print("LHEScaleWeight  vector has length ", len(events.nLHEScaleWeight[0]))
            
    return [nom, renorm_up, renorm_down, factr_up, factr_down]


def add_HiggsEW_kFactors(genHiggs, dataset):

    hew_kfactors = correctionlib.CorrectionSet.from_file("data/EWHiggsCorrection/EWHiggsCorrections.json")
    def get_hpt():
        boson = ak.firsts(genHiggs[
            (genHiggs.pdgId == 25)
            & genHiggs.hasFlags(["fromHardProcess", "isLastCopy"])
        ])
        return np.array(ak.fill_none(boson.pt, 0.))

    if "VBF" in dataset:
        hpt = get_hpt()
        ewkcorr = hew_kfactors["VBF_EW"]
        ewknom = ewkcorr.evaluate(hpt)
        return  "VBF_EW", ewknom

    elif "WH" in dataset or "ZH" in dataset:
        hpt = get_hpt()
        ewkcorr = hew_kfactors["VH_EW"]
        ewknom = ewkcorr.evaluate(hpt)
        return "VH_EW", ewknom
    

    elif "ttH" in dataset:
        hpt = get_hpt()
        ewkcorr = hew_kfactors["ttH_EW"]
        ewknom = ewkcorr.evaluate(hpt)
        return "ttH_EW", ewknom
    else :
        return None


def get_JER_and_JES(events, FatJets, year, shift_syst=""):
    #UL2018 -> (19UL18_V5 , 19UL18_JRV2) / UL17 -> (19UL17_V5, 19UL17_JRV2) / UL2016APV -> (19UL16APV_V7, 20UL16APV_JRV3) / UL2016 -> (19UL16_V7, 20UL16_JRV3)  
    #UL17 https://cms-talk.web.cern.ch/t/ak8-jets-jec-for-summer19ul17-mc/23154/8
    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC#Recommended_for_MC
    #https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    FatJets["pt_raw"], FatJets["mass_raw"] = (1 - FatJets.rawFactor) * FatJets.pt, (1 - FatJets.rawFactor) * FatJets.mass
    FatJets['pt_gen'] = ak.values_astype(ak.fill_none(FatJets.matched_gen.pt, 0), np.float32)
    FatJets['rho'] = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, FatJets.pt)[0]
    events_cache = events.caches[0]

    Jetext = extractor()
    Jetext.add_weight_sets([
        f"* * data/JERS/{year}UL_V_MC_L1FastJet_AK8PFPuppi.jec.txt",
        f"* * data/JERS/{year}UL_V_MC_L2Relative_AK8PFPuppi.jec.txt",
        f"* * data/JERS/{year}UL_V_MC_Uncertainty_AK8PFPuppi.junc.txt",
        f"* * data/JERS/{year}UL_JR_MC_PtResolution_AK8PFPuppi.jr.txt",
        f"* * data/JERS/{year}UL_JR_MC_SF_AK8PFPuppi.jersf.txt",
    ])
    Jetext.finalize()
    Jetevaluator = Jetext.make_evaluator()

    jec_names = [f"{year}UL_V_MC_L1FastJet_AK8PFPuppi", f"{year}UL_V_MC_L2Relative_AK8PFPuppi",
                 f"{year}UL_V_MC_Uncertainty_AK8PFPuppi", f"{year}UL_JR_MC_PtResolution_AK8PFPuppi",
                 f"{year}UL_JR_MC_SF_AK8PFPuppi"]
    jec_stack = JECStack({name: Jetevaluator[name] for name in jec_names})
    
    name_map = jec_stack.blank_name_map
    name_map.update({"JetPt": "pt", "JetMass": "mass", "JetEta": "eta", "JetA": "area",
                     "ptGenJet": "pt_gen", "ptRaw": "pt_raw", "massRaw": "mass_raw", "Rho": "rho"})
    
    corrected_jets = CorrectedJetsFactory(name_map, jec_stack).build(FatJets, lazy_cache=events.caches[0])
    
    if shift_syst == "JERUp":
        FatJets = corrected_jets.JER.up
    elif shift_syst == "JERDown":
        FatJets = corrected_jets.JER.down
    elif shift_syst == "JESUp":
        FatJets = corrected_jets.JES_jes.up
    elif shift_syst == "JESDown":
        FatJets = corrected_jets.JES_jes.down
    else:
        # either nominal or some shift systematic unrelated to jets
        FatJets = corrected_jets

    return FatJets

def get_JMR_JMS(Jet, year, shift_syst=""):
    # jet mass https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsDP23044

    substr_cset = correctionlib.CorrectionSet.from_file("data/jms/Substructure_jmssf.json")

    jet_pt   = Jet.pt #np.array(ak.fill_none(Jet.pt, 0.))

    jms_nom  = substr_cset[f"jmssf_{year}"].evaluate(jet_pt,"")
    jms_up   = substr_cset[f"jmssf_{year}"].evaluate(jet_pt,"up")
    jms_down = substr_cset[f"jmssf_{year}"].evaluate(jet_pt,"down")

    mass = Jet.msoftdrop

    corrected_mass_up   = mass * jms_up
    corrected_mass_down = mass * jms_down
    corrected_mass_nomi = mass * jms_nom

    for index, value in enumerate(corrected_mass_up):
        if value < corrected_mass_nomi[index] and value > 50 :
            print("corrected_mass_up is less than corrected_mass_nomi = ", corrected_mass_nomi[index], " corrected_mass_up = ", corrected_mass_up[index])
    for index, value in enumerate(corrected_mass_nomi):
        if value < corrected_mass_down[index] and value > 50 :
            print("corrected_mass_nomi is less than corrected_mass_down = ", corrected_mass_nomi[index], " corrected_mass_down = ", corrected_mass_down[index])

    if shift_syst == "JMSUp":
        corrected_mass = mass * jms_up
    elif shift_syst == "JMSDown":
        corrected_mass = mass * jms_down
    else:
        corrected_mass = mass * jms_nom
    return corrected_mass

#smearing = np.random.normal(mass[:,])
# scale to JMR nom, down, up (minimum at 0)
#jmr_central, jmr_down, jmr_up = (
    #((smearing * max(jmrValues[year][i] - 1, 0)) + 1) for i in range(3)
    #)



def get_jetTriggerSF(events, year, selection):

    leadingjet = ak.firsts(events.FatJet)
    jet_triggerSF = correctionlib.CorrectionSet.from_file("data/correction/mc/TrgEffSF/fatjet_triggerSF.json") # correctionlib.CorrectionSet.from_file("data/trigger/fatjet_triggerSF.json")

    def mask(w):
        return np.where(selection.all('JetID'), w, 1.)

    # Same for 2016 and 2016APV
    if '2016' in year:
        year = '2016'

    jet_pt   = np.array(ak.fill_none(leadingjet.pt, 0.))
    jet_msd  = np.array(ak.fill_none(leadingjet.msoftdrop, 0.))  # note: uncorrected
    nom_trg  = mask(jet_triggerSF[f'fatjet_triggerSF{year}'].evaluate("nominal", jet_pt, jet_msd))
    up_trg   = mask(jet_triggerSF[f'fatjet_triggerSF{year}'].evaluate("stat_up", jet_pt, jet_msd))
    down_trg = mask(jet_triggerSF[f'fatjet_triggerSF{year}'].evaluate("stat_dn", jet_pt, jet_msd))

    return [nom_trg, up_trg, down_trg]


    
def selGenPartsWithStatusFlag(GenPart_StatusFlags_list, statusFlag_toSelect):  
    # Check if statusFlag_toSelect th bit is 1 in binary version of GenPart_StatusFlags
    return ( GenPart_StatusFlags_list & (2 ** int(statusFlag_toSelect)) ) > 0  


def DfColLabel_convert_bytes_to_string(df):
    cols_rename = {}
    for col in df.columns:
        if isinstance(col, (bytes, bytearray)):
            cols_rename[col] = col.decode()
    print("DfColLabel_convert_bytes_to_string:: cols_rename: {}".format(cols_rename))
    df.rename(columns=cols_rename, inplace=True)
    return df


def cut_ObjectMultiplicity(nObjects, nObjects_min=None, nObjects_max=None):
    '''
    Check if "nObjects are with nObjects_min and nObjects_max", if they are specified.
    If either of nObjects_min and nObjects_max are not specified, then corresponding condition is not checked.

    Return:
        True: if nObjects passes the condition.
        False: if nObjects fails the condition
    '''
    mask = mask_low = mask_up  = None
    if nObjects_min is not None: mask_low = (nObjects >= nObjects_min)    
    if nObjects_max is not None: mask_up  = (nObjects <= nObjects_max)

    if (nObjects_min is not None) and (nObjects_max is not None):
        mask = (mask_low and mask_up)
    elif (nObjects_min is not None):
        mask = mask_low
    else:
        mask = mask_up

    return mask



def cut_ObjectPt(objects_Pt, PtThrsh_Lead=None, PtThrsh_Sublead=None, PtThrsh_Third=None, PtThrsh_Fourth=None, PtThrsh_Fifth=None):
    '''
    Check Objects Pt is above their resepctive thresholds.
    If PtThrsh_<rank> is not set, then their Pt condition is not checked.

    Return:
        True: All objects' Pt is about respective threshold
        False: Else false
    '''
    print("objects_Pt ({}) : {}".format(type(objects_Pt),  objects_Pt))
    condition = True
    if                                    objects_Pt[0] < PtThrsh_Lead:     condition = False    
    if (PtThrsh_Sublead is not None) and (objects_Pt[1] < PtThrsh_Sublead): condition = False
    if (PtThrsh_Third   is not None) and (objects_Pt[2] < PtThrsh_Third):   condition = False
    if (PtThrsh_Fourth  is not None) and (objects_Pt[3] < PtThrsh_Fourth):  condition = False
    if (PtThrsh_Fifth   is not None) and (objects_Pt[4] < PtThrsh_Fifth):   condition = False

    return condition
    
def cut_ObjectPt_1(objects_Pt, PtThrshs):
    print("objects_Pt ({}) : {}, \t\t PtThrshs ({}) : {}".format(type(objects_Pt),  objects_Pt, type(PtThrshs), PtThrshs))

    return True

def cut_ObjectEta(objects_Eta, EtaThrsh, nObjects):
    '''
    Check Objects abs(Eta) is greater than thresholds set in EtaThrshs list.

    Return:
        True: All objects' Eta is about respective threshold
        False: Else false
    '''    
    condition = True
    for iObject in range(nObjects):
        if abs(objects_Eta[iObject]) > EtaThrsh:
            condition = False
            break
    return condition


## General python functions ----------------------------------------------------------------------------------

def GetDictFromJsonFile(filePath):
    # Lines starting with '#' are not read out, and also content between '/* .... */' are not read.
    # Content between " '''   ....  ''' " are not read
    # Source: https://stackoverflow.com/questions/29959191/how-to-parse-json-file-with-c-style-comments
    
    contents = ""
    fh = open(filePath)
    for line in fh:
        #cleanedLine = line.split("//", 1)[0]
        cleanedLine = line.split("#", 1)[0]
        if len(cleanedLine) > 0 and line.endswith("\n") and "\n" not in cleanedLine:
            cleanedLine += "\n"
        contents += cleanedLine
    fh.close
    
    #while "/*" in contents:
    #    preComment, postComment = contents.split("/*", 1)ß
    #    contents = preComment + postComment.split("*/", 1)[1]
    while "'''" in contents:
        preComment, postComment = contents.split("'''", 1)
        contents = preComment + postComment.split("'''", 1)[1]

    dictionary =  json.loads( contents )
    return dictionary


def executeBashCommand(sCmd1):
    cmd1 = None
    if isinstance(sCmd1, str):
        cmd1 = shlex.split(sCmd1)
    elif isinstance(sCmd1, list):
        cmd1 = sCmd1
    #print(f"executeBashCommand():: {sCmd1 = } \t {cmd1 = }")

    result = subprocess.run(
        cmd1,
        capture_output=True,
        text = True
    )
    print(f"{sCmd1}: ")
    print("result.stdout: %s" % (result.stdout))
    if result.stderr:
        #print(f"{result.stderr = }")
        print("result.stderr: %s" % (result.stderr))
        
    return result.stdout


def fillCoffeaHist(
        h = coffea_hist.Hist('tmp'),
        dataset = '',
        syst = None, 
        xValue = None,
        yValue = None,
        zValue = None,        
        wgt = None
):
    mask_ = ~ ak.is_none(xValue, axis=0) # do not fill 'None'
    kwargs = {
        'dataset'              : dataset,
        'systematic'           : syst,
        'weight'               : wgt[mask_],
        h.dense_axes()[0].name : xValue[mask_]
    } 
    if yValue:
        kwargs[h.dense_axes()[1].name] = yValue[mask_]
    if zValue:
        kwargs[h.dense_axes()[2].name] = zValue[mask_]
    h.fill( **kwargs )



def rebinTH1(h1_, nRebins):
    #print(f"rebinTH1():: histogram type {type(h1_) = },  {isinstance(h1_, hist.Hist) = }  ")
    if not (isinstance(h1_, hist.Hist) or  isinstance(h1_, coffea_hist.Hist)):
        print(f"rebinTH1():: histogram type {type(h1_)} not implemented... so could not rebin histogram ")
        return h1_
    
    if len(h1_.axes) != 1:
        print(f"rebinTH1:: histogram is not 1D")
        return h1_

    h1Rebin_ = None
    if   nRebins == 1:
        h1Rebin_ = h1_
    elif   nRebins == 2:
        h1Rebin_ = h1_[::2j]
    elif nRebins == 3:
        h1Rebin_ = h1_[::3j]
    elif nRebins == 4:
        h1Rebin_ = h1_[::4j]
    elif nRebins == 5:
        h1Rebin_ = h1_[::5j]
    elif nRebins == 6:
        h1Rebin_ = h1_[::6j]
    elif nRebins == 8:
        h1Rebin_ = h1_[::8j]
    elif nRebins == 10:
        h1Rebin_ = h1_[::10j]
    elif nRebins == 12:
        h1Rebin_ = h1_[::12j]        
    elif nRebins == 20:
        h1Rebin_ = h1_[::20j]
    elif nRebins == 40:
        h1Rebin_ = h1_[::40j]
    elif nRebins == 50:
        h1Rebin_ = h1_[::50j]
    elif nRebins == 100:
        h1Rebin_ = h1_[::100j]
        print("Rebin 100 <<<")
    else:
        print(f"nRebins={nRebins} is not yet implemented... Implement it \t\t **** ERROR ****")        
        
    #print(f"h1_ values ({len(h1_.values())}): {h1_.values()} \n variances ({len(h1_.variances())}): {h1_.variances()}")
    #print(f"h1Rebin_ values ({len(h1Rebin_.values())}): {h1Rebin_.values()} \n variances ({len(h1Rebin_.variances())}): {h1Rebin_.variances()}")
    if   nRebins > 1:    
        h1_ = h1Rebin_

    return h1_


def rebinTH2(h1_, nRebinX, nRebinY):
    #print(f"rebinTH1():: histogram type {type(h1_) = },  {isinstance(h1_, hist.Hist) = }  ")
    if not (isinstance(h1_, hist.Hist) or  isinstance(h1_, coffea_hist.Hist)):
        print(f"rebinTH1():: histogram type {type(h1_)} not implemented... so could not rebin histogram ")
        return h1_
    
    if len(h1_.axes) != 2:
        print(f"rebinTH1:: histogram is not 2D")
        return h1_

    h1Rebin_ = None
    if   nRebinX == 1:
        if   nRebinY == 1:
            h1Rebin_ = h1_
        elif nRebinY == 2:
            h1Rebin_ = h1_[::1j, ::2j]
        elif nRebinY == 3:
            h1Rebin_ = h1_[::1j, ::3j]
        elif nRebinY == 4:
            h1Rebin_ = h1_[::1j, ::4j]
        elif nRebinY == 5:
            h1Rebin_ = h1_[::1j, ::5j]
        elif nRebinY == 6:
            h1Rebin_ = h1_[::1j, ::6j]
        elif nRebinY == 10:
            h1Rebin_ = h1_[::1j, ::10j]
        elif nRebinY == 20:
            h1Rebin_ = h1_[::1j, ::20j]
        elif nRebinY == 40:
            h1Rebin_ = h1_[::1j, ::40j]
        elif nRebinY == 50:
            h1Rebin_ = h1_[::1j, ::50j]
        elif nRebinY == 100:
            h1Rebin_ = h1_[::1j, ::100j]
        else:
            print(f"{nRebinX = }, {nRebinY = } is not yet implemented... Implement it \t\t **** ERROR ****")        

    elif   nRebinX == 2:
        if   nRebinY == 1:
            h1Rebin_ = h1_
        elif nRebinY == 2:
            h1Rebin_ = h1_[::2j, ::2j]
        elif nRebinY == 3:
            h1Rebin_ = h1_[::2j, ::3j]
        elif nRebinY == 4:
            h1Rebin_ = h1_[::2j, ::4j]
        elif nRebinY == 5:
            h1Rebin_ = h1_[::2j, ::5j]
        elif nRebinY == 6:
            h1Rebin_ = h1_[::2j, ::6j]
        elif nRebinY == 10:
            h1Rebin_ = h1_[::2j, ::10j]
        elif nRebinY == 20:
            h1Rebin_ = h1_[::2j, ::20j]
        elif nRebinY == 40:
            h1Rebin_ = h1_[::2j, ::40j]
        elif nRebinY == 50:
            h1Rebin_ = h1_[::2j, ::50j]
        elif nRebinY == 100:
            h1Rebin_ = h1_[::2j, ::100j]
        else:
            print(f"{nRebinX = }, {nRebinY = } is not yet implemented... Implement it \t\t **** ERROR ****")        

    elif   nRebinX == 3:
        if   nRebinY == 1:
            h1Rebin_ = h1_
        elif nRebinY == 2:
            h1Rebin_ = h1_[::3j, ::2j]
        elif nRebinY == 3:
            h1Rebin_ = h1_[::3j, ::3j]
        elif nRebinY == 4:
            h1Rebin_ = h1_[::3j, ::4j]
        elif nRebinY == 5:
            h1Rebin_ = h1_[::3j, ::5j]
        elif nRebinY == 6:
            h1Rebin_ = h1_[::3j, ::6j]
        elif nRebinY == 10:
            h1Rebin_ = h1_[::3j, ::10j]
        elif nRebinY == 20:
            h1Rebin_ = h1_[::3j, ::20j]
        elif nRebinY == 40:
            h1Rebin_ = h1_[::3j, ::40j]
        elif nRebinY == 50:
            h1Rebin_ = h1_[::3j, ::50j]
        elif nRebinY == 100:
            h1Rebin_ = h1_[::3j, ::100j]
        else:
            print(f"{nRebinX = }, {nRebinY = } is not yet implemented... Implement it \t\t **** ERROR ****")        

    elif   nRebinX == 4:
        if   nRebinY == 1:
            h1Rebin_ = h1_
        elif nRebinY == 2:
            h1Rebin_ = h1_[::4j, ::2j]
        elif nRebinY == 3:
            h1Rebin_ = h1_[::4j, ::3j]
        elif nRebinY == 4:
            h1Rebin_ = h1_[::4j, ::4j]
        elif nRebinY == 5:
            h1Rebin_ = h1_[::4j, ::5j]
        elif nRebinY == 6:
            h1Rebin_ = h1_[::4j, ::6j]
        elif nRebinY == 10:
            h1Rebin_ = h1_[::4j, ::10j]
        elif nRebinY == 20:
            h1Rebin_ = h1_[::4j, ::20j]
        elif nRebinY == 40:
            h1Rebin_ = h1_[::4j, ::40j]
        elif nRebinY == 50:
            h1Rebin_ = h1_[::4j, ::50j]
        elif nRebinY == 100:
            h1Rebin_ = h1_[::4j, ::100j]
        else:
            print(f"{nRebinX = }, {nRebinY = } is not yet implemented... Implement it \t\t **** ERROR ****")        

    elif   nRebinX == 5:
        if   nRebinY == 1:
            h1Rebin_ = h1_
        elif nRebinY == 2:
            h1Rebin_ = h1_[::5j, ::2j]
        elif nRebinY == 3:
            h1Rebin_ = h1_[::5j, ::3j]
        elif nRebinY == 4:
            h1Rebin_ = h1_[::5j, ::4j]
        elif nRebinY == 5:
            h1Rebin_ = h1_[::5j, ::5j]
        elif nRebinY == 6:
            h1Rebin_ = h1_[::5j, ::6j]
        elif nRebinY == 10:
            h1Rebin_ = h1_[::5j, ::10j]
        elif nRebinY == 20:
            h1Rebin_ = h1_[::5j, ::20j]
        elif nRebinY == 40:
            h1Rebin_ = h1_[::5j, ::40j]
        elif nRebinY == 50:
            h1Rebin_ = h1_[::5j, ::50j]
        elif nRebinY == 100:
            h1Rebin_ = h1_[::5j, ::100j]
        else:
            print(f"{nRebinX = }, {nRebinY = } is not yet implemented... Implement it \t\t **** ERROR ****")        

    elif   nRebinX == 6:
        if   nRebinY == 1:
            h1Rebin_ = h1_
        elif nRebinY == 2:
            h1Rebin_ = h1_[::6j, ::2j]
        elif nRebinY == 3:
            h1Rebin_ = h1_[::6j, ::3j]
        elif nRebinY == 4:
            h1Rebin_ = h1_[::6j, ::4j]
        elif nRebinY == 5:
            h1Rebin_ = h1_[::6j, ::5j]
        elif nRebinY == 6:
            h1Rebin_ = h1_[::6j, ::6j]
        elif nRebinY == 10:
            h1Rebin_ = h1_[::6j, ::10j]
        elif nRebinY == 20:
            h1Rebin_ = h1_[::6j, ::20j]
        elif nRebinY == 40:
            h1Rebin_ = h1_[::6j, ::40j]
        elif nRebinY == 50:
            h1Rebin_ = h1_[::6j, ::50j]
        elif nRebinY == 100:
            h1Rebin_ = h1_[::6j, ::100j]
        else:
            print(f"{nRebinX = }, {nRebinY = } is not yet implemented... Implement it \t\t **** ERROR ****")        

    elif   nRebinX == 10:
        if   nRebinY == 1:
            h1Rebin_ = h1_
        elif nRebinY == 2:
            h1Rebin_ = h1_[::10j, ::2j]
        elif nRebinY == 3:
            h1Rebin_ = h1_[::10j, ::3j]
        elif nRebinY == 4:
            h1Rebin_ = h1_[::10j, ::4j]
        elif nRebinY == 5:
            h1Rebin_ = h1_[::10j, ::5j]
        elif nRebinY == 6:
            h1Rebin_ = h1_[::10j, ::6j]
        elif nRebinY == 10:
            h1Rebin_ = h1_[::10j, ::10j]
        elif nRebinY == 20:
            h1Rebin_ = h1_[::10j, ::20j]
        elif nRebinY == 40:
            h1Rebin_ = h1_[::10j, ::40j]
        elif nRebinY == 50:
            h1Rebin_ = h1_[::10j, ::50j]
        elif nRebinY == 100:
            h1Rebin_ = h1_[::10j, ::100j]
        else:
            print(f"{nRebinX = }, {nRebinY = } is not yet implemented... Implement it \t\t **** ERROR ****")        

    else:
        print(f"{nRebinX = }, {nRebinY = } is not yet implemented... Implement it \t\t **** ERROR ****")   

       
    #print(f"h1_ values ({len(h1_.values())}): {h1_.values()} \n variances ({len(h1_.variances())}): {h1_.variances()}")
    #print(f"h1Rebin_ values ({len(h1Rebin_.values())}): {h1Rebin_.values()} \n variances ({len(h1Rebin_.variances())}): {h1Rebin_.variances()}")
    if   nRebinX > 1 or nRebinY > 1 :    
        h1_ = h1Rebin_

    return h1_


def variableRebinTH1(h1_, xNewEdges):
    #print(f"rebinTH1():: histogram type {type(h1_) = },  {isinstance(h1_, hist.Hist) = }  ")
    if not (isinstance(h1_, hist.Hist) or  isinstance(h1_, coffea_hist.Hist)):
        print(f"rebinTH1():: histogram type {type(h1_)} not implemented... so could not rebin histogram ")
        return h1_
    
    if len(h1_.axes) != 1:
        print(f"rebinTH1:: histogram is not 1D")
        return h1_

    if not ( isinstance(xNewEdges, list)  or isinstance(xNewEdges, (np.ndarray, np.generic)) ):
        print(f"xNewEdges ({type(xNewEdges)}) needs to be list type")
        return h1_


    
    xOldEdges = h1_.axes[0].edges
    # bin numbers along the Xold axis that correspond to bin-edges of Xnew axis
    xOldIdx_pointing_xNewEdges = np.digitize(xNewEdges, xOldEdges) - 1 # xOldBinNumber (starting from 0) corresponds to xNewBinEdges
    print(f"h1_ x-axis ({type(xOldEdges)}) ({len(xOldEdges)}): {xOldEdges}")
    print(f"new x_axis ({type(xNewEdges)}) ({len(xNewEdges)}): {xNewEdges}")
    print(f"xOldIdx_pointing_xNewEdges ({type(xOldIdx_pointing_xNewEdges)}) ({len(xOldIdx_pointing_xNewEdges)}): {xOldIdx_pointing_xNewEdges}")

    h1Rebin_ = hist.Hist(hist.axis.Variable(xNewEdges, name=h1_.axes[0].name, label=h1_.axes[0].label), storage=hist.storage.Double())

    print(f"h1Rebin_.axes[0].centers ({type(h1Rebin_.axes[0].centers)}) ({len(h1Rebin_.axes[0].centers)}) ({h1Rebin_.axes[0].centers.shape[0]}) {h1Rebin_.axes[0].centers}")    
    print(f"h1_.values() ({h1_.values().shape[0]}): {h1_.values()}")
    for iBinXnew in range(h1Rebin_.axes[0].centers.shape[0]):
        firstBinInRangeXold        = xOldIdx_pointing_xNewEdges[iBinXnew    ]
        lasttBinInRangeXold_plus_1 = xOldIdx_pointing_xNewEdges[iBinXnew + 1]
        nEvents_iBinXnew = h1_.values()[firstBinInRangeXold : lasttBinInRangeXold_plus_1].sum()
        dX = (xNewEdges[iBinXnew+1] - xNewEdges[iBinXnew])
        print(f"iBinXnew: {iBinXnew}, ({xNewEdges[iBinXnew]}, {xNewEdges[iBinXnew+1]}) binXold: ({firstBinInRangeXold}, {lasttBinInRangeXold_plus_1}), ({xOldEdges[firstBinInRangeXold]}, {xOldEdges[lasttBinInRangeXold_plus_1]})", )
        print(f"{h1_.values()[firstBinInRangeXold : lasttBinInRangeXold_plus_1] = }, {nEvents_iBinXnew = }, {nEvents_iBinXnew/dX = } ")

    #print(f"h1Rebin_: {h1Rebin_}")
    print(f"h1Rebin_.axes[0].centers ({type(h1Rebin_.axes[0].centers)}) ({len(h1Rebin_.axes[0].centers)}) ({h1Rebin_.axes[0].centers.shape[0]}) {h1Rebin_.axes[0].centers}")
    print(f"h1Rebin_.values() ({type(h1Rebin_.values())}) ({len(h1Rebin_.values())}): {h1Rebin_.values()}")
    print(f"{h1Rebin_.variances() = }")






    return h1_

def calculateAverageOfArrays(array_list):
    #printVariable('\n ')
    #print(f"{len(array_list)}")
    #print(f"{type(array_list) = }")
    sTmp_ = 'len(array_list): %d \n' %(len(array_list))
    printNow_ = False
    for a_ in array_list:
        #print(f"{len(a_) = }: {ak.sum(ak.is_none(a_, axis=0), axis=0) = }")
        sTmp_ += 'len(a_): %d,  ak.sum(ak.is_none(a_, axis=0), axis=0): %g \n' % (len(a_), ak.sum(ak.is_none(a_, axis=0), axis=0) )
        if ak.sum(ak.is_none(a_, axis=0), axis=0)>0: printNow_ = True
    array_list = ak.fill_none(array_list, 0, axis=-1)
    for a_ in array_list:
        #print(f"{len(a_) = }: {ak.sum(ak.is_none(a_, axis=0), axis=0) = }")
        sTmp_ += 'len(a_): %d,  ak.sum(ak.is_none(a_, axis=0), axis=0): %g \t after update\n' % (len(a_), ak.sum(ak.is_none(a_, axis=0), axis=0) )
        #if ak.sum(ak.is_none(a_, axis=0), axis=0)>0: printNow_ = True
    if printNow_:
        print(sTmp_)
        #printVariable('\n array_list', array_list); sys.stdout.flush()
    a = np.vstack(array_list) # same as np.concatenate(array_list, axis=0)
    #avg_ = np.sum(a, axis=0) / len(array_list)
    return np.sum(a, axis=0) / len(array_list) 

def calculateMaxOfTwoArrays(array_a, array_b):
    a_new = ak.where(
        (array_a > array_b),
        array_a,
        array_b
    )
    return a_new

def calculateMaxOfArrays(array_list):


    a_max = array_list[0]
    for i in range(1, len(array_list)):
        a_max = calculateMaxOfTwoArrays(a_max, array_list[i])
    return a_max

def calculateMinOfTwoArrays(array_a, array_b):
    a_new = ak.where(
        (array_a < array_b),
        array_a,
        array_b
    )
    return a_new

def calculateMinOfArrays(array_list):
    a_min = array_list[0]
    for i in range(1, len(array_list)):
        a_min = calculateMinOfTwoArrays(a_min, array_list[i])
    return a_min

def array_PutLowerBound(array_list, k):
    a_new = ak.where(
        (array_list > k),
        array_list,
        ak.full_like(array_list, k)
    )
    return a_new

def array_PutUpperBound(array_list, k):
    a_new = ak.where(
        (array_list < k),
        array_list,
        ak.full_like(array_list, k)
    )
    return a_new


def stringHasSubstring(string, substringList):
    hasSubstring = False
    for s_ in substringList:
        if s_ in string: hasSubstring = True
    return hasSubstring


def printVariable(sName, var):
    printInDetail=True
    #if nEventsToAnalyze == -1: printInDetail = False
    if str(type(var)) in ['numpy.ndarray', "<class 'numpy.ndarray'>", "<class 'numpy.ma.core.MaskedArray'>"]: printInDetail = False # as gave error
    #print(f"printInDetail: {printInDetail} {sName} ({type(var)}) ({len(var)}): {var}")
    if not printInDetail:
        #print(f"{sName} ({type(var)}) ({len(var)}): {var}")
        try:
            print(f"{sName} ({type(var)}) ({len(var)}): {var.tolist()}")
        except:
            print(f"{sName} ({type(var)}) ({len(var)}): {var}")
    else:
        try:
            print(f"{sName} ({type(var)}) ({len(var)}): {var.to_list()}")
        except:
            print(f"{sName} ({type(var)}) ({len(var)}): {var}")


def printVariablePtEtaPhi(sName, var):
    var_PtEtaPhi = ak.zip([
        var.pt,
        var.eta,
        var.phi
    ])
    printVariable('%s PtEtaPhi' % sName, var_PtEtaPhi)


def akArray_isin(testArray, referenceArray):
    '''
    Compare each element of testArray (along axis=1) with referenceArray (along axis=1), and return boolean array (with shape of testArray) 
    '''
    return ak.from_iter( [ np.isin(testArray[idx_], referenceArray[idx_]) for idx_ in range(len(testArray)) ] )


def insertInListBeforeThisElement(list1, sConditionToAdd, addBeforeThisCondition):
    for idx_ in range(len(list1)):
        if list1[idx_] == addBeforeThisCondition:
            list1.insert(idx_, sConditionToAdd)
            break
    return list1
