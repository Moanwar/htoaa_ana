import os
import sys
from datetime import datetime
print(f"htoaa_CommonTools:: here1 {datetime.now() = }"); sys.stdout.flush()
import subprocess
import shlex
import logging
import json
print(f"htoaa_CommonTools:: here2 {datetime.now() = }"); sys.stdout.flush()
import numpy as np
import math
import awkward as ak
print(f"htoaa_CommonTools:: here3 {datetime.now() = }"); sys.stdout.flush()
#import uproot
import uproot3 as uproot
print(f"htoaa_CommonTools:: here4 {datetime.now() = }"); sys.stdout.flush()
#import ROOT as R
print(f"htoaa_CommonTools:: here5 {datetime.now() = }"); sys.stdout.flush()
from parse import *
import logging

print(f"htoaa_CommonTools:: here6 {datetime.now() = }"); sys.stdout.flush()
from htoaa_Settings import * 
print(f"htoaa_CommonTools:: here7 {datetime.now() = }"); sys.stdout.flush()
from htoaa_Samples import (
    kData, kQCD_bEnrich, kQCD_bGen, kQCDIncl, kZJets, kWJets
)
print(f"htoaa_CommonTools:: here8 {datetime.now() = }"); sys.stdout.flush()
#from numba import jit
print(f"htoaa_CommonTools:: here9 {datetime.now() = }"); sys.stdout.flush()

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
    print(f"htoaa_CommonTools::getNanoAODFile() here1 {datetime.now() = }"); sys.stdout.flush()

    if useLocalFileIfExists and fileName.startswith("/store/") and server in ['lxplus']: # 'eos cp' works on lxplus
        # check if NanoAOD exists in /eos area
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
        print(f"Checking for eos file: {fileName_EOS = }") 
        print(f"htoaa_CommonTools::getNanoAODFile() here2 {datetime.now() = }"); sys.stdout.flush()

        print(f"{fileName_EOS = }: {os.path.exists(fileName_EOS) = } ")
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
