import os
import sys
import subprocess
import shlex
import logging
import json
import numpy as np
#import uproot
import uproot3 as uproot
import ROOT as R
from parse import *

from htoaa_Settings import * 
from htoaa_Samples import (
    kData, kQCD_bEnrich, kQCD_bGen, kQCDIncl, kZJets, kWJets
)


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
    
    

def setXRootDRedirector(fileName):
    if not fileName.startswith("/store/"):
        return fileName
    
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

def xrdcpFile(sFileName, sFileNameLocal, nTry = 3):
    command_ = "time xrdcp %s %s" % (sFileName, sFileNameLocal)
    command_list_ = command_.split(" ")
    print(f"{command_ = }")
    for iTry in range(nTry):
        process = subprocess.Popen(command_list_,
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True
                                   )
        stdout, stderr = process.communicate()
        print(f"  {iTry = } {stdout = }, {stderr = }");  sys.stdout.flush()
        if 'FATAL' not in stderr and 'ERROR' not in stderr : # download was successful
            return True

    return False
    

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
    #    preComment, postComment = contents.split("/*", 1)
    #    contents = preComment + postComment.split("*/", 1)[1]
    while "'''" in contents:
        preComment, postComment = contents.split("'''", 1)
        contents = preComment + postComment.split("'''", 1)[1]

    dictionary =  json.loads( contents )
    return dictionary



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
