
import os
import sys
from pathlib import Path
import json
import glob
import argparse
from collections import OrderedDict as OD
import numpy as np
import time
from datetime import datetime
import copy

import ROOT as R

from htoaa_Settings import *


printLevel = 10

def th1GenBinContent(hist, x):
    #print(f"hist ({type(hist)}): {hist}")
    iBin = hist.FindBin(x)
    return hist.GetBinContent(iBin)
    

if __name__ == "__main__":
    
    print("samplesList_update.py:: main: {}".format(sys.argv)); sys.stdout.flush()

    parser = argparse.ArgumentParser(description='samplesList_update')
    parser.add_argument('-era', dest='era',   type=str, default=Era_2018,                    choices=[Era_2016, Era_2017, Era_2018], required=False)
    args=parser.parse_args()
    print("args: {}".format(args))

    era              = args.era


    
    sFIn_analysis_stage1_2018 = "data/countSumEventsInSamples/analyze_htoaa_stage1_2018.root"

    sHistogramName = 'evt/$DATASETNAME/hCutFlowWeighted_central'

    

    samplesInfo               = None
    sFIn_analysis_stage1      = None
    sFileSamplesInfo_updated  = None
    if era == Era_2018:
        sFIn_analysis_stage1 = sFIn_analysis_stage1_2018
        with open(sFileSamplesInfo[era]) as fSamplesInfo:
            samplesInfo = json.load(fSamplesInfo) # Samples_Era.json

        sFileSamplesInfo_updated = sFileSamplesInfo[era]
        sFileSamplesInfo_updated = sFileSamplesInfo_updated.replace('.json', '_updated.json')

    fIn_analysis_stage1 = R.TFile(sFIn_analysis_stage1_2018)
    print(f"{sFIn_analysis_stage1 = } \n")

    if printLevel >= 10:
        print(f"samplesInfo ({type(samplesInfo)}) ({len(samplesInfo)}): {json.dumps(samplesInfo, indent=4)}")

    for sampleName in samplesInfo:
        datasetName           = samplesInfo[sampleName]["dataset"]
        
        datasetName_trim      = datasetName[1:] if datasetName.startswith('/') else datasetName
        datasetName_parts     = datasetName.split('/')
        isMC                  = True if datasetName_parts[-1] == 'NANOAODSIM' else False
        #print(f"sampleName: {sampleName} \t\t\t dataset: {datasetName} ")

        if not isMC: continue

        nEvents   = samplesInfo[sampleName]["nEvents"]
        sumEvents = samplesInfo[sampleName]["sumEvents"]
        
        #sHistoName = 'evt/%s/hCutFlowWeighted_central' % (datasetName_trim)
        sHistoName_toUse = sHistogramName
        sHistoName_toUse = sHistoName_toUse.replace('$DATASETNAME', datasetName_trim)
        
        h = None
        h = fIn_analysis_stage1.Get(sHistoName_toUse)
        if not h: continue
        
        nEvents_cal           = th1GenBinContent(h, 0)
        sumEvents_cal         = th1GenBinContent(h, 1)
        nEvents_posGenWgt_cal = th1GenBinContent(h, 2)
        nEvents_negGenWgt_cal = th1GenBinContent(h, 3)

        samplesInfo[sampleName]["sumEvents"] = sumEvents_cal

        if printLevel >= 1:
            print(" %-80s: %12d / %12d (%6d),  %12d,  %12d, %12d \t %g" % (sampleName, nEvents_cal,nEvents,(nEvents_cal-nEvents), sumEvents_cal, nEvents_posGenWgt_cal, nEvents_negGenWgt_cal, abs(nEvents_negGenWgt_cal/nEvents_cal)))


    with open(sFileSamplesInfo_updated, "w") as fSampleInfo:
        json.dump(samplesInfo, fSampleInfo, indent=4)

        print(f"\n\n Updated sample list wrote to {sFileSamplesInfo_updated}")
    
