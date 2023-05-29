
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

sys.path.insert(1, '../') # to import file from other directory (../ in this case)

from htoaa_Settings import *


def cal_HistogramBins_wMaxRelUncertainty(h, maxRelUnc=10):
    ds
    



if __name__ == "__main__":
    
    print("cal_HTRewght_QCDbGen.py:: main: {}".format(sys.argv)); sys.stdout.flush()

    parser = argparse.ArgumentParser(description='cal_HTRewght_QCDbGen')
    parser.add_argument('-era', dest='era',   type=str, default=Era_2018,                    choices=[Era_2016, Era_2017, Era_2018], required=False)
    args=parser.parse_args()
    print("args: {}".format(args))

    era              = args.era


    nRebin = 10
    sFIn_dict = {
        Era_2016: "",
        Era_2017: "",        
        Era_2018: "../data/correction/mc/HTSamplesStitch/LHE_HT_2018.root"
    }

    sFOut_dict = {
        Era_2016: "",
        Era_2017: "",        
        Era_2018: "../data/correction/mc/HTSamplesStitch/HTSamplesStitchSF_2018.root"
    }

    sHistogramNameForSF = {
        "QCD_bGen": {
            "N": "evt/QCD_Incl_PSWeight/hGenLHE_HT_SelQCDbGen_central", # SF numerator histogram
            "D": "evt/QCD_bGen/hGenLHE_HT_SelQCDbGen_central", # SF denominator histogram
            "Range": [100, 1800], # axis (HT) range for SF computation
        }
    }

    
    sFIn  = sFIn_dict[era]
    sFOut = sFOut_dict[era]


        
    # Read input file --------------------------------
    print(f"Input file: {sFIn}")
    fIn = R.TFile(sFIn)
    if not fIn.IsOpen():
        print(f"Could not open input file: {sFIn}")
        exit(0)


    hTmp = R.TH1D("hTMp", "", 1,0,1)
    hTmp.SetDefaultSumw2()
    
    fOut = R.TFile(sFOut, 'RECREATE')

    for sample_category in sHistogramNameForSF:
        sHisto_N               = sHistogramNameForSF[sample_category]["N"]
        sHisto_D               = sHistogramNameForSF[sample_category]["D"]
        xRangeForSFComputation = sHistogramNameForSF[sample_category]["Range"]

        hN       = fIn.Get(sHisto_N)
        hD       = fIn.Get(sHisto_D)

        hN.Rebin(nRebin)
        hD.Rebin(nRebin)
        
        sample_N = sHisto_N.split('/')[1]
        sample_D = sHisto_D.split('/')[1]

        histoNameLastPart_N = sHisto_N.split('/')[-1]
        histoNameLastPart_D = sHisto_D.split('/')[-1]

        hN.SetName("%s_%s" % (histoNameLastPart_N, sample_N))
        hD.SetName("%s_%s" % (histoNameLastPart_D, sample_D))

        hSF = hD.Clone("%s_SF" % (histoNameLastPart_D))
        hSF.Divide(hN, hD)
        hSF.GetYaxis().SetTitle('SF')

        hSF_0 = hSF.Clone("%s_Ratio" % (histoNameLastPart_D))

        for iBin in range(1, hSF.GetNbinsX()+1):
            y        = hSF.GetBinContent(iBin)
            xBinMin = hSF.GetXaxis().GetBinLowEdge(iBin)
            xBinMax = hSF.GetXaxis().GetBinUpEdge(iBin)
            #if abs(y - 0) < 1e-3: continue
            if xBinMin >= xRangeForSFComputation[0] and xBinMax <= xRangeForSFComputation[1]: continue

            #print(f"{iBin = }, {y = }, {xBinMin = }, {xBinMax = }, { = }, { = }")
            hSF.SetBinContent(iBin, 1)
        
        fOut.cd()
        fOut.mkdir(sample_category)
        fOut.cd(sample_category)

        hN.Write()
        hD.Write()
        hSF.Write()
        hSF_0.Write()
        
        
    fOut.Close()
    print(f"Wrote SFs in {sFOut}")
