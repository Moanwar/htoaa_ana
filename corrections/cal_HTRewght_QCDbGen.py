
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
            "HTBins": [
                [ 100,  200],
                [ 200,  300],
                [ 300,  500],
                [ 500,  700],
                [ 700, 1000],
                [1000, 1500],
                [1500, 2000],
                [2000, 3000],
            ],
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

    HTRewgt = OD()
    for sample_category in sHistogramNameForSF:
        sHisto_N               = sHistogramNameForSF[sample_category]["N"]
        sHisto_D               = sHistogramNameForSF[sample_category]["D"]
        xRangeForSFComputation = sHistogramNameForSF[sample_category]["Range"]
        HTBins                 = sHistogramNameForSF[sample_category]["HTBins"]
        HTRewgt[sample_category] = OD()
        HTRewgt[sample_category][era] = OD()

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

        '''
        for iBin in range(1, hSF.GetNbinsX()+1):
            y        = hSF.GetBinContent(iBin)
            xBinMin = hSF.GetXaxis().GetBinLowEdge(iBin)
            xBinMax = hSF.GetXaxis().GetBinUpEdge(iBin)
            #if abs(y - 0) < 1e-3: continue
            if xBinMin >= xRangeForSFComputation[0] and xBinMax <= xRangeForSFComputation[1]: continue

            #print(f"{iBin = }, {y = }, {xBinMin = }, {xBinMax = }, { = }, { = }")
            hSF.SetBinContent(iBin, 1)
        '''

        c1 = R.TCanvas("c1_%s" % (sample_category),"c1_%s" % (sample_category), 600,450)
        c1.cd()
        legend = R.TLegend(0.1,0.1, 0.99,0.6)
        legend.AddEntry(hSF, hSF.GetName(), "lep")

        #fSFs = []
        HTRewgt[sample_category][era]["FitFunctionFormat"] = "{p0} + ({p1} * (x - {HTBinMin}))" 
        for HTBin in HTBins:
            HTBinMin = HTBin[0]
            HTBinMax = HTBin[1]

            #sFitFunction = "[0] + [1]*(x - %d)" % (HTBinMin)
            sFitFunction = f"[0] + ([1] * (x - {HTBinMin}))"
            fSF = R.TF1("fSF_%s_HT%dto%d" % (sample_category, HTBinMin,HTBinMax), sFitFunction, HTBinMin,HTBinMax-5)
            fSF.SetLineColor(2)
            fSF.SetLineWidth(2)
            hSF.Fit(fSF, "R+")

            sFitFunction_afterFit = sFitFunction
            sFitFunction_afterFit = sFitFunction_afterFit.replace('[0]', '%.6f' % fSF.GetParameter(0) )
            sFitFunction_afterFit = sFitFunction_afterFit.replace('[1]', '%.6f' % fSF.GetParameter(1) )
            legend.AddEntry(fSF, "%s: %s" % (fSF.GetName(), sFitFunction_afterFit), "l")
            HTRewgt[sample_category][era]['HT%dto%d' % (HTBinMin,HTBinMax)] = sFitFunction_afterFit
            
        
        hSF.Draw()
        legend.Draw()

        #c1.Draw()
        
        fOut.cd()
        fOut.mkdir(sample_category)
        fOut.cd(sample_category)

        hN.Write()
        hD.Write()
        hSF.Write()
        hSF_0.Write()
        c1.Write()
        
        
        
    fOut.Close()
    print(f"Wrote SFs in {sFOut}")
    print(f'\n"HTRewgt" : {json.dumps(HTRewgt, indent = 4)} ')
    

    
