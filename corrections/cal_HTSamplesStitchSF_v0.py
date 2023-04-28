
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



if __name__ == "__main__":
    
    print("cal_HTSamplesStitchSF.py:: main: {}".format(sys.argv)); sys.stdout.flush()

    parser = argparse.ArgumentParser(description='cal_HTSamplesStitchSF')
    parser.add_argument('-era', dest='era',   type=str, default=Era_2018,                    choices=[Era_2016, Era_2017, Era_2018], required=False)
    args=parser.parse_args()
    print("args: {}".format(args))

    era              = args.era



    sFIn_dict = {
        Era_2016: "",
        Era_2017: "",        
        Era_2018: "../data/correction/mc/HTSamplesStitch/LHE_HT_2018.root"
    }
    sHistogramName = 'evt/$SAMPLECATEGORYNAME/hGenLHE_HT_all_central'
    sFOut_dict = {
        Era_2016: "",
        Era_2017: "",        
        Era_2018: "../data/correction/mc/HTSamplesStitch/HTSamplesStitchSF_2018.root"
    }


    kInf = 999999.0
    sHTRange     = 'HTRange'
    sHTBins      = 'HTBins'
    sHTBinEdges  = 'HTBinEdges'
    sHTBinAddiSF = 'HTBinAddiSF'

    HTBinsDetails_forSFCal = OD([
        ("QCD_bGen", {
            sHTRange:    [0, 2500],
            sHTBins: OD([
                ('100to200', {
                    sHTRange: [ 100,  200],
                    sHTBinAddiSF:    1266000   / 1275000,
                }),
                ('200to300', {
                    sHTRange: [ 200,  300],
                    sHTBinAddiSF:     109900   /  111700,
                }),
                ('300to500', {
                    sHTRange: [ 300,  500],
                    sHTBinAddiSF:      27360   /   27960,
                }),
                ('500to700', {
                    sHTRange: [ 500,  700],
                    sHTBinAddiSF:       2991   /    3078,
                }),
                ('700to1000', {
                    sHTRange: [ 700, 1000],
                    sHTBinAddiSF:      731.8  /    721.8,
                }),
                ('1000to1500', {
                    sHTRange: [1000, 1500],
                    sHTBinAddiSF:     139.3  /    138.2,
                }),
                ('1500to2000', {
                    sHTRange: [1500, 2000],
                    sHTBinAddiSF:      14.74 /     13.61,
                }),
                ('2000toInf', {
                    sHTRange: [2000, kInf],
                    sHTBinAddiSF:        3.09 /      2.92,
                }),
                
            ]),
        }),

        
    ])

    nBins_forHTExtrapolation = 5
    calHTCorrSFFromHighToLowHT = True 
    
    sFIn  = sFIn_dict[era]
    sFOut = sFOut_dict[era]


        
    # Read input file --------------------------------
    print(f"Input file: {sFIn}")
    fIn = R.TFile(sFIn)
    if not fIn.IsOpen():
        print(f"Could not open input file: {sFIn}")
        exit(0)

    fOut = R.TFile(sFOut, 'RECREATE')
    
    for sampleName, HTBinsDetails in HTBinsDetails_forSFCal.items():
        HTRangeFull_min = HTBinsDetails[sHTRange][0]
        HTRangeFull_max = HTBinsDetails[sHTRange][1]
        HTBins_dict     = HTBinsDetails[sHTBins]

        nBins_HT = int(HTRangeFull_max - HTRangeFull_min)
        hSF_step0 = R.TH1D("hHTSamplesStitchSF_step0", "", nBins_HT, HTRangeFull_min, HTRangeFull_max)
        for iBin in range(1, nBins_HT+1):
            hSF_step0.SetBinContent(iBin, 1.0)

        sHistoName_toUse = sHistogramName
        sHistoName_toUse = sHistoName_toUse.replace('$SAMPLECATEGORYNAME', sampleName)
        hGenLHE_HT = fIn.Get(sHistoName_toUse)
        hGenLHE_HT_wHTSamplesStitch = hGenLHE_HT.Clone('%s_wHTSamplesStitch' % (hGenLHE_HT.GetName()))
        print(f"GEN LHE_HT: {sHistoName_toUse}")


        # Update GenLHE_HT spectrum using cross-section derived using UL samples
        hGenLHE_HT_wAddiSF = hGenLHE_HT.Clone('%s_wAdditionalSF' % (hGenLHE_HT.GetName()))        
        for HTBinName in HTBinsDetails[sHTBins]:
            HTBin_min    = HTBinsDetails[sHTBins][HTBinName][sHTRange][0]
            HTBin_max    = HTBinsDetails[sHTBins][HTBinName][sHTRange][1]
            SFAdditional = HTBinsDetails[sHTBins][HTBinName][sHTBinAddiSF]

            bin_HTBin_min = hGenLHE_HT_wAddiSF.FindBin(HTBin_min)
            bin_HTBin_max = hGenLHE_HT_wAddiSF.FindBin(HTBin_max)
            for iBin in range(bin_HTBin_min, bin_HTBin_max):
                nEvents    = hGenLHE_HT.GetBinContent(iBin)
                errNEvents = hGenLHE_HT.GetBinError(iBin)
                hGenLHE_HT_wAddiSF.SetBinContent(iBin, nEvents    * SFAdditional)
                hGenLHE_HT_wAddiSF.SetBinError(  iBin, errNEvents * SFAdditional)


        HTBinEdges = []
        for HTBinName in HTBinsDetails[sHTBins]:
            for i in range(2):
                if HTBinsDetails[sHTBins][HTBinName][sHTRange][i] not in HTBinEdges:
                    HTBinEdges.append( HTBinsDetails[sHTBins][HTBinName][sHTRange][i] )

        print(f"HTBinEdges: {HTBinEdges}")
        print(f" {HTBinEdges[::-1][1:-1] = }")

        HTBins_forSFComputation = None
        if calHTCorrSFFromHighToLowHT:
            HTBins_forSFComputation  = HTBinEdges[::-1][1:-1] # HT bins edges excluding end-points and in reverse order in HT
        else:
            HTBins_forSFComputation  = HTBinEdges[1:-1]       # HT bins edges excluding end-points
            
        for HTBinEdge in HTBins_forSFComputation:
            bin_HTBinEdge = hGenLHE_HT.FindBin(HTBinEdge)
            print(f"HTBinEdge: {HTBinEdge}")
        
            
            HTFitHigherRange_min = hGenLHE_HT.GetBinCenter(bin_HTBinEdge)
            HTFitHigherRange_max = hGenLHE_HT.GetBinCenter(bin_HTBinEdge + nBins_forHTExtrapolation -1 )
            print(f"{HTFitHigherRange_min = }, {HTFitHigherRange_max = }")
            fit_higher = R.TF1("lHigherAt%g"%(HTBinEdge), "pol1", HTFitHigherRange_min, HTFitHigherRange_max)
            
            hGenLHE_HT_wHTSamplesStitch.Fit(fit_higher, "R0Q")

            print(f" nEvents: {hGenLHE_HT.GetBinContent(bin_HTBinEdge)} +- {hGenLHE_HT.GetBinError(bin_HTBinEdge)},   nEvents_fit: {fit_higher.Eval(HTBinEdge)},   diff: {hGenLHE_HT.GetBinContent(bin_HTBinEdge) - fit_higher.Eval(HTBinEdge)},  {(hGenLHE_HT.GetBinContent(bin_HTBinEdge) - fit_higher.Eval(HTBinEdge))/hGenLHE_HT.GetBinContent(bin_HTBinEdge)*100} %")

            
            HTFitLowerRange_min = hGenLHE_HT.GetBinCenter(bin_HTBinEdge - nBins_forHTExtrapolation)
            HTFitLowerRange_max = hGenLHE_HT.GetBinCenter(bin_HTBinEdge - 1 )
            print(f"{HTFitLowerRange_min = }, {HTFitLowerRange_max = }")
            fit_lower = R.TF1("lLowerAt%g"%(HTBinEdge), "pol1", HTFitLowerRange_min, HTFitLowerRange_max)
            
            hGenLHE_HT_wHTSamplesStitch.Fit(fit_lower, "R0Q")

            kFactor = None
            if calHTCorrSFFromHighToLowHT:
                kFactor = fit_higher.Eval(HTBinEdge) / fit_lower.Eval(HTBinEdge)
                print(f" kFactor: {kFactor} =  {fit_higher.Eval(HTBinEdge)} / {fit_lower.Eval(HTBinEdge)} ")
            else:
                kFactor = fit_lower.Eval(HTBinEdge)  / fit_higher.Eval(HTBinEdge)
                print(f" kFactor: {kFactor} =  {fit_lower.Eval(HTBinEdge)}  / {fit_higher.Eval(HTBinEdge)} ")
            

            binsToCorrect = None
            if calHTCorrSFFromHighToLowHT:
                binsToCorrect = range(1, bin_HTBinEdge)
            else:
                binsToCorrect = range(bin_HTBinEdge, hGenLHE_HT.GetNbinsX()+1)
                
            for iBin in binsToCorrect:
                SF_ = hSF_step0.GetBinContent(iBin) * kFactor
                
                hSF_step0.SetBinContent(iBin, SF_)

                nEvents_    = hGenLHE_HT.GetBinContent(iBin)
                errNEvents_ = hGenLHE_HT.GetBinError(iBin)
                hGenLHE_HT_wHTSamplesStitch.SetBinContent(iBin, nEvents_    * SF_)
                hGenLHE_HT_wHTSamplesStitch.SetBinError(  iBin, errNEvents_ * SF_)
            

        fOut.cd()
        fOut.mkdir(sampleName)
        fOut.cd(sampleName)
        hSF_step0.Write()
        hGenLHE_HT.Write()
        hGenLHE_HT_wHTSamplesStitch.Write()
        hGenLHE_HT_wAddiSF.Write()

    fOut.Close()
    print(f"Wrote SFs in {sFOut}")
