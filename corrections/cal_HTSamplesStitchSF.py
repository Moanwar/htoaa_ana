
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





def calculate_HTStitchSFs(
        hGenLHE_HT,
        HTRangeFull_min,
        HTRangeFull_max,
        HTBinEdges,
        calHTCorrSFFromHighToLowHT,
        sVersion
):
    print(f"\n\ncalculate_HTStitchSFs():: hGenLHE_HT: {hGenLHE_HT.GetName()}, HTBinEdges: {HTBinEdges}, calHTCorrSFFromHighToLowHT: {calHTCorrSFFromHighToLowHT}, sVersion: {sVersion}")
        
    sHTCorrSFDirection = 'HTCorrSFFromHighToLowHT' if calHTCorrSFFromHighToLowHT else 'HTCorrSFFromLowToHighHT'    
    
    nBins_HT = int(HTRangeFull_max - HTRangeFull_min)
    hSF_step0 = R.TH1D("hHTSamplesStitchSF_step0_%s%s" % (sHTCorrSFDirection, sVersion), "", nBins_HT, HTRangeFull_min, HTRangeFull_max)
    for iBin in range(1, nBins_HT+1):
        hSF_step0.SetBinContent(iBin, 1.0)

    hGenLHE_HT_wHTSamplesStitch = hGenLHE_HT.Clone('%s_wHTSamplesStitch_%s' % (hGenLHE_HT.GetName(), sHTCorrSFDirection))

    HTBins_forSFComputation = None
    if calHTCorrSFFromHighToLowHT:
        HTBins_forSFComputation  = HTBinEdges[::-1][1:-1] # HT bins edges excluding end-points and in reverse order in HT
    else:
        HTBins_forSFComputation  = HTBinEdges[1:-1]       # HT bins edges excluding end-points
    
    print(f"calculate_HTStitchSFs():: HTBinEdges: {HTBinEdges}")
    print(f"calculate_HTStitchSFs():: HTBins_forSFComputation: {HTBins_forSFComputation}")
    for HTBinEdge in HTBins_forSFComputation:
        bin_HTBinEdge = hGenLHE_HT.FindBin(HTBinEdge)
        print(f"calculate_HTStitchSFs():: HTBinEdge: {HTBinEdge}")


        HTFitHigherRange_min = hGenLHE_HT.GetBinCenter(bin_HTBinEdge)
        HTFitHigherRange_max = hGenLHE_HT.GetBinCenter(bin_HTBinEdge + nBins_forHTExtrapolation -1 )
        print(f"{HTFitHigherRange_min = }, {HTFitHigherRange_max = }")
        fit_higher = R.TF1("%s_%s_lHigherAt%g"%(hGenLHE_HT.GetName(), sHTCorrSFDirection, HTBinEdge), "pol1", HTFitHigherRange_min, HTFitHigherRange_max)

        hGenLHE_HT_wHTSamplesStitch.Fit(fit_higher, "R0Q")
        #hGenLHE_HT_wHTSamplesStitch.Fit(fit_higher, "RQ")

        print(f"calculate_HTStitchSFs()::  nEvents: {hGenLHE_HT.GetBinContent(bin_HTBinEdge)} +- {hGenLHE_HT.GetBinError(bin_HTBinEdge)},   nEvents_fit: {fit_higher.Eval(HTBinEdge)},   diff: {hGenLHE_HT.GetBinContent(bin_HTBinEdge) - fit_higher.Eval(HTBinEdge)},  {(hGenLHE_HT.GetBinContent(bin_HTBinEdge) - fit_higher.Eval(HTBinEdge))/hGenLHE_HT.GetBinContent(bin_HTBinEdge)*100} %")


        HTFitLowerRange_min = hGenLHE_HT.GetBinCenter(bin_HTBinEdge - nBins_forHTExtrapolation)
        HTFitLowerRange_max = hGenLHE_HT.GetBinCenter(bin_HTBinEdge - 1 )
        print(f"{HTFitLowerRange_min = }, {HTFitLowerRange_max = }")
        fit_lower = R.TF1("%s_%s_lLowerAt%g"%(hGenLHE_HT.GetName(), sHTCorrSFDirection, HTBinEdge), "pol1", HTFitLowerRange_min, HTFitLowerRange_max)

        hGenLHE_HT_wHTSamplesStitch.Fit(fit_lower, "R0Q")
        #hGenLHE_HT_wHTSamplesStitch.Fit(fit_lower, "RQ")

        eNEvents = hGenLHE_HT.GetBinError(bin_HTBinEdge) if hGenLHE_HT.GetBinError(bin_HTBinEdge) > hGenLHE_HT.GetBinError(bin_HTBinEdge - 1) else hGenLHE_HT.GetBinError(bin_HTBinEdge - 1)
        dY = fit_higher.Eval(HTBinEdge) - fit_lower.Eval(HTBinEdge)

        # ignore jump if statistically insignificant
        if abs(dY) < 3 * eNEvents: continue
        

        kFactor = None
        if calHTCorrSFFromHighToLowHT:
            kFactor = fit_higher.Eval(HTBinEdge) / fit_lower.Eval(HTBinEdge)
            print(f"calculate_HTStitchSFs()::  kFactor: {kFactor} =  {fit_higher.Eval(HTBinEdge)} / {fit_lower.Eval(HTBinEdge)} ")
        else:
            kFactor = fit_lower.Eval(HTBinEdge)  / fit_higher.Eval(HTBinEdge)
            print(f"calculate_HTStitchSFs()::  kFactor: {kFactor} =  {fit_lower.Eval(HTBinEdge)}  / {fit_higher.Eval(HTBinEdge)} ")


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


    # scale hSF to have totoal cross-section of the process over full HT range before and after applying HTSamplesStitch remains same
    histoIntegralRatio = hGenLHE_HT.Integral() / hGenLHE_HT_wHTSamplesStitch.Integral() 
    hSF_step1 = hSF_step0.Clone(hSF_step0.GetName().replace('step0', 'step1'))
    hSF_step1.Scale( histoIntegralRatio )
    print(f"{hGenLHE_HT.Integral() = }, {hGenLHE_HT_wHTSamplesStitch.Integral() = },  IntegralRatio: {histoIntegralRatio}")
    
    hGenLHE_HT_wHTSamplesStitch1 = hGenLHE_HT.Clone('%s_wHTSamplesStitch1_%s' % (hGenLHE_HT.GetName(), sHTCorrSFDirection))
    for iBin in range(1, hGenLHE_HT.GetNbinsX()+1):
        SF_         = hSF_step1.GetBinContent(iBin)
        nEvents_    = hGenLHE_HT.GetBinContent(iBin)
        errNEvents_ = hGenLHE_HT.GetBinError(iBin)
        hGenLHE_HT_wHTSamplesStitch1.SetBinContent(iBin, nEvents_    * SF_)
        hGenLHE_HT_wHTSamplesStitch1.SetBinError(  iBin, errNEvents_ * SF_)
        hSF_step1.SetBinError(iBin, 0.0)
        
    print(f"{hGenLHE_HT_wHTSamplesStitch1.Integral() = }")
    
    return hSF_step0, hGenLHE_HT_wHTSamplesStitch, hSF_step1, hGenLHE_HT_wHTSamplesStitch1




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

        ("ZJets", {
            sHTRange:    [0, 1000],
            sHTBins: OD([
                ('200to400', {
                    sHTRange: [ 200,  400],
                }),
                ('400to600', {
                    sHTRange: [ 400,  600],
                }),
                ('600to800', {
                    sHTRange: [ 600,  800],
                }),
                ('800toInt', {
                    sHTRange: [ 800,  kInf],
                }),
            ]),
        }),

         ("WJets", {
            sHTRange:    [0, 1000],
            sHTBins: OD([
                ('200to400', {
                    sHTRange: [ 200,  400],
                }),
                ('400to600', {
                    sHTRange: [ 400,  600],
                }),
                ('600to800', {
                    sHTRange: [ 600,  800],
                }),
                ('800toInt', {
                    sHTRange: [ 800,  kInf],
                }),
            ]),
        }),
       
    ])

    nBins_forHTExtrapolation = 7
    
    
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


        sHistoName_toUse = sHistogramName
        sHistoName_toUse = sHistoName_toUse.replace('$SAMPLECATEGORYNAME', sampleName)
        hGenLHE_HT = fIn.Get(sHistoName_toUse)
        print(f"GEN LHE_HT: {sHistoName_toUse}")


        HTBinEdges = []
        for HTBinName in HTBinsDetails[sHTBins]:
            for i in range(2):
                if HTBinsDetails[sHTBins][HTBinName][sHTRange][i] not in HTBinEdges:
                    HTBinEdges.append( HTBinsDetails[sHTBins][HTBinName][sHTRange][i] )

        print(f"HTBinEdges: {HTBinEdges}")
        print(f" {HTBinEdges[::-1][1:-1] = }")

        # SF computation from high to low HT
        hSF_step0_1, hGenLHE_HT_wHTSamplesStitch_1, hSF_step1_1, hGenLHE_HT_wHTSamplesStitch1_1 = calculate_HTStitchSFs(
            hGenLHE_HT = hGenLHE_HT,
            HTRangeFull_min = HTRangeFull_min,
            HTRangeFull_max = HTRangeFull_max,
            HTBinEdges = HTBinEdges,
            calHTCorrSFFromHighToLowHT = True,
            sVersion=""
        )
        # SF computation from low to high HT
        hSF_step0_2, hGenLHE_HT_wHTSamplesStitch_2, hSF_step1_2, hGenLHE_HT_wHTSamplesStitch1_2 = calculate_HTStitchSFs(
            hGenLHE_HT = hGenLHE_HT,
            HTRangeFull_min = HTRangeFull_min,
            HTRangeFull_max = HTRangeFull_max,
            HTBinEdges = HTBinEdges,
            calHTCorrSFFromHighToLowHT = False,
            sVersion=""
        )

        fOut.cd()
        fOut.mkdir(sampleName)
        fOut.cd(sampleName)

        hGenLHE_HT.Write()

        hSF_step0_1.Write()
        hGenLHE_HT_wHTSamplesStitch_1.Write()
        hSF_step1_1.Write()
        hGenLHE_HT_wHTSamplesStitch1_1.Write()
        
        hSF_step0_2.Write()
        hGenLHE_HT_wHTSamplesStitch_2.Write()
        hSF_step1_2.Write()
        hGenLHE_HT_wHTSamplesStitch1_2.Write()

        
        # check using additional SF if provided
        sHTBin_first = list(HTBinsDetails[sHTBins].keys())[0]
        if sHTBinAddiSF in HTBinsDetails[sHTBins][sHTBin_first]:

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


            # SF computation from high to low HT
            hSF_step0_3, hGenLHE_HT_wHTSamplesStitch_3, hSF_step1_3, hGenLHE_HT_wHTSamplesStitch1_3 = calculate_HTStitchSFs(
                hGenLHE_HT = hGenLHE_HT_wAddiSF,
                HTRangeFull_min = HTRangeFull_min,
                HTRangeFull_max = HTRangeFull_max,
                HTBinEdges = HTBinEdges,
                calHTCorrSFFromHighToLowHT = True,
                sVersion="_wAddiSF"
            )
            # SF computation from low to high HT
            hSF_step0_4, hGenLHE_HT_wHTSamplesStitch_4, hSF_step1_4, hGenLHE_HT_wHTSamplesStitch1_4 = calculate_HTStitchSFs(
                hGenLHE_HT = hGenLHE_HT_wAddiSF,
                HTRangeFull_min = HTRangeFull_min,
                HTRangeFull_max = HTRangeFull_max,
                HTBinEdges = HTBinEdges,
                calHTCorrSFFromHighToLowHT = False,
                sVersion="_wAddiSF"
            )


            fOut.cd(sampleName)
        
            hSF_step0_3.Write()
            hGenLHE_HT_wHTSamplesStitch_3.Write()
            hSF_step1_3.Write()
            hGenLHE_HT_wHTSamplesStitch1_3.Write()

            hSF_step0_4.Write()
            hGenLHE_HT_wHTSamplesStitch_4.Write()
            hSF_step1_4.Write()
            hGenLHE_HT_wHTSamplesStitch1_4.Write()


    fOut.Close()
    print(f"Wrote SFs in {sFOut}")
