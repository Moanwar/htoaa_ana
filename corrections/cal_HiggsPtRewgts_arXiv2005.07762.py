import os
import sys
#sys.argv.append( '-b-' )
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1D, TH2D, TH1, TF1, TEfficiency, TLegend
from ROOT import gROOT, gBenchmark, gRandom, gSystem, gStyle
import ctypes
import re
import argparse
import enum
import copy
import math
import numpy as np
import random

ROOT.gROOT.SetBatch(True)

sIpFileName = 'ipFileName'
sHistName   = 'histogramName'


def readHistFromFile(sFile, sHistNameFull, nRebinX=1, nRebinY=1):
    h = None
    f = TFile(sFile)
    if not f.IsOpen():
        print(f"{sFile} could not open")
        exit(0)
    
    h = copy.deepcopy( f.Get(sHistNameFull) )
    if h == None:
        print(f"Could not read {sHistNameFull} histogram from {sFile} file. \t\t\t *** ERROR ***")
        exit(0)
    f.Close()
    h.Rebin(nRebinX)
    return h

def readAndAddHistsFromFile(histograms_list, hTotName, nRebinX=1, nRebinY=1):
    hAdded = None
    for hist_i_dict in histograms_list:
        sFIn_      = hist_i_dict[sIpFileName]
        sHistName_ = hist_i_dict[sHistName]
        hTmp_ = readHistFromFile(sFIn_, sHistName_, nRebinX=nRebinX)    

        if hAdded == None: hAdded = hTmp_.Clone(hTotName)
        else:              hAdded.Add(hTmp_)
        #print(f"{sHistName_}: {hTmp_.Integral()}, {hAdded.Integral()}, ")

    return hAdded

def histIntegralAbvX(h, x1):
    integral = h.Integral(h.FindBin(x1), h.GetNbinsX())
    return integral

def histIntegralAbvX_wUncrt(h, x1):
    xLowBin = h.FindBin(x1)    
    return histIntegralAndError(h, xLowBin)

def histIntegralAndError(h, xLowBin=1, xHighBin=-1):
    if xHighBin<=0: xHighBin = h.GetNbinsX()
    eN = ctypes.c_double()
    N  = h.IntegralAndError(xLowBin, xHighBin, eN)
    return N, eN





if __name__ == "__main__":
    print("Running plot_HiggsPtRewgts", flush=True)
    gStyle.SetOptStat(0)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)

    Era = '2018'
    HiggsProdModes = [
        'ggHtoaato4b',
        'VBFHtoaato4b',
        'WHtoaato4b',
        'ZHtoaato4b',
        'ttHtoaato4b',        
    ]
    HiggsExclCrossSections = {
        'ggHtoaato4b': 48.61   * 0.057,
        'VBFHtoaato4b': 3.766  * 0.173,
        'WHtoaato4b':   1.358  * 0.132,
        'ZHtoaato4b':   0.880  * 0.129,
        'ttHtoaato4b':  0.5071 * 0.285,        
    }
    mAs = ['12','15','20','25','30','35','40','45','50','55','60']
    #mAs = ['12']

    pTStepsForCumulative = [400, 450, 500, 550, 600, 650, 700, 750, 800]


    sIpFileAllHist = '/eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts_arXiv2005.07762/%s/analyze_htoaa_stage1.root' % (Era)
    sHistNamesShort = ['hGenHiggsPt_Nom', 'hGenHiggsPt_wHiggsPtRewgt_Nom'] # 'hGenHiggsPt_Nom' 


    LuminosityPetYear = {
        '2018': 59.83,
    }

    sOpDir = '/eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts_arXiv2005.07762/%s/plots' % (Era)
    os.makedirs(sOpDir, exist_ok=True) 


    print(f"{Era = }")
    for HiggsProdMode in HiggsProdModes:
        print(f"\n\n{HiggsProdMode = }")

        for sHistNameShort in sHistNamesShort:
            print(f"\n{sHistNameShort = }")

            hHToAATo4B_list = []

            for mA in mAs:
                hHToAATo4B_list.append({
                    sIpFileName: sIpFileAllHist,
                    sHistName:  'evt/%s_mA_%s/%s'     % (HiggsProdMode,mA, sHistNameShort)
                })

            hHToAATo4B_mAAll_avg = readAndAddHistsFromFile(hHToAATo4B_list, '%s_%s_mAAll_avg'%(sHistNameShort, HiggsProdMode), nRebinX=1)
            print("Full integral (cross-section in fb-1) for %s before normalization_0: %f " % (HiggsProdMode, hHToAATo4B_mAAll_avg.Integral()))

            hHToAATo4B_mAAll_avg.Scale( 1. / len(mAs)) # Average distribution over all mAs

            print("Full integral (cross-section in fb-1) for %s before normalization_1: %f " % (HiggsProdMode, hHToAATo4B_mAAll_avg.Integral()))
            hHToAATo4B_mAAll_avg.Scale( 1. / LuminosityPetYear[Era]) # Normalize to 1 fb-1 luminosity

            print("Full integral (cross-section in fb-1) for %s before normalization: %f " % (HiggsProdMode, hHToAATo4B_mAAll_avg.Integral()))

            # Normalized distribution to Higgs cross-section for pT(Higgs)>150 GeV samples
            #hHToAATo4B_mAAll_avg.Scale( HiggsExclCrossSections[HiggsProdMode] / hHToAATo4B_mAAll_avg.Integral() )

            print("Full integral (cross-section in fb-1) for %s: %f " % (HiggsProdMode, hHToAATo4B_mAAll_avg.Integral()))
            for pT_thsh in pTStepsForCumulative:
                integralCumlt, eIntegralCumlt = histIntegralAbvX_wUncrt(hHToAATo4B_mAAll_avg, pT_thsh)
                #print("%g \t %f \t %f" % (pT_thsh, integralCumlt, eIntegralCumlt.value))
                print("%g \t %f " % (pT_thsh, integralCumlt))















