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


def readHistFromFile(sFile, sHistNameFull, nRebinX=1, nRebinY=1, maintainScale=0):
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

    integral_0 = h.Integral()
    h.Rebin(nRebinX)
    if maintainScale: # divide h by nRebinX to nullify effect of histogram rebining to maintain scale of the histogram
        h.Scale( 1/nRebinX )
    return h


def plotHistograms(histogram_dict, sCanvasName, sLegendHeader='', xLable='', yLable='', xRange=[], setLogY=1, saveAs='c.png'):
    colors_list = [1, 2, 4, 6, 28, 46, 7, 3]

    rnd = str(random.randint(0,10000))

    c1 = TCanvas(sCanvasName+rnd, sCanvasName+rnd, 600,500)
    c1.SetObjectStat(0)
    c1.SetLogy(setLogY)
    #c1.SetGrid()
    c1.cd()

    i = 0
    if  len(histogram_dict.keys()) == 1:
        leg = TLegend(0.45,0.79,0.85,0.85)
    elif len(histogram_dict.keys()) == 2:
        leg = TLegend(0.45,0.70,0.85,0.85)
    else:
        leg = TLegend(0.45,0.65,0.85,0.85)
    if sLegendHeader:  leg.SetHeader(sLegendHeader, 'C')
    for sHistName, h_ in histogram_dict.items():
        h_.SetMarkerStyle(20)
        h_.SetMarkerSize(0.5)
        h_.SetMarkerColor(colors_list[i])
        h_.SetLineColor(colors_list[i])     
        if  xLable != '': h_.GetXaxis().SetTitle(xLable)
        if  yLable != '': h_.GetYaxis().SetTitle(yLable) 
        if xRange:        h_.GetXaxis().SetRangeUser(xRange[0], xRange[1])     
        h_.GetXaxis().SetTitleSize(0.045) 
        h_.GetYaxis().SetTitleSize(0.045) 
        h_.GetYaxis().SetTitleOffset(1.1) 
        h_.SetTitle('')

        if i == 0: h_.Draw('')
        else:      h_.Draw('same')  
        leg.AddEntry(h_, sHistName, 'lep') # lep
        i += 1

    c1.cd()
    leg.Draw()
    c1.Update()
    c1.SaveAs(saveAs)

    #return copy.deepcopy( c1 )
    #return c1
    return







if __name__ == "__main__":
    print("Running plot_HiggsPtRewgts", flush=True)
    gStyle.SetOptStat(0)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)

    sOpDir = '/eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts/2018/plots'
    os.makedirs(sOpDir, exist_ok=True) 

    
    ## GGH ----------------------------------------------------------------------------------------------
    sFIn = '/eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts/2018/ggHHiggsPtRewgt_HToAATo4B.root'
    sHistName_Sig  = 'hGenHiggsPt_Nom_HToAATo4B_Stitch'
    sHistName_Ref  = 'hGenHiggsPt_Nom_HqtStitched'
    sHistName_kFct = 'hGenHiggsPt_Nom_Wgt_Hqt'
    sOpPlot_pt     = '%s/HiggsPt_ggH.png' % (sOpDir)
    sOpPlot_kFct   = '%s/kFactor_ggH.png' % (sOpDir)
    nRebinsX = 2
    xRange  = [0, 800]
    
    # Pt
    h_dict = {
        r'gg #rightarrow H signal (LO)': readHistFromFile(sFile=sFIn, sHistNameFull=sHistName_Sig, nRebinX=nRebinsX, maintainScale=1),
        r'Hqt v2.0':                        readHistFromFile(sFile=sFIn, sHistNameFull=sHistName_Ref, nRebinX=nRebinsX, maintainScale=1),
    }
    plotHistograms(
        histogram_dict = h_dict, 
        sCanvasName = r'cHiggsPt_ggH',
        sLegendHeader = r'gg #rightarrow H', 
        xLable='Higgs p_{T} [GeV]', yLable=r'1/d#sigma d#sigma / dp_{T} [GeV^{-1}]', 
        xRange = xRange,
        setLogY=1,
        saveAs=sOpPlot_pt
    )    

    

    # kFactor
    nRebinsX = 5
    h_dict = {
        r'gg #rightarrow H signal (LO)': readHistFromFile(sFile=sFIn, sHistNameFull=sHistName_kFct, nRebinX=nRebinsX, maintainScale=1)
    }
    plotHistograms(
        histogram_dict = h_dict, 
        sCanvasName = r'ckFactor_ggH',
        sLegendHeader = '', 
        xLable='Higgs p_{T} [GeV]', yLable=r'k-factor', 
        xRange = xRange,
        setLogY=0,
        saveAs=sOpPlot_kFct
    ) 
       



    ## VBFH ----------------------------------------------------------------------------------------------
    sFIn = '/eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts/2018/VBFHHiggsPtRewgt_HToAATo4B.root'
    sHistName_Sig  = 'hGenHiggsPt_Nom_HToAATo4B_Stitch'
    sHistName_Ref  = 'hGenHiggsPt_Nom_HiggsNLO'
    sHistName_kFct = 'hGenHiggsPt_Nom_Wgt_NLO'
    sOpPlot_pt     = '%s/HiggsPt_VBFH.png' % (sOpDir)
    sOpPlot_kFct   = '%s/kFactor_VBFH.png' % (sOpDir)
    nRebinsX = 2
    xRange  = [0, 650]
    
    # Pt
    h_dict = {
        r'VBF signal (LO)': readHistFromFile(sFile=sFIn, sHistNameFull=sHistName_Sig, nRebinX=nRebinsX, maintainScale=1),
        r'SM VBF (NLO)':    readHistFromFile(sFile=sFIn, sHistNameFull=sHistName_Ref, nRebinX=nRebinsX, maintainScale=1),
    }
    plotHistograms(
        histogram_dict = h_dict, 
        sCanvasName = r'cHiggsPt_VBFH',
        sLegendHeader = r'VBF', 
        xLable='Higgs p_{T} [GeV]', yLable=r'1/d#sigma d#sigma / dp_{T} [GeV^{-1}]', 
        xRange = xRange,
        setLogY=1,
        saveAs=sOpPlot_pt
    )    

    

    # kFactor
    nRebinsX = 5
    h_dict = {
        r'VBF signal (LO)': readHistFromFile(sFile=sFIn, sHistNameFull=sHistName_kFct, nRebinX=nRebinsX, maintainScale=1)
    }
    plotHistograms(
        histogram_dict = h_dict, 
        sCanvasName = r'ckFactor_VBFH',
        sLegendHeader = '', 
        xLable='Higgs p_{T} [GeV]', yLable=r'k-factor', 
        xRange = xRange,
        setLogY=0,
        saveAs=sOpPlot_kFct
    )    

    
    ## ttH ----------------------------------------------------------------------------------------------
    sFIn = '/eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts/2018/ttHHiggsPtRewgt_HToAATo4B.root'
    sHistName_Sig  = 'hGenHiggsPt_Nom_HToAATo4B_Stitch'
    sHistName_Ref  = 'hGenHiggsPt_Nom_HiggsNLO'
    sHistName_kFct = 'hGenHiggsPt_Nom_Wgt_NLO'
    sOpPlot_pt     = '%s/HiggsPt_ttH.png' % (sOpDir)
    sOpPlot_kFct   = '%s/kFactor_ttH.png' % (sOpDir)
    nRebinsX = 2
    xRange  = [0, 650]
    
    # Pt
    h_dict = {
        r't#bar{t}H signal (LO)': readHistFromFile(sFile=sFIn, sHistNameFull=sHistName_Sig, nRebinX=nRebinsX, maintainScale=1),
        r'SM t#bar{t}H (NLO)':    readHistFromFile(sFile=sFIn, sHistNameFull=sHistName_Ref, nRebinX=nRebinsX, maintainScale=1),
    }
    plotHistograms(
        histogram_dict = h_dict, 
        sCanvasName = r'cHiggsPt_ttH',
        sLegendHeader = r't#bar{t}H', 
        xLable='Higgs p_{T} [GeV]', yLable=r'1/d#sigma d#sigma / dp_{T} [GeV^{-1}]', 
        xRange = xRange,
        setLogY=1,
        saveAs=sOpPlot_pt
    )    

    

    # kFactor
    nRebinsX = 5
    h_dict = {
        r't#bar{t}H signal (LO)': readHistFromFile(sFile=sFIn, sHistNameFull=sHistName_kFct, nRebinX=nRebinsX, maintainScale=1)
    }
    plotHistograms(
        histogram_dict = h_dict, 
        sCanvasName = r'ckFactor_VBFH',
        sLegendHeader = '', 
        xLable='Higgs p_{T} [GeV]', yLable=r'k-factor', 
        xRange = xRange,
        setLogY=0,
        saveAs=sOpPlot_kFct
    )    
    