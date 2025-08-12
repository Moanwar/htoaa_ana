'''
Python environment with ROOT and cmsstyle: https://cms-analysis.docs.cern.ch/guidelines/plotting/#installation, https://cms-sw.github.io/venv.html
    cd /afs/cern.ch/work/s/ssawant/private/htoaa/cmsplots/CMSSW_14_1_0_pre4/src
    cmsenv

    python3 Plots
'''

import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1D, TH2D, TH1, TF1, TEfficiency, TLegend
from ROOT import gROOT, gBenchmark, gRandom, gSystem, gStyle
import cmsstyle
import copy

ROOT.gROOT.SetBatch(ROOT.kTRUE)

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

    return h


if __name__ == '__main__':

    # set the global ROOT style to CMS style
    cmsstyle.setCMSStyle()
    cmsstyle1 = cmsstyle.getCMSStyle()
    cmsstyle1.SetPadRightMargin(0.10)

    # set the luminosity, the COM energy, the Run period to show in the canvases
    cmsstyle.SetLumi(34.8)
    cmsstyle.SetEnergy(13)
    # default extra text is "Preliminary", set it to an empty string to remove it
    cmsstyle.SetExtraText('Preliminary')

    sInFile = "data/correction/mc/BtagSF/2018/jetBtagEfficiency.root"
    sHisto = "hJetBtagEffi_b_TT_Presel"

    h = readHistFromFile(sInFile, sHisto)

    square = True
    c = cmsstyle.cmsCanvas( #https://github.com/cms-cat/cmsstyle/blob/master/src/cmsstyle/cmsstyle.py#L962
        canvName="c",
        x_min=h.GetXaxis().GetXmin(),
        x_max=h.GetXaxis().GetXmax(),
        y_min=h.GetYaxis().GetXmin(),
        y_max=h.GetYaxis().GetXmax(),
        nameXaxis="[T]",
        nameYaxis="eta",
        square=False, #cmsstyle.kSquare,
        iPos=0,
        extraSpace=0.01,
        with_z_axis=True,
        #scaleLumi=1,
        #yTitOffset=None,     
    )   
    #h.GetZaxis().SetTitle("Efficiency")
    #h.GetZaxis().SetTitleOffset(1.4 if square else 0.8)
    c.SetLogx(1)
    #c.SetRightMargin(0.25)

    cmsstyle.SetCMSPalette()

    gStyle.SetPaintTextFormat("4.2f");
    h.SetMarkerSize(1.5)

    h.Draw("same colz TEXT")

    # Set a new palette
    #cmsstyle.SetAlternative2DColor(h)

    # Allow to adjust palette position
    cmsstyle.UpdatePalettePosition(h, c)
    #cmsstyle.UpdatePalettePosition(h, c, 0.85,0.99,0.1,0.9)
    #cmsstyle.UpdatePalettePosition(hist=h, canv=None, X1=0.8, X2=0.9, Y1=0.13, Y2=0.95, isNDC=True)   

    cmsstyle.UpdatePad() 

    cmsstyle.SaveCanvas(c, './test.pdf')

        
        
  