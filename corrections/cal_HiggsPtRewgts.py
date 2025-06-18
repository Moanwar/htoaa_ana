
import os
import sys
#sys.argv.append( '-b-' )
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TH1, TF1, TEfficiency, TLegend
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import ctypes
import re
import argparse
import enum
import copy
import math
import numpy as np

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

def histIntegralAndError(h, xLowBin=1, xHighBin=-1):
    if xHighBin<=0: xHighBin = h.GetNbinsX()
    eN = ctypes.c_double()
    N  = h.IntegralAndError(xLowBin, xHighBin, eN)
    return N, eN

def makeRatioHist(hNume, hDenom, sRatioName):
    hRatio = hNume.Clone(sRatioName)
    hRatio.Divide(hNume, hDenom)
    return hRatio

def makeCumulativeHist(h):
    hCumulative = h.Clone('%s_cumulative' % (h.GetName()))
    nBinsX = h.GetNbinsX()
    for iBin in range(1, nBinsX):
        eN = ctypes.c_double()
        N = h.IntegralAndError(iBin, nBinsX, eN)
        hCumulative.SetBinContent(iBin, N)
        hCumulative.SetBinError(iBin, eN)
    return hCumulative

def fitHistogram(h, sFitFuncLocal, FitRangeLocal, sFitFull, FitRangeFull, sCanvasName=''):
    c1 = TCanvas(sCanvasName, sCanvasName, 600,500)
    #c1.SetLogy(setLogY)
    c1.SetGrid()
    c1.cd()

    h.Draw()

    fLocal = TF1('fLocal_%s'%(h.GetName()), sFitFuncLocal, FitRangeLocal[0], FitRangeLocal[1])
    fLocal.SetLineColor(2)
    fLocal.SetLineWidth(2)
    fResultsPtr = h.Fit(fLocal, "RVS")
    fitParams = fLocal.GetParameters()
    print(f"{fResultsPtr = }")
    print(f'{fResultsPtr.Print("V") = }')
    print(f"{fLocal.GetName() = }, {fitParams = }, {type(fitParams) = }  ")
    print(f"{len(fitParams) = }")
    print(f"{fitParams[0] = }")
    print(f"{fitParams[1] = }")
    print(f"{fitParams[2] = }")
    print(f"{fitParams[3] = }")
    #print(f"{np.frombuffer(fitParams, dtype=np.float32, count=4) = }")
    
    

    return copy.deepcopy( c1 )
    



def plotHistograms(histogram_dict, sCanvasName, xLable='', yLable='', setLogY=1):
    colors_list = [1, 2, 4, 6, 28, 46, 7, 3]

    c1 = TCanvas(sCanvasName, sCanvasName, 600,500)
    c1.SetLogy(setLogY)
    c1.SetGrid()
    c1.cd()

    i = 0
    leg = TLegend(0.5,0.85,0.99,0.99)
    for sHistName, h_ in histogram_dict.items():
        h_.SetMarkerStyle(20)
        h_.SetMarkerSize(0.5)
        h_.SetMarkerColor(colors_list[i])
        h_.SetLineColor(colors_list[i])     
        if  xLable != '': h_.GetXaxis().SetTitle(xLable)
        if  yLable != '': h_.GetYaxis().SetTitle(yLable)        

        if i == 0: h_.Draw()
        else:      h_.Draw('same')  
        leg.AddEntry(h_, sHistName, 'lep')
        i += 1

    c1.cd()
    leg.Draw()

    return copy.deepcopy( c1 )



def stitchInclAndExclHistogramsAlongXaxis(hInclusive, hExclusive, X_toStitchHistograms):
    # Normalize hInclusive to unit area
    print(f"Before: {histIntegralAndError(hInclusive) = }")
    hInclusive.Scale( 1./hInclusive.Integral() )
    print(f"After:  {histIntegralAndError(hInclusive) = }")

    # Scale hHighPt (= h_Excl) to match hInclusive area for Pt > X_toStitchHistograms (400 GeV)
    integral1_Incl     = histIntegralAbvX(hInclusive,     X_toStitchHistograms)
    integral1_Excl = histIntegralAbvX(hExclusive, X_toStitchHistograms)
    hExclusive.Scale( integral1_Incl / integral1_Excl )

    # hStitch =  
    #           hInclusive (when Pt <  X_toStitchHistograms)
    #           +
    #           hHighPt    (when Pt >= X_toStitchHistograms)
    hHToAATo4B_Stitch = hInclusive.Clone(hInclusive.GetName().replace('Incl','Stitch'))
    print(f"Before: {histIntegralAndError(hHToAATo4B_Stitch) = }")
    for iBin in range(hHToAATo4B_Stitch.FindBin(X_toStitchHistograms), hHToAATo4B_Stitch.GetNbinsX()+1):
        hHToAATo4B_Stitch.SetBinContent(iBin, hExclusive.GetBinContent(iBin))
        hHToAATo4B_Stitch.SetBinError(iBin,   hExclusive.GetBinError(iBin))
    print(f"After:  {histIntegralAndError(hHToAATo4B_Stitch) = }")

    return hInclusive, hExclusive, hHToAATo4B_Stitch


def makeHqt_HiggsPt_Hist(sFInHqtHiggsSpectrum, sFOutHqtHiggsHist):
    hHiggsCrosssecStitchedInPt     = TH1F('hGenHiggsPt_Nom_HqtStitched',     'Calculated with Hqt - stitched',  2000,            0,               2000)
    hHiggsCrosssecStitchedInLog2Pt = TH1F('hGenHiggsLog2Pt_Nom_HqtStitched', 'Calculated with Hqt- stitched',   200,  math.log2(1),    math.log2(2000))
    hHiggsCrosssecMatchedInPt      = TH1F('hGenHiggsPt_Nom_HqtMatched',      'Calculated with Hqt - matched',  2000,            0,               2000)
    hHiggsCrosssecMatchedInLog2Pt  = TH1F('hGenHiggsLog2Pt_Nom_HqtMatched',  'Calculated with Hqt - matched',   200,  math.log2(1),    math.log2(2000))
    
    ## Read Hqt Higgs Pt spectrum and store it into histogram
    HiggsCrosssection_dict = {}
    with open(sFInHqtHiggsSpectrum) as fInHqtHiggsSpectrum:
        for line in fInHqtHiggsSpectrum:
            data_list = re.split(r'\s+', line.rstrip()) # split line into list of words. 
            if '' in data_list:  data_list.remove('')   # remove empy string if any
            print(f"{data_list = }")

            if data_list[0] == '(': continue  # comment line starts with '(' charecter

            # data_list columns: 'qt', 'res', 'asym', 'fixor', 'matched', 'switched'
            # Read 'qt': Higgs pT,      'switched': Higgs differential cross-section NNLO+LO
            # https://www.physik.uzh.ch/~grazzini/codes/note20.pdf#page=4
            Pt           = float(data_list[0])
            Crosssection_stitched = float(data_list[-1])
            Crosssection_matched  = float(data_list[-2])
            print(f"{Pt}, \t {Crosssection_matched = }, \t {Crosssection_stitched = }")

            # skip -ve cross-section readings
            if (Crosssection_matched < 0) and (Crosssection_stitched < 0): continue

            HiggsCrosssection_dict[Pt] = {
                'matched':  Crosssection_matched,
                'switched': Crosssection_stitched
            }
            
    xsErrorFraction = 0.05

    # Set cross-section for hHiggsCrosssecInPt bins
    hTmp_dict = {
        'matched':  hHiggsCrosssecMatchedInPt,
        'switched': hHiggsCrosssecStitchedInPt
    }
    for XSType, hTmp_ in hTmp_dict.items():
        for iBin in range(1, hTmp_.GetNbinsX()+1):
            # calculate average of cross-sections over pT points fall in the histogram pT bin
            PtPoints_list = []
            for Pt in HiggsCrosssection_dict.keys():
                if Pt >= hTmp_.GetXaxis().GetBinLowEdge(iBin) and \
                   Pt <  hTmp_.GetXaxis().GetBinUpEdge(iBin):
                    PtPoints_list.append(Pt)

            #hTmp_.SetBinContent(iBin, avgXS)
            nPtPoints = len(PtPoints_list)
            if nPtPoints == 0: continue

            # Pick up middle point in Pt that lie within bin of hitogram
            iPtPointSel = int( np.median( range(nPtPoints) ) )
            PtPointSel  = PtPoints_list[iPtPointSel]
            XS = HiggsCrosssection_dict[PtPointSel][XSType]
            if XS < 0: continue
            hTmp_.SetBinContent(iBin, XS)
            hTmp_.SetBinError(  iBin, XS * xsErrorFraction) # error: 5%

    
    # Set cross-section for hHiggsCrosssecInLog2Pt bins
    hTmp_dict = {
        'matched':  hHiggsCrosssecMatchedInLog2Pt,
        'switched': hHiggsCrosssecStitchedInLog2Pt
    }
    for XSType, hTmp_ in hTmp_dict.items():
        for iBin in range(1, hTmp_.GetNbinsX()+1):
            # calculate average of cross-sections over pT points fall in the histogram pT bin
            PtPoints_list = []
            for Pt in HiggsCrosssection_dict.keys():
                if math.log2(Pt) >= hTmp_.GetXaxis().GetBinLowEdge(iBin) and \
                   math.log2(Pt) <  hTmp_.GetXaxis().GetBinUpEdge(iBin):
                    PtPoints_list.append(Pt)

            nPtPoints = len(PtPoints_list)
            if nPtPoints == 0: continue

            # Pick up middle point in Pt that lie within bin of hitogram
            iPtPointSel = int( np.median( range(nPtPoints) ) )
            PtPointSel  = PtPoints_list[iPtPointSel]
            XS = HiggsCrosssection_dict[PtPointSel][XSType]
            if XS < 0: continue
            hTmp_.SetBinContent(iBin, XS)
            hTmp_.SetBinError(  iBin, XS * xsErrorFraction) # error: 5%
    



    ## Write Higgs Pt histograms into output root file
    sDir_ = os.path.dirname( os.path.realpath(sFOutHqtHiggsHist) )
    if not os.path.exists(sDir_): 
        os.makedirs( sDir_ )
    fOutHqtHiggsHist = TFile(sFOutHqtHiggsHist, 'recreate')
    fOutHqtHiggsHist.cd()
    hHiggsCrosssecStitchedInPt.Write()
    hHiggsCrosssecStitchedInLog2Pt.Write()   
    hHiggsCrosssecMatchedInPt.Write()
    hHiggsCrosssecMatchedInLog2Pt.Write()   
    fOutHqtHiggsHist.Close()

    print(f"\nWrote Higgs pT cross-section (Hqt) histograms into {sFOutHqtHiggsHist}.")   


#def calculateGGFHiggsPtRewgt(
#    hHToAATo4B_Incl_list, hHToAATo4B_Excl_list,  hGGFHTo2B_Incl_list, 
#):
            





    

if __name__ == "__main__":
    print("Running cal_HiggsPtRewgts", flush=True)

    parser = argparse.ArgumentParser(description='cal_HiggsPtRewgts')
    parser.add_argument('-era', dest='era',   type=str, default='2018',                    choices=['2016','2017','2018'], required=False)
    parser.add_argument('-prodMode', dest='prodMode',   type=str, default='ggH',           choices=['ggH','VBFH', 'WH', 'WplusH','WminusH', 'ZH', 'ttH'], required=True)
    parser.add_argument('-useHToAATo4TauSignal',  action='store_true', default=False, help='Use HToAATo4Tau signal for cross-checks')
    args=parser.parse_args()
    print("args: {}".format(args), flush=True)
    era                     = args.era
    productionMode          = args.prodMode   
    useHToAATo4TauSignal    = args.useHToAATo4TauSignal # Default: False


    # Hqt Higgs cross-section in pT
    # /afs/cern.ch/work/s/ssawant/private/htoaa/HqT_HiggsPtCode/HqT2.0/HqTspectrum_13TeV.out
    sFInHqtHiggsSpectrum_mTopInfinite = "/afs/cern.ch/work/s/ssawant/private/htoaa/HqT_HiggsPtCode/HqT2.0/HqTspectrum_13TeV_pt1To2kGeVBin1GeV.out"
    sFInHqtHiggsSpectrum_mTopFinite   = "/afs/cern.ch/work/s/ssawant/private/htoaa/HqT_HiggsPtCode/HqT2.0/HqTspectrum_13TeV_mTopFinite_pt1To2kGeVBin1GeV.out"
    sFInHiggsSpectrum_GGH_NNLO        = "/eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts/2018/GenHiggsPt_GGH_NNLO.root"
    sOutDir                           = "/eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts/%s" % (era)
    sFOutHqtHiggsHist_mTopInfinite    = "%s/Hqt_HiggsPtHist_mTopInfinite.root" % (sOutDir)
    sFOutHqtHiggsHist_mTopFinite      = "%s/Hqt_HiggsPtHist_mTopFinite.root" % (sOutDir)

    sIpFileAllHist = '/eos/cms/store/user/ssawant/htoaa/analysis/20250603_CalHiggsPtRewgt_1/2018/analyze_htoaa_stage1.root'
    sHistNameShort_list = ['hGenHiggsPt_Nom'] # ['hGenHiggsPt_Nom', 'hGenHiggsLog2Pt_Nom']  'hGenHiggsPt_Nom', 'hGenHiggsPt_wHiggsPtRewgt_Nom'
    nRebinX_dict = {'hGenHiggsPt_Nom': 10, 'hGenHiggsLog2Pt_Nom': 4}
    HiggsPtPoint_toStitchHToAATo4BSamples = 300 # GeV
    HiggsPtPoint_toStitchHTo2BSamples     = 300 # GeV
    HiggsPt_fitRange                  = [150, 1100] # GeV
    HiggsPt_fullRange                 = [  1, 2000] # GeV
    sFOutHiggsPtRewgt    = "%s/%sHiggsPtRewgt_%s.root" % (sOutDir, productionMode, 'HToAATo4Tau' if useHToAATo4TauSignal else 'HToAATo4B')
    xsHiggs_dict = {
        'ggH': 48.61 * 1000, # 33.8 * 1000, # 48.61 * 1000, # fb
    }
    


    sProcessesHiggsLO = productionMode
    if productionMode == 'WplusH' or productionMode == 'WminusH':  sProcessesHiggsLO = 'WH'
    
    sProcessesHiggsNLO = []
    if productionMode == 'VBFH':    sProcessesHiggsNLO = ['VBFHToBB_powheg'] 
    if productionMode == 'WH':      sProcessesHiggsNLO = ['WplusHToBBQQ', 'WplusHToBBLNu', 'WminusHToBBQQ', 'WminusHToBBLNu'] 
    if productionMode == 'WplusH':  sProcessesHiggsNLO = ['WplusHToBBQQ', 'WplusHToBBLNu'] 
    if productionMode == 'WminusH': sProcessesHiggsNLO = ['WminusHToBBQQ', 'WminusHToBBLNu'] 
    #if productionMode == 'ZH':      sProcessesHiggsNLO = ['ZHToBBX'] # powheg
    if productionMode == 'ZH':      sProcessesHiggsNLO = ['ZHToMuMuG'] # dalitz amcatnloFXFX
    if productionMode == 'ttH':     sProcessesHiggsNLO = ['ttH'] 
    

    
    os.makedirs( os.path.dirname( os.path.realpath(sFOutHiggsPtRewgt) ), exist_ok=True )
    if os.path.exists(sFOutHiggsPtRewgt): os.remove(sFOutHiggsPtRewgt) # Delete output file if it exists

    if productionMode == 'ggH' and 1==1:
        # Read Hqt Higgs Pt spectrum and store it into histogram
        makeHqt_HiggsPt_Hist(sFInHqtHiggsSpectrum_mTopInfinite, sFOutHqtHiggsHist_mTopInfinite)
        makeHqt_HiggsPt_Hist(sFInHqtHiggsSpectrum_mTopFinite,   sFOutHqtHiggsHist_mTopFinite)
    
    for sHistNameShort in sHistNameShort_list:
        hHToAATo4B_Excl_list = []
        hHToAATo4B_Incl_list     = []
        mAs = [ '12', '15', '20', '25', '30', '35', '40', '45', '50', '55', '60'  ]
        sHistName_HToAATo4B_Excl_0 = 'evt/%stoaato4b_mA_%s/%s'
        sHistName_HToAATo4B_Incl_0 = 'evt/%stoaato4b_Incl_mA_%s/%s'
        if useHToAATo4TauSignal:
            mAs = ['All']
            sHistName_HToAATo4B_Incl_0 = 'evt/%stoaato4tau_mA_%s/%s'    
            if productionMode in ['WH', 'WplusH', 'WminusH', 'ZH']:
                sProcessesHiggsLO = 'VH'        
        for mA in mAs:
            hHToAATo4B_Excl_list.append({
                sIpFileName: sIpFileAllHist,
                sHistName:  sHistName_HToAATo4B_Excl_0     % (sProcessesHiggsLO,mA, sHistNameShort)
            })
            hHToAATo4B_Incl_list.append({
                sIpFileName: sIpFileAllHist,
                sHistName:  sHistName_HToAATo4B_Incl_0     % (sProcessesHiggsLO,mA, sHistNameShort)
            })

        hGGFHTo2B_Excl_list = [
            {
                sIpFileName: sIpFileAllHist,
                sHistName: 'evt/GluGluHToBB_Pt-200ToInf/%s' % (sHistNameShort)
            }
        ]
        hGGFHTo2B_Incl_list  = [
            {
                sIpFileName: sIpFileAllHist,
                sHistName: 'evt/GluGluHToBB_Incl/%s' % (sHistNameShort)
            }
        ]

        hHiggsNLO_list = []
        for sProcessHiggsNLO in sProcessesHiggsNLO:
            hHiggsNLO_list.append({
                sIpFileName: sIpFileAllHist,
                sHistName: 'evt/%s/%s'      % (sProcessHiggsNLO, sHistNameShort)
            })



        nRebinX = nRebinX_dict[sHistNameShort]
        HiggsPtPoint_toStitchHToAATo4BSamples_toUse = math.log2(HiggsPtPoint_toStitchHToAATo4BSamples) if 'Log2Pt' in sHistNameShort else HiggsPtPoint_toStitchHToAATo4BSamples
        HiggsPtPoint_toStitchHTo2BSamples_toUse     = math.log2(HiggsPtPoint_toStitchHTo2BSamples)     if 'Log2Pt' in sHistNameShort else HiggsPtPoint_toStitchHTo2BSamples
        HiggsPtPoint_toStitchHToAATo4BSamples_toUse = HiggsPtPoint_toStitchHToAATo4BSamples
        HiggsPtPoint_toStitchHTo2BSamples_toUse     = HiggsPtPoint_toStitchHTo2BSamples 
        HiggsPt_fitRange_toUse                      = HiggsPt_fitRange
        HiggsPt_fullRange_toUse                     = HiggsPt_fullRange
        if 'Log2Pt' in sHistNameShort:
            HiggsPtPoint_toStitchHToAATo4BSamples_toUse = math.log2(HiggsPtPoint_toStitchHToAATo4BSamples)
            HiggsPtPoint_toStitchHTo2BSamples_toUse     = math.log2(HiggsPtPoint_toStitchHTo2BSamples) 
            HiggsPt_fitRange_toUse                      = np.log2(HiggsPt_fitRange_toUse)
            HiggsPt_fullRange_toUse                     = np.log2(HiggsPt_fullRange)


        ## HToAATo4B histogram ------------
        print(f"\n\nWorking with HToAATo4B")
        if not useHToAATo4TauSignal:
            hHToAATo4B_Excl = readAndAddHistsFromFile(hHToAATo4B_Excl_list, '%s_HToAATo4B_Excl'%(sHistNameShort), nRebinX=nRebinX)
            hHToAATo4B_Incl = readAndAddHistsFromFile(hHToAATo4B_Incl_list, '%s_HToAATo4B_Incl'%(sHistNameShort), nRebinX=nRebinX)
            # Stitch hInclusive and hExclusive histograms
            hHToAATo4B_Incl, hHToAATo4B_Excl, hHToAATo4B_Stitch = stitchInclAndExclHistogramsAlongXaxis(
                hInclusive = hHToAATo4B_Incl, 
                hExclusive = hHToAATo4B_Excl, 
                X_toStitchHistograms = HiggsPtPoint_toStitchHToAATo4BSamples_toUse
                )
            print(f"After: {histIntegralAndError(hHToAATo4B_Incl) = }, {histIntegralAndError(hHToAATo4B_Excl) = }, {histIntegralAndError(hHToAATo4B_Stitch) = }, ")
        else:
            hHToAATo4B_Stitch = readAndAddHistsFromFile(hHToAATo4B_Incl_list, '%s_HToAATo4Tau_mA_All'%(sHistNameShort), nRebinX=nRebinX)
            hHToAATo4B_Stitch.Scale( 1./hHToAATo4B_Stitch.Integral() )
            hHToAATo4B_Incl = hHToAATo4B_Stitch.Clone('%s_cloneIncl'%hHToAATo4B_Stitch.GetName())
            hHToAATo4B_Excl = hHToAATo4B_Stitch.Clone('%s_cloneExcl'%hHToAATo4B_Stitch.GetName())

        hHiggsNLO = None
        if productionMode == 'ggH':
            ## GGFHTo2B histogram ------------
            print(f"\n\nWorking with GGFHTo2B")
            hGGFHTo2B_Excl = readAndAddHistsFromFile(hGGFHTo2B_Excl_list, '%s_GGFHTo2B_Excl'%(sHistNameShort), nRebinX=nRebinX)
            hGGFHTo2B_Incl = readAndAddHistsFromFile(hGGFHTo2B_Incl_list, '%s_GGFHTo2B_Incl'%(sHistNameShort), nRebinX=nRebinX)
            # Stitch hInclusive and hExclusive histograms
            hGGFHTo2B_Incl, hGGFHTo2B_Excl, hGGFHTo2B_Stitch = stitchInclAndExclHistogramsAlongXaxis(
                hInclusive = hGGFHTo2B_Incl, 
                hExclusive = hGGFHTo2B_Excl, 
                X_toStitchHistograms = HiggsPtPoint_toStitchHTo2BSamples_toUse
                )
            print(f"After: {histIntegralAndError(hGGFHTo2B_Incl) = }, {histIntegralAndError(hGGFHTo2B_Excl) = }, {histIntegralAndError(hGGFHTo2B_Stitch) = }, ")

            ## Hqt Higgs Pt histogram
            hGGFH_Hqt_stitched_mTopInfinite = readHistFromFile(sFOutHqtHiggsHist_mTopInfinite, '%s_HqtStitched'%(sHistNameShort), nRebinX=nRebinX) 
            hGGFH_Hqt_stitched_mTopInfinite.Scale( 1./hGGFH_Hqt_stitched_mTopInfinite.Integral() ) # Normalize Hqt histogram to unit area
            hGGFH_Hqt_stitched = readHistFromFile(sFOutHqtHiggsHist_mTopFinite, '%s_HqtStitched'%(sHistNameShort), nRebinX=nRebinX)            
            hGGFH_Hqt_stitched.Scale( 1./hGGFH_Hqt_stitched.Integral() ) # Normalize Hqt histogram to unit area

            ## GGFH NNLO MC (Yihui) histogram
            hGGFH_NNLO = readHistFromFile(sFInHiggsSpectrum_GGH_NNLO, sHistNameShort, nRebinX=nRebinX)
            hGGFH_NNLO.Scale( 1./hGGFH_NNLO.Integral() ) # Normalize Hqt histogram to unit area


            
            hHiggsNLO = hGGFHTo2B_Stitch

            hHToAATo4B_Stitch_cloneXSNorm               = hHToAATo4B_Stitch.Clone( '%s_XSNorm'%(hHToAATo4B_Stitch.GetName()))
            hHiggsNLO_cloneXSNorm                       = hHiggsNLO.Clone(         '%s_XSNorm'%(hHiggsNLO.GetName()))
            hGGFH_Hqt_stitched_cloneXSNorm              = hGGFH_Hqt_stitched.Clone('%s_XSNorm'%(hGGFH_Hqt_stitched.GetName()))
            hGGFH_Hqt_stitched_mTopInfinite_cloneXSNorm = hGGFH_Hqt_stitched_mTopInfinite.Clone('%s_mTopInfinite_XSNorm'%(hGGFH_Hqt_stitched_mTopInfinite.GetName()))
            hGGFH_NNLO_cloneXSNorm                      = hGGFH_NNLO.Clone('%s_XSNorm'%(hGGFH_NNLO.GetName()))
            for h_ in [hHToAATo4B_Stitch_cloneXSNorm, hHiggsNLO_cloneXSNorm, hGGFH_Hqt_stitched_cloneXSNorm, hGGFH_Hqt_stitched_mTopInfinite_cloneXSNorm, hGGFH_NNLO_cloneXSNorm]:
                h_.Scale( xsHiggs_dict['ggH'] )
            histograms_dict_ = {
                r"$gg\to H\to aa \to 4b, \sigma = %.2f pb$" %(xsHiggs_dict['ggH']/1000):            hHToAATo4B_Stitch_cloneXSNorm,
                r"$gg\to H\to 2b, NLO, \sigma = %.2f pb$" %(xsHiggs_dict['ggH']/1000):              hHiggsNLO_cloneXSNorm,
                r"$gg\to H$ HqT2.0$, mTop infinite, $\sigma = %.2f pb" %(xsHiggs_dict['ggH']/1000): hGGFH_Hqt_stitched_mTopInfinite_cloneXSNorm,
                r"$gg\to H$ HqT2.0$, mTop finite, $\sigma = %.2f pb" %(xsHiggs_dict['ggH']/1000):   hGGFH_Hqt_stitched_cloneXSNorm,
                r"$gg\to H\to 2b, NNLO, \sigma = %.2f pb$" %(xsHiggs_dict['ggH']/1000):             hGGFH_NNLO_cloneXSNorm,
            }
            cCompareSamplesXSNorm = plotHistograms(histograms_dict_, 'c%s_CompareSamplesXSNorm'%(sHistNameShort), yLable=r'$\frac{d\sigma}{dpT}$')

            # Plot cumulative differential cross-section in 10 GeV bins
            hHToAATo4B_Stitch_cloneXSCumul  = makeCumulativeHist(hHToAATo4B_Stitch_cloneXSNorm)
            hHiggsNLO_cloneXSCumul          = makeCumulativeHist(hHiggsNLO_cloneXSNorm)
            hGGFH_Hqt_stitched_cloneXSCumul = makeCumulativeHist(hGGFH_Hqt_stitched_cloneXSNorm)
            hGGFH_Hqt_stitched_mTopInfinite_cloneXSCumul = makeCumulativeHist(hGGFH_Hqt_stitched_mTopInfinite_cloneXSNorm)
            histograms_dict_ = {
                r"$gg\to H\to aa \to 4b, \sigma = %.2f pb$" %(xsHiggs_dict['ggH']/1000):            hHToAATo4B_Stitch_cloneXSCumul,
                r"$gg\to H\to 2b, \sigma = %.2f pb$" %(xsHiggs_dict['ggH']/1000):                   hHiggsNLO_cloneXSCumul,
                r"$gg\to H$ HqT2.0$, mTop finite, $\sigma = %.2f pb" %(xsHiggs_dict['ggH']/1000):                  hGGFH_Hqt_stitched_cloneXSCumul,
                r"$gg\to H$ HqT2.0$, mTop infinite, $\sigma = %.2f pb" %(xsHiggs_dict['ggH']/1000):                  hGGFH_Hqt_stitched_mTopInfinite_cloneXSCumul,
            }
            cCompareSamplesXSCumul = plotHistograms(histograms_dict_, 'c%s_CompareSamplesXSCumulative'%(sHistNameShort), yLable=r'$\frac{d\sigma}{dpT}  [fb/GeV]$')

            '''
            # Plot cumulative differential cross-section in 50 GeV bins
            hHToAATo4B_Stitch_cloneXSCumul_1  = hHToAATo4B_Stitch_cloneXSNorm.Clone('%s_50GeVBin'%(hHToAATo4B_Stitch_cloneXSNorm.GetName()))
            hHiggsNLO_cloneXSCumul_1          = hHiggsNLO_cloneXSNorm.Clone('%s_50GeVBin'%(hHiggsNLO_cloneXSNorm.GetName()))
            hGGFH_Hqt_stitched_cloneXSCumul_1 = hGGFH_Hqt_stitched_cloneXSNorm.Clone('%s_50GeVBin'%(hGGFH_Hqt_stitched_cloneXSNorm.GetName()))
            nRebinsX_1 = int( 50 / nRebinX )
            hHToAATo4B_Stitch_cloneXSCumul_1.Rebin( nRebinsX_1 )
            hHiggsNLO_cloneXSCumul_1.Rebin( nRebinsX_1 )
            hGGFH_Hqt_stitched_cloneXSCumul_1.Rebin( nRebinsX_1 ) 
            hHToAATo4B_Stitch_cloneXSCumul_1  = makeCumulativeHist(hHToAATo4B_Stitch_cloneXSCumul_1)
            hHiggsNLO_cloneXSCumul_1          = makeCumulativeHist(hHiggsNLO_cloneXSCumul_1)
            hGGFH_Hqt_stitched_cloneXSCumul_1 = makeCumulativeHist(hGGFH_Hqt_stitched_cloneXSCumul_1)
            histograms_dict_ = {
                r"$gg\to H\to aa \to 4b, \sigma = %.2f pb$" %(xsHiggs_dict['ggH']/1000):            hHToAATo4B_Stitch_cloneXSCumul_1,
                r"$gg\to H\to 2b, \sigma = %.2f pb$" %(xsHiggs_dict['ggH']/1000):                   hHiggsNLO_cloneXSCumul_1,
                r"$gg\to H$ HqT2.0, \sigma = %.2f pb" %(xsHiggs_dict['ggH']/1000):                  hGGFH_Hqt_stitched_cloneXSCumul_1,
            }
            cCompareSamplesXSCumul_1 = plotHistograms(histograms_dict_, 'c%s_CompareSamplesXSCumulative_50GeVBin'%(sHistNameShort), yLable=r'$\frac{d\sigma}{dpT}$')
            '''

            
            histograms_dict_ = {
                r"$gg\to H\to aa \to 4b$":            hHToAATo4B_Stitch,
                r"$gg\to H\to 2b$":                   hHiggsNLO,
                r"$gg\to H$ HqT2.0":                  hGGFH_Hqt_stitched,
            }
        else: # Production modes other than ggH
            hHiggsNLO = readAndAddHistsFromFile(hHiggsNLO_list, '%s_HiggsNLO'%(sHistNameShort), nRebinX=nRebinX)
            hHiggsNLO.Scale( 1./hHiggsNLO.Integral() )
            histograms_dict_ = {
                r"$%s\to aa \to 4b$" % (productionMode):            hHToAATo4B_Stitch,
                r"%s (NLO)" % (productionMode):       hHiggsNLO,
            }


        ## Plot compare diffent samples
        cCompareSamples = plotHistograms(histograms_dict_, 'c%s_CompareSamples'%(sHistNameShort), yLable='A.U.')

        ## Higgs Pt reweights: w/ NLO
        hGGFHiggsPtReweights_NLO = makeRatioHist(
            hNume      = hHiggsNLO, 
            hDenom     = hHToAATo4B_Stitch, 
            sRatioName = '%s_Wgt_NLO'%(sHistNameShort)
            )
        
        if productionMode == 'ggH':
            ## Higgs Pt reweights: w/ Hqt
            hGGFHiggsPtReweights_Hqt_stitched = makeRatioHist(
                hNume      = hGGFH_Hqt_stitched, 
                hDenom     = hHToAATo4B_Stitch, 
                sRatioName = '%s_Wgt_Hqt'%(sHistNameShort)
                )
        

        
        
        cFitWeights = fitHistogram(
            h = hGGFHiggsPtReweights_NLO, 
            sFitFuncLocal = '[0] + [1]*x + [2]*pow(x,2) + [3]*pow(x,3)', 
            FitRangeLocal = HiggsPt_fitRange_toUse, 
            sFitFull = '[0] + [1]*x + [2]*pow(x,2) + [3]*pow(x,3)', 
            FitRangeFull = HiggsPt_fullRange_toUse, 
            sCanvasName = 'c%s_FitWeights'%(sHistNameShort))
        

        
        

        
        
        fOutHiggsPtRewgt_GGF = TFile(sFOutHiggsPtRewgt, 'update')
        fOutHiggsPtRewgt_GGF.cd()
        hHToAATo4B_Incl.Write()
        hHToAATo4B_Excl.Write()
        hHToAATo4B_Stitch.Write()

        if productionMode == 'ggH':
            hGGFHTo2B_Incl.Write()
            hGGFHTo2B_Excl.Write()
            #hGGFHTo2B_Stitch.Write()          
            
        hHiggsNLO.Write()

        hGGFHiggsPtReweights_NLO.Write()
        if productionMode == 'ggH': 
            hGGFH_Hqt_stitched_mTopInfinite.Write()
            hGGFH_Hqt_stitched.Write()
            hGGFHiggsPtReweights_Hqt_stitched.Write()
            cCompareSamplesXSNorm.Write()
            cCompareSamplesXSCumul.Write()
            #cCompareSamplesXSCumul_1.Write()
        
        cCompareSamples.Write()

        fOutHiggsPtRewgt_GGF.Close()
        print(f"\nWrote Higgs pT reweights histograms ({sHistNameShort}) into {sFOutHiggsPtRewgt}.") 




























    

    





        
    

    


