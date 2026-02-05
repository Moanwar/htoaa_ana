import os
import sys
import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TH1, TEfficiency, TLegend, TPad
from ROOT import gROOT, gBenchmark, gRandom, gSystem, gStyle
from ROOT import kBlack, kBlue, kRed
import math
import ctypes
import copy
import cmsstyle as cmsstyle

ROOT.gROOT.SetBatch(True)

'''
Slack message from Andrew on 22/09/2025: https://baylorhep.slack.com/archives/DUJFTPGTY/p1758562760457979
for the Data-MC plots in the next round of the AN, I think it would be nice to correct the Higgs AK8 pT spectrum for the QCD background.  We should be able to use the data vs. MC for pT > 500 GeV directly, as in this region the trigger efficiency is close to 100%, so we should be mostly free of trigger sculpting effects.  As you can see, for > 500 GeV, the data/MC ratio drops in a quasi-linear way.
https://ssawant.web.cern.ch/ssawant/HToAA/DatavsMC/20250805_DataMC_AN/gg0l/Run2/gg0lIncl/?match=hLeadingFatJetPt*
Could you try:
Normalize the full MC to data
Fit the data/MC ratio (including uncertainties, i.e. sqrt[ 1/data + errMC^2 / MC^2 ]) in the [500, 1000] GeV range with:
A + B*x
A * exp(B*x)
Apply the ratios 2a and 2b to each MC bin to re-make the stack plots in the full [250, 1000] GeV range.  (This should not require event reprocessing, since there is just one weight for each histogram bin.)
Measure the overall data/MC normalization after applying these weights.
If fit 2a or 2b gives improved data/MC agreement, we can use this fit (with the norm factor from #4) to re-weight the QCD background, event by event, in all hadronic categories in the next round of data/MC plots. 

Input histogram file: /eos/cms/store/user/ssawant/htoaa/analysis/20250805_DataMC_1/Run2/gg0l/analyze_htoaa_stage1.root
'''

'''
A + B*x
TFitEditor::DoFit - using function PrevFitTMP  0x163a9bf90
 FCN=33.0076 FROM MIGRAD    STATUS=CONVERGED      41 CALLS          42 TOTAL
                     EDM=1.34914e-19    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           1.06111e+00   5.34558e-03   2.93359e-06   4.78364e-07
   2  p1          -2.60981e-04   8.09487e-06   1.26918e-06   3.30251e-04

A * exp(B*x):    
TFitEditor::DoFit - using function PrevFitTMP  0x163a9bf90
 FCN=28.0616 FROM MIGRAD    STATUS=CONVERGED      75 CALLS          76 TOTAL
                     EDM=1.6471e-11    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           8.43691e-02   7.00406e-03   2.71188e-06   5.00981e-03
   2  p1          -3.31077e-04   1.33464e-05   5.16750e-09   2.78861e+00

'''

colors_list = [1, 4, 6, 28, 46, 7, 3, 2]
Luminosity = 138
CMEnergy = 13

def createRatio(h1, h2):
    h3 = h1.Clone("h3")
    #h3.SetLineColor(kBlack)
    #h3.SetMarkerStyle(21)
    h3.SetTitle("")
    h3.SetMinimum(0.8)
    h3.SetMaximum(1.35)
    # Set up plot for markers and errors
    h3.Sumw2()
    h3.SetStats(0)
    h3.Divide(h2)
    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetTitle("ratio h1/h2 ")
    y.SetNdivisions(505)
    y.SetTitleSize(20)
    y.SetTitleFont(43)
    y.SetTitleOffset(1.55)
    y.SetLabelFont(43)
    y.SetLabelSize(15)
    # Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetTitleSize(20)
    x.SetTitleFont(43)
    x.SetTitleOffset(4.0)
    x.SetLabelFont(43)
    x.SetLabelSize(15)
    return h3


def plotHistogramsAndRatioPlots(histo_list, sSaveAs):
    c = TCanvas("c", "canvas", 800, 800)
    # Upper histogram plot is pad1
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)  # joins upper and lower plot
    pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.Draw()

    for iH in range(len(histo_list)):
        pad1.cd()
        drawOptions = '' if iH==0 else 'same'
        histo_list[iH].Draw(drawOptions)

        if iH > 0:
            hRatio = createRatio(histo_list[iH], histo_list[0])
            pad2.cd()
            hRatio.Draw("ep")

    c.SaveAs(sSaveAs)


def draw_stack_with_ratio (
        histos, colors, labels, htotal_prediction,
        label,
        xLable = '', yLable = '', ratioLable = 'Data / MC',        
        xRange = [], yRange = [], 
        ):
    """Make the plot with the histograms on the fly with a ratio panel.
        Ref: https://cms-analysis.docs.cern.ch/guidelines/plotting/examples/#stack-plot-with-cmsstyle
    """

    # set the global ROOT style to CMS style
    cmsstyle.setCMSStyle()

    # tweak the number of divisions on the X axis
    # since the number of jets is an integer, disable secondary tick marks
    cmsstyle.getCMSStyle().SetNdivisions(5, "X")

    # set the luminosity, the COM energy, the Run period to show in the canvases
    cmsstyle.SetLumi(Luminosity) #cmsstyle.SetLumi(Luminosity, run=None)
    cmsstyle.SetEnergy(CMEnergy)
    # default extra text is "Preliminary", set it to an empty string to remove it
    cmsstyle.SetExtraText('Preliminary')

    # Build the stack
    hs = cmsstyle.buildTHStack(histos[1:], colors[1:], LineColor=-1, FillColor=-1)

    # prepare the canvas, specifying a name, min and maxes on X and Y axes and axis titles.
    xRange = xRange if len(xRange) > 0 else [histos[0].GetXaxis().GetXmin(), histos[0].GetXaxis().GetXmax()]
    yRange = yRange if len(yRange) > 0 else [histos[0].GetYaxis().GetXmin(), 2*cmsstyle.cmsReturnMaxY(histos[0])]
    xLable = xLable if xLable else histos[0].GetXaxis().GetTitle()
    yLable = yLable if yLable else histos[0].GetYaxis().GetTitle()
    yLable = yLable if yLable else 'Events'
    print(f"{xLable = }, {yLable = }, {ratioLable = }, ")
    
    c = cmsstyle.cmsDiCanvas(histos[0].GetName(),
                             xRange[0], xRange[1],
                            yRange[0], yRange[1],
                            0, 2,
                            xLable, yLable, ratioLable)
    
    # prepare a legend and fill it
    # note that this plot has a slightly different size, so the legend bounds and the font size are different
    plotlegend = cmsstyle.cmsLeg(0.42,0.50,0.92,0.88, textSize=0.045, columns=2)  # The legend!
    cmsstyle.addToLegend(plotlegend, *[(histos[i], labels[i], 'lpe' if i==0 else 'f') for i in range(len(histos))])
    cmsstyle.addToLegend(plotlegend, (htotal_prediction, 'Uncertainty', 'f'))

    pad1 = c.cd(1)
    # draw the stack
    cmsstyle.cmsObjectDraw(hs,"HIST")
    pad1.SetLogy()

    # tweak the use of scientific notation for the y axis.
    cmsstyle.GetCmsCanvasHist(ROOT.gPad).GetYaxis().SetMaxDigits(3)

    # draw the uncertainty band 
    cmsstyle.cmsObjectDraw(htotal_prediction, "E2", FillStyle=3345, LineWidth=0, FillColor=12, MarkerSize=0)

    # draw the data histogram 
    cmsstyle.cmsObjectDraw(histos[0], "E", MarkerStyle=ROOT.kFullCircle)

    # add an extra label on the plot a position specified in normalized coordinates
    extraLabel = ROOT.TLatex(0.216, 0.70, label)
    extraLabel.SetNDC()
    # make sure to set the proper font
    cmsstyle.cmsObjectDraw(extraLabel, TextFont=cmsstyle.additionalInfoFont)

    # force redraw axis
    cmsstyle.UpdatePad()

    c.cd(2)
    data_ratio = histos[0].Clone()
    data_ratio.Divide(htotal_prediction)

    #prediction_ratio = htotal_prediction.Clone()
    #prediction_ratio.Divide(htotal_prediction)
    #cmsstyle.cmsObjectDraw(prediction_ratio, "E2", FillStyle=3345, LineWidth=0, FillColor=12, MarkerSize=0)
    cmsstyle.cmsObjectDraw(data_ratio, "E", MarkerStyle=ROOT.kFullCircle)

    # Saving the result!
    cmsstyle.UpdatePad(c)

    c.SaveAs(f"{label}_ratio.png".replace(" ", "_"))



def plotHistogramsAndRatio_cmsstyle (
        histos, colors, labels, htotal_prediction,
        label, canvasSaveAs, 
        xLable = '', yLable = '', ratioLable = 'Data / MC',        
        xRange = [], yRange = [], #ratioPlotRange = [0.4, 1.6]
        fitFunction_RatioPlot = '', fitRange_RatioPlot = []
        ):
    """Make the plot with the histograms on the fly with a ratio panel.
        Ref: https://cms-analysis.docs.cern.ch/guidelines/plotting/examples/#stack-plot-with-cmsstyle
    """

    # set the global ROOT style to CMS style
    cmsstyle.setCMSStyle()

    # tweak the number of divisions on the X axis
    # since the number of jets is an integer, disable secondary tick marks
    #cmsstyle.getCMSStyle().SetNdivisions(5, "X")

    # set the luminosity, the COM energy, the Run period to show in the canvases
    cmsstyle.SetLumi(Luminosity) #cmsstyle.SetLumi(Luminosity, run=None)
    cmsstyle.SetEnergy(CMEnergy)
    # default extra text is "Preliminary", set it to an empty string to remove it
    cmsstyle.SetExtraText('Preliminary')

    # Build the stack
    #hs = cmsstyle.buildTHStack(histos[1:], colors[1:], LineColor=-1, FillColor=-1)

    # prepare the canvas, specifying a name, min and maxes on X and Y axes and axis titles.
    xRange = xRange if len(xRange) > 0 else [histos[0].GetXaxis().GetXmin(), histos[0].GetXaxis().GetXmax()]
    yRange = yRange if len(yRange) > 0 else [histos[0].GetYaxis().GetXmin(), 2*cmsstyle.cmsReturnMaxY(histos[0])]
    xLable = xLable if xLable else histos[0].GetXaxis().GetTitle()
    yLable = yLable if yLable else histos[0].GetYaxis().GetTitle()
    yLable = yLable if yLable else 'Events'
    print(f"{xLable = }, {yLable = }, {ratioLable = }, ")
    
    c = cmsstyle.cmsDiCanvas(histos[0].GetName()+label,
                             xRange[0], xRange[1],
                            yRange[0], yRange[1],
                            0.5, 1.5,
                            xLable, yLable, ratioLable,
                            extraSpace=0.01, iPos=0)
    
    # prepare a legend and fill it
    # note that this plot has a slightly different size, so the legend bounds and the font size are different
    plotlegend = cmsstyle.cmsLeg(0.7,0.7,0.92,0.90, textSize=0.045, columns=1)  # The legend!
    #cmsstyle.addToLegend(plotlegend, *[(histos[i], labels[i], 'lpe' if i==0 else 'f') for i in range(len(histos))])
    #cmsstyle.addToLegend(plotlegend, (htotal_prediction, 'Uncertainty', 'f'))
    cmsstyle.addToLegend(plotlegend, *[(histos[i], labels[i], 'lpe') for i in range(len(histos))])

    c.cd(1)
    # draw the stack
    #cmsstyle.cmsObjectDraw(hs,"HIST")
    ROOT.gPad.SetLogy()

    # tweak the use of scientific notation for the y axis.
    cmsstyle.GetCmsCanvasHist(ROOT.gPad).GetYaxis().SetMaxDigits(3)
    ROOT.TGaxis.SetExponentOffset(-0.10, 0.01, "Y")

    # draw the uncertainty band 
    #cmsstyle.cmsObjectDraw(htotal_prediction, "E2", FillStyle=3345, LineWidth=0, FillColor=12, MarkerSize=0)

    # draw the data histogram 
    #cmsstyle.cmsObjectDraw(histos[0], "E", MarkerStyle=ROOT.kFullCircle)
    for i in range(len(histos)):
        #drawOption = 'EP'
        #if i > 0: drawOption += ' same'
        #cmsstyle.cmsObjectDraw(histos[i], drawOption, MarkerStyle=ROOT.kFullCircle, MarkerColor=colors_list[i], LineColor=colors_list[i])
        cmsstyle.cmsObjectDraw(histos[i], "E", MarkerStyle=20, MarkerSize=0.9, MarkerColor=colors_list[i], LineColor=colors_list[i])

    # add an extra label on the plot a position specified in normalized coordinates
    extraLabel = ROOT.TLatex(0.4, 0.80, label)
    extraLabel.SetNDC()
    # make sure to set the proper font
    cmsstyle.cmsObjectDraw(extraLabel, TextFont=cmsstyle.additionalInfoFont)

    # force redraw axis
    cmsstyle.UpdatePad()

    c.cd(2)
    prediction_ratio = htotal_prediction.Clone()
    prediction_ratio.Divide(htotal_prediction)
    cmsstyle.cmsObjectDraw(prediction_ratio, "E2", FillStyle=3345, LineWidth=0, FillColor=12, MarkerSize=0)
    data_ratio_list = {}
    for i in range(1, len(histos)):
        data_ratio_list[i] = histos[0].Clone()
        data_ratio_list[i].Divide(histos[i])
        cmsstyle.cmsObjectDraw(data_ratio_list[i], "E", MarkerStyle=20, MarkerSize=0.9, MarkerColor=colors_list[i], LineColor=colors_list[i])

    if fitFunction_RatioPlot:
        nParameters = fitFunction_RatioPlot.count('[') 
        fRatioPlot = ROOT.TF1(histos[0].GetName()+label+'_ratioPlot', fitFunction_RatioPlot, fitRange_RatioPlot[0], fitRange_RatioPlot[1])
        fRatioPlot.SetLineColor(colors_list[-1])
        data_ratio_list[1].Fit(fRatioPlot, "R0")
        cmsstyle.cmsObjectDraw(fRatioPlot, "l", LineColor=ROOT.kRed, LineWidth=2)
        print(f"{nParameters = }, {fRatioPlot.GetParameters() = }, {fRatioPlot.GetParameters()[0] = }")
        sFitResult = fitFunction_RatioPlot #'Fit: '+fitFunction_RatioPlot
        for i in range(nParameters): 
            sFitResult = sFitResult.replace('[%d]'%(i), '%.3E'%(fRatioPlot.GetParameter(i)) )
        print(f"{sFitResult = }")
        '''
        extraLabel = ROOT.TLatex(0.7, 0.85, sFitResult)
        extraLabel.SetNDC()
        # make sure to set the proper font
        cmsstyle.cmsObjectDraw(extraLabel, TextFont=cmsstyle.additionalInfoFont+6)
        '''
        plotlegend1 = cmsstyle.cmsLeg(0.4,0.75,0.92,0.90, textSize=0.09, columns=1)  # The legend!
        cmsstyle.addToLegend(plotlegend1, *[(fRatioPlot, sFitResult, 'l')])
        

    # Saving the result!
    cmsstyle.UpdatePad(c)

    c.SaveAs(canvasSaveAs)
    
    return copy.deepcopy(c)







if __name__ == "__main__":

    sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250805_DataMC_1/Run2/gg0l/analyze_htoaa_stage1.root'
    sOpDir = '/eos/cms/store/user/ssawant/htoaa/analysis/20250805_DataMC_1/Run2/gg0l/plots_QCDScale'
    sHistogramName0 = 'evt/$PROCESS/hLeadingFatJetPt_gg0lIncl_$SYSTEMATIC'
    dataProcessList = [
        'JetHT_Run2016B-ver2_HIPM', 'JetHT_Run2016C-HIPM', 'JetHT_Run2016D-HIPM', 'JetHT_Run2016E-HIPM', 'JetHT_Run2016F-HIPM',
        'BTagCSV_Run2016B-ver2_HIPM', 'BTagCSV_Run2016C-HIPM', 'BTagCSV_Run2016D-HIPM', 'BTagCSV_Run2016E-HIPM', 'BTagCSV_Run2016F-HIPM',
        'JetHT_Run2016F', 'JetHT_Run2016G', 'JetHT_Run2016H',
        'BTagCSV_Run2016F', 'BTagCSV_Run2016G', 'BTagCSV_Run2016H',
        'JetHT_Run2017B', 'JetHT_Run2017C', 'JetHT_Run2017D', 'JetHT_Run2017E', 'JetHT_Run2017F',
        'BTagCSV_Run2017B', 'BTagCSV_Run2017C', 'BTagCSV_Run2017D', 'BTagCSV_Run2017E', 'BTagCSV_Run2017F',    
        'JetHT_Run2018A', 'JetHT_Run2018B', 'JetHT_Run2018C', 'JetHT_Run2018D', 
    ]
    mcBkgProcessList = [
        "QCD_bEnr", "QCD_BGen", "QCD_Incl",
        "TT0l", "TT1l", "TT2l",
        "STop_t", "STbar_t", "ST_s_0l", "ST_s_1l", "STop_tW_Incl", "STbar_tW_Incl",
        "Zqq", "Zvv", "Zll", "Wqq", "Wlv",
        "ZZ", "WZ", "WW",
        "ttZ", "ttW", "tZq",   "ZZZ", "WZZ", "WWZ", "WWW",
        'GluGluHToBB_Pt-200ToInf', 'VBFH_dipoleRecoilOn', 'WplusHToBBQQ', 'WplusHToBBLNu', 'WminusHToBBQQ', 'WminusHToBBLNu', 'ZHToBBX', 'ttHToBB'
    ]
    nRebinX = 4
    

    ## Open i/p file
    fIn = TFile(sIpFile)
    print(f"Reading ip file {sIpFile}: {fIn}")

    if not os.path.exists(sOpDir): 
        os.makedirs( sOpDir )
    gStyle.SetOptStat(0)


    ## Read Data histogram
    hData = None
    for iProcess in dataProcessList:
        sHistoName = sHistogramName0
        sHistoName = sHistoName.replace('$PROCESS', iProcess)
        sHistoName = sHistoName.replace('$SYSTEMATIC', 'noweight')
        h_ = fIn.Get(sHistoName)
        if nRebinX != 1: h_.Rebin(nRebinX)
        if hData == None: hData = h_
        else:             hData.Add(h_)
        #print(f"Data histogram: {sHistoName}, {hData.GetEntries() = }, {h_.GetEntries() = }")
        

    ## Read MCBkg histogram
    hMCBkgTot = None
    for iProcess in mcBkgProcessList:
        sHistoName = sHistogramName0
        sHistoName = sHistoName.replace('$PROCESS', iProcess)
        sHistoName = sHistoName.replace('$SYSTEMATIC', 'Nom')
        h_ = fIn.Get(sHistoName)
        if nRebinX != 1: h_.Rebin(nRebinX)
        if hMCBkgTot == None: hMCBkgTot = h_
        else:                 hMCBkgTot.Add(h_)
        #print(f"Data histogram: {sHistoName}, {hMCBkgTot.GetEntries() = }, {h_.GetEntries() = }")


    hMCBkgTot_original = hMCBkgTot.Clone('hMCBkg_0')

    '''
    c1 = TCanvas('c1', 'c1', 600,500)
    c1.SetLogy(1)
    c1.SetGrid()
    c1.cd()
    leg1 = TLegend(0.5,0.85,0.99,0.99)
    hData.SetMarkerStyle(20)
    hData.SetMarkerSize(0.5)
    hData.SetMarkerColor(colors_list[0])
    hData.SetLineColor(colors_list[0]) 
    hMCBkgTot_original.SetMarkerStyle(20)
    hMCBkgTot_original.SetMarkerSize(0.5)
    hMCBkgTot_original.SetMarkerColor(colors_list[1])
    hMCBkgTot_original.SetLineColor(colors_list[1]) 
    hData.GetXaxis().SetTitle('pT(Higgs candidate AK8 jet) [GeV]')
    hData.GetYaxis().SetTitle('Events')
    hData.Draw()
    hMCBkgTot_original.Draw('same')
    leg1.AddEntry(hData, 'Data', 'lep')
    leg1.AddEntry(hMCBkgTot_original, 'MC', 'lep')
    leg1.Draw()
    #c1.Update()
    c1.SaveAs('%s/HiggsPt_original.png' % (sOpDir))
    '''

    '''
    plotHistogramsAndRatioPlots(
        histo_list = [hData, hMCBkgTot_original], 
        sSaveAs = '%s/HiggsPt_original.png' % (sOpDir)
    )
    '''
    '''
    draw_stack_with_ratio(
        histos = [hData, hMCBkgTot_original],
        colors = colors_list[:2],
        labels = ['Data',  "MC"],
        htotal_prediction = hMCBkgTot_original, 
        label = "HiggsPt_original",

    )
    '''

    histos = [hData, hMCBkgTot_original]
    cDataMC_original = plotHistogramsAndRatio_cmsstyle(
        histos = histos,
        colors = colors_list,
        labels = ['Data',  "MC"],
        htotal_prediction = histos[1], 
        label = "Original",
        canvasSaveAs = '%s/HiggsPt_original.png' % (sOpDir),
        xLable = 'pT(Higgs candidate AK8 jet) [GeV]', yLable = 'Events', ratioLable = 'Data / MC',
        xRange = [200, 1000], yRange=[1e3, 8e6]
    )

    ## Scale MC to Data
    kDataMCRatio_original = hData.Integral() / hMCBkgTot_original.Integral()
    print(f"{hData.Integral() = }, {hMCBkgTot_original.Integral() = }, {kDataMCRatio_original = }")
    hMCBkgTot_Normalized = hMCBkgTot_original.Clone('hMCBkg_1')
    hMCBkgTot_Normalized.Scale(kDataMCRatio_original)
    print(f"{hData.Integral() = }, {hMCBkgTot_Normalized.Integral() = }, {hData.Integral() / hMCBkgTot_Normalized.Integral() = }")

    
    histos = [hData, hMCBkgTot_Normalized]
    print(f"Normalized MC to Data. Fit ratio with linear function", flush=True)
    cDataMC_normalized_linearFit = plotHistogramsAndRatio_cmsstyle(
        histos = histos,
        colors = colors_list,
        labels = ['Data',  "MC"],
        htotal_prediction = histos[1], 
        label = "Normalized",
        canvasSaveAs = '%s/HiggsPt_Normalized_RatioPlotFitLinear.png' % (sOpDir),
        xLable = 'pT(Higgs candidate AK8 jet) [GeV]', yLable = 'Events', ratioLable = 'Data / MC',
        xRange = [200, 1000], yRange=[1e3, 8e6],
        fitFunction_RatioPlot = '[0] + [1]*x', fitRange_RatioPlot = [500, 1000]
    )
    '''
    Fit result of ratio plot fitted to '[0] + [1]*x' in [500, 1000] range
 FCN=33.0076 FROM MIGRAD    STATUS=CONVERGED      41 CALLS          42 TOTAL
                     EDM=1.22089e-19    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           1.06111e+00   5.34558e-03   2.93359e-06   4.72308e-07
   2  p1          -2.60981e-04   8.09487e-06   1.26918e-06   3.16679e-04
    '''

    print(f"Normalized MC to Data. Fit ratio with exponential function", flush=True)    
    cDataMC_normalized_exponentialFit = plotHistogramsAndRatio_cmsstyle(
        histos = histos,
        colors = colors_list,
        labels = ['Data',  "MC"],
        htotal_prediction = histos[1], 
        label = "Normalized",
        canvasSaveAs = '%s/HiggsPt_Normalized_RatioPlotExponential.png' % (sOpDir),
        xLable = 'pT(Higgs candidate AK8 jet) [GeV]', yLable = 'Events', ratioLable = 'Data / MC',
        xRange = [200, 1000], yRange=[1e3, 8e6],
        fitFunction_RatioPlot = '[0] + exp([1]*x)', fitRange_RatioPlot = [500, 1000]
    )
    '''
    Fit result of ratio plot fitted to '[0] + [1]*x' in [500, 1000] range
 FCN=28.0616 FROM MIGRAD    STATUS=CONVERGED      75 CALLS          76 TOTAL
                     EDM=1.6472e-11    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           8.43691e-02   7.00406e-03   2.71188e-06   5.00983e-03
   2  p1          -3.31077e-04   1.33464e-05   5.16750e-09   2.78865e+00
    '''

    ## Reweight MC using linear
    hMCBkgTot_Rewgted_LinearFit = hMCBkgTot_original.Clone('%s_rewgted_LinearFit'%(hMCBkgTot_original.GetName())) 
    xRange_Reweighting = [250, 1000]
    print(f"{hMCBkgTot_original.FindBin(xRange_Reweighting[0]) = }, {hMCBkgTot_original.FindBin(xRange_Reweighting[1]) = }")
    for iBin in range(
        hMCBkgTot_original.FindBin(xRange_Reweighting[0]),
        hMCBkgTot_original.FindBin(xRange_Reweighting[1])
        ):
        x = hMCBkgTot_original.GetBinCenter(iBin)
        wgt = 1.061 - 2.610E-04*x
        hMCBkgTot_Rewgted_LinearFit.SetBinContent(iBin, hMCBkgTot_original.GetBinContent(iBin) * wgt)
        hMCBkgTot_Rewgted_LinearFit.SetBinError(iBin,   hMCBkgTot_original.GetBinError(iBin) * wgt)
        #print(f"{iBin = }, {x = }, {wgt = }, {hMCBkgTot_original.GetBinContent(iBin) = }, {hMCBkgTot_Rewgted_LinearFit.GetBinContent(iBin) = } ")
    print(f"{hData.Integral() = }, {hMCBkgTot_Rewgted_LinearFit.Integral() = }, {hData.Integral() / hMCBkgTot_Rewgted_LinearFit.Integral() = }, {hMCBkgTot_original.Integral() / hMCBkgTot_Rewgted_LinearFit.Integral() = }")
    histos = [hData, hMCBkgTot_original, hMCBkgTot_Rewgted_LinearFit]
    cDataMC_Rewgted_LinearFit = plotHistogramsAndRatio_cmsstyle(
        histos = histos,
        colors = colors_list,
        labels = ['Data',  "MC", "MC reweighted"],
        htotal_prediction = histos[1], 
        label = "Rewgted_LinearFit",
        canvasSaveAs = '%s/HiggsPt_Rewgted_LinearFit.png' % (sOpDir),
        xLable = 'pT(Higgs candidate AK8 jet) [GeV]', yLable = 'Events', ratioLable = 'Data / MC',
        xRange = [200, 1000], yRange=[1e3, 8e6]
    )        

        
    ## Reweight MC using exponential fit
    hMCBkgTot_Rewgted_ExponentialFit = hMCBkgTot_original.Clone('%s_rewgted_ExponentialFit'%(hMCBkgTot_original.GetName())) 
    for iBin in range(
        hMCBkgTot_original.FindBin(xRange_Reweighting[0]),
        hMCBkgTot_original.FindBin(xRange_Reweighting[1])
        ):
        x = hMCBkgTot_original.GetBinCenter(iBin)
        wgt = 8.437E-02 + math.exp(-3.311E-04*x)
        hMCBkgTot_Rewgted_ExponentialFit.SetBinContent(iBin, hMCBkgTot_original.GetBinContent(iBin) * wgt)
        hMCBkgTot_Rewgted_ExponentialFit.SetBinError(iBin,   hMCBkgTot_original.GetBinError(iBin) * wgt)
        #print(f"{iBin = }, {x = }, {wgt = }, {hMCBkgTot_original.GetBinContent(iBin) = }, {hMCBkgTot_Rewgted_ExponentialFit.GetBinContent(iBin) = } ")
    print(f"{hData.Integral() = }, {hMCBkgTot_Rewgted_ExponentialFit.Integral() = }, {hData.Integral() / hMCBkgTot_Rewgted_ExponentialFit.Integral() = }, {hMCBkgTot_original.Integral() / hMCBkgTot_Rewgted_ExponentialFit.Integral() = }, ")
    histos = [hData, hMCBkgTot_original, hMCBkgTot_Rewgted_ExponentialFit]
    cDataMC_Rewgted_ExponentialFit = plotHistogramsAndRatio_cmsstyle(
        histos = histos,
        colors = colors_list,
        labels = ['Data',  "MC", "MC reweighted"],
        htotal_prediction = histos[1], 
        label = "Rewgted_ExponentialFit",
        canvasSaveAs = '%s/HiggsPt_Rewgted_ExponentialFit.png' % (sOpDir),
        xLable = 'pT(Higgs candidate AK8 jet) [GeV]', yLable = 'Events', ratioLable = 'Data / MC',
        xRange = [200, 1000], yRange=[1e3, 8e6]
    )  

    # Andrew: Data/MC ratio plot with the pT reweighting using exponential fit is slightly more flat. Use exponential fit for pT reweighting
    # Total number of events should not change after QCD pT reweighting
    # Total number of events before and after QCD pT reweighting out-of-box: hMCBkgTot_original.Integral() / hMCBkgTot_Rewgted_ExponentialFit.Integral() = 1.0326494383386389. 

    # Renormalize MC_QCDPtRewgted to keep total yield unchanged: Way 1: Andrew: Add 3.265% into A of 'A + exp(Bx)' fit
    hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way1 = hMCBkgTot_original.Clone('%s_rewgted_ExponentialFit'%(hMCBkgTot_original.GetName())) 
    for iBin in range(
        hMCBkgTot_original.FindBin(xRange_Reweighting[0]),
        hMCBkgTot_original.FindBin(xRange_Reweighting[1])
        ):
        x = hMCBkgTot_original.GetBinCenter(iBin)
        wgt = 0.11697 + math.exp(-3.311E-04*x) # (8.437E-02 + 0.0326) + math.exp(-3.311E-04*x)
        hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way1.SetBinContent(iBin, hMCBkgTot_original.GetBinContent(iBin) * wgt)
        hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way1.SetBinError(iBin,   hMCBkgTot_original.GetBinError(iBin) * wgt)
        #print(f"{iBin = }, {x = }, {wgt = }, {hMCBkgTot_original.GetBinContent(iBin) = }, {hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way1.GetBinContent(iBin) = } ")
    print(f"{hData.Integral() = }, {hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way1.Integral() = }, {hData.Integral() / hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way1.Integral() = }, {hMCBkgTot_original.Integral() / hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way1.Integral() = }, ")
    histos = [hData, hMCBkgTot_original, hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way1]
    cDataMC_Rewgted_ExponentialFit = plotHistogramsAndRatio_cmsstyle(
        histos = histos,
        colors = colors_list,
        labels = ['Data',  "MC", "MC reweighted"],
        htotal_prediction = histos[1], 
        label = "Rewgted_ExponentialFit",
        canvasSaveAs = '%s/HiggsPt_Rewgted_ExponentialFit_Renormalized_way1AddResidue.png' % (sOpDir),
        xLable = 'pT(Higgs candidate AK8 jet) [GeV]', yLable = 'Events', ratioLable = 'Data / MC',
        xRange = [200, 1000], yRange=[1e3, 8e6]
    )    

    # Renormalize MC_QCDPtRewgted to keep total yield unchanged: Way 2: Siddhesh: Scale 'A + exp(Bx)' fit by 1.03265    
    hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way2 = hMCBkgTot_original.Clone('%s_rewgted_ExponentialFit'%(hMCBkgTot_original.GetName())) 
    for iBin in range(
        hMCBkgTot_original.FindBin(xRange_Reweighting[0]),
        hMCBkgTot_original.FindBin(xRange_Reweighting[1])
        ):
        x = hMCBkgTot_original.GetBinCenter(iBin)
        wgt = (8.437E-02 + math.exp(-3.311E-04*x)) * 1.03265
        hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way2.SetBinContent(iBin, hMCBkgTot_original.GetBinContent(iBin) * wgt)
        hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way2.SetBinError(iBin,   hMCBkgTot_original.GetBinError(iBin) * wgt)
        #print(f"{iBin = }, {x = }, {wgt = }, {hMCBkgTot_original.GetBinContent(iBin) = }, {hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way2.GetBinContent(iBin) = } ")
    print(f"{hData.Integral() = }, {hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way2.Integral() = }, {hData.Integral() / hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way2.Integral() = }, {hMCBkgTot_original.Integral() / hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way2.Integral() = }, ")
    histos = [hData, hMCBkgTot_original, hMCBkgTot_Rewgted_ExponentialFit_Renormalized_way2]
    cDataMC_Rewgted_ExponentialFit = plotHistogramsAndRatio_cmsstyle(
        histos = histos,
        colors = colors_list,
        labels = ['Data',  "MC", "MC reweighted"],
        htotal_prediction = histos[1], 
        label = "Rewgted_ExponentialFit",
        canvasSaveAs = '%s/HiggsPt_Rewgted_ExponentialFit_Renormalized_way2ScaleResidue.png' % (sOpDir),
        xLable = 'pT(Higgs candidate AK8 jet) [GeV]', yLable = 'Events', ratioLable = 'Data / MC',
        xRange = [200, 1000], yRange=[1e3, 8e6]
    ) 



    hRatio_1 = hData.Clone('%s_ratio_1'%(hData.GetName()))
    hRatio_1.Divide(hMCBkgTot_Normalized)
    

    
    sOpFile = '%s/QCDScale_histograms.root'%(sOpDir)
    fOp = TFile(sOpFile, 'recreate')
    fOp.cd();
    #cDataMC_0.Write()
    #cDataMC_1.Write()
    hRatio_1.Write()
    cDataMC_original.Write()
    cDataMC_normalized_linearFit.Write()
    cDataMC_normalized_exponentialFit.Write()

    fOp.Close()
    print(f"Wrote {sOpFile}")
    
    







