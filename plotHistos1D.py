
import os
import sys
from collections import OrderedDict as OD

import ROOT as R
R.PyConfig.IgnoreCommandLineOptions = True
R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing




if __name__ == '__main__':

    sIpFile  = "~/Work/CMS/htoaa/analysis/20221020/analyze_htoaa_tmp.root"
    sOutDir  = "./plots/20221020"

    sDir0 = "evt"

    sHistoScale = 'histoScale'
    sSampleCategories = OD([
        ('QCD', {sHistoScale: 1}),
        ('SUSY_GluGluH_01J_HToAATo4B', {sHistoScale: 1e4}),
    ])

    sXaxisLabel = 'xAxisLabel'
    sXaxisRange = 'xAxisRange'
    sHistograms = OD([
        ('hLeadingFatJetPt', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetEta', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetPhi', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetMass', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetMSoftDrop', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetBtagDeepB', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetBtagDDBvLV2', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetBtagDDCvBV2', {sXaxisLabel: '', sXaxisRange: []} ),

        ('hLeadingFatJetBtagHbb', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetDeepTagMD_H4qvsQCD', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetDeepTagMD_HbbvsQCD', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetDeepTagMD_ZHbbvsQCD', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetDeepTagMD_ZHccvsQCD', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetDeepTagMD_ZbbvsQCD', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetDeepTagMD_ZvsQCD', {sXaxisLabel: '', sXaxisRange: []} ),
        
        ('hLeadingFatJetDeepTagMD_bbvsLight', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetDeepTagMD_ccvsLight', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetDeepTag_H', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetDeepTag_QCDothers', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetN2b1', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetN3b1', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetTau1', {sXaxisLabel: '', sXaxisRange: []} ),
        
        ('hLeadingFatJetTau2', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetTau3', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetTau4', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetNBHadrons', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetNCHadrons', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetNConstituents', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetParticleNetMD_QCD', {sXaxisLabel: '', sXaxisRange: []} ),
        
        ('hLeadingFatJetParticleNetMD_Xbb', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetParticleNetMD_Xcc', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetParticleNetMD_Xqq', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetParticleNet_H4qvsQCD', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetParticleNet_HbbvsQCD', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetParticleNet_HccvsQCD', {sXaxisLabel: '', sXaxisRange: []} ),
        ('hLeadingFatJetParticleNet_QCD', {sXaxisLabel: '', sXaxisRange: []} ),
        
        ('hLeadingFatJetParticleNet_mass', {sXaxisLabel: '', sXaxisRange: []} ),
        #(, {sXaxisLabel: '', sXaxisRange: []} ),
       
    ])
    sHistogram_template = 'evt/$SAMPLECAT/$HISTONAME_central'




    R.gStyle.SetPadTopMargin(0.055);
    R.gStyle.SetOptTitle(0)
    
    R.gStyle.SetPadTopMargin(0.055);
    R.gStyle.SetPadRightMargin(0.05);
    R.gStyle.SetPadBottomMargin(0.12);
    R.gStyle.SetPadLeftMargin(0.12);
    
    #// use large Times-Roman fonts
    R.gStyle.SetTextFont(132);
    R.gStyle.SetTextSize(0.05);
    
    R.gStyle.SetLabelFont(132,"x");
    R.gStyle.SetLabelFont(132,"y");
    R.gStyle.SetLabelFont(132,"z");
    
    R.gStyle.SetLabelSize(0.045,"x");
    R.gStyle.SetLabelSize(0.045,"y");
    R.gStyle.SetLabelSize(0.045,"z");
    
    R.gStyle.SetTitleSize(0.05,"x");
    R.gStyle.SetTitleSize(0.05,"y");
    R.gStyle.SetTitleSize(0.05,"z");
    
    R.gStyle.SetTitleOffset(0.85,"x");
    R.gStyle.SetTitleOffset(1.0,"y");
    
    #R.gStyle.SetNdivisions(520, "x");
    R.gStyle.SetNdivisions(505, "y");
    
    #// legend attributes
    R.gStyle.SetLegendFillColor(-1);
    R.gStyle.SetLegendBorderSize(1);
    
    R.gStyle.SetOptTitle(0);
    R.gStyle.SetOptStat(0);

    R.gStyle.SetLegendFillColor(0);


    colors   = [R.kRed, R.kBlue, R.kMagenta, R.kCyan+1, R.kOrange+2, R.kGreen+2, R.kGray+22, R.kTeal, R.kViolet, R.kSpring, R.kBlack, ]
    markers  = [20, 24, 23, 20]
    colors2  = [R.kBlack, R.kRed, R.kBlue, R.kMagenta, R.kCyan+1, R.kOrange+2, R.kGray+2, R.kGreen+2, R.kBlue-10, R.kViolet, R.kTeal, R.kSpring, R.kBlack, ]
    markers2 = [20, 22, 24, 23, 21]




    if not os.path.exists(sOutDir):
        #os.mkdir(sOutDir)
        os.makedirs(sOutDir)    

    fIp = R.TFile(sIpFile)
    if not fIp.IsOpen():
        print(f"{sIpFile} could not open")
        exit(0)
        
    for sHistoName, histoDetails in sHistograms.items():
        c1 = R.TCanvas("c1", "c1", 500, 400)
        c1.cd()
        legend = R.TLegend(0.4, 0.9, 1,1)
        
        for iSampleCat, sampleCat in enumerate(list(sSampleCategories.keys())):
            sampleCatDetails = sSampleCategories[sampleCat]
            hScale           = sampleCatDetails[sHistoScale]
            
            sHisto = sHistogram_template
            sHisto = sHisto.replace('$SAMPLECAT', sampleCat)
            sHisto = sHisto.replace('$HISTONAME', sHistoName)
            print(f"sHisto: {sHisto}")
            h = fIp.Get(sHisto)
            if h is None:
                print(f"sHisto: {sHisto} could not fetch")

            h.SetLineColor(colors[iSampleCat])
            h.SetMarkerColor(colors[iSampleCat])
            h.SetMarkerSize(1.3)

            #print(f"histoDetails: {histoDetails}")
            
            if histoDetails[sXaxisLabel] != '':
                h.GetXaxis().SetTitle(histoDetails[sXaxisLabel])
            if len(histoDetails[sXaxisRange]) > 0:
                h.GetXaxis().SetRangeUser(histoDetails[sXaxisRange][0], histoDetails[sXaxisRange][1])

            h.Scale(hScale)
                
            if iSampleCat == 0: h.Draw('P')
            else:               h.Draw('P same')

            sLegend = sampleCat
            if hScale != 1: sLegend += ' #times %.1e' % hScale
            legend.AddEntry(h, sLegend, "lep")

        legend.Draw()
        
        c1.Update()
        c1.SaveAs('%s/%s_yLin.png' % (sOutDir, sHistoName))

        c1.SetLogy()
        c1.Update()
        c1.SaveAs('%s/%s_yLog.png' % (sOutDir, sHistoName))
        
