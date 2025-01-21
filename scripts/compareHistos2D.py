from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TH1, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem, gPad, gStyle
import ctypes
import copy

'''
histoInfoDict = {
    '0': {
        'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250120_JetBtagEffi/2018/jetBtagEfficiency.root',
        'histoName': 'hJetBtagEffi_b_QCD_TT_Presel'
    },
    '1': {
        'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250120_JetBtagEffi/2018/jetBtagEfficiency.root',
        'histoName': 'hJetBtagEffi_b_QCD_TT_WP60SRplusSB'
    }
}
'''
histoInfoDict = {
    '0': {
        'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250120_JetBtagEffi/2018/jetBtagEfficiency.root',
        'histoName': 'hJetBtagEffi_b_TT_Presel'
    },
    '1': {
        'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250120_JetBtagEffi/2018/jetBtagEfficiency.root',
        'histoName': 'hJetBtagEffi_b_QCD_Presel'
    }
}

histograms = {}
for iHisto, histoInfo in histoInfoDict.items():
    f = TFile(histoInfo['ipFile'])
    histograms[iHisto] = copy.deepcopy( f.Get(histoInfo['histoName']) )
    f.Close()

print(f"{histograms = }")

hRatio = histograms['0'].Clone('hRatio')
for iBinX in range(1, hRatio.GetNbinsX()+1):
    for iBinY in range(1, hRatio.GetNbinsY()+1):
        N0  = histograms['0'].GetBinContent(iBinX, iBinY)
        eN0 = histograms['0'].GetBinError(iBinX, iBinY)
        N1  = histograms['1'].GetBinContent(iBinX, iBinY)
        eN1 = histograms['0'].GetBinError(iBinX, iBinY)
        errMax = max(eN0, eN1)
        pull = (N1 / N0 - 1) / errMax
        hRatio.SetBinContent(iBinX, iBinY, pull)
        
        
c1 = TCanvas( 'c1', 'Dynamic Filling Example', 500, 400 )
gPad.SetLogx()
gStyle.SetOptStat(0)
hRatio.Draw('colz')
c1.Update()
#c1.SaveAs()
input("Type anything")