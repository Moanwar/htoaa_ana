import os, sys
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TH1, TEfficiency, TLegend
from ROOT import gROOT, gBenchmark, gRandom, gSystem, gPad, gStyle
import ctypes
import copy

histoInfoDict0 = {
    'b': { # b jetFlavour
        'b_TT': {
            'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/jetBtagEfficiency.root',
            'histoName': 'hJetBtagEffi_b_TT_Presel'
        },
        'b_QCD': {
            'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/jetBtagEfficiency.root',
            'histoName': 'hJetBtagEffi_b_QCD_Presel'
        },
        'b_QCD_TT': {
            'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/jetBtagEfficiency.root',
            'histoName': 'hJetBtagEffi_b_QCD_TT_Presel'
        }   
    },

    'c': {# c jetFlavour
        'c_TT': {
            'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/jetBtagEfficiency.root',
            'histoName': 'hJetBtagEffi_c_TT_Presel'
        },
        'c_QCD': {
            'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/jetBtagEfficiency.root',
            'histoName': 'hJetBtagEffi_c_QCD_Presel'
        },
        'c_QCD_TT': {
            'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/jetBtagEfficiency.root',
            'histoName': 'hJetBtagEffi_c_QCD_TT_Presel'
        }   
    },

    'l': {# light jetFlavour
        'l_TT': {
            'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/jetBtagEfficiency.root',
            'histoName': 'hJetBtagEffi_l_TT_Presel'
        },
        'l_QCD': {
            'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/jetBtagEfficiency.root',
            'histoName': 'hJetBtagEffi_l_QCD_Presel'
        },
        'l_QCD_TT': {
            'ipFile': '/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/jetBtagEfficiency.root',
            'histoName': 'hJetBtagEffi_l_QCD_TT_Presel'
        }   
    },
}


sOutDir='/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/plots'
if not os.path.exists(sOutDir):
    os.makedirs(sOutDir)

for jetFlavour, histoInfoDict in histoInfoDict0.items():
    #if jetFlavour != 'b': continue

    histograms = {}
    for iHisto, histoInfo in histoInfoDict.items():
        print(f"{iHisto = }, {histoInfo['ipFile'] = }, {histoInfo['histoName'] = }")
        f = TFile(histoInfo['ipFile'])
        histograms[iHisto] = copy.deepcopy( f.Get(histoInfo['histoName']) )
        f.Close()

    key0 = list(histograms.keys())[0]
    print(f"{histograms = }, {histograms[key0].GetYaxis().GetNbins() = }")

    for ibinY in range(1,histograms[key0].GetYaxis().GetNbins()+1):
    #for ibinY in range(1,2):
        h0 = histograms[key0]
        ibinY_min = h0.GetYaxis().GetBinLowEdge(ibinY)
        ibinY_max = h0.GetYaxis().GetBinUpEdge(ibinY)
        print(f"{ibinY = }, {ibinY_min = }, {ibinY_max = }")

        c1 = TCanvas( 'c1', 'Dynamic Filling Example', 500, 400 )
        gPad.SetLogx()
        gStyle.SetOptStat(0)

        leg = TLegend(0.8,0.8,1,1)
        leg.SetHeader('%g < |eta| < %g'%(ibinY_min,ibinY_max), 'C')
        
        i=0
        for iHisto, h in histograms.items():
            hProjX = h.ProjectionX('Px_%s_%d'%(iHisto,ibinY), ibinY,ibinY)

            c1.cd()
            hProjX.SetLineColor(i+1)
            hProjX.SetMarkerColor(i+1)
            #hProjX.GetXaxis().SetRangeUser(10,3000)
            hProjX.GetYaxis().SetRangeUser(hProjX.GetMinimum()*0.8,hProjX.GetMaximum()*1.3)
            hProjX.GetYaxis().SetTitle('Efficiency')
            hProjX.Draw("same")

            leg.AddEntry(hProjX, iHisto, 'lep')

            i+=1

            '''
            c2 = TCanvas( 'c2', 'Dynamic Filling Example', 1000, 800 )
            gPad.SetLogx()
            gStyle.SetOptStat(0)  
            h.Draw('colz text') 
            c2.Update()         

            print(f"{jetFlavour = }, {ibinY = }, {ibinY_min = }, {ibinY_max = }, {iHisto = }. Print hProjX::")
            for ibinX in range(1,hProjX.GetXaxis().GetNbins()+1):
                print(f"\t {ibinX = }, {hProjX.GetBinContent(ibinX) = }")
            input('Type anything')
            '''




        leg.Draw()
        #c1.Update()
        c1.SaveAs('%s/hJetBtagEffi_%s_Eta%gto%g.png' % (sOutDir, jetFlavour, ibinY_min,ibinY_max))
        #input('Type anything')

print(f"\n {sOutDir = }")