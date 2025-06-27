from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TH1, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import ctypes
 
# 2018
#sIpFile = "/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/analyze_htoaa_stage1.root"
#sOpFile = "/eos/cms/store/user/ssawant/htoaa/analysis/20250402_BtagEffi_QCDTT/2018/jetBtagEfficiency.root"
# 2017
#sIpFile = "/eos/cms/store/user/ssawant/htoaa/analysis/20250612_BtagEffi_QCDTT/2017/analyze_htoaa_stage1.root"
#sOpFile = "/eos/cms/store/user/ssawant/htoaa/analysis/20250612_BtagEffi_QCDTT/2017/jetBtagEfficiency.root"
# 2016preVFP
#sIpFile = "/eos/cms/store/user/ssawant/htoaa/analysis/20250625_BtagEffi/2016preVFP/analyze_htoaa_stage1.root"
#sOpFile = "/eos/cms/store/user/ssawant/htoaa/analysis/20250625_BtagEffi/2016preVFP/jetBtagEfficiency.root"
# 2016postVFP
sIpFile = "/eos/cms/store/user/ssawant/htoaa/analysis/20250625_BtagEffi/2016postVFP/analyze_htoaa_stage1.root"
sOpFile = "/eos/cms/store/user/ssawant/htoaa/analysis/20250625_BtagEffi/2016postVFP/jetBtagEfficiency.root"

samples = {
    'TT':     ['TT0l','TT1l'],
    'QCD':    ['QCD_bEnr', 'QCD_BGen', 'QCD_Incl'],
    'QCD_TT': ['QCD_bEnr', 'QCD_BGen', 'QCD_Incl','TT0l','TT1l'],
    
}
WPs = {
    'Presel':       ['Presel'], 
    'SRWP60':       ['gg0lIncl_SRWP60'], 
    'SBWP60':       ['gg0lIncl_SBWP60'],
    'WP60SRplusSB': ['gg0lIncl_SRWP60', 'gg0lIncl_SBWP60']
}
JetFlavours = ['b', 'c', 'l'] #['b', 'c', 'l']
EffiCompoments = ['deom', 'nume'] # ['denominator', 'numerator']
HistogramName = 'hJetBtagEffi'


fIn  = TFile(sIpFile)
fOut = TFile(sOpFile, 'recreate')


hDummy = TH1F('hDummy','',1,0,1)
hDummy.SetDefaultSumw2()


for jetFlavour in JetFlavours:
    print(f"{jetFlavour = }")

    for sampleName, sampleList in samples.items():
        print(f"{sampleName = }, {sampleList = }")

        for WPName, WPList in WPs.items():
            print(f"{WPName = }")

            hEffiComponents = {} # dict of hDenominator and hNumerator
            for effiComponent in EffiCompoments:

                for sample in sampleList:
                    for WP in WPList:
                        sHistoNameFull = 'evt/%s/%s_%s_%s_%s_Nom' % (sample,  HistogramName, jetFlavour, effiComponent, WP, ) # hJetBtagEffi_b_deom_Presel_Nom
                        h_ = fIn.Get(sHistoNameFull)
                        #print(f"h_ ({type(h_)}){h_}")
                        if effiComponent not in hEffiComponents:
                            hEffiComponents[effiComponent] = h_
                        else: 
                            hEffiComponents[effiComponent].Add(h_)

                        eIntegral = ctypes.c_double(0.0)
                        integral = hEffiComponents[effiComponent].IntegralAndError(1, h_.GetNbinsX(), 1, h_.GetNbinsY(), eIntegral)

                        print(f"h_ ({type(h_)}){h_}, {h_.Integral() = }, {integral = }, {eIntegral = }")

            sEffiName = "%s_%s_%s_%s" % (HistogramName, jetFlavour, sampleName, WPName)
            hDenominator = hEffiComponents[EffiCompoments[0]]
            hNumerator   = hEffiComponents[EffiCompoments[1]]
            #hEfficiency = TEfficiency(hNumerator, hDenominator)
            hEfficiency = hNumerator.Clone(sEffiName); hEfficiency.Divide(hNumerator, hDenominator);
            hEfficiency.SetName(sEffiName)
            hEfficiency.SetTitle(sEffiName)
            hDenominator.SetName(sEffiName+'_denominator')
            hDenominator.SetTitle(sEffiName+'_denominator')
            hNumerator.SetName(sEffiName+'_numerator')
            hNumerator.SetTitle(sEffiName+'_numerator')
            


            fOut.cd()
            hEfficiency.Write()
            hDenominator.Write()
            hNumerator.Write()

            

                    


fIn.Close()
fOut.Close()
                


