import os
import numpy as np
from collections import OrderedDict as OD

sXRange = "xAxisRange"; sYRange = "yAxisRange";
sXLabel = 'xAxisLabel'; sYLabel = 'yAxisLabel';
sXScale = 'xAxisScale';
sNRebin = "nRebin"
sHistosToOverlay = 'histosToOverlay'
sHistosToHadd = 'histosToHadd'
sIpFileNameNice = 'ipFileNameNice'
sHistName   = 'histogramName'
sScaleFactors = 'sScaleFactors'
sMakeRatioPlot = 'sMakeRatioPlot'

'''
sIpFiles = OD([
    # (<file name to refer>, <file path+name>)
    ('2018',        '/eos/cms/store/user/ssawant/htoaa/analysis/20250507_gg0l_DataMCPlots_NoSyst/2018/analyze_htoaa_stage1.root'),
    ('2017',        '/eos/cms/store/user/ssawant/htoaa/analysis/20250628_gg0lDatMC_JetHTplusBTagCSV/2017/analyze_htoaa_stage1.root'),
    ('2016preVFP',  '/eos/cms/store/user/ssawant/htoaa/analysis/20250630_DataMC/2016preVFP/gg0l/analyze_htoaa_stage1.root'),
    ('2016postVFP', '/eos/cms/store/user/ssawant/htoaa/analysis/20250630_DataMC/2016postVFP/gg0l/analyze_htoaa_stage1.root'),
      
])
sAnaVersion = 'CompareGGHSignal'
print(f"sAnaVersion: {sAnaVersion}")

#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/plots'
sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250630_DataMC/Run2/plots/%s' % (sAnaVersion)

histograms_dict = OD([
    #("hLeadingPtGenBquark_pt_all", {sXLabel: 'Leading FatJet mass [GeV]', sYLabel: 'Events', sXRange: [0, 200]}),
   
    ("CompareRun2GGHSignal", {
        sXLabel: 'pT(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events',
        sXRange: [180, 1000], #sXScale: 'log_10',
        sNRebin: 4, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("2018 GGH signal", [
                {sIpFileNameNice: '2018', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},
            ]),
            ("2017 GGH signal", [
                {sIpFileNameNice: '2017', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},                
            ]),     
            ("2016preVFP GGH signal", [
                {sIpFileNameNice: '2016preVFP', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},                
            ]),    
            ("2016postVFP GGH signal", [
                {sIpFileNameNice: '2016postVFP', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},
            ]),                             
        ])
    }),    
    
    ("CompareRun2GGHSignalNorm1fb", {
        sXLabel: 'pT(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events/1 fb-1',
        sXRange: [180, 1000], #sXScale: 'log_10',
        sNRebin: 4, 
        sMakeRatioPlot: True,
        sScaleFactors: {"2018 GGH signal": 1/59.8, "2017 GGH signal": 1/41.5, "2016preVFP GGH signal": 1/19.5, "2016postVFP GGH signal": 1/16.8},
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("2018 GGH signal", [
                {sIpFileNameNice: '2018', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},
            ]),
            ("2017 GGH signal", [
                {sIpFileNameNice: '2017', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},                
            ]),     
            ("2016preVFP GGH signal", [
                {sIpFileNameNice: '2016preVFP', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},                
            ]),    
            ("2016postVFP GGH signal", [
                {sIpFileNameNice: '2016postVFP', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},
            ]),                              
        ])
    }),      
])
'''


sIpFiles = OD([
    # (<file name to refer>, <file path+name>)
    ('2018',        '/eos/cms/store/user/ssawant/htoaa/analysis/20250625_gg0lDataMC_1/2018/analyze_htoaa_stage1.root'),
    ('2017',        '/eos/cms/store/user/ssawant/htoaa/analysis/20250625_gg0lDataMC_1/2017/analyze_htoaa_stage1.root'),
    ('2016preVFP',  '/eos/cms/store/user/ssawant/htoaa/analysis/20250627_gg0lDataMC/2016preVFP/analyze_htoaa_stage1.root'),
    ('2016postVFP', '/eos/cms/store/user/ssawant/htoaa/analysis/20250627_gg0lDataMC/2016postVFP/analyze_htoaa_stage1.root'),
      
])
sAnaVersion = 'CompareGGHSignal'
print(f"sAnaVersion: {sAnaVersion}")

sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250627_gg0lDataMC/Run2/plots/%s' % (sAnaVersion)

histograms_dict = OD([
    #("hLeadingPtGenBquark_pt_all", {sXLabel: 'Leading FatJet mass [GeV]', sYLabel: 'Events', sXRange: [0, 200]}),
   
    ("CompareRun2GGHSignal", {
        sXLabel: 'pT(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events',
        sXRange: [180, 1000], #sXScale: 'log_10',
        sNRebin: 4, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("2018 GGH signal", [
                {sIpFileNameNice: '2018', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},
            ]),
            ("2017 GGH signal", [
                {sIpFileNameNice: '2017', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},                
            ]),     
            ("2016preVFP GGH signal", [
                {sIpFileNameNice: '2016preVFP', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},                
            ]),    
            ("2016postVFP GGH signal", [
                {sIpFileNameNice: '2016postVFP', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},
            ]),                             
        ])
    }),    
    
    ("CompareRun2GGHSignalNorm1fb", {
        sXLabel: 'pT(Higgs candidate AK8 jet) [GeV]', sYLabel: 'Events/1 fb-1',
        sXRange: [180, 1000], #sXScale: 'log_10',
        sNRebin: 4, 
        sMakeRatioPlot: True,
        sScaleFactors: {"2018 GGH signal": 1/59.8, "2017 GGH signal": 1/41.5, "2016preVFP GGH signal": 1/19.5, "2016postVFP GGH signal": 1/16.8},
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("2018 GGH signal", [
                {sIpFileNameNice: '2018', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},
            ]),
            ("2017 GGH signal", [
                {sIpFileNameNice: '2017', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},                
            ]),     
            ("2016preVFP GGH signal", [
                {sIpFileNameNice: '2016preVFP', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},                
            ]),    
            ("2016postVFP GGH signal", [
                {sIpFileNameNice: '2016postVFP', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPt_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom'},
            ]),                              
        ])
    }),      
])

