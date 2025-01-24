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


sIpFiles = OD([
    # (<file name to refer>, <file path+name>)
    ('X4bTaggerV1', '/eos/cms/store/user/ssawant/htoaa/analysis/20240627_gg0l_1/2018/analyze_htoaa_stage1.root'),   
    ('X4bTaggerV2', '/eos/cms/store/user/ssawant/htoaa/analysis/20250122_gg0l_DataMC_1/2018/analyze_htoaa_stage1.root')        
])
#sAnaVersion = list(sIpFiles.keys())[0]
sAnaVersion = "gg->H->aa->4b sample mA = 30 GeV"
print(f"sAnaVersion: {sAnaVersion}")

#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/plots'
#sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250122_gg0l_DataMC_1/2018/plots_extra/%s' % (sAnaVersion)
sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250122_gg0l_DataMC_1/2018/plots_extra'

if not os.path.exists(sOpDir):
    os.makedirs(sOpDir)


histograms_dict = OD([
    #("hLeadingPtGenBquark_pt_all", {sXLabel: 'Leading FatJet mass [GeV]', sYLabel: 'Events', sXRange: [0, 200]}),
   
    ("mA_ggHtoaato4b_mA_30", {
        sXLabel: 'mA [GeV]', sYLabel: 'Events',
        sXRange: [0, 70], #sXScale: 'log_10',
        sNRebin: 40, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("mA_v1", [
                {sIpFileNameNice: 'X4bTaggerV1', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetParticleNet_massA_Hto4b_avg_v013_gg0lIncl_SBWP95to60_central'},
            ]),
            ("mA_v2", [
                {sIpFileNameNice: 'X4bTaggerV2', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetPNet_massAa_gg0lIncl_SBWP60_Nom'},
            ]),
        ])
    }),    
    
    
])


'''
histograms_dict = OD([
    #("hLeadingPtGenBquark_pt_all", {sXLabel: 'Leading FatJet mass [GeV]', sYLabel: 'Events', sXRange: [0, 200]}),
   
    ("mH_ggHtoaato4b_mA_30", {
        sXLabel: 'mH [GeV]', sYLabel: 'Events',
        sXRange: [0, 200], #sXScale: 'log_10',
        sNRebin: 5, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("mH_v1", [
                {sIpFileNameNice: 'X4bTaggerV1', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetParticleNet_massH_Hto4b_avg_v0123_gg0lIncl_SBWP95to60_central'},
            ]),
            ("mH_v2", [
                {sIpFileNameNice: 'X4bTaggerV2', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetMassH_v2b_gg0lIncl_SBWP60_Nom'},
            ]),
            ("mass", [
                {sIpFileNameNice: 'X4bTaggerV2', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetMass_gg0lIncl_SBWP60_Nom'},
            ]),
            ("m soft-drop", [
                {sIpFileNameNice: 'X4bTaggerV2', sHistName: 'evt/ggHtoaato4b_mA_30/hLeadingFatJetMSoftDrop_gg0lIncl_SBWP60_Nom'},
            ]),
            
        ])
    }),    
    
    
])
'''