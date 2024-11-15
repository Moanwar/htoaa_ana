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


sIpFiles = OD()

sIpFiles = OD([
    # (<file name to refer>, <file path+name>)
    #('SystVariations', '/eos/cms/store/user/ssawant/htoaa/analysis/20240809_gg0l_FullSyst/2018/analyze_htoaa_stage1.root')   
    ('SystVariations', '/eos/cms/store/user/ssawant/htoaa/analysis/20240809_gg0l_FullSyst/2018/analyze_htoaa_SUSY_GluGluH_01J_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8_0_0.root')    
])
sAnaVersion = list(sIpFiles.keys())[0]
print(f"sAnaVersion: {sAnaVersion}")

sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240809_gg0l_FullSyst/2018/plots_20240823_1/%s' % (sAnaVersion)

CATAGORIES = ["gg0lIncl"]
SYSTMATICS = [
    "PU", 
    "ISR", "FSR", "QCDRenorm", "QCDFactr", "PDF", 
    "GGHPtRewgt",
    "JES", "JER",
]

histoNamesShorts_dict = {
    #'hLeadingFatJetMass_vs_massA_Hto4b_avg':                        'mass',
    #'hLeadingFatJetMSoftDrop_vs_massA_Hto4b_avg':                   'msoft',
    #'hLeadingFatJetParticleNet_massH_Hto4b_avg_vs_massA_Hto4b_avg': 'pnet',  
    'hLeadingFatJetMSoftDrop':                   'msoft',  
    'hLeadingFatJetPt': 'FatJetPt'
}


if not os.path.exists(sOpDir):
    os.makedirs(sOpDir)


histograms_dict = OD()

print(f"{CATAGORIES = }, \n{SYSTMATICS = }, \n{histograms_dict = }")

for category in CATAGORIES:
    print(f"{category = }")
    
    for syst in SYSTMATICS:
        print(f"{syst = }")    

        for histo_name, histo_name_toSave in histoNamesShorts_dict.items():
            histo_name_toUse = '%s_%s_%s' % (category, histo_name, syst)
            print(f"{histo_name = }")    
            xRange_toUse = None
            sNRebin_toUse = 1
            if histo_name == 'hLeadingFatJetMSoftDrop': 
                xRange_toUse = [50, 170]
                sNRebin_toUse = 4
            if histo_name == 'hLeadingFatJetPt': 
                xRange_toUse = [250, 800]
                sNRebin_toUse = 4


            histograms_dict[histo_name_toUse] = {
                sXLabel: '%s [GeV]'%(histo_name_toSave), sYLabel: 'A. U.',
                sXRange: xRange_toUse, #sXScale: 'log_10',
                sNRebin: sNRebin_toUse, 
                sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
                    ("Nom", [
                        {sIpFileNameNice: sAnaVersion, sHistName: 'evt/ggHtoaato4b_mA_30/%s_%s_SRWP40_Nom' % (histo_name, category)},
                    ]),
                    ("%sUp" % (syst), [
                        {sIpFileNameNice: sAnaVersion, sHistName: 'evt/ggHtoaato4b_mA_30/%s_%s_SRWP40_%sUp' % (histo_name, category, syst)},
                    ]),
                    ("%sDown" % (syst), [
                        {sIpFileNameNice: sAnaVersion, sHistName: 'evt/ggHtoaato4b_mA_30/%s_%s_SRWP40_%sDown' % (histo_name, category, syst)},
                    ]),
                ])
            }

