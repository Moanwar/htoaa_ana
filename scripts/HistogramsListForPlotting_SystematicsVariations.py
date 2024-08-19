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
    ('SystVariations', '/eos/cms/store/user/ssawant/htoaa/analysis/20240809_gg0l_FullSyst/2018/analyze_htoaa_stage1.root')    
])
sAnaVersion = list(sIpFiles.keys())[0]
print(f"sAnaVersion: {sAnaVersion}")

sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240809_gg0l_FullSyst/2018/plots/%s' % (sAnaVersion)

CATAGORIES = ["gg0lIncl"]
SYSTMATICS = ["PU"]

histoNamesShorts_dict = {
    #'hLeadingFatJetMass_vs_massA_Hto4b_avg':                        'mass',
    'hLeadingFatJetMSoftDrop_vs_massA_Hto4b_avg':                   'msoft',
    #'hLeadingFatJetParticleNet_massH_Hto4b_avg_vs_massA_Hto4b_avg': 'pnet',    
}


if not os.path.exists(sOpDir):
    os.makedirs(sOpDir)


histograms_dict = OD()

print(f"{CATAGORIES = }, \n{SYSTMATICS = }, \n{}")

for category in CATAGORIES:
    print(f"{category = }")
    
    for syst in SYSTMATICS:
        print(f"{syst = }")    

        for histo_name, histo_name_toSave in histoNamesShorts_dict.items():
            print(f"{histo_name = }")    

            histograms_dict[syst] = {
                sXLabel: 'LHE HT [GeV]', sYLabel: 'A. U.',
                sXRange: [50, 2500], #sXScale: 'log_10',
                sNRebin: 10, 
                sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
                    ("Nom", [
                        {sIpFileNameNice: sAnaVersion, sHistName: 'evt/ggHtoaato4b_mA_20/%s_SRWP40_Nom' % (histo_name)},
                    ]),
                    ("%sUp" % (syst), [
                        {sIpFileNameNice: sAnaVersion, sHistName: 'evt/ggHtoaato4b_mA_20/%s_SRWP40_%sUp' % (histo_name, syst)},
                    ]),
                    ("%sDown" % (syst), [
                        {sIpFileNameNice: sAnaVersion, sHistName: 'evt/ggHtoaato4b_mA_20/%s_SRWP40_%sDown' % (histo_name, syst)},
                    ]),
                ])
            }

