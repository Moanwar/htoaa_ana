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
    ('2018_woHTSamplesStitch', '/home/siddhesh/Work/CMS/htoaa/analysis/20230504_v2_QCD_woHTSamplesStitch/2018/analyze_htoaa_stage1.root')    
    #('2018_HTSamplesStitchSF_step1', '/home/siddhesh/Work/CMS/htoaa/analysis/20230504_v2_QCD_wHTSamplesStitchSF_step1/2018/analyze_htoaa_stage1.root')    
])
sAnaVersion = list(sIpFiles.keys())[0]
print(f"sAnaVersion: {sAnaVersion}")

#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/plots'
sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230504_v2_QCD_woHTSamplesStitch/2018/plots/%s' % (sAnaVersion)
#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230504_v2_QCD_wHTSamplesStitchSF_step1/2018/plots/%s' % (sAnaVersion)

HTRangesForPlotting = [
    [50, 150],
    [150, 250],
    [250, 350],
    [450, 550],
    [650, 750],
    [950, 1050],
    [1450, 1550],
    [1950, 2050],
    
    #150, 350],
    #[175, 325],
]

print(f"{HTRangesForPlotting = }")


histograms_dict = OD()

for HTRange in HTRangesForPlotting:
    HTMin = HTRange[0]
    HTMax = HTRange[1]


    
    histograms_dict["LHE_HT_QCD_bEnrich_cut_2_HT%dto%d_all" % (HTMin, HTMax)] =  {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [HTMin, HTMax], #sXScale: 'log_10',
        #sNRebin: 2, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl + QCD_bEnrich selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenLHE_HT_SelQCDbEnrich_central'},
            ]), 
            ("QCD_bEnrich_all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_all_central'},
            ]),
            ("QCD_bEnrich + QCD_bEnrich selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_SelQCDbEnrich_central'},
            ]), 
        ])
    }

    histograms_dict["LHE_HT_QCD_bGen_cut_2_HT%dto%d_all" % (HTMin, HTMax)] =  {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [HTMin, HTMax], #sXScale: 'log_10',
        #sNRebin: 2, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl + QCD_bGen selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenLHE_HT_SelQCDbGen_central'},
            ]), 
            ("QCD_bGen_all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_all_central'},
            ]),
            ("QCD_bGen + QCD_bGen selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_SelQCDbGen_central'},
            ]), 
        ])
    }   

