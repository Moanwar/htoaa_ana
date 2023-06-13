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
    ('LHE_HT_2018', '../data/correction/mc/HTSamplesStitch/LHE_HT_2018.root')    
])
sAnaVersion = list(sIpFiles.keys())[0]
print(f"sAnaVersion: {sAnaVersion}")

#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/plots'
sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230504_v2_QCD_woHTSamplesStitch/2018/plots/%s' % (sAnaVersion)

histograms_dict = OD([
    #("hLeadingPtGenBquark_pt_all", {sXLabel: 'Leading FatJet mass [GeV]', sYLabel: 'Events', sXRange: [0, 200]}),
   
    ("LHE_HT_QCDIncl_vs_QCDInclPSWeight", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [50, 2500], #sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl_PSWeight", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl_PSWeight/hGenLHE_HT_all_central'},
            ]),
            ("QCD_Incl", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenLHE_HT_all_central'},
            ]),
        ])
    }),    
    
    
])

