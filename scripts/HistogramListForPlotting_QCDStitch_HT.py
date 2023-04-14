import os
import numpy as np
from collections import OrderedDict as OD

sXRange = "xAxisRange"; sYRange = "yAxisRange";
sXLabel = 'xAxisLabel'; sYLabel = 'yAxisLabel';
sNRebin = "nRebin"
sHistosToOverlay = 'histosToOverlay'
sHistosToHadd = 'histosToHadd'
sIpFileNameNice = 'ipFileNameNice'
sHistName   = 'histogramName'


sIpFiles = OD([
    # (<file name to refer>, <file path+name>)
    #(sAnaVersion, '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/analyze_hadded_QCD_HT100to200.root')
    #(sAnaVersion, '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD/analyze_hadded_QCD.root')
    #(sAnaVersion, '')
    ('QCD_fullHT', '/home/siddhesh/Work/CMS/htoaa/analysis/20230407_QCD/analyze_hadded_QCD.root')
    #('QCD_HT100to200', '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD/analyze_hadded_QCD_HT100to200.root')
    #('QCD_HT200to300', '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD/analyze_hadded_QCD_HT200to300.root')
    #('QCD_HT300to500', '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD/analyze_hadded_QCD_HT300to500.root')
    #('QCD_HT500to700', '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD/analyze_hadded_QCD_HT500to700.root')
    #('QCD_HT700to1000', '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD/analyze_hadded_QCD_HT700to1000.root')
    #('QCD_HT1000to1500', '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD/analyze_hadded_QCD_HT1000to1500.root')
    #('QCD_HT1500to2000', '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD/analyze_hadded_QCD_HT1500to2000.root')
    #('QCD_HT2000toInf', '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD/analyze_hadded_QCD_HT2000toInf.root')
    
])
sAnaVersion = list(sIpFiles.keys())[0]
print(f"HistogramListForPlotting_QCDStitch_HT:: sAnaVersion: {sAnaVersion}")

#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/plots'
sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230407_QCD/plots_HT/%s' % (sAnaVersion)

HTRangesForPlotting = [
    #[50, 150],
    #[150, 250],
    #[250, 350],
    #[450, 550],
    #[650, 750],
    #[950, 1050],
    #[1450, 1550],
    #[1950, 2050],
    
    #150, 350],
    [175, 325],
]

print(f"{HTRangesForPlotting = }")

histograms_dict = OD()

for HTRange in HTRangesForPlotting:
    HTMin = HTRange[0]
    HTMax = HTRange[1]
    
    histograms_dict["hGenLHE_HT%dto%d_all" % (HTMin, HTMax)] = {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [HTMin, HTMax],
        #sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenLHE_HT_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenLHE_HT_all_central'}
            ]),
        ])
    }
    
    histograms_dict["hGenLHE_HT%dto%d_QCDStitch" % (HTMin, HTMax)] = {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [HTMin, HTMax],
        #sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenLHE_HT_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]),
        ])
    }

    histograms_dict["hGenLHE_HT%dto%d_QCD_bEnrich" % (HTMin, HTMax)] = {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [HTMin, HTMax],
        #sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_all_central'},
            ]),
        ])
    }
    

    histograms_dict["hGenLHE_HT%dto%d_QCD_bGen" % (HTMin, HTMax)] = {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [HTMin, HTMax],
        #sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_all_central'}
            ]),
        ])
    }
    
    histograms_dict["hGenLHE_HT%dto%d_QCD_Incl" % (HTMin, HTMax)] = {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [HTMin, HTMax],
        #sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenLHE_HT_all_central'}
            ]),
        ])
    }
    
