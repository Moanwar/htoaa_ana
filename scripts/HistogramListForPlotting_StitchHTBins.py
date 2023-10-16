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



sAnaVersion = "DYJets" # 'QCD', "ZJetsToQQ", "WJetsToQQ", "WJetsToLNu", "DYJets"
sIpFiles = OD([
    # (<file name to refer>, <file path+name>)
    #(sAnaVersion, '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/analyze_hadded_QCD_HT100to200.root')
    #(sAnaVersion, '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD/analyze_hadded_QCD.root')
    #(sAnaVersion, '')
    #('QCD_fullHT', '/home/siddhesh/Work/CMS/htoaa/analysis/20230407_QCD/analyze_hadded_QCD.root')
    #(sAnaVersion, '/home/siddhesh/Work/CMS/htoaa/analysis/20230519_StitchOverlapPhSpRemoval/2018/analyze_htoaa_stage1.root')
    #(sAnaVersion, '/Users/siddhesh/Work/CMS/htoaa/analysis/20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt/2018/analyze_htoaa_stage1.root')
    (sAnaVersion, '/eos/cms/store/user/ssawant/htoaa/analysis/20231011_DY/2018/analyze_htoaa_stage1.root')
])
#sAnaVersion = list(sIpFiles.keys())[0]
print(f"HistogramListForPlotting_QCDStitch_HT:: sAnaVersion: {sAnaVersion}")

#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/plots'
sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20231011_DY/2018/plots/%s' % (sAnaVersion)

HTRangesForPlotting_QCD = [
    [50, 150],
    [150, 250],
    [250, 350],
    [450, 550],
    [650, 750],
    [950, 1050],
    [1450, 1550],
    [1950, 2050],
    [2900, 3100],
    
    #150, 350],
    #[175, 325],
]

HTRangesForPlotting_ZJetsToQQ = [
    [350, 450],
    [550, 650],
    [750, 850],    
]

HTRangesForPlotting_WJetsToQQ = [
    [350, 450],
    [550, 650],
    [750, 850],   
]

HTRangesForPlotting_WJetsToLNu = [
    [70, 150],   
    [150, 250],   
    [350, 450],   
    [550, 650],   
    [750, 850],   
    [1150, 1250],   
    [2450, 2550],   
]

HTRangesForPlotting_DYJets = [
    [85, 115],   
    [150, 250],   
    [350, 450],   
    [550, 650],   
    [750, 850],   
    [1150, 1250],   
    [2450, 2550],   
]


if 'QCD' in sAnaVersion:
    HTRangesForPlotting = HTRangesForPlotting_QCD
elif 'ZJetsToQQ' in sAnaVersion:
    HTRangesForPlotting = HTRangesForPlotting_ZJetsToQQ
elif 'WJetsToQQ' in sAnaVersion:
    HTRangesForPlotting = HTRangesForPlotting_WJetsToQQ
elif 'WJetsToLNu' in sAnaVersion:
    HTRangesForPlotting = HTRangesForPlotting_WJetsToLNu
elif 'DYJets' in sAnaVersion:
    HTRangesForPlotting = HTRangesForPlotting_DYJets


print(f"{HTRangesForPlotting = }")

histograms_dict = OD()

for HTRange in HTRangesForPlotting:
    HTMin = HTRange[0]
    HTMax = HTRange[1]
    
    if 'QCD' in sAnaVersion:
        '''
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
        '''

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
                ("QCD_Incl-madgraph_PSWeights all", [
                    {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_all_central'}
                ]),
                ("QCD_Incl-madgraphMLM all", [
                    {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl/hGenLHE_HT_all_central'}
                ]),
            ])
        }
        '''
        histograms_dict["hGenLHE_HT%dto%d_QCD_Incl_PSWeight" % (HTMin, HTMax)] = {
            sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
            sXRange: [HTMin, HTMax],
            #sNRebin: 2,
            sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
                ("QCD_Incl-madgraph_PSWeights all", [
                    {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_all_central'}
                ]),
            ])
        }
        '''
        
    elif 'ZJetsToQQ' in sAnaVersion:
        histograms_dict["hGenLHE_HT%dto%d_ZJetsToQQ_HT" % (HTMin, HTMax)] = {
            sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
            sXRange: [HTMin, HTMax],
            #sNRebin: 2,
            sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
                ("ZJetsToQQ_HT all", [
                    {sIpFileNameNice: sAnaVersion, sHistName: 'evt/ZJetsToQQ_HT/hGenLHE_HT_all_central'}
                ]),
            ])
        }

    elif 'WJetsToQQ' in sAnaVersion:
        histograms_dict["hGenLHE_HT%dto%d_WJetsToQQ_HT" % (HTMin, HTMax)] = {
            sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
            sXRange: [HTMin, HTMax],
            #sNRebin: 2,
            sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
                ("WJetsToQQ_HT all", [
                    {sIpFileNameNice: sAnaVersion, sHistName: 'evt/WJetsToQQ_HT/hGenLHE_HT_all_central'}
                ]),
            ])
        }
        

    elif 'WJetsToLNu' in sAnaVersion:
        histograms_dict["hGenLHE_HT%dto%d_WJetsToLNu_HT" % (HTMin, HTMax)] = {
            sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
            sXRange: [HTMin, HTMax],
            #sNRebin: 2,
            sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
                ("WJetsToLNu_HT all", [
                    {sIpFileNameNice: sAnaVersion, sHistName: 'evt/WJetsToLNu_HT/hGenLHE_HT_all_central'}
                ]),
            ])
        }

    elif 'DYJets' in sAnaVersion:
        histograms_dict["hGenLHE_HT%dto%d_DYJets_HT" % (HTMin, HTMax)] = {
            sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
            sXRange: [HTMin, HTMax],
            #sNRebin: 2,
            sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
                ("DYJets_HT all", [
                    {sIpFileNameNice: sAnaVersion, sHistName: 'evt/DYJets_HT_LO/hGenLHE_HT_all_central'}
                ]),
            ])
        }
