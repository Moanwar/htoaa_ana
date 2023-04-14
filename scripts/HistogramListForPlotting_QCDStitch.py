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
print(f"sAnaVersion: {sAnaVersion}")

#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/plots'
sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230407_QCD/plots_test/%s' % (sAnaVersion)

histograms_dict = OD([
    #("hLeadingPtGenBquark_pt_all", {sXLabel: 'Leading FatJet mass [GeV]', sYLabel: 'Events', sXRange: [0, 200]}),
   
    ("hLeadingPtGenBquark_pt", {
        sXLabel: 'pT(pT leading b-quark) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquark_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquark_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquark_pt_all_central'}
            ]),
        ])
    }),    
    
    
    ("hLeadingPtGenBquark_eta", {
        sXLabel: 'eta(pT leading b-quark)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquark_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquark_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquark_eta_all_central'}
            ]),
        ])
    }),

    
   
    ("hLeadingPtGenBquarkHardSctred_pt", {
        sXLabel: 'pT(pT leading b-quark, pythia status 23) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquarkHardSctred_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquarkHardSctred_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquarkHardSctred_pt_all_central'}
            ]),
        ])
    }),        
    
    ("hLeadingPtGenBquarkHardSctred_eta", {
        sXLabel: 'eta(pT leading b-quark, pythia status 23)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBquarkHardSctred_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBquarkHardSctred_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBquarkHardSctred_eta_all_central'}
            ]),
        ])
    }),

    
   
    ("hLeadingPtGenBHadron_pt", {
        sXLabel: 'pT(pT leading b-hadron) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadron_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadron_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadron_pt_all_central'}
            ]),
        ])
    }),    
        
    ("hLeadingPtGenBHadron_eta", {
        sXLabel: 'eta(pT leading b-hadron)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadron_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadron_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadron_eta_all_central'}
            ]),
        ])
    }),
    

   
    ("hLeadingPtGenBHadronStatus2_pt", {
        sXLabel: 'pT(pT leading b-hadron with status=2) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadronStatus2_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadronStatus2_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadronStatus2_pt_all_central'}
            ]),
        ])
    }),    
    
    
    ("hLeadingPtGenBHadronStatus2_eta", {
        sXLabel: 'eta(pT leading b-hadron with status=2)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hLeadingPtGenBHadronStatus2_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hLeadingPtGenBHadronStatus2_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hLeadingPtGenBHadronStatus2_eta_all_central'}
            ]),
        ])
    }),
    
    
    
    ("hGenBquark_leadingPt", {
        sXLabel: 'pT(leading pT b-quark) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_leadingPt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_leadingPt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_leadingPt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_leadingPt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_leadingPt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_leadingPt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_leadingPt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_leadingPt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_leadingPt_all_central'}
            ]),
        ])
    }),    
   
    ("hGenBquark_subleadingPt", {
        sXLabel: 'pT(subleading pT b-quark) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_subleadingPt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_subleadingPt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_subleadingPt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_subleadingPt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_subleadingPt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_subleadingPt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_subleadingPt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_subleadingPt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_subleadingPt_all_central'}
            ]),
        ])
    }),    
   
    ("hGenBquark_thirdLeadingPt", {
        sXLabel: 'pT(3rd leading pT b-quark) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_thirdLeadingPt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_thirdLeadingPt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_thirdLeadingPt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_thirdLeadingPt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_thirdLeadingPt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_thirdLeadingPt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_thirdLeadingPt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_thirdLeadingPt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_thirdLeadingPt_all_central'}
            ]),
        ])
    }),    
    
    ("hGenBquark_forthLeadingPt", {
        sXLabel: 'pT(4th leading pT b-quark) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_forthLeadingPt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_forthLeadingPt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_forthLeadingPt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_forthLeadingPt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_forthLeadingPt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_forthLeadingPt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenBquark_forthLeadingPt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenBquark_forthLeadingPt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hGenBquark_forthLeadingPt_all_central'}
            ]),
        ])
    }),    
    
   
    
    ("hGenLHE_HT", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [0, 2500],
        sNRebin: 2,
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
    }),    
    

## Subleading ---------
   
    ("hSubleadingPtGenBquark_pt", {
        sXLabel: 'pT(pT sub-leading b-quark) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquark_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquark_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquark_pt_all_central'}
            ]),
        ])
    }),    
    
    
    ("hSubleadingPtGenBquark_eta", {
        sXLabel: 'eta(pT sub-leading b-quark)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquark_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquark_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquark_eta_all_central'}
            ]),
        ])
    }),

    
   
    ("hSubleadingPtGenBquarkHardSctred_pt", {
        sXLabel: 'pT(pT sub-leading b-quark, pythia status 23) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquarkHardSctred_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquarkHardSctred_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquarkHardSctred_pt_all_central'}
            ]),
        ])
    }),        
    
    ("hSubleadingPtGenBquarkHardSctred_eta", {
        sXLabel: 'eta(pT sub-leading b-quark, pythia status 23)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBquarkHardSctred_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBquarkHardSctred_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBquarkHardSctred_eta_all_central'}
            ]),
        ])
    }),

    
   
    ("hSubleadingPtGenBHadron_pt", {
        sXLabel: 'pT(pT sub-leading b-hadron) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadron_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadron_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadron_pt_all_central'}
            ]),
        ])
    }),    
        
    ("hSubleadingPtGenBHadron_eta", {
        sXLabel: 'eta(pT sub-leading b-hadron)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadron_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadron_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadron_eta_all_central'}
            ]),
        ])
    }),
    

   
    ("hSubleadingPtGenBHadronStatus2_pt", {
        sXLabel: 'pT(pT sub-leading b-hadron with status=2) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadronStatus2_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadronStatus2_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadronStatus2_pt_all_central'}
            ]),
        ])
    }),    
    
    
    ("hSubleadingPtGenBHadronStatus2_eta", {
        sXLabel: 'eta(pT sub-leading b-hadron with status=2)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hSubleadingPtGenBHadronStatus2_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hSubleadingPtGenBHadronStatus2_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hSubleadingPtGenBHadronStatus2_eta_all_central'}
            ]),
        ])
    }),
    

    
    
## Third-leading -----------------------------------------------------------------------------------------------------------------------------
   
    ("hThirdLeadingPtGenBquark_pt", {
        sXLabel: 'pT(pT third-leading b-quark) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquark_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquark_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquark_pt_all_central'}
            ]),
        ])
    }),    
    
    
    ("hThirdLeadingPtGenBquark_eta", {
        sXLabel: 'eta(pT third-leading b-quark)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquark_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquark_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquark_eta_all_central'}
            ]),
        ])
    }),

    
   
    ("hThirdLeadingPtGenBquarkHardSctred_pt", {
        sXLabel: 'pT(pT third-leading b-quark, pythia status 23) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquarkHardSctred_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquarkHardSctred_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquarkHardSctred_pt_all_central'}
            ]),
        ])
    }),        
    
    ("hThirdLeadingPtGenBquarkHardSctred_eta", {
        sXLabel: 'eta(pT third-leading b-quark, pythia status 23)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBquarkHardSctred_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBquarkHardSctred_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBquarkHardSctred_eta_all_central'}
            ]),
        ])
    }),

    
   
    ("hThirdLeadingPtGenBHadron_pt", {
        sXLabel: 'pT(pT third-leading b-hadron) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadron_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadron_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadron_pt_all_central'}
            ]),
        ])
    }),    
        
    ("hThirdLeadingPtGenBHadron_eta", {
        sXLabel: 'eta(pT third-leading b-hadron)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadron_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadron_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadron_eta_all_central'}
            ]),
        ])
    }),
    

   
    ("hThirdLeadingPtGenBHadronStatus2_pt", {
        sXLabel: 'pT(pT third-leading b-hadron with status=2) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadronStatus2_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadronStatus2_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadronStatus2_pt_all_central'}
            ]),
        ])
    }),    
    
    
    ("hThirdLeadingPtGenBHadronStatus2_eta", {
        sXLabel: 'eta(pT third-leading b-hadron with status=2)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hThirdLeadingPtGenBHadronStatus2_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hThirdLeadingPtGenBHadronStatus2_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hThirdLeadingPtGenBHadronStatus2_eta_all_central'}
            ]),
        ])
    }),
    
    

    
## Fourth-leading ------------------------------------------------------------------------------------------------------------------------------------------------
   
    ("hFourthLeadingPtGenBquark_pt", {
        sXLabel: 'pT(pT fourth-leading b-quark) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquark_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquark_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquark_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquark_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquark_pt_all_central'}
            ]),
        ])
    }),    
    
    
    ("hFourthLeadingPtGenBquark_eta", {
        sXLabel: 'eta(pT fourth-leading b-quark)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquark_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquark_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquark_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquark_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquark_eta_all_central'}
            ]),
        ])
    }),

    
   
    ("hFourthLeadingPtGenBquarkHardSctred_pt", {
        sXLabel: 'pT(pT fourth-leading b-quark, pythia status 23) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquarkHardSctred_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquarkHardSctred_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquarkHardSctred_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquarkHardSctred_pt_all_central'}
            ]),
        ])
    }),        
    
    ("hFourthLeadingPtGenBquarkHardSctred_eta", {
        sXLabel: 'eta(pT fourth-leading b-quark, pythia status 23)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquarkHardSctred_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBquarkHardSctred_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBquarkHardSctred_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBquarkHardSctred_eta_all_central'}
            ]),
        ])
    }),

    
   
    ("hFourthLeadingPtGenBHadron_pt", {
        sXLabel: 'pT(pT fourth-leading b-hadron) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadron_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadron_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadron_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadron_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadron_pt_all_central'}
            ]),
        ])
    }),    
        
    ("hFourthLeadingPtGenBHadron_eta", {
        sXLabel: 'eta(pT fourth-leading b-hadron)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadron_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadron_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadron_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadron_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadron_eta_all_central'}
            ]),
        ])
    }),
    

   
    ("hFourthLeadingPtGenBHadronStatus2_pt", {
        sXLabel: 'pT(pT fourth-leading b-hadron with status=2) [GeV]', sYLabel: 'Events',
        sXRange: [0, 125],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadronStatus2_pt_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadronStatus2_pt_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadronStatus2_pt_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadronStatus2_pt_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadronStatus2_pt_all_central'}
            ]),
        ])
    }),    
    
    
    ("hFourthLeadingPtGenBHadronStatus2_eta", {
        sXLabel: 'eta(pT fourth-leading b-hadron with status=2)', sYLabel: 'Events',
        sXRange: [-6, 6],
        sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD Stitch bQuark pT (option 1)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadronStatus2_eta_QCDStitchCutBQuarkPt_central'},
            ]),
            ("QCD Stitch bHadron (option 2)", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadronStatus2_eta_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hFourthLeadingPtGenBHadronStatus2_eta_all_central'},
            ]),
            ("QCD_bGen all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hFourthLeadingPtGenBHadronStatus2_eta_all_central'}
            ]),
            ("QCDIncl all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCDIncl/hFourthLeadingPtGenBHadronStatus2_eta_all_central'}
            ]),
        ])
    }),
    
    
    
])

