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
    # (<file name to refer>, <file path+name>
    ('20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt','/Users/siddhesh/Work/CMS/htoaa/analysis/20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt/2018/analyze_htoaa_stage1.root'), 
    ('20230519_StitchOverlapPhSpRemoval',               '/Users/siddhesh/Work/CMS/htoaa/analysis/20230519_StitchOverlapPhSpRemoval/2018/analyze_htoaa_stage1.root')
])
sAnaVersion = list(sIpFiles.keys())[0]
print(f"sAnaVersion: {sAnaVersion}")

#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/plots'
#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230519_StitchOverlapPhSpRemoval/2018/plots/%s' % (sAnaVersion)
#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230504_v2_QCD_wHTSamplesStitchSF_step1/2018/plots/%s' % (sAnaVersion)
sOpDir  = '/Users/siddhesh/Work/CMS/htoaa/analysis/20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt/2018/plots/%s' % (sAnaVersion)

histograms_dict = OD([
    #("hLeadingPtGenBquark_pt_all", {sXLabel: 'Leading FatJet mass [GeV]', sYLabel: 'Events', sXRange: [0, 200]}),
   
    ("LHE_HT_QCD_bEnrich_cut_0", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [50, 3000], #sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_bEnrich_all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_all_central'},
            ]),
            ("QCD_bEnrich + QCD_bEnrich selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_SelQCDbEnrich_central'},
            ]), 
        ])
    }),    

     ("LHE_HT_QCD_bGen_cut_0", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [50, 3000], #sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_bGen_all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_all_central'},
            ]),
            ("QCD_bGen + QCD_bGen selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_SelQCDbGen_central'},
            ]), 
        ])
    }),    
   
   
    ("LHE_HT_QCD_bEnrich_cut_2", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [50, 3000], #sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl-madgraphMLM + QCD_bEnrich selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl/hGenLHE_HT_SelQCDbEnrich_central'},
            ]), 
            ("QCD_bEnrich_all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_all_central'},
            ]),
            ("QCD_bEnrich + QCD_bEnrich selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_SelQCDbEnrich_central'},
            ]), 
        ])
    }),    

     ("LHE_HT_QCD_bGen_cut_2", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [50, 3000], #sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl-madgraphMLM + QCD_bGen selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl/hGenLHE_HT_SelQCDbGen_central'},
            ]), 
            ("QCD_bGen_all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_all_central'},
            ]),
            ("QCD_bGen + QCD_bGen selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_SelQCDbGen_central'},
            ]), 
        ])
    }),    
 

   
    ("LHE_HT_QCD_bEnrich_cut_3", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [50, 3000], #sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl-madgraph_PSWeights + QCD_bEnrich selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_SelQCDbEnrich_central'},
            ]), 
            ("QCD_bEnrich_all", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_all_central'},
            ]),
            ("QCD_bEnrich + QCD_bEnrich selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_SelQCDbEnrich_central'},
            ]), 
        ])
    }),    

     ("LHE_HT_QCD_bGen_cut_3", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [50, 4000], #sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl-madgraph_PSWeights + QCD_bGen selection", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_SelQCDbGen_central'},
            ]), 
            
            
#            ("QCD_bGen_all, w/o HT reweight", [
#                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_bGen/hGenLHE_HT_all_central'},
#            ]),
#            ("QCD_bGen + QCD_bGen selection, w/o HT reweight", [
#                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_bGen/hGenLHE_HT_SelQCDbGen_central'},
#            ]), 
            
            
            ("QCD_bGen_all, w/ HT reweight", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_all_central'},
            ]),
            ("QCD_bGen + QCD_bGen selection, w/ HT reweight", [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bGen/hGenLHE_HT_SelQCDbGen_central'},
            ]), 
            
        ])
    }),    
    
])

