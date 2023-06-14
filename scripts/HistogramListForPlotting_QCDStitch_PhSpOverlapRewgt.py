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
    ('20230519_StitchOverlapPhSpRemoval',                   '/eos/cms/store/user/ssawant/htoaa/analysis/20230519_StitchOverlapPhSpRemoval/2018/analyze_htoaa_stage1.root'), 
    ('20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt',    '/eos/cms/store/user/ssawant/htoaa/analysis/20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt/2018/analyze_htoaa_stage1.root'), 
    ('20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt',  '/eos/cms/store/user/ssawant/htoaa/analysis/20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt/2018/analyze_htoaa_stage1.root'),
    ('20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg',  '/eos/cms/store/user/ssawant/htoaa/analysis/20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg/2018/analyze_htoaa_stage1.root')
])
sAnaVersion = 'QCDStitch' #list(sIpFiles.keys())[0]
print(f"sAnaVersion: {sAnaVersion}")

#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/plots'
#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230519_StitchOverlapPhSpRemoval/2018/plots/%s' % (sAnaVersion)
#sOpDir  = '/home/siddhesh/Work/CMS/htoaa/analysis/20230504_v2_QCD_wHTSamplesStitchSF_step1/2018/plots/%s' % (sAnaVersion)
#sOpDir  = '/Users/siddhesh/Work/CMS/htoaa/analysis/20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt/2018/plots/' 
sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt/2018/plots/' 

histograms_dict = OD([
    #("hLeadingPtGenBquark_pt_all", {sXLabel: 'Leading FatJet mass [GeV]', sYLabel: 'Events', sXRange: [0, 200]}),
   
    ("LHE_HT_QCD_StitchPhSpOverlapRemove_woQCDbGenHTRewgt", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [100, 3000], sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl_all", [
                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_all_central'},
            ]),
            ("QCD stitch w/ ph.sp. overlap removed, w/o QCD_bGen HT reweight", [
                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]), 
            ("QCD_Incl stitch w/ ph.sp. overlap removed", [
                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]), 
            ("QCD_bGen stitch w/ ph.sp. overlap removed, w/o QCD_bGen HT reweight", [
                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]), 
            ("QCD_bEnrich stitch w/ ph.sp. overlap removed", [
                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]), 
        ])
    }),    
   
    ("LHE_HT_QCD_StitchPhSpOverlapRemove_wQCDbGenHTRewgt", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [100, 4000], sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl_all", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_all_central'},
            ]),
            ("QCD stitch w/ ph.sp. overlap removed, w/ QCD_bGen HT reweight ", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]), 
            ("QCD_Incl stitch w/ ph.sp. overlap removed ", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]), 
            ("QCD_bGen stitch w/ ph.sp. overlap removed, w/ QCD_bGen HT reweight ", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]), 
            ("QCD_bEnrich stitch w/ ph.sp. overlap removed ", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]), 
        ])
    }),    

    ("LHE_HT_QCD_StitchPhSpOverlapRewgt_wQCDbGenHTRewgt", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [100, 4000], sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl_all", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_all_central'},
            ]),
            ("QCD stitch w/ ph.sp. overlap reweighted, w/ QCD_bGen HT reweight ", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCDStitch_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitch_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitch_central'},
            ]), 

            ("QCD_Incl_Remnant PhSp w/ ph.sp. overlap reweighted stitch, w/ QCD_bGen HT reweight ", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
            ]), 
            ("QCD_bGen PhSp w/ ph.sp. overlap reweighted stitch, w/ QCD_bGen HT reweight ", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_bGen_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_bGen_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_bGen_PhSp_central'},
            ]), 
            ("QCD_bEnrich PhSp w/ ph.sp. overlap reweighted stitch, w/ QCD_bGen HT reweight ", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
            ]), 


        ])
    }),    

 
    ("LHE_HT_QCD_StitchPhSpOverlapRewgt_woQCDbGenHTRewgt", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [100, 4000], sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl_all", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_all_central'},
            ]),
            ("QCD stitch w/ ph.sp. overlap reweighted, w/o QCD_bGen HT reweight ", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCDStitch_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitch_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitch_central'},
            ]), 

            ("QCD_Incl_Remnant PhSp w/ ph.sp. overlap reweighted stitch", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
            ]), 
            ("QCD_bGen PhSp w/ ph.sp. overlap reweighted stitch", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_bGen_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_bGen_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_bGen_PhSp_central'},
            ]), 
            ("QCD_bEnrich PhSp w/ ph.sp. overlap reweighted stitch", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
            ]), 

        ])
    }),    

    
    ("LHE_HT_QCD_StitchCompare", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [50, 4000], #sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl_all", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_all_central'},
            ]),
#            ("QCD stitch w/ ph.sp. overlap removed, w/o QCD_bGen HT reweight", [
#                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_QCDStitchCutBHadron_central'},
#                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBHadron_central'},
#                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBHadron_central'},
#            ]), 
            ("QCD stitch w/ ph.sp. overlap removed, w/ QCD_bGen HT reweight ", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]), 
            ("QCD stitch w/ ph.sp. overlap reweighted, w/ QCD_bGen HT reweight ", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCDStitch_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitch_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitch_central'},
            ]), 
#            ("QCD stitch w/ ph.sp. overlap reweighted, w/o QCD_bGen HT reweight ", [
#                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCDStitch_central'},
#                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitch_central'},
#                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitch_central'},
#            ]), 
       ])
    }),    
    
     ("LHE_HT_QCD_StitchCompare_logX", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [50, 4000], sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl_all", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_all_central'},
            ]),
#            ("QCD stitch w/ ph.sp. overlap removed, w/o QCD_bGen HT reweight", [
#                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_QCDStitchCutBHadron_central'},
#                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBHadron_central'},
#                {sIpFileNameNice: '20230519_StitchOverlapPhSpRemoval', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBHadron_central'},
#            ]), 
            ("QCD stitch w/ ph.sp. overlap removed, w/ QCD_bGen HT reweight ", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBHadron_central'},
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]), 
            ("QCD stitch w/ ph.sp. overlap reweighted, w/ QCD_bGen HT reweight ", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCDStitch_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitch_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitch_central'},
            ]), 
#            ("QCD stitch w/ ph.sp. overlap reweighted, w/o QCD_bGen HT reweight ", [
#                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCDStitch_central'},
#                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitch_central'},
#                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitch_central'},
#            ]), 
       ])
    }),    
  

 
    ("LHE_HT_QCD_Incl_Remnant_Compare", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [100, 4000], sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_Incl_Remnant w/ ph.sp. overlap removed", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_Incl_Remnant w/ ph.sp. overlap reweighted, w/o QCD_bGen HT reweight", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
            ]),
            ("QCD_Incl_Remnant w/ ph.sp. overlap reweighted, w/ QCD_bGen HT reweight", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_Incl_Remnant_PhSp_central'},
            ]),

        ])
    }),    

    ("LHE_HT_QCD_bGen_Compare", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [100, 4000], sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_bGen w/ ph.sp. overlap removed", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bGen w/ ph.sp. overlap removed, v2", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_SelQCDbGen_central'},
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_SelQCDbGen_central'},
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_SelQCDbGen_central'},
            ]),

            ("QCD_bGen w/ ph.sp. overlap reweighted, w/o QCD_bGen HT reweight", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_bGen_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_bGen_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_bGen_PhSp_central'},
            ]),
            ("QCD_bGen w/ ph.sp. overlap reweighted, w/ QCD_bGen HT reweight", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_bGen_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_bGen_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_bGen_PhSp_central'},
            ]),

        ])
    }),    


    ("LHE_HT_QCD_bEnrich_Compare", {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [100, 4000], sXScale: 'log_10',
        sNRebin: 10, 
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("QCD_bEnrich w/ ph.sp. overlap removed", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCDStitchCutBHadron_central'},
            ]),
            ("QCD_bEnrich w/ ph.sp. overlap removed, v2", [
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_Incl_PSWeight/hGenLHE_HT_SelQCDbEnrich_central'},
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_SelQCDbEnrich_central'},
                {sIpFileNameNice: '20230602_StitchOverlapPhSpRemoval_QCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_SelQCDbEnrich_central'},
            ]),

            ("QCD_bEnrich w/ ph.sp. overlap reweighted, w/o QCD_bGen HT reweight", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_woQCDbGenHTRewg', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
            ]),
            ("QCD_bEnrich w/ ph.sp. overlap reweighted, w/ QCD_bGen HT reweight", [
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_Incl/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bGen/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
                {sIpFileNameNice: '20230613_QCDStitchPhSpOverlapRewgt_wQCDbGenHTRewgt', sHistName: 'evt/QCD_bEnrich/hGenLHE_HT_QCD_bEnrich_PhSp_central'},
            ]),

        ])
    }),    


])

