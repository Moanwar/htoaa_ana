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
    #(sAnaVersion, '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD_HT100to200/analyze_hadded_QCD_HT100to200.root')
    #(sAnaVersion, '/home/siddhesh/Work/CMS/htoaa/analysis/20230324_QCD/analyze_hadded_QCD.root')
    #(sAnaVersion, '')
    ('QCD_fullHT', '/eos/cms/store/user/ssawant/htoaa/analysis/20250314_CR_QCD4b_QCD/2018/analyze_htoaa_stage1.root')
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
sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250314_CR_QCD4b_QCD/2018/plots_1/%s' % (sAnaVersion)

histograms_dict = OD()

for cat in ["CR4b_3M2T", "CR4b_3M3T", "CR4b_4M3T", "CR4b_4M4T"]:
    histograms_dict["hGenLHE_HT_%s_QCD_bEnr" % (cat)] = {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [100, 2600],
        sNRebin: 1,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            (r'QCD $leq$2b', [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnr_0bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnr_1bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnr_2bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},                
            ]),
            (r'QCD 3b', [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnr_3bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},               
            ]),
            (r'QCD $geq$4b', [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_bEnr_4bAndMoreCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},               
            ]),
            
        ])
    }

    histograms_dict["hGenLHE_HT_%s_QCD_BGen" % (cat)] = {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [100, 2600],
        sNRebin: 1,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            (r'QCD $leq$2b', [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_BGen_0bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_BGen_1bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_BGen_2bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},                
            ]),
            (r'QCD 3b', [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_BGen_3bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},               
            ]),
            (r'QCD $geq$4b', [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_BGen_4bAndMoreCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},               
            ]),
            
        ])
    }

    histograms_dict["hGenLHE_HT_%s_QCD_Incl" % (cat)] = {
        sXLabel: 'LHE HT [GeV]', sYLabel: 'Events',
        sXRange: [100, 2600],
        sNRebin: 1,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            (r'QCD $leq$2b', [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl_0bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl_1bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl_2bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},                
            ]),
            (r'QCD 3b', [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl_3bCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},               
            ]),
            (r'QCD $geq$4b', [
                {sIpFileNameNice: sAnaVersion, sHistName: 'evt/QCD_Incl_4bAndMoreCat/hGenLHE_HT_SelQCD_%s_Nom'%(cat)},               
            ]),
            
        ])
    }

    
