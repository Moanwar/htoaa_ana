# %%
import os, sys
import numpy as np
from collections import OrderedDict as OD
import math
#import uproot3
import uproot as uproot
import hist
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import enum
import mplhep as hep
from parse import *
import copy
import json

sys.path.append( os.path.abspath('../') )
print(f"{os.path.abspath('../') = }")

from htoaa_CommonTools import (
    rebinTH1, rebinTH2, variableRebinTH1,
)

#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20240809_gg0l_FullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240809_gg0l_FullSyst/2018/2DAlphabet_inputFiles'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250123_gg0l_NoSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250123_gg0l_NoSyst/2018/2DAlphabet_inputFiles'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250123_gg0l_PNetX4b_versions/2018/analyze_htoaa_stage1.root'
#sOpDir0  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250123_gg0l_PNetX4b_versions/2018/2DAlphabet_inputFiles'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250207_gg0l_NoSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250207_gg0l_NoSyst/2018/2DAlphabet_inputFiles'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250207_gg0l_NoSyst_1/2018/analyze_htoaa_stage1.root'
#sOpDir0  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250207_gg0l_NoSyst_1/2018/2DAlphabet_inputFiles'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_gg0l_FullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_gg0l_FullSyst/2018/2DAlphabet_inputFiles'
#CAT0 = 'gg0l'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad', 'Zvv'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_ttHad_FullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_ttHad_FullSyst/2018/2DAlphabet_inputFiles'
#CAT0 = 'ttHad'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad', 'Zvv'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_Vjj_FullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_Vjj_FullSyst/2018/2DAlphabet_inputFiles'
#CAT0 = 'Vjj'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad', 'Zvv'
sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250502_CR_QCD4b_FullSyst/2018/analyze_htoaa_stage1.root'
sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250502_CR_QCD4b_FullSyst/2018/CombineFit_inputFiles'
CAT0 = 'CR_QCD4b'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad', 'Zvv', 'CR_QCD4b' 


Era = '2018' # '2018, 'Run2
# Original categories name: 'gg0l', 'VBFjj', 'Wlv', 'Zll', 'Zvv',  'Vjj'. 'ZvvIncl','ZvvLo', 'ZvvHi', 'gg0lIncl', 'gg0lLo', 'gg0lHi',

CATAGORIES_gg0l = {
    "gg0lIncl" : "gg0lIncl_Xto4bv2",
    "gg0lHi" :   "gg0lHi_Xto4bv2",
    "gg0lLo" :   "gg0lLo_Xto4bv2",    
}
CATAGORIES_ttHad = {
    "ttHad" : "ttHad_Xto4bv2",   
}
CATAGORIES_Vjj = {
    "Vjj" : "Vjj_Xto4bv2",   
}
CATAGORIES_Zvv = {
    "ZvvIncl" : "ZvvIncl_Xto4bv2",
    "ZvvHi" :   "ZvvHi_Xto4bv2",
    "ZvvLo" :   "ZvvLo_Xto4bv2",    
}
CATAGORIES_CR_QCD4b = {
    "CR4b_3M2T" : "CR4b_3M2T",
    "CR4b_3M3T" : "CR4b_3M3T",
    "CR4b_4M3T" : "CR4b_4M3T",
    "CR4b_4M4T" : "CR4b_4M4T",        
}

WPs_perCategory = {
    'gg0l':  ['WP40'],
    'VBFjj': ['WP40'],
    'Vjj':   ['WP60'],
    'ttHad': ['WP60'],
    'Zvv':   ['WP60'],     
    'CR_QCD4b':   ['NoWP'],     
}

fIpFile = uproot.open(sIpFile)

# %%

processes_dict = {
    'DataJetHT': ['JetHT_Run2018A', 'JetHT_Run2018B', 'JetHT_Run2018C', 'JetHT_Run2018D'],
    'DataMET':   ['MET_Run2018A', 'MET_Run2018B', 'MET_Run2018C', 'MET_Run2018D'],
}

if CAT0 not in ['CR_QCD4b']:
    processes_dict.update( {
        'ggHtoaato4b_mA_12': ['ggHtoaato4b_mA_12'],
        'ggHtoaato4b_mA_15': ['ggHtoaato4b_mA_15'],
        'ggHtoaato4b_mA_20': ['ggHtoaato4b_mA_20'],
        'ggHtoaato4b_mA_25': ['ggHtoaato4b_mA_25'],
        'ggHtoaato4b_mA_30': ['ggHtoaato4b_mA_30'],
        'ggHtoaato4b_mA_35': ['ggHtoaato4b_mA_35'],
        'ggHtoaato4b_mA_40': ['ggHtoaato4b_mA_40'],
        'ggHtoaato4b_mA_45': ['ggHtoaato4b_mA_45'],
        'ggHtoaato4b_mA_50': ['ggHtoaato4b_mA_50'],
        'ggHtoaato4b_mA_55': ['ggHtoaato4b_mA_55'],
        'ggHtoaato4b_mA_60': ['ggHtoaato4b_mA_60'],

        'VBFHtoaato4b_mA_12': ['VBFHtoaato4b_mA_12'],
        'VBFHtoaato4b_mA_15': ['VBFHtoaato4b_mA_15'],
        'VBFHtoaato4b_mA_20': ['VBFHtoaato4b_mA_20'],
        'VBFHtoaato4b_mA_25': ['VBFHtoaato4b_mA_25'],
        'VBFHtoaato4b_mA_30': ['VBFHtoaato4b_mA_30'],
        'VBFHtoaato4b_mA_35': ['VBFHtoaato4b_mA_35'],
        'VBFHtoaato4b_mA_40': ['VBFHtoaato4b_mA_40'],
        'VBFHtoaato4b_mA_45': ['VBFHtoaato4b_mA_45'],
        'VBFHtoaato4b_mA_50': ['VBFHtoaato4b_mA_50'],
        'VBFHtoaato4b_mA_55': ['VBFHtoaato4b_mA_55'],
        'VBFHtoaato4b_mA_60': ['VBFHtoaato4b_mA_60'],

        'WHtoaato4b_mA_12': ['WHtoaato4b_mA_12'],
        'WHtoaato4b_mA_15': ['WHtoaato4b_mA_15'],
        'WHtoaato4b_mA_20': ['WHtoaato4b_mA_20'],
        'WHtoaato4b_mA_25': ['WHtoaato4b_mA_25'],
        'WHtoaato4b_mA_30': ['WHtoaato4b_mA_30'],
        'WHtoaato4b_mA_35': ['WHtoaato4b_mA_35'],
        'WHtoaato4b_mA_40': ['WHtoaato4b_mA_40'],
        'WHtoaato4b_mA_45': ['WHtoaato4b_mA_45'],
        'WHtoaato4b_mA_50': ['WHtoaato4b_mA_50'],
        'WHtoaato4b_mA_55': ['WHtoaato4b_mA_55'],
        'WHtoaato4b_mA_60': ['WHtoaato4b_mA_60'],

        'ZHtoaato4b_mA_12': ['ZHtoaato4b_mA_12'],
        'ZHtoaato4b_mA_15': ['ZHtoaato4b_mA_15'],
        'ZHtoaato4b_mA_20': ['ZHtoaato4b_mA_20'],
        'ZHtoaato4b_mA_25': ['ZHtoaato4b_mA_25'],
        'ZHtoaato4b_mA_30': ['ZHtoaato4b_mA_30'],
        'ZHtoaato4b_mA_35': ['ZHtoaato4b_mA_35'],
        'ZHtoaato4b_mA_40': ['ZHtoaato4b_mA_40'],
        'ZHtoaato4b_mA_45': ['ZHtoaato4b_mA_45'],
        'ZHtoaato4b_mA_50': ['ZHtoaato4b_mA_50'],
        'ZHtoaato4b_mA_55': ['ZHtoaato4b_mA_55'],
        'ZHtoaato4b_mA_60': ['ZHtoaato4b_mA_60'],

        'ttHtoaato4b_mA_12': ['ttHtoaato4b_mA_12'],
        'ttHtoaato4b_mA_15': ['ttHtoaato4b_mA_15'],
        'ttHtoaato4b_mA_20': ['ttHtoaato4b_mA_20'],
        'ttHtoaato4b_mA_25': ['ttHtoaato4b_mA_25'],
        'ttHtoaato4b_mA_30': ['ttHtoaato4b_mA_30'],
        'ttHtoaato4b_mA_35': ['ttHtoaato4b_mA_35'],
        'ttHtoaato4b_mA_40': ['ttHtoaato4b_mA_40'],
        'ttHtoaato4b_mA_45': ['ttHtoaato4b_mA_45'],
        'ttHtoaato4b_mA_50': ['ttHtoaato4b_mA_50'],
        'ttHtoaato4b_mA_55': ['ttHtoaato4b_mA_55'],
        'ttHtoaato4b_mA_60': ['ttHtoaato4b_mA_60'],
    } )

if CAT0 not in ['CR_QCD4b']:
    processes_dict.update( {
        'QCD_bEnr': ['QCD_bEnr'],
        'QCD_BGen': ['QCD_BGen'],
        'QCD_Incl': ['QCD_Incl'],
    } )
else: 
    processes_dict.update( {
        'QCD_012bCat': ['QCD_0bCat', 'QCD_1bCat', 'QCD_2bCat', ],
        'QCD_3bGen': ['QCD_3bCat'],
        'QCD_4bAndMoreCat': ['QCD_4bAndMoreCat'],
    } )

processes_dict.update( {
    'TT0l': ['TT0l'],
    'TT1l': ['TT1l'],
    'TT2l': ['TT2l'],
    'STop_t': ['STop_t'],
    'STbar_t': ['STbar_t'],
    'ST_s_0l': ['ST_s_0l'],
    'ST_s_1l': ['ST_s_1l'],
    'STop_tW_Incl': ['STop_tW_Incl'],
    'STbar_tW_Incl': ['STbar_tW_Incl'],
    'STop_tW_12l': ['STop_tW_12l'],
    'STbar_tW_12l': ['STbar_tW_12l'],
    #'ttZ': ['ttZ'],    
    #'ttW': ['ttW'],
    #'tZq': ['tZq'],
    'Zqq': ['Zqq'],
    'Zvv': ['Zvv'],
    
    'Zll': ['Zll'],
    'Wqq': ['Wqq'],
    'Wlv': ['Wlv'],
    
    'ZZ': ['ZZ'],
    'WZ': ['WZ'],
    'WW': ['WW'],
    
    'ZZZ': ['ZZZ'],
    'WZZ': ['WZZ'],
    'WWZ': ['WWZ'],
    'WWW': ['WWW'],
    'ggH': ['ggH'],

} )
systematics_forData_dict = {'noweight': 'Nom'}
systematics_dict = {
    'Nom': 'Nom',
    'PUUp': 'PUUp',
    'PUDown': 'PUDown',  
    'LPRewgtUp': 'LPRewgtUp',
    'LPRewgtDown': 'LPRewgtDown',
    'GGHPtRewgtUp': 'GGHPtRewgtUp',
    'GGHPtRewgtDown': 'GGHPtRewgtDown',
    'TopPtReWeightUp': 'TopPtReWeightUp',
    'TopPtReWeightDown': 'TopPtReWeightDown',
    'ISRUp':'ISRUp',
    'ISRDown':'ISRDown',
    'FSRUp':'FSRUp',
    'FSRDown':'FSRDown',
    'QCDRenormUp':'QCDRenormUp',
    'QCDRenormDown':'QCDRenormDown',
    'QCDFactrUp':'QCDFactrUp',
    'QCDFactrDown':'QCDFactrDown',
    'PDFUp':'PDFUp',
    'PDFDown':'PDFDown',
    'BtagUp': 'BtagUp',
    'BtagDown': 'BtagDown',
    'BtagUncorrUp': 'BtagUncorrUp',
    'BtagUncorrDown': 'BtagUncorrDown',
    'BtagCorrUp': 'BtagCorrUp',
    'BtagCorrDown': 'BtagCorrDown', 

    'JERUp': 'JERUp',
    'JERDown': 'JERDown',
    'JESUp': 'JESUp',
    'JESDown': 'JESDown',
    'JESHEMIssueUp': 'JESHEMIssueUp',
    'JESHEMIssueDown': 'JESHEMIssueDown',   
     
}
#systematics_dict = {'Nom': 'Nom'} #{'central': 'Nom'}

syst_MCAll       = ['Nom', 'PUUp', 'PUDown',              
                    'ISRUp','ISRDown','FSRUp','FSRDown','QCDRenormUp','QCDRenormDown','QCDFactrUp','QCDFactrDown','PDFUp','PDFDown',              
                    'JERUp','JERDown','JESUp','JESDown','JESHEMIssueUp','JESHEMIssueDown', 
                    ]
syst_MCSingalH    = ['LPRewgtUp', 'LPRewgtDown',]
syst_MCSingalGGH  = ['GGHPtRewgtUp','GGHPtRewgtDown',]
syst_MCTT         = ['TopPtReWeightUp','TopPtReWeightDown',]
syst_MCBtagYerly  = ['BtagUp','BtagDown',]
syst_MCBtagRun2   = ['BtagUncorrUp','BtagCorrDown', ]
syst_MCBtag = syst_MCBtagRun2 if 'run' in Era.lower() else syst_MCBtagYerly
syst_MCAll += syst_MCBtag

systematics_perProcess = {
    'Data': ['Nom'],

    'ggHtoaato4b_mA_12': syst_MCAll + syst_MCSingalH + syst_MCSingalGGH,
    'ggHtoaato4b_mA_15': syst_MCAll + syst_MCSingalH + syst_MCSingalGGH,
    'ggHtoaato4b_mA_20': syst_MCAll + syst_MCSingalH + syst_MCSingalGGH,
    'ggHtoaato4b_mA_25': syst_MCAll + syst_MCSingalH + syst_MCSingalGGH,
    'ggHtoaato4b_mA_30': syst_MCAll + syst_MCSingalH + syst_MCSingalGGH,
    'ggHtoaato4b_mA_35': syst_MCAll + syst_MCSingalH + syst_MCSingalGGH,
    'ggHtoaato4b_mA_40': syst_MCAll + syst_MCSingalH + syst_MCSingalGGH,
    'ggHtoaato4b_mA_45': syst_MCAll + syst_MCSingalH + syst_MCSingalGGH,
    'ggHtoaato4b_mA_50': syst_MCAll + syst_MCSingalH + syst_MCSingalGGH,
    'ggHtoaato4b_mA_55': syst_MCAll + syst_MCSingalH + syst_MCSingalGGH,
    'ggHtoaato4b_mA_60': syst_MCAll + syst_MCSingalH + syst_MCSingalGGH,

    'VBFHtoaato4b_mA_12': syst_MCAll + syst_MCSingalH,
    'VBFHtoaato4b_mA_15': syst_MCAll + syst_MCSingalH,
    'VBFHtoaato4b_mA_20': syst_MCAll + syst_MCSingalH,
    'VBFHtoaato4b_mA_25': syst_MCAll + syst_MCSingalH,
    'VBFHtoaato4b_mA_30': syst_MCAll + syst_MCSingalH,
    'VBFHtoaato4b_mA_35': syst_MCAll + syst_MCSingalH,
    'VBFHtoaato4b_mA_40': syst_MCAll + syst_MCSingalH,
    'VBFHtoaato4b_mA_45': syst_MCAll + syst_MCSingalH,
    'VBFHtoaato4b_mA_50': syst_MCAll + syst_MCSingalH,
    'VBFHtoaato4b_mA_55': syst_MCAll + syst_MCSingalH,
    'VBFHtoaato4b_mA_60': syst_MCAll + syst_MCSingalH,

    'WHtoaato4b_mA_12': syst_MCAll + syst_MCSingalH,
    'WHtoaato4b_mA_15': syst_MCAll + syst_MCSingalH,
    'WHtoaato4b_mA_20': syst_MCAll + syst_MCSingalH,
    'WHtoaato4b_mA_25': syst_MCAll + syst_MCSingalH,
    'WHtoaato4b_mA_30': syst_MCAll + syst_MCSingalH,
    'WHtoaato4b_mA_35': syst_MCAll + syst_MCSingalH,
    'WHtoaato4b_mA_40': syst_MCAll + syst_MCSingalH,
    'WHtoaato4b_mA_45': syst_MCAll + syst_MCSingalH,
    'WHtoaato4b_mA_50': syst_MCAll + syst_MCSingalH,
    'WHtoaato4b_mA_55': syst_MCAll + syst_MCSingalH,
    'WHtoaato4b_mA_60': syst_MCAll + syst_MCSingalH,

    'ZHtoaato4b_mA_12': syst_MCAll + syst_MCSingalH,
    'ZHtoaato4b_mA_15': syst_MCAll + syst_MCSingalH,
    'ZHtoaato4b_mA_20': syst_MCAll + syst_MCSingalH,
    'ZHtoaato4b_mA_25': syst_MCAll + syst_MCSingalH,
    'ZHtoaato4b_mA_30': syst_MCAll + syst_MCSingalH,
    'ZHtoaato4b_mA_35': syst_MCAll + syst_MCSingalH,
    'ZHtoaato4b_mA_40': syst_MCAll + syst_MCSingalH,
    'ZHtoaato4b_mA_45': syst_MCAll + syst_MCSingalH,
    'ZHtoaato4b_mA_50': syst_MCAll + syst_MCSingalH,
    'ZHtoaato4b_mA_55': syst_MCAll + syst_MCSingalH,
    'ZHtoaato4b_mA_60': syst_MCAll + syst_MCSingalH,

    'ttHtoaato4b_mA_12': syst_MCAll + syst_MCSingalH,
    'ttHtoaato4b_mA_15': syst_MCAll + syst_MCSingalH,
    'ttHtoaato4b_mA_20': syst_MCAll + syst_MCSingalH,
    'ttHtoaato4b_mA_25': syst_MCAll + syst_MCSingalH,
    'ttHtoaato4b_mA_30': syst_MCAll + syst_MCSingalH,
    'ttHtoaato4b_mA_35': syst_MCAll + syst_MCSingalH,
    'ttHtoaato4b_mA_40': syst_MCAll + syst_MCSingalH,
    'ttHtoaato4b_mA_45': syst_MCAll + syst_MCSingalH,
    'ttHtoaato4b_mA_50': syst_MCAll + syst_MCSingalH,
    'ttHtoaato4b_mA_55': syst_MCAll + syst_MCSingalH,
    'ttHtoaato4b_mA_60': syst_MCAll + syst_MCSingalH,
}
if CAT0 not in ['CR_QCD4b']:
    systematics_perProcess.update( {
        'QCD_bEnr': syst_MCAll,
        'QCD_BGen': syst_MCAll,
        'QCD_Incl': syst_MCAll,
    } )
else:
    systematics_perProcess.update( {
        'QCD_012bCat': syst_MCAll,
        'QCD_3bGen': syst_MCAll,
        'QCD_4bAndMoreCat': syst_MCAll,
    } )

systematics_perProcess.update( {
    'TT0l': syst_MCAll + syst_MCTT,
    'TT1l': syst_MCAll + syst_MCTT,
    'TT2l': syst_MCAll + syst_MCTT,
    'STop_t': syst_MCAll,
    'STbar_t': syst_MCAll,
    'ST_s_0l': syst_MCAll,
    'ST_s_1l': syst_MCAll,
    'STop_tW_Incl': syst_MCAll,
    'STbar_tW_Incl': syst_MCAll,
    'STop_tW_12l': syst_MCAll,
    'STbar_tW_12l': syst_MCAll,
    #'ttZ': syst_MCAll,    
    #'ttW': syst_MCAll,
    #'tZq': syst_MCAll,
    'Zqq': syst_MCAll,
    'Zvv': syst_MCAll,
    
    'Zll': syst_MCAll,
    'Wqq': syst_MCAll,
    'Wlv': syst_MCAll,
    
    'ZZ': syst_MCAll,
    'WZ': syst_MCAll,
    'WW': syst_MCAll,
    
    'ZZZ': syst_MCAll,
    'WZZ': syst_MCAll,
    'WWZ': syst_MCAll,
    'WWW': syst_MCAll,
    'ggH': syst_MCAll,

} )



#histograms_list = [
#    'hLeadingFatJetMass_vs_massA_Hto4b_avg',
#    'hLeadingFatJetMSoftDrop_vs_massA_Hto4b_avg',
#    'hLeadingFatJetParticleNet_massH_Hto4b_avg_vs_massA_Hto4b_avg',
#]
histograms_dict_v1 = {
    'hLeadingFatJetMass_vs_massA_Hto4b_avg':                        'mass',
    'hLeadingFatJetMSoftDrop_vs_massA_Hto4b_avg':                   'msoft',
    'hLeadingFatJetParticleNet_massH_Hto4b_avg_vs_massA_Hto4b_avg': 'pnet',    
}
histograms_dict_v0 = {
    'hLeadingFatJetMass_vs_massAa':                        'mass',
    'hLeadingFatJetMSoftDrop_vs_massAa':                   'msoft',
    'hLeadingFatJetPNet_massH_v2b_vs_massAa': 'pnet',    
}
histograms_dict = {
    #'hLeadingFatJetMass_vs_massAa':                        'mass_vs_massAa',
    #'hLeadingFatJetMSoftDrop_vs_massAa':                   'msoft_vs_massAa',
    'hLeadingFatJetPNet_massH_v2b_vs_massAa':              'pnet_vs_massAa',    

    #'hLeadingFatJetMass_vs_massA34a':                        'mass_vs_massA34a',
    #'hLeadingFatJetMSoftDrop_vs_massA34a':                   'msoft_vs_massA34a',
    'hLeadingFatJetPNet_massH_v2b_vs_massA34a':              'pnet_vs_massA34a',  ##  

    #'hLeadingFatJetMass_vs_massA34b':                        'mass_vs_massA34b',

    #'hLeadingFatJetMass_vs_massA34d':                        'mass_vs_massA34d',
    #'hLeadingFatJetMSoftDrop_vs_massA34d':                   'msoft_vs_massA34d',
    'hLeadingFatJetPNet_massH_v2b_vs_massA34d':              'pnet_vs_massA34d',  ##   
    
}

nRebinsX = 10
nRebinsY = 4

if CAT0 in ['CR_QCD4b']:
    histograms_dict = {
        'hLeadingFatJetPNet_X4b_v1_Haa4b_vs_QCD': 'X4b_v1_Haa4b_vs_QCD',
        'hLeadingFatJetPNet_X4b_v1_Haa4b_score': 'X4b_v1_Haa4b_score',
        'hLeadingFatJetPNet_X4b_v2a_Haa4b_score': 'X4b_v2a_Haa4b_score',
        'hLeadingFatJetPNet_X4b_v2b_Haa4b_score': 'X4b_v2b_Haa4b_score',
        'hLeadingFatJetPNet_X4b_v2ab_Haa4b_score': 'X4b_v2ab_Haa4b_score',
        'hLeadingFatJetPNet_X4b_v2a_Haa34b_score': 'X4b_v2a_Haa34b_score',
        'hLeadingFatJetPNet_X4b_v2b_Haa34b_score': 'X4b_v2b_Haa34b_score',
        'hLeadingFatJetPNet_X4b_v2ab_Haa34b_score': 'X4b_v2ab_Haa34b_score',
    }
    nRebinsX = 1
    nRebinsY = 1


selectionTags_dict_v1 = {
    'WP40': {
        'Pass': 'SRWP40',
        'Fail': 'SBWP80to40'
    },
    'WP60': {
        'Pass': 'SRWP60',
        'Fail': 'SBWP95to60'
    }, 
    'WP80': {
        'Pass': 'SRWP80',
        'Fail': 'SBWP99to80'
    },        
}  
selectionTags_dict = {
    'WP40': {
        'Pass': 'SRWP40',
        'Fail': 'SBWP40'
    },
    'WP45a': {
        'Pass': 'SRWP45a',
        'Fail': 'SBWP45a'
    },
    'WP45b': {
        'Pass': 'SRWP45b',
        'Fail': 'SBWP45b'
    },
    'WP50': {
        'Pass': 'SRWP50',
        'Fail': 'SBWP50'
    },
    'WP60': {
        'Pass': 'SRWP60',
        'Fail': 'SBWP60'
    },
    'WP65': {
        'Pass': 'SRWP65',
        'Fail': 'SBWP65'
    },
    'WP70': {
        'Pass': 'SRWP70',
        'Fail': 'SBWP70'
    },
    'WP80': {
        'Pass': 'SRWP80',
        'Fail': 'SBWP80'
    },
    'NoWP': {
        'NoWP': 'NoWP'
    }
            
}  


if 'gg0l' in CAT0:
    CATAGORIES = CATAGORIES_gg0l
elif 'ttHad' in CAT0:
    CATAGORIES = CATAGORIES_ttHad
elif 'Vjj' in CAT0:
    CATAGORIES = CATAGORIES_Vjj
elif 'Zvv' in CAT0:
    CATAGORIES = CATAGORIES_Zvv   
elif 'CR_QCD4b' in CAT0:
    CATAGORIES = CATAGORIES_CR_QCD4b 


# %%
for CAT, CAT_original in CATAGORIES.items():
    sOpDir = '%s/%s' % (sOpDir0, CAT)
    if not os.path.exists(sOpDir):   os.makedirs(sOpDir)

    # Add 'Data' according to category and later skip 'DataJetHT' and 'DataMET' from processes_dict
    if 'Zvv' in CAT: processes_dict['Data'] = processes_dict['DataMET']
    else:            processes_dict['Data'] = processes_dict['DataJetHT']
  

    sHistoNames_NotRead   = []
    sProcessNames_NotRead = {}
    for processNameToUse, processNameList in processes_dict.items():
        if processNameToUse in ['DataJetHT', 'DataMET']: continue
        
        PROC = processNameToUse
        YEAR = Era
        
        sOpFile = '%s/%s_%s_%s.root' % (sOpDir, CAT,PROC,YEAR)
        fOpFile = uproot.recreate(sOpFile)

        #for histo_name in histograms_list:
        for histo_name, histo_name_toSave in histograms_dict.items():
            for selectionWP, selectionRegions_dict in selectionTags_dict.items():
                # Skip WPs not listed in WPs_perCategory for current CAT0
                if selectionWP not in WPs_perCategory[CAT0]: continue

                for selectionRegion, selectionRegionNameOriginal0 in selectionRegions_dict.items():
                    #selectionRegionNameOriginal = selectionRegionNameOriginal0
                    #if ('Zvv' in CAT) or ('gg0l' in CAT):
                    #    selectionRegionNameOriginal = '%s_%s' %(CAT, selectionRegionNameOriginal0)
                    selectionRegionNameOriginal = '%s_%s' %(CAT_original, selectionRegionNameOriginal0)
                    if 'NoWP' in selectionWP:
                        selectionRegionNameOriginal = CAT_original
                        
                    for systematic, systematic_toSave in systematics_dict.items():
                        # Read systematics relavant to a given process
                        if systematic not in systematics_perProcess[processNameToUse]: continue
                    
                        hAdded = None
                        for processName in processNameList:
                            systematicNameToUse = systematic 
                            if 'Run' in processName or 'ata' in processNameToUse: # for Data
                                systematics_forData_original = list(systematics_forData_dict.keys())[0]
                                systematic_toSave            = systematics_forData_dict[systematics_forData_original]
                                systematicNameToUse          = systematics_forData_original
                            histo_name_toUse_full = 'evt/%s/%s_%s_%s' % (processName, histo_name, selectionRegionNameOriginal, systematicNameToUse)
                            #print(f"{histo_name_toUse_full = }")
                            try:
                                h = fIpFile[histo_name_toUse_full].to_hist()
                            except:
                                # histogram could not read
                                print(f"{histo_name_toUse_full = } could not read")
                                sHistoNames_NotRead.append(histo_name_toUse_full)
                                if processNameToUse not in sProcessNames_NotRead.keys():
                                    sProcessNames_NotRead[processNameToUse] = []
                                if processName not in sProcessNames_NotRead[processNameToUse]:
                                    sProcessNames_NotRead[processNameToUse].append( processName )

                                continue
                            #print(f"Before {h.axes = }, {len(h.axes) = }, {h.axes[0] =  }") 
                            if   len(h.axes) == 1: h = h[::hist.rebin(nRebinsX)]
                            elif len(h.axes) == 2: h = h[::hist.rebin(nRebinsX), ::hist.rebin(nRebinsY)]
                            else: 
                                print(f"hRebin for {len(h.axes) = } not implemented")
                            #print(f"After {h.axes = }, {h.axes[0] =  }")
                            if hAdded == None: hAdded = h
                            else:              hAdded = hAdded + h

                        if not hAdded: continue

                        ## Set bins with nEvents < 0 to nEvents = 0
                        nEvts          = hAdded.values()
                        varNEvts       = hAdded.variances()                
                        nEvts_modified = np.where(
                            (nEvts < 0),
                            np.full_like(nEvts, 1e-6),
                            nEvts
                        )
                        #print(f"Before: {hAdded = }")
                        #print(f"Before: {hAdded[:] = }")
                        #print(f"Before: {hAdded.values() = }")
                        #print(f"Before: {hAdded.variances() = }")
                        if   len(hAdded.axes) == 1: hAdded[:] = np.stack((nEvts_modified, varNEvts), axis=-1)
                        elif len(hAdded.axes) == 2: hAdded[:, :] = np.stack((nEvts_modified, varNEvts), axis=-1)
                        #hAdded[:, :] = np.stack((nEvts_modified, varNEvts), axis=-1)
                        #print(f"After: {hAdded = }")
                        #print(f"After: {hAdded[:] = }")
                        #print(f"After: {hAdded.values() = }")
                        #print(f"After: {hAdded.variances() = }")
        
                        '''
                        #histoNameToSave = '%s/%s_%s_%s' % (processNameToUse, histo_name,selectionTagNameToUse,systematic)
                        histoNameToSave = '%s/%s_%s' % (processNameToUse, histo_name,selectionTagNameToUse)
                        print(f"{histoNameToSave = }")
                        fOpFile[histoNameToSave] = hAdded
                        '''

                        #histoNameToSave = '%s_%s' % (histo_name,selectionTagNameToUse)
                        MASS = histo_name_toSave
                        WP = selectionWP
                        PassOrFail = selectionRegion
                        selRegion1 = '_%s_%s' % (WP, PassOrFail) if 'NoWP' not in selectionWP else ''                        
                        SYST = systematic_toSave
                        histoNameToSave = '%s_%s_%s_%s%s_%s' % (CAT, PROC, YEAR, MASS, selRegion1, SYST)
                        fOpFile[histoNameToSave] = hAdded

        fOpFile.close()
    




# %%
print(f"sHistoNames_NotRead: ")
for sHistoName_NotRead in sHistoNames_NotRead:
    print(f"\t {sHistoName_NotRead}")

# %%
#print(f"{json.dumps(sProcessNames_NotRead, indent=4) = }")
print("sProcessNames_NotRead: ")
print(json.dumps(sProcessNames_NotRead, indent=4))


