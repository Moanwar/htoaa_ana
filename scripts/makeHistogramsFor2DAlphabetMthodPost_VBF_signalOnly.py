'''
To run:
    python3 makeHistogramsFor2DAlphabetMthod.py <sAnaVersion> <Year> <CAT0>
        sAnaVersion: version name of analysis folder. E.g. 20250713_DatacardsFullSyst
        Year: 2016preVFP, 2016postVFP, 2017, 2018
        CAT0: 'gg0l', 'VBFjj', 'Vjj', 'tt0l', 'Zvv'
    e.g. time python3 makeHistogramsFor2DAlphabetMthod.py 20250713_DatacardsFullSyst 2018 gg0l 2>&1 | tee cout_makeHistogramsFor2DAlphabetMthod_20250713_DatacardsFullSyst_2018_gg0l.txt
 time python3 makeHistogramsFor2DAlphabetMthod.py 20251021_DataMC 2016preVPF VBFjj
'''

# %%
import os, sys
import psutil
import time
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
print(f"\n{sys.argv = }")

from htoaa_Settings import *
from htoaa_CommonTools import (
    rebinTH1, rebinTH2, variableRebinTH1,
)

def getAllocatedMomory():
    return "Memory: %g MB" % ((psutil.Process().memory_info().rss / (1024 * 1024))) 

global Year;
#sAnaVersion = '20250713_DatacardsFullSyst';    Year         = '2016preVFP';
#CAT0 = 'gg0l'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad'/'tt0l', 'Zvv'
sAnaVersion = sys.argv[1]
Year        = sys.argv[2]
CAT0        = sys.argv[3]
print(f"\n{sAnaVersion = }, {Year = }, {CAT0} \n")
#print(f"{psutil.Process().memory_info().rss / (1024 * 1024) = }", flush=True)
print(f"here1 {getAllocatedMomory()}", flush=True)

printLevel = 0


anaSuperCat = ''
if 'gg0l' in CAT0:  anaSuperCat = 'gg0l'
if 'VBF'  in CAT0:  anaSuperCat = 'VBFjj'
if 'Vjj'  in CAT0:  anaSuperCat = 'Vjj'
if 'Zvv'  in CAT0:  anaSuperCat = 'Zvv'
if 'tt0l' in CAT0:  anaSuperCat = 'tt0l'

Year_4Letters = Year[:4] #'2018' # '2018, 'Run2
Era = Year
# Original categories name: 'gg0l', 'VBFjj', 'Wlv', 'Zll', 'Zvv',  'Vjj'. 'ZvvIncl','ZvvLo', 'ZvvHi', 'gg0lIncl', 'gg0lLo', 'gg0lHi',
#time python3 makeHistogramsFor2DAlphabetMthod.py 20251021_DataMC 2018 VBF
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
#CAT0 = 'gg0l'  # 'gg0l', 'VBFjj', 'Vjj', 'tt0l', 'Zvv'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_ttHad_FullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_ttHad_FullSyst/2018/2DAlphabet_inputFiles'
#CAT0 = 'ttHad'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad', 'Zvv'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_Vjj_FullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_Vjj_FullSyst/2018/2DAlphabet_inputFiles'
#CAT0 = 'Vjj'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad', 'Zvv'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_Zvv_FullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250305_Zvv_FullSyst/2018/2DAlphabet_inputFiles'
#CAT0 = 'Zvv'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad', 'Zvv'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250502_gg0l_FullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250502_gg0l_FullSyst/2018/2DAlphabet_inputFiles'
#CAT0 = 'gg0l'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad', 'Zvv'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250603_gg0lDatacardsFullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250603_gg0lDatacardsFullSyst/2018/2DAlphabet_inputFiles'
#CAT0 = 'gg0l'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad', 'Zvv'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250603_VjjDatacardsFullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250603_VjjDatacardsFullSyst/2018/2DAlphabet_inputFiles'
#CAT0 = 'Vjj'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad', 'Zvv'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250603_ZvvDatacardsFullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250603_ZvvDatacardsFullSyst/2018/2DAlphabet_inputFiles'
#CAT0 = 'Zvv'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad', 'Zvv'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250603_tt0lDatacardsFullSyst/2018/analyze_htoaa_stage1.root'
#sOpDir0 = '/eos/cms/store/user/ssawant/htoaa/analysis/20250603_tt0lDatacardsFullSyst/2018/2DAlphabet_inputFiles'
#CAT0 = 'tt0l'  # 'gg0l', 'VBFjj', 'Vjj', 'ttHad'/'tt0l', 'Zvv'
#sIpFile = '/eos/user/m/moanwar/htoaa/%s/%s/%s/analyze_htoaa_stage1.root' % (sAnaVersion,Year, anaSuperCat)
sIpFile = f'/eos/user/m/moanwar/htoaa/{sAnaVersion}/{Year}/{anaSuperCat}/analyze_htoaa_stage1.root'

sOpDir0 = '/eos/user/m/moanwar/htoaa/analysis/VBF_channel/2DAlphabetfiles_VBF_sys_13Oct/%s/%s/%s/2DAlphabet_inputFiles' % (sAnaVersion, Year, anaSuperCat)



print(f"{sIpFile = }")
if 'Zvv'       in CAT0:
    ExpDatasetNames = ['MET']
else:
    ExpDatasetNames = ['JetHT']
    if Year != '2018':
        ExpDatasetNames.append( 'BTagCSV' )
ExpData_dict = {
    'Data': ['%s_Run%s%s' % (ExpDatasetName, Year[:4],EraInYear) for EraInYear in YearsAndEras_dict[Year] for ExpDatasetName in ExpDatasetNames]
}
print(f"{ExpData_dict = }", flush=True)

#exit(0)


CATAGORIES_gg0l = {
    #"gg0lIncl" : "gg0lIncl_Xto4bv2",
    "gg0lHi" :   "gg0lHi_Xto4bv2",
    "gg0lLo" :   "gg0lLo_Xto4bv2",    
}
#CATAGORIES_VBFjj = {
#    "VBFTight" :   "VBFTight_Xto4bv2",
#    "VBFLoose" :   "VBFLoose_Xto4bv2",    
#}
#"VBFIncl" : "VBFIncl_Xto4bv2",

CATAGORIES_VBFjj = {
    #"VBFHi"   : "VBFHi_Xto4bv2",
    #"VBFLo"   : "VBFLo_Xto4bv2",
    "VBFLoPTLo" : "VBFLoPTLo_Xto4bv2",
    "VBFLoPTHi" : "VBFLoPTHi_Xto4bv2",
    "VBFHiPTLo"	: "VBFHiPTLo_Xto4bv2",
    "VBFHiPTHi" : "VBFHiPTHi_Xto4bv2",
}

CATAGORIES_Vjj = {
    "VjjIncl" : "VjjIncl_Xto4bv2", 
    "VjjHi"   : "VjjHi_Xto4bv2", 
    "VjjLo"   : "VjjLo_Xto4bv2",       
}
CATAGORIES_Zvv = {
    "ZvvIncl" : "ZvvIncl_Xto4bv2",
    "ZvvHi" :   "ZvvHi_Xto4bv2",
    "ZvvLo" :   "ZvvLo_Xto4bv2",    
}
CATAGORIES_tt0l = {
    #"tt0l" : "tt0l_1TFJ_ge0BOutsideSelFJ_Xto4bv2",   
    "tt0lIncl" : "tt0l_1TFJ_ge0BOutsideSelFJ_Xto4bv2",   
    "tt0l0b" : "tt0l_1TFJ_0BOutsideSelFJ_Xto4bv2",   
    "tt0l1b" : "tt0l_1TFJ_ge1BOutsideSelFJ_Xto4bv2",   
}

WPs_perCategory = {
    'gg0l':  ['WP40', 'WP60', ], # 'WP60',
    'VBFjj': ['WP40','WP60'],
    'Vjj':   ['WP60'],
    'tt0l':  ['WP60'],
    'Zvv':   ['WP60'],        
}

if printLevel >= 6:
    print(f"here2 {getAllocatedMomory()}", flush=True)

fIpFile = uproot.open(sIpFile)
list_fIpFile_keys = list(fIpFile.keys())
if printLevel >= 6:
    print(f"here3 ipfile {getAllocatedMomory()}", flush=True)


# %%

processes_dict = {
    #'DataJetHT': ['JetHT_Run2018A', 'JetHT_Run2018B', 'JetHT_Run2018C', 'JetHT_Run2018D'],
    #'DataMET':   ['MET_Run2018A', 'MET_Run2018B', 'MET_Run2018C', 'MET_Run2018D'],
    #'Data': ExpData_dict['Data'],

    "ggHtoaato4b_mA_12p0": ["ggHtoaato4b_mA_12p0"],
    "ggHtoaato4b_mA_15p0": ["ggHtoaato4b_mA_15p0"],
    "ggHtoaato4b_mA_20p0": ["ggHtoaato4b_mA_20p0"],
    "ggHtoaato4b_mA_25p0": ["ggHtoaato4b_mA_25p0"],
    "ggHtoaato4b_mA_30p0": ["ggHtoaato4b_mA_30p0"],
    "ggHtoaato4b_mA_35p0": ["ggHtoaato4b_mA_35p0"],
    "ggHtoaato4b_mA_40p0": ["ggHtoaato4b_mA_40p0"],
    "ggHtoaato4b_mA_45p0": ["ggHtoaato4b_mA_45p0"],
    "ggHtoaato4b_mA_50p0": ["ggHtoaato4b_mA_50p0"],
    "ggHtoaato4b_mA_55p0": ["ggHtoaato4b_mA_55p0"],
    "ggHtoaato4b_mA_60p0": ["ggHtoaato4b_mA_60p0"],
    "ggHtoaato4b_mA_11p0": ["ggHtoaato4b_mA_11p0"],
    "ggHtoaato4b_mA_11p5": ["ggHtoaato4b_mA_11p5"],
    "ggHtoaato4b_mA_12p5": ["ggHtoaato4b_mA_12p5"],
    "ggHtoaato4b_mA_13p0": ["ggHtoaato4b_mA_13p0"],
    "ggHtoaato4b_mA_13p5": ["ggHtoaato4b_mA_13p5"],
    "ggHtoaato4b_mA_14p0": ["ggHtoaato4b_mA_14p0"],
    "ggHtoaato4b_mA_16p0": ["ggHtoaato4b_mA_16p0"],
    "ggHtoaato4b_mA_17p0": ["ggHtoaato4b_mA_17p0"],
    "ggHtoaato4b_mA_18p5": ["ggHtoaato4b_mA_18p5"],
    "ggHtoaato4b_mA_21p5": ["ggHtoaato4b_mA_21p5"],
    "ggHtoaato4b_mA_23p0": ["ggHtoaato4b_mA_23p0"],
    "ggHtoaato4b_mA_27p5": ["ggHtoaato4b_mA_27p5"],
    "ggHtoaato4b_mA_32p5": ["ggHtoaato4b_mA_32p5"],
    "ggHtoaato4b_mA_37p5": ["ggHtoaato4b_mA_37p5"],
    "ggHtoaato4b_mA_42p5": ["ggHtoaato4b_mA_42p5"],
    "ggHtoaato4b_mA_47p5": ["ggHtoaato4b_mA_47p5"],
    "ggHtoaato4b_mA_52p5": ["ggHtoaato4b_mA_52p5"],
    "ggHtoaato4b_mA_57p5": ["ggHtoaato4b_mA_57p5"],
    "ggHtoaato4b_mA_62p5": ["ggHtoaato4b_mA_62p5"],

    "VBFHtoaato4b_mA_12p0": ["VBFHtoaato4b_mA_12p0"],
    "VBFHtoaato4b_mA_15p0": ["VBFHtoaato4b_mA_15p0"],
    "VBFHtoaato4b_mA_20p0": ["VBFHtoaato4b_mA_20p0"],
    "VBFHtoaato4b_mA_25p0": ["VBFHtoaato4b_mA_25p0"],
    "VBFHtoaato4b_mA_30p0": ["VBFHtoaato4b_mA_30p0"],
    "VBFHtoaato4b_mA_35p0": ["VBFHtoaato4b_mA_35p0"],
    "VBFHtoaato4b_mA_40p0": ["VBFHtoaato4b_mA_40p0"],
    "VBFHtoaato4b_mA_45p0": ["VBFHtoaato4b_mA_45p0"],
    "VBFHtoaato4b_mA_50p0": ["VBFHtoaato4b_mA_50p0"],
    "VBFHtoaato4b_mA_55p0": ["VBFHtoaato4b_mA_55p0"],
    "VBFHtoaato4b_mA_60p0": ["VBFHtoaato4b_mA_60p0"],
    "VBFHtoaato4b_mA_11p0": ["VBFHtoaato4b_mA_11p0"],
    "VBFHtoaato4b_mA_11p5": ["VBFHtoaato4b_mA_11p5"],
    "VBFHtoaato4b_mA_12p5": ["VBFHtoaato4b_mA_12p5"],
    "VBFHtoaato4b_mA_13p0": ["VBFHtoaato4b_mA_13p0"],
    "VBFHtoaato4b_mA_13p5": ["VBFHtoaato4b_mA_13p5"],
    "VBFHtoaato4b_mA_14p0": ["VBFHtoaato4b_mA_14p0"],
    "VBFHtoaato4b_mA_16p0": ["VBFHtoaato4b_mA_16p0"],
    "VBFHtoaato4b_mA_17p0": ["VBFHtoaato4b_mA_17p0"],
    "VBFHtoaato4b_mA_18p5": ["VBFHtoaato4b_mA_18p5"],
    "VBFHtoaato4b_mA_21p5": ["VBFHtoaato4b_mA_21p5"],
    "VBFHtoaato4b_mA_23p0": ["VBFHtoaato4b_mA_23p0"],
    "VBFHtoaato4b_mA_27p5": ["VBFHtoaato4b_mA_27p5"],
    "VBFHtoaato4b_mA_32p5": ["VBFHtoaato4b_mA_32p5"],
    "VBFHtoaato4b_mA_37p5": ["VBFHtoaato4b_mA_37p5"],
    "VBFHtoaato4b_mA_42p5": ["VBFHtoaato4b_mA_42p5"],
    "VBFHtoaato4b_mA_47p5": ["VBFHtoaato4b_mA_47p5"],
    "VBFHtoaato4b_mA_52p5": ["VBFHtoaato4b_mA_52p5"],
    "VBFHtoaato4b_mA_57p5": ["VBFHtoaato4b_mA_57p5"],
    "VBFHtoaato4b_mA_62p5": ["VBFHtoaato4b_mA_62p5"],

   "WHtoaato4b_mA_12p0": ["WHtoaato4b_mA_12p0"],
    "WHtoaato4b_mA_15p0": ["WHtoaato4b_mA_15p0"],
    "WHtoaato4b_mA_20p0": ["WHtoaato4b_mA_20p0"],
    "WHtoaato4b_mA_25p0": ["WHtoaato4b_mA_25p0"],
    "WHtoaato4b_mA_30p0": ["WHtoaato4b_mA_30p0"],
    "WHtoaato4b_mA_35p0": ["WHtoaato4b_mA_35p0"],
    "WHtoaato4b_mA_40p0": ["WHtoaato4b_mA_40p0"],
    "WHtoaato4b_mA_45p0": ["WHtoaato4b_mA_45p0"],
    "WHtoaato4b_mA_50p0": ["WHtoaato4b_mA_50p0"],
    "WHtoaato4b_mA_55p0": ["WHtoaato4b_mA_55p0"],
    "WHtoaato4b_mA_60p0": ["WHtoaato4b_mA_60p0"],
    "WHtoaato4b_mA_11p0": ["WHtoaato4b_mA_11p0"],
    "WHtoaato4b_mA_11p5": ["WHtoaato4b_mA_11p5"],
    "WHtoaato4b_mA_12p5": ["WHtoaato4b_mA_12p5"],
    "WHtoaato4b_mA_13p0": ["WHtoaato4b_mA_13p0"],
    "WHtoaato4b_mA_13p5": ["WHtoaato4b_mA_13p5"],
    "WHtoaato4b_mA_14p0": ["WHtoaato4b_mA_14p0"],
    "WHtoaato4b_mA_16p0": ["WHtoaato4b_mA_16p0"],
    "WHtoaato4b_mA_17p0": ["WHtoaato4b_mA_17p0"],
    "WHtoaato4b_mA_18p5": ["WHtoaato4b_mA_18p5"],
    "WHtoaato4b_mA_21p5": ["WHtoaato4b_mA_21p5"],
    "WHtoaato4b_mA_23p0": ["WHtoaato4b_mA_23p0"],
    "WHtoaato4b_mA_27p5": ["WHtoaato4b_mA_27p5"],
    "WHtoaato4b_mA_32p5": ["WHtoaato4b_mA_32p5"],
    "WHtoaato4b_mA_37p5": ["WHtoaato4b_mA_37p5"],
    "WHtoaato4b_mA_42p5": ["WHtoaato4b_mA_42p5"],
    "WHtoaato4b_mA_47p5": ["WHtoaato4b_mA_47p5"],
    "WHtoaato4b_mA_52p5": ["WHtoaato4b_mA_52p5"],
    "WHtoaato4b_mA_57p5": ["WHtoaato4b_mA_57p5"],
    "WHtoaato4b_mA_62p5": ["WHtoaato4b_mA_62p5"],

    "ZHtoaato4b_mA_12p0": ["ZHtoaato4b_mA_12p0"],
    "ZHtoaato4b_mA_15p0": ["ZHtoaato4b_mA_15p0"],
    "ZHtoaato4b_mA_20p0": ["ZHtoaato4b_mA_20p0"],
    "ZHtoaato4b_mA_25p0": ["ZHtoaato4b_mA_25p0"],
    "ZHtoaato4b_mA_30p0": ["ZHtoaato4b_mA_30p0"],
    "ZHtoaato4b_mA_35p0": ["ZHtoaato4b_mA_35p0"],
    "ZHtoaato4b_mA_40p0": ["ZHtoaato4b_mA_40p0"],
    "ZHtoaato4b_mA_45p0": ["ZHtoaato4b_mA_45p0"],
    "ZHtoaato4b_mA_50p0": ["ZHtoaato4b_mA_50p0"],
    "ZHtoaato4b_mA_55p0": ["ZHtoaato4b_mA_55p0"],
    "ZHtoaato4b_mA_60p0": ["ZHtoaato4b_mA_60p0"],
    "ZHtoaato4b_mA_11p0": ["ZHtoaato4b_mA_11p0"],
    "ZHtoaato4b_mA_11p5": ["ZHtoaato4b_mA_11p5"],
    "ZHtoaato4b_mA_12p5": ["ZHtoaato4b_mA_12p5"],
    "ZHtoaato4b_mA_13p0": ["ZHtoaato4b_mA_13p0"],
    "ZHtoaato4b_mA_13p5": ["ZHtoaato4b_mA_13p5"],
    "ZHtoaato4b_mA_14p0": ["ZHtoaato4b_mA_14p0"],
    "ZHtoaato4b_mA_16p0": ["ZHtoaato4b_mA_16p0"],
    "ZHtoaato4b_mA_17p0": ["ZHtoaato4b_mA_17p0"],
    "ZHtoaato4b_mA_18p5": ["ZHtoaato4b_mA_18p5"],
    "ZHtoaato4b_mA_21p5": ["ZHtoaato4b_mA_21p5"],
    "ZHtoaato4b_mA_23p0": ["ZHtoaato4b_mA_23p0"],
    "ZHtoaato4b_mA_27p5": ["ZHtoaato4b_mA_27p5"],
    "ZHtoaato4b_mA_32p5": ["ZHtoaato4b_mA_32p5"],
    "ZHtoaato4b_mA_37p5": ["ZHtoaato4b_mA_37p5"],
    "ZHtoaato4b_mA_42p5": ["ZHtoaato4b_mA_42p5"],
    "ZHtoaato4b_mA_47p5": ["ZHtoaato4b_mA_47p5"],
    "ZHtoaato4b_mA_52p5": ["ZHtoaato4b_mA_52p5"],
    "ZHtoaato4b_mA_57p5": ["ZHtoaato4b_mA_57p5"],
    "ZHtoaato4b_mA_62p5": ["ZHtoaato4b_mA_62p5"],

    "ttHtoaato4b_mA_12p0": ["ttHtoaato4b_mA_12p0"],
    "ttHtoaato4b_mA_15p0": ["ttHtoaato4b_mA_15p0"],
    "ttHtoaato4b_mA_20p0": ["ttHtoaato4b_mA_20p0"],
    "ttHtoaato4b_mA_25p0": ["ttHtoaato4b_mA_25p0"],
    "ttHtoaato4b_mA_30p0": ["ttHtoaato4b_mA_30p0"],
    "ttHtoaato4b_mA_35p0": ["ttHtoaato4b_mA_35p0"],
    "ttHtoaato4b_mA_40p0": ["ttHtoaato4b_mA_40p0"],
    "ttHtoaato4b_mA_45p0": ["ttHtoaato4b_mA_45p0"],
    "ttHtoaato4b_mA_50p0": ["ttHtoaato4b_mA_50p0"],
    "ttHtoaato4b_mA_55p0": ["ttHtoaato4b_mA_55p0"],
    "ttHtoaato4b_mA_60p0": ["ttHtoaato4b_mA_60p0"],
    "ttHtoaato4b_mA_11p0": ["ttHtoaato4b_mA_11p0"],
    "ttHtoaato4b_mA_11p5": ["ttHtoaato4b_mA_11p5"],
    "ttHtoaato4b_mA_12p5": ["ttHtoaato4b_mA_12p5"],
    "ttHtoaato4b_mA_13p0": ["ttHtoaato4b_mA_13p0"],
    "ttHtoaato4b_mA_13p5": ["ttHtoaato4b_mA_13p5"],
    "ttHtoaato4b_mA_14p0": ["ttHtoaato4b_mA_14p0"],
    "ttHtoaato4b_mA_16p0": ["ttHtoaato4b_mA_16p0"],
    "ttHtoaato4b_mA_17p0": ["ttHtoaato4b_mA_17p0"],
    "ttHtoaato4b_mA_18p5": ["ttHtoaato4b_mA_18p5"],
    "ttHtoaato4b_mA_21p5": ["ttHtoaato4b_mA_21p5"],
    "ttHtoaato4b_mA_23p0": ["ttHtoaato4b_mA_23p0"],
    "ttHtoaato4b_mA_27p5": ["ttHtoaato4b_mA_27p5"],
    "ttHtoaato4b_mA_32p5": ["ttHtoaato4b_mA_32p5"],
    "ttHtoaato4b_mA_37p5": ["ttHtoaato4b_mA_37p5"],
    "ttHtoaato4b_mA_42p5": ["ttHtoaato4b_mA_42p5"],
    "ttHtoaato4b_mA_47p5": ["ttHtoaato4b_mA_47p5"],
    "ttHtoaato4b_mA_52p5": ["ttHtoaato4b_mA_52p5"],
    "ttHtoaato4b_mA_57p5": ["ttHtoaato4b_mA_57p5"],
    "ttHtoaato4b_mA_62p5": ["ttHtoaato4b_mA_62p5"],

    
}
systematics_forData_dict = {'noweight': 'Nom'}
systematics_dict = {'Nom': 'Nom',}
for systNameShort, systName0  in SystNameConvs.items(): 
    #systName_ = systName0.replace('$YEAR', Year_4Letters)
    systName_ = systName0.replace('$YEAR', Era)
    if systNameShort == 'BtagCorr':
        systName_ = systName0
    systematics_dict[systNameShort+SystNameConvUp  ] = systName_+SystNameConvUp
    systematics_dict[systNameShort+SystNameConvDown] = systName_+SystNameConvDown    
print(f"{systematics_dict = }")


systNameShort_MCNom = ['Nom']
systNameShort_MCAll = [
    'PU', 'AK8JetJES', 'AK8JetJER', 'AK4JetJES', 'AK4JetJER', 'METUnclE'
    'ISR', 'FSR', 'QCDScale', 'PDF', 'massScaleH','massResolH','massScaleA','massResolA'
    ]
if 'Zvv'       in CAT0:
    systNameShort_MCAll.extend( ['MetTrigEffi',  'METUnclE' ])
else:
    systNameShort_MCAll.extend( ['JetTrigEffi'])
if Year != '2018':
    #systNameShort_MCAll.extend( ['2018HEM1516Issue'])
    systNameShort_MCAll.extend( ['L1Prefire'])
if kDatasetToAnalyze == DatasetToAnalyze.SingleYear:
    systNameShort_MCAll.extend( ['Btag'])
else:
    systNameShort_MCAll.extend( ['BtagCorr', 'BtagUncorr'])

systNameShort_MCTT = ['TopPtReWeight']

systNameShort_MCSignalH    = ['LPRewgt']
systNameShort_MCSignalGGH  = ['ggHPtRewgt']
systNameShort_MCSignalVBFH = ['VBFHPtRewgt']
systNameShort_MCSignalWH   = ['WHPtRewgt']
systNameShort_MCSignalZH   = ['ZHPtRewgt']
systNameShort_MCSignalTTH  = ['ttHPtRewgt']


systematics_perProcess = {
    #'Data': systNameShort_MCNom,

    "ggHtoaato4b_mA_12p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_15p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_20p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_25p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_30p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_35p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_40p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_45p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_50p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_55p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_60p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_11p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_11p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_12p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_13p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_13p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_14p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_16p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_17p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_18p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_21p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_23p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_27p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_32p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_37p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_42p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_47p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_52p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_57p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,
    "ggHtoaato4b_mA_62p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalGGH,

    "VBFHtoaato4b_mA_12p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_15p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_20p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_25p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_30p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_35p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_40p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_45p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_50p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_55p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_60p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_11p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_11p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_12p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_13p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_13p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_14p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_16p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_17p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_18p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_21p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_23p0":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_27p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_32p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_37p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_42p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_47p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_52p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_57p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,
    "VBFHtoaato4b_mA_62p5":  systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalVBFH,

    "WHtoaato4b_mA_12p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_15p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_20p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_25p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_30p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_35p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_40p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_45p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_50p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_55p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_60p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_11p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_11p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_12p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_13p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_13p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_14p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_16p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_17p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_18p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_21p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_23p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_27p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_32p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_37p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_42p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_47p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_52p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_57p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    "WHtoaato4b_mA_62p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalWH,
    
    "ZHtoaato4b_mA_12p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_15p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_20p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_25p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_30p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_35p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_40p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_45p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_50p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_55p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_60p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_11p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_11p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_12p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_13p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_13p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_14p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_16p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_17p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_18p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_21p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_23p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_27p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_32p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_37p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_42p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_47p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_52p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_57p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ZHtoaato4b_mA_62p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalZH,
    "ttHtoaato4b_mA_12p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_15p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_20p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_25p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_30p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_35p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_40p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_45p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_50p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_55p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_60p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_11p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_11p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_12p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_13p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_13p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_14p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_16p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_17p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_18p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_21p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_23p0": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_27p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_32p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_37p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_42p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_47p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_52p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_57p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,
    "ttHtoaato4b_mA_62p5": systNameShort_MCNom + systNameShort_MCAll + systNameShort_MCSignalH + systNameShort_MCSignalTTH,

}
print(f"{systematics_perProcess = }")



#histograms_list = [
#    'hLeadingFatJetMass_vs_massA_Hto4b_avg',
#    'hLeadingFatJetMSoftDrop_vs_massA_Hto4b_avg',
#    'hLeadingFatJetParticleNet_massH_Hto4b_avg_vs_massA_Hto4b_avg',
#]
'''
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
'''
histograms_dict = {
    #'hLeadingFatJetMass_vs_massAa':                        'mass_4a',
    #'hLeadingFatJetMSoftDrop_vs_massAa':                   'msoft_4a',
    'hLeadingFatJetPNet_massH_v2b_vs_massAa':              'pnet_4a', # 'pnet_vs_massAa',   
    #'hLeadingFatJetMass_vs_massA34a':                        'mass_34a',
    #'hLeadingFatJetMSoftDrop_vs_massA34a':                   'msoft_34a',
    'hLeadingFatJetPNet_massH_v2b_vs_massA34a':              'pnet_34a',  ##  'pnet_vs_massA34a'
    #'hLeadingFatJetMass_vs_massA34b':                        'mass_vs_massA34b',
    #'hLeadingFatJetMass_vs_massA34d':                        'mass_34d',
    #'hLeadingFatJetMSoftDrop_vs_massA34d':                   'msoft_34d',
    #'hLeadingFatJetPNet_massH_v2b_vs_massA34d':              'pnet_34d',  ##  'pnet_vs_massA34d' 
    
}


nRebinsX = 1 # 10
nRebinsY = 1 #  4

'''
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
'''  
selectionTags_dict = {
    'WP40': {
        'Pass': 'SRWP40',
        'Fail': 'SBWP40'
    },
    #'WP45a': {
    #    'Pass': 'SRWP45a',
    #    'Fail': 'SBWP45a'
    #},
    #'WP45b': {
    #    'Pass': 'SRWP45b',
    #    'Fail': 'SBWP45b'
    #},
    #'WP50': {
    #    'Pass': 'SRWP50',
    #    'Fail': 'SBWP50'
    #},
    'WP60': {
        'Pass': 'SRWP60',
        'Fail': 'SBWP60'
    },
    #'WP65': {
    #    'Pass': 'SRWP65',
    #    'Fail': 'SBWP65'
    #},
    #'WP70': {
    #    'Pass': 'SRWP70',
    #    'Fail': 'SBWP70'
    #},
    'WP80': {
        'Pass': 'SRWP80',
        'Fail': 'SBWP80'
    },
            
}  


if 'gg0l' in CAT0:
    CATAGORIES = CATAGORIES_gg0l
elif 'VBF' in CAT0:
    CATAGORIES = CATAGORIES_VBFjj
elif 'tt0l' in CAT0:
    CATAGORIES = CATAGORIES_tt0l
elif 'Vjj' in CAT0:
    CATAGORIES = CATAGORIES_Vjj
elif 'Zvv' in CAT0:
    CATAGORIES = CATAGORIES_Zvv   

if printLevel >= 6:
    print(f"here4 setting all dict {getAllocatedMomory()}", flush=True)
    #print(f"\n\n\n {list_fIpFile_keys = } \n\n\n")
    lastTimeStamp = time.time()


# %%
for CAT, CAT_original in CATAGORIES.items():
    
    sOpDir = '%s/%s' % (sOpDir0, CAT)
    if not os.path.exists(sOpDir):   os.makedirs(sOpDir)

    # Add 'Data' according to category and later skip 'DataJetHT' and 'DataMET' from processes_dict
    #if 'Zvv' in CAT: processes_dict['Data'] = processes_dict['DataMET']
    #else:            processes_dict['Data'] = processes_dict['DataJetHT']
  

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
                    selectionRegionNameOriginal = selectionRegionNameOriginal0
                    #if ('Zvv' in CAT) or ('gg0l' in CAT):
                    #    selectionRegionNameOriginal = '%s_%s' %(CAT, selectionRegionNameOriginal0)
                    selectionRegionNameOriginal = '%s_%s' %(CAT_original, selectionRegionNameOriginal0)
                        
                    for systematic, systematic_toSave in systematics_dict.items():
                        systematic_woUpDown = systematic.replace(SystNameConvUp, '')
                        systematic_woUpDown = systematic_woUpDown.replace(SystNameConvDown, '')
                        # Read systematics relavant to a given process
                        if systematic_woUpDown not in systematics_perProcess[processNameToUse]: continue

                        if printLevel >= 6:
                            lapTime       = time.time() - lastTimeStamp
                            lastTimeStamp = time.time()
                            print(f"here10 {getAllocatedMomory()}, {lapTime} sec, {CAT}, {processNameToUse}, {histo_name}, {selectionWP}, {selectionRegion}, {systematic}", flush=True)

                    
                        hAdded = None
                        for processName in processNameList:
                            systematicNameToUse = systematic_toSave # systematic 
                            if 'Run' in processName or 'ata' in processNameToUse: # for Data
                                systematics_forData_original = list(systematics_forData_dict.keys())[0]
                                systematic_toSave            = systematics_forData_dict[systematics_forData_original]
                                systematicNameToUse          = systematics_forData_original
                            histo_name_toUse_full = 'evt/%s/%s_%s_%s' % (processName, histo_name, selectionRegionNameOriginal, systematicNameToUse)
                            #print(f"{histo_name_toUse_full = }")
                            '''
                            # try except 
                            print(f"\t\t\t\t {histo_name_toUse_full}:  {(histo_name_toUse_full in list_fIpFile_keys) = }, \t {('%s;1'%(histo_name_toUse_full) in list_fIpFile_keys) = }")
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
                            '''
                            if ((histo_name_toUse_full in list_fIpFile_keys) or ('%s;1'%(histo_name_toUse_full) in list_fIpFile_keys)):
                                h = fIpFile[histo_name_toUse_full].to_hist()
                            else:
                                # histogram could not read
                                print(f"{histo_name_toUse_full = } could not read")
                                sHistoNames_NotRead.append(histo_name_toUse_full)
                                if processNameToUse not in sProcessNames_NotRead.keys():
                                    sProcessNames_NotRead[processNameToUse] = []
                                if processName not in sProcessNames_NotRead[processNameToUse]:
                                    sProcessNames_NotRead[processNameToUse].append( processName )
                                continue

                            if ((nRebinsX != 1) or (nRebinsY != 1)):
                                h = h[::hist.rebin(nRebinsX), ::hist.rebin(nRebinsY)]
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
                        hAdded[:, :] = np.stack((nEvts_modified, varNEvts), axis=-1)
        
                        '''
                        #histoNameToSave = '%s/%s_%s_%s' % (processNameToUse, histo_name,selectionTagNameToUse,systematic)
                        histoNameToSave = '%s/%s_%s' % (processNameToUse, histo_name,selectionTagNameToUse)
                        print(f"{histoNameToSave = }")
                        fOpFile[histoNameToSave] = hAdded
                        '''

                        #histoNameToSave = '%s_%s' % (histo_name,selectionTagNameToUse)
                        MASS       = histo_name_toSave
                        WP         = selectionWP
                        PassOrFail = selectionRegion
                        SYST       = systematic_toSave
                        histoNameToSave = '%s_%s_%s_%s_%s_%s_%s' % (CAT, PROC, YEAR, MASS, WP, PassOrFail, SYST)
                        fOpFile[histoNameToSave] = hAdded

        fOpFile.close()
        del fOpFile
    
if printLevel >= 6:
    print(f"here100 Saved histograms {getAllocatedMomory()}", flush=True)






# %%
print(f"sHistoNames_NotRead: ")
for sHistoName_NotRead in sHistoNames_NotRead:
    print(f"\t {sHistoName_NotRead}")

# %%
#print(f"{json.dumps(sProcessNames_NotRead, indent=4) = }")
print("sProcessNames_NotRead: ")
print(json.dumps(sProcessNames_NotRead, indent=4))


if printLevel >= 6:
    print(f"here1000 End {getAllocatedMomory()}", flush=True)
