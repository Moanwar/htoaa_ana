'''
To run:
    python3 getEventYield_DataVsMC.py <sAnaDir> <Dataset> <CAT>
        sAnaDir: full path of analysis directory where output histograms are stored. E.g. /eos/cms/store/user/ssawant/htoaa/analysis/20250713_DatacardsFullSyst
        Datasets: comman separated dataset strings . Options: '2016preVFP', '2016postVFP', '2017', '2018', 'Run2', 'All'
        CAT0: 'gg0l', 'VBFjj', 'Wlv', 'Zll', 'Zvv',  'Vjj'. 'ZvvIncl','ZvvLo', 'ZvvHi', 'gg0lIncl', 'gg0lLo', 'gg0lHi', VjjLo, VjjHi, VjjIncl, 'tt0l', 'tt0l_1TFJ_ge0BOutsideSelFJ', 'CR_QCD4b'
            'tt0l_ge1NonHFatJet_0BExtra', 'tt0l_ge1NonHFatJet_1BExtra', 'tt0l_ge1NonHFatJet_ge2BExtra', 'tt0l_0NonHFatJet_ge2B'
            tt0l_1TFJ_0BOutsideSelFJ, tt0l_1TFJ_ge1BOutsideSelFJ, tt0l_1TFJ_ge0BOutsideSelFJ
            'trigEffi', 'CR_QCD4b'
    e.g. time python3 getEventYield_DataVsMC.py /eos/cms/store/user/ssawant/htoaa/analysis/20250713_DatacardsFullSyst 2018 gg0l 
'''

import os, sys
import numpy as np
from collections import OrderedDict as OD
#import uproot3 as uproot
#import uproot3 
import uproot as uproot
import hist
import matplotlib.pyplot as plt
import mplhep as hep
import json

sys.path.insert(1, '../') # to import file from other directory (../ in this case)

from htoaa_Settings import *
from htoaa_CommonTools import(
    getTH1BinContent
)

sAnaDir     = sys.argv[1]
Datasets     = sys.argv[2]
CAT         = sys.argv[3]

#global Year;
#sAnaVersion = '20250712_DataMCwoTrgSFCorrectVjjVeto';    Year         = '2018';

runMode = 'ExcelCompatible' #'LatexCompatible', 'ExcelCompatible' 

#CAT = 'gg0lIncl'  # 'gg0l', 'VBFjj', 'Wlv', 'Zll', 'Zvv',  'Vjj'. 'ZvvIncl','ZvvLo', 'ZvvHi', 'gg0lIncl', 'gg0lLo', 'gg0lHi', 'tt0l', 'VBFjj', 'tt0l_1TFJ_ge0BOutsideSelFJ' 
# 'tt0l_ge1NonHFatJet_0BExtra', 'tt0l_ge1NonHFatJet_1BExtra', 'tt0l_ge1NonHFatJet_ge2BExtra', 'tt0l_0NonHFatJet_ge2B'

anaSuperCat = ''
if 'gg0l'     in CAT:  anaSuperCat = 'gg0l'
if 'VBF'      in CAT:  anaSuperCat = 'VBFjj'
if 'Vjj'      in CAT:  anaSuperCat = 'Vjj'
if 'Zvv'      in CAT:  anaSuperCat = 'Zvv'
if 'tt0l'     in CAT:  anaSuperCat = 'tt0l'
if 'trigEffi' in CAT:  anaSuperCat = 'trigEffi'
if 'CR_QCD4b' in CAT:  anaSuperCat = 'CR_QCD4b'

subCats = []
if   'gg0l'     in CAT:
    subCats = ["gg0lIncl", "gg0lHi", "gg0lLo", "gg0l0bLo", "gg0l0bHi", "gg0l1bLo", "gg0l1bHi"]
elif 'VBF'      in CAT:    
    subCats = ["VBFHi", "VBFLo"]
elif 'Vjj'      in CAT:    
    subCats = ["VjjIncl", "VjjHi", "VjjLo", "VjjHi350", "VjjLo350", "VjjHi400", "VjjLo400"]
elif 'Zvv'      in CAT:    
    subCats = ["ZvvIncl", "ZvvHi", "ZvvLo"]
elif 'tt0l'      in CAT:    
    subCats = ["tt0l_1TFJ_ge0BOutsideSelFJ", "tt0l_1TFJ_0BOutsideSelFJ", "tt0l_1TFJ_ge1BOutsideSelFJ"]

## Set selection tags
selectionTags = []
for subCat_ in subCats:
    selectionTags.extend( [subCat_,] ) # '%sMsdLt50' % (CAT),'%sMsdGt50' % (CAT)]
    if 'gg0l' in CAT: selectionTags.extend([ '%s_Xto4bv2_SBplusSRWP40' % (subCat_),] )
    else:             selectionTags.extend([ '%s_Xto4bv2_SBplusSRWP60' % (subCat_), ] )
if 'trigEffi' in CAT: 
    selectionTags = ['JetTrgEffiDenom', 'JetTrgEffiNume_Trg_Combo_AK4AK8Jet_HT_VBF']
elif 'CR_QCD4b'      in CAT: 
    selectionTags = ["CR4b_3M2T", "CR4b_3M3T", "CR4b_4M3T", "CR4b_4M4T"]


# Datasets to use
Datasets_list = Datasets.split(',')
print(f"{Datasets_list = }")

# Year, Era are set internally to one of the following: '2016preVFP', '2016postVFP', '2017', '2018'
YearsAll_list = [Era_2016preVFP, Era_2016postVFP, Era_2017, Era_2018]

YearsToRun_dict = {}
Years = []
for Dataset  in Datasets_list:
    if Dataset==Era_Run2: # make Run2 data-mc plots
        YearsToRun_dict[Era_Run2] = YearsAll_list
        Years = YearsAll_list
    elif Dataset=='All': # make Run2 and individual 4 years data-mc plots
        YearsToRun_dict[Era_Run2] = YearsAll_list
        for Year_ in YearsAll_list:
            YearsToRun_dict[Year_] = [Year_]
        Years = YearsAll_list
    else: # make individual year's data-mc plots
        for Year_ in YearsAll_list:
            if Year_ != Dataset: continue
            YearsToRun_dict[Year_] = [Year_]   
            Years.append(Year_)

print(f"{YearsToRun_dict = }, \n{Years = }")



## Read input files
sIpFiles = {}
for Era in Years:
    sIpFiles[Era] = '%s/%s/%s/analyze_htoaa_stage1.root' % (sAnaDir, Era, anaSuperCat) # 20250612_gg0lDataMC_1, 20250613_gg0lDataMC_1, 20250617_gg0lDataMC, 20250617_gg0lDataMC_1
sOpDirNameShort = 'plots'



DataObs_DirName_dict = {}
luminosity_total_dict = {}
for DatasetName_, YearsToRun_list_ in YearsToRun_dict.items():
    luminosity_total_ = 0.0
    for Year_ in YearsToRun_list_:
        if 'Zvv'       in CAT:
            ExpDatasetNames = ['MET']
            HLT_toUse       = 'Trg_Combo_MET'
        elif 'trigEffi'       in CAT:
            ExpDatasetNames = ['SingleMuon']
            HLT_toUse       = 'Trg_Combo_Mu'
        else:
            ExpDatasetNames = ['JetHT']
            HLT_toUse       = 'Trg_Combo_AK4AK8Jet_HT_VBF'
            if Year_ != '2018':
                ExpDatasetNames.append( 'BTagCSV' )
        DataObs_DirName_list_i_ = ['%s_Run%s%s' % (ExpDatasetName, Year_[:4],EraInYear) for ExpDatasetName in ExpDatasetNames for EraInYear in YearsAndEras_dict[Year_]]
        if Year_ not in DataObs_DirName_dict: DataObs_DirName_dict[Year_] = DataObs_DirName_list_i_
        luminosity_total_ += Luminosities_TotalPerYear[Year_][HLT_toUse][0]
    luminosity_total_dict[DatasetName_] = luminosity_total_



RunMode = '' # '', 'test'
printLevel = 0 #


print(f"{selectionTags = }")
print(f"{sIpFiles = } ")

#print(f"{DataObs_DirName_dict = }")
print("DataObs_DirName_dict: ", json.dumps(DataObs_DirName_dict, indent=4))
print(f"{luminosity_total_dict = }")

fIpFiles = {} 
sIpFileKeys_dict = {}
for Year_ in Years:   
    fIpFiles[Year_] = uproot.open(sIpFiles[Year_])
    sIpFileKeys_dict[Year_] = list(fIpFiles[Year_].keys())

MCBkg_dict = {}
if CAT not in ['CR_QCD4b']:
    MCBkg_dict.update({
        #'QCD': ["QCD_bEnr", "QCD_BGen", "QCD_Incl"],
        "QCD_bEnr": ["QCD_bEnr"],
        "QCD_BGen": ["QCD_BGen"],
        "QCD_Incl": ["QCD_Incl"],  
    })
else: 
    MCBkg_dict.update( {
        'QCD_012b':      ['QCD_0bCat', 'QCD_1bCat', 'QCD_2bCat', ],
        'QCD_3b':        ['QCD_3bCat'],
        'QCD_4bAndMore': ['QCD_4bAndMoreCat'],
    } )
MCBkg_dict.update( {  
    #r't$\bar{t}$+X': ["TT0l", "TT1l", "TT2l"],
    'TT0l': ['TT0l'],
    'TT1l': ['TT1l'],
    'TT2l': ['TT2l'],
    'Single top': ["STop_t", "STbar_t", "ST_s_0l", "ST_s_1l", "STop_tW_Incl", "STbar_tW_Incl"], # "STop_tW_12l", "STbar_tW_12l"],
    #'V+X': ["Zqq", "Zvv", "Zll", "Wqq", "Wlv"],
    'ttZ': ['ttZ'],    
    'ttW': ['ttW'],
    'tZq': ['tZq'],
    'Zqq': ['Zqq'],
    'Zll': ['Zll'],
    'Zvv': ['Zvv'],
    'Wqq': ['Wqq'],
    'Wlv': ['Wlv'],
    'Diboson': ["ZZ", "WZ", "WW"],
    'VVV': ['ZZZ', 'WZZ', 'WWZ', 'WWW'],

    'ggH':  ['GluGluHToBB_Pt-200ToInf'], #['ggH'],
    'VBFH': ['VBFH_dipoleRecoilOn'], #['VBFHToBB_herwig'], #['VBFHToBB_powheg'],
    'WH':   ['WplusHToBBQQ', 'WplusHToBBLNu', 'WminusHToBBQQ', 'WminusHToBBLNu'],
    'ZH':   ['ZHToBBX'],
    'ttH':  ['ttHToBB'], # ['ttHToBB', 'ttHToNonBB'],

    #'': [],
} )

MCSig_dict = {
    'ggHtoaato4b_mA_15': ['ggHtoaato4b_mA_15'],
    'ggHtoaato4b_mA_30': ['ggHtoaato4b_mA_30'],
    'ggHtoaato4b_mA_55': ['ggHtoaato4b_mA_55'],

    'VBFHtoaato4b_mA_15': ['VBFHtoaato4b_mA_15'],
    'VBFHtoaato4b_mA_30': ['VBFHtoaato4b_mA_30'],
    'VBFHtoaato4b_mA_55': ['VBFHtoaato4b_mA_55'],
    #'VBFHtoaato4b_mA_60': ['VBFHtoaato4b_mA_60'],

    'WHtoaato4b_mA_15': ['WHtoaato4b_mA_15'],
    'WHtoaato4b_mA_30': ['WHtoaato4b_mA_30'],
    'WHtoaato4b_mA_55': ['WHtoaato4b_mA_55'],
    
    'ZHtoaato4b_mA_15': ['ZHtoaato4b_mA_15'],
    'ZHtoaato4b_mA_30': ['ZHtoaato4b_mA_30'],
    'ZHtoaato4b_mA_55': ['ZHtoaato4b_mA_55'],
    
    'ttHtoaato4b_mA_15': ['ttHtoaato4b_mA_15'],
    'ttHtoaato4b_mA_30': ['ttHtoaato4b_mA_30'],
    'ttHtoaato4b_mA_55': ['ttHtoaato4b_mA_55'],    
}
if CAT in ['CR_QCD4b']:
    MCSig_dict = {}

systematics_list = ['Nom'] # ['central']
systematics_forData = 'noweight'
xValueForEventYield = 4

sSpace     = ' & '     if runMode == 'LatexCompatible' else ' \t '
sPlusMinus = ' $\pm$ ' if runMode == 'LatexCompatible' else ' \t '
sEndline   = ' \\ \n ' if runMode == 'LatexCompatible' else ' \n'

sProcessNames_NotRead = {}

evtYields_dict = {} # dict levels: evtYields_dict[<dataset name>][<lumi used: Total or 1>][<selectionTag_syst>][<process>][<nEvt or errN or nEvtUnwgt>]
for sDatasetName, YearsToRun_list in YearsToRun_dict.items():
    luminosity_total         = luminosity_total_dict[sDatasetName]
    if sDatasetName not in evtYields_dict:
        evtYields_dict[sDatasetName] = {}

    for luminosity_toUse in [luminosity_total, 1.0]: # make event yield table twice 1) total luminosity, 2) 1 fb-1
        sLuminosity_toUse = 'Full' if luminosity_toUse == luminosity_total else str(luminosity_toUse)

        luminosity_Scaling_toUse = luminosity_toUse / luminosity_total
        if luminosity_toUse not in evtYields_dict[sDatasetName]:
            evtYields_dict[sDatasetName][sLuminosity_toUse] = {}

        sOpDir  = '%s/%s/%s/%s' % (sAnaDir, sDatasetName, anaSuperCat, sOpDirNameShort)
        if not os.path.exists(sOpDir):
            os.makedirs(sOpDir)

        if printLevel >= 0: 
            print(f"{sDatasetName}: {luminosity_toUse = }, {luminosity_total = },  {luminosity_Scaling_toUse = }, \n{YearsToRun_list = }, \n{sOpDir = }\n", flush=True)


        ExpData_dict = {} #{ 'Year': DataObs_DirName_dict[sDatasetName] }
        for Year_ in YearsToRun_list:
            ExpData_dict[Year_] = DataObs_DirName_dict[Year_]
        samples_dict = {**MCSig_dict, **MCBkg_dict, **ExpData_dict}
        #print(f"{samples_dict = }")

        eventYields = OD()
        sEventYield = "%s %s %s\n" % (sAnaDir, Dataset, CAT)
        MCBkg_total = 0
        MCSig_total = 0
        for selectionTag in selectionTags:   
            for systematics_name in systematics_list:
                sSelTag_Syst = '%s_%s' % (selectionTag, systematics_name)
                if sSelTag_Syst not in evtYields_dict[sDatasetName][sLuminosity_toUse]:
                    evtYields_dict[sDatasetName][sLuminosity_toUse][sSelTag_Syst] = {}

                sEventYield += "Selection%s %s %s %s" % (sSpace, selectionTag, systematics_name, sEndline)
                sEventYield += "%-45s %s  %10s %s %10s %s %12s %s" % ('Samples', sSpace,'nEvt', sPlusMinus, 'Uncert.', sSpace, 'nEvt unwgt', sEndline)
                eventYields[systematics_name] = OD()
                MCBkg_total = 0;   MCBkg_total_variance = 0
                MCSig_total = 0;   MCSig_total_variance = 0
                Data_total  = 0;   Data_total_variance  = 0

                if printLevel >= 5:
                    print(f"\t\t{selectionTag}, {systematics_name}")

                #for sampleNameShort, sample_category_list in (MCBkg_dict+MCSig_dict+ExpData_dict).items():
                for sampleNameShort, sample_category_list in samples_dict.items():
                    if printLevel >= 5:
                        print(f"\t\t\t{sampleNameShort}")
                    nEvts_sum_    = 0; nEvts_unwgt_sum_ = 0
                    variance_sum_ = 0
                    if sampleNameShort not in evtYields_dict[sDatasetName][sLuminosity_toUse][sSelTag_Syst]:
                        evtYields_dict[sDatasetName][sLuminosity_toUse][sSelTag_Syst][sampleNameShort] = {}
                    for sample_category in sample_category_list:
                        #systematics_name_toUse = systematics_name if sampleNameShort not in ExpData_dict.keys() else systematics_forData 
                        #systematics_name_toUse = '_%s' % (systematics_name_toUse) # Previous version had systematics name for 'hCutFlow' histograqms
                        systematics_name_toUse = ''                               # New version dropped systematics name for 'hCutFlow' histograqms
                        if printLevel >= 5:
                            print(f"\t\t\t\t{sample_category}")

                        h = None; iYear_ = 0
                        for Year_ in YearsToRun_list:
                            if sampleNameShort in ExpData_dict.keys(): # data
                                if sampleNameShort != Year_: # <data_dir> (sampleNameShort for data) is found only in fIpFiles of the corresponding year
                                    continue
                            # "hCutFlowWeighted"
                            histo_name_full = 'evt/%s/%s_%s%s' % (sample_category, "hCutFlowWeighted", selectionTag, systematics_name_toUse)
                            if printLevel >= 6:
                                print(f"\t\t\t\t\t{histo_name_full} {Year_}")
                            if ((histo_name_full in sIpFileKeys_dict[Year_]) or ('%s;1'%(histo_name_full) in sIpFileKeys_dict[Year_])):
                                h_i = fIpFiles[Year_][histo_name_full].to_hist()
                            else:
                                if Year_ not in sProcessNames_NotRead:
                                    sProcessNames_NotRead[Year_] = []
                                if sample_category not in sProcessNames_NotRead[Year_]:
                                    sProcessNames_NotRead[Year_].append( sample_category )
                                continue
                            if iYear_ == 0:  h = h_i
                            else:            h = h + h_i
                            iYear_ += 1
                        if h == None: continue
                        h = h * luminosity_Scaling_toUse # scale nEvts_Wgt with luminosity if requires
                        nEvt_and_var = getTH1BinContent(h, xValue=xValueForEventYield)

                        h = None; iYear_ = 0
                        for Year_ in YearsToRun_list:
                            if sampleNameShort in ExpData_dict.keys(): # data
                                if sampleNameShort != Year_: # <data_dir> (sampleNameShort for data) is found only in fIpFiles of the corresponding year
                                    continue
                            # "hCutFlow"
                            histo_name_full = 'evt/%s/%s_%s%s' % (sample_category, "hCutFlow", selectionTag, systematics_name_toUse)
                            if ((histo_name_full in sIpFileKeys_dict[Year_]) or ('%s;1'%(histo_name_full) in sIpFileKeys_dict[Year_])):
                                h_i = fIpFiles[Year_][histo_name_full].to_hist()
                            else:
                                if Year_ not in sProcessNames_NotRead:
                                    sProcessNames_NotRead[Year_] = []
                                if sample_category not in sProcessNames_NotRead[Year_]:
                                    sProcessNames_NotRead[Year_].append( sample_category )
                                continue
                            if iYear_ == 0:  h = h_i
                            else:            h = h + h_i
                            iYear_ += 1
                        nEvt_and_var_unwgt = getTH1BinContent(h, xValue=xValueForEventYield)

                        nEvts_sum_       += nEvt_and_var.value
                        variance_sum_    += nEvt_and_var.variance
                        nEvts_unwgt_sum_ += nEvt_and_var_unwgt.value

                        #sEventYield += "\t %-45s:  %15g +- %15g  (%15d) \n" % (sample_category, nEvt_and_var.value, np.sqrt(nEvt_and_var.variance), nEvt_and_var_unwgt.value)

                    if sampleNameShort not in ExpData_dict.keys():
                        sEventYield += "%-45s %s  %10.1f %s %10.1f %s %12d %s" % (sampleNameShort, sSpace, nEvts_sum_, sPlusMinus, np.sqrt(variance_sum_), sSpace, nEvts_unwgt_sum_, sEndline)

                    if sampleNameShort in MCBkg_dict.keys():
                        MCBkg_total          += nEvts_sum_
                        MCBkg_total_variance += variance_sum_
                    if sampleNameShort in ExpData_dict.keys():
                        Data_total          += nEvts_sum_
                        Data_total_variance += variance_sum_

                    if sampleNameShort not in ExpData_dict.keys():
                        evtYields_dict[sDatasetName][sLuminosity_toUse][sSelTag_Syst][sampleNameShort] = {
                            'nEvents': nEvts_sum_,
                            'Uncertainty': np.sqrt(variance_sum_),
                            'nEventsUnwgt': nEvts_unwgt_sum_
                        }

                    eventYields[systematics_name][sampleNameShort] = OD([
                        ('nEvents', nEvts_sum_),
                        ('Uncertainty', np.sqrt(variance_sum_)),
                        ('nEventsUnwgt', nEvts_unwgt_sum_),
                        ('isMCBkg', sampleNameShort in MCBkg_dict.keys())
                    ]) 
                        
                    
                sEventYield += "\n"
                sEventYield += "%-45s %s  %10.1f %s %10.1f  %s" % ('Bkg_total ', sSpace, MCBkg_total, sPlusMinus, np.sqrt(MCBkg_total_variance), sEndline)
                sEventYield += "%-45s %s  %10.1f %s %10.1f  %s" % ('Data_total', sSpace, Data_total, sPlusMinus, np.sqrt(Data_total_variance), sEndline)
                sEventYield += "%-45s %s  %10g  %s" % ('Data/MCBkg', sSpace, Data_total / MCBkg_total if MCBkg_total > 0 else 0, sEndline)

                sEventYield += "\n\n" 

                evtYields_dict[sDatasetName][sLuminosity_toUse][sSelTag_Syst]['Bkg_total '] = {
                    'nEvents':      MCBkg_total,
                    'Uncertainty':  np.sqrt(MCBkg_total_variance),
                    'nEventsUnwgt': -1
                }
                evtYields_dict[sDatasetName][sLuminosity_toUse][sSelTag_Syst]['Data_total'] = {
                    'nEvents':      Data_total,
                    'Uncertainty':  np.sqrt(Data_total_variance),
                    'nEventsUnwgt': -1
                }
                evtYields_dict[sDatasetName][sLuminosity_toUse][sSelTag_Syst]['Data/MCBkg'] = {
                    'nEvents':      Data_total / MCBkg_total if MCBkg_total > 0 else 0,
                    'Uncertainty':  -1,
                    'nEventsUnwgt': -1
                }

                continue

                for idx_, sample_category in enumerate(MCSig_dict.keys()):
                    if abs(MCBkg_total - 0) < 1e-6: continue
                    sEventYield += "%-45s:  %15g  \n" % ('S/sqrt(B) @ %s' % (sample_category), eventYields[systematics_name][sample_category]['nEvents'] / np.sqrt(MCBkg_total))
                sEventYield += "\n\n" 
            

                ## print sorted event yield table
                #bkgs_list = list(eventYields[systematics_name].keys())
                bkgs_list = [ process for process in eventYields[systematics_name].keys() if eventYields[systematics_name][process]['isMCBkg']  ]
                nEvent_bkgs_list = [ eventYields[systematics_name][bkg]['nEvents'] for bkg in bkgs_list]
                bkgRanks  = np.argsort(np.array(nEvent_bkgs_list))[::-1]
                #print(f"{bkgs_list = }")
                #print(f"{nEvent_bkgs_list = }")
                #print(f"{np.argsort(np.array(nEvent_bkgs_list)) = }")
                #print(f"{bkgRanks = }")

                sEventYield += "Sorted background processes::\n"
                for bkgRank in bkgRanks:
                    processName_ = bkgs_list[bkgRank]
                    sEventYield += "%-45s:  %15.1f +- %15.1f  (%15d) \n" % (
                        processName_, 
                        eventYields[systematics_name][processName_]['nEvents'],
                        eventYields[systematics_name][processName_]['Uncertainty'],
                        eventYields[systematics_name][processName_]['nEventsUnwgt'])
                sEventYield += "\n\n" 

                

        print("%s" % (sEventYield))
        #print(f"{json.dumps(eventYields, indent=4) = }")
        #json.dumps(eventYields, indent=4)

        print("\n\nsProcessNames_NotRead: ")
        print(json.dumps(sProcessNames_NotRead, indent=4))

        with open('%s/eventYieldTable_lumi%g.txt' % (sOpDir, round(luminosity_toUse,1)), 'w') as fOut:
            fOut.write(sEventYield)
            fOut.write(json.dumps(sProcessNames_NotRead, indent=4))

if len(Years) == 1: exit(0)


## Print event yields for all datasets side-by-side --> Innermost loop over Datasets 
sDatasetName_0_ = list(evtYields_dict.keys())[0]
sDatasetName_last_ = list(evtYields_dict.keys())[-1]
for sLuminosity_toUse in evtYields_dict[sDatasetName_0_]: # Luminosity: 'Full', '1.0'
    sOpDir  = '%s/%s/%s/%s' % (sAnaDir, 'Run2', anaSuperCat, sOpDirNameShort)
    if not os.path.exists(sOpDir):
        os.makedirs(sOpDir)
    
    sEventYield = "%s %s %s %s\n" % (sAnaDir, 'Run2', CAT, sLuminosity_toUse)
    for selectionTag in selectionTags:   
        for systematics_name in systematics_list:
            sSelTag_Syst = '%s_%s' % (selectionTag, systematics_name)
            sEventYield += "Selection%s %s %s %s" % (sSpace, selectionTag, systematics_name, sEndline)
            sEventYield += "%-45s %s" % ('       ', sSpace )
            for sDatasetName in YearsToRun_dict:
                sEndline_toUse = sEndline if sDatasetName == sDatasetName_last_ else sSpace
                sEventYield += "  %10s %s %10s %s %12s %s" % (sDatasetName, sPlusMinus, ' ', sSpace, ' ', sEndline_toUse)
            sEventYield += "%-45s %s" % ('Samples', sSpace)
            for sDatasetName in YearsToRun_dict:
                sEndline_toUse = sEndline if sDatasetName == sDatasetName_last_ else sSpace
                sEventYield += "  %10s %s %10s %s %12s %s" % ('nEvt', sPlusMinus, 'Uncert.', sSpace, 'nEvt unwgt', sEndline_toUse)
                

            for sampleNameShort in evtYields_dict[sDatasetName_0_][sLuminosity_toUse][sSelTag_Syst]:
                if sampleNameShort in YearsToRun_dict: continue # skip sampleName==Year as 'Data' is stored separately
                sEventYield += "%-45s %s" % (sampleNameShort, sSpace)

                for sDatasetName in YearsToRun_dict:
                    nEvents_      = evtYields_dict[sDatasetName][sLuminosity_toUse][sSelTag_Syst][sampleNameShort]['nEvents']
                    Uncertainty_  = evtYields_dict[sDatasetName][sLuminosity_toUse][sSelTag_Syst][sampleNameShort]['Uncertainty']
                    nEventsUnwgt_ = evtYields_dict[sDatasetName][sLuminosity_toUse][sSelTag_Syst][sampleNameShort]['nEventsUnwgt']
                    sEndline_toUse = sEndline if sDatasetName == sDatasetName_last_ else sSpace
                    sTemplate_ = "  %10.1f %s %10.1f %s %12d %s"
                    if sampleNameShort == 'Data/MCBkg':
                        sTemplate_ = "  %10g %s %10.1f %s %12d %s"
                    sEventYield += sTemplate_ % (nEvents_, sPlusMinus, Uncertainty_, sSpace, nEventsUnwgt_, sEndline_toUse)
            sEventYield += "\n\n" 

    print("%s" % (sEventYield))

    with open('%s/eventYieldTable_compare_lumi%s.txt' % (sOpDir, sLuminosity_toUse), 'w') as fOut:
        fOut.write(sEventYield)
        fOut.write(json.dumps(sProcessNames_NotRead, indent=4))    
