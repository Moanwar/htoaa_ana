
from collections import OrderedDict as OD
import enum

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookXrootdService
# cms-xrd-global.cern.ch "global redirector"
# xrootd-cms.infn.it for Europe and Asia
xrootd_redirectorName  = "root://xrootd-cms.infn.it//"
# try / or // at the end
xrootd_redirectorNames = [
    "root://xrootd-cms.infn.it/",
    "root://cms-xrd-global.cern.ch/",
    "root://cmsxrootd.fnal.gov/",
    "root://cms-xrd-global.cern.ch:1094/",
]
sampleFormat = "nanoAOD" 

### Miscellaneous constants
kPt_Max     = 99999.0
kLHE_HT_Max = 99999.0
NanoAODFileSize_Min = 0.3 # in MB
#------------------------------------

Era_2016 = '2016'
Era_2017 = '2017'
Era_2018 = '2018'

sFileSamplesInfo = {
    Era_2016: "Samples_2016UL.json",
    Era_2017: "Samples_2017UL.json",
    Era_2018: "Samples_2018UL.json"
}

# Refer https://docs.google.com/spreadsheets/d/1xDLsr3ikLJxuMPNiSRs79YjTzbN64RetXL3A-tL6-hY/edit?usp=sharing
# /eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/QCD*/r1/PNet_*.root
sPathSkimmedNanoAODs = {
    Era_2018: {
        'unskimmed': {
            'Data': '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/PNet_v1_2023_10_06/$SAMPLETAG/$SAMPLENAME/r*/PNet_*.root',
            'MC':   '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/$SAMPLENAME/r1/PNet_*.root' 
        },
        'skim_Hto4b_0p8': {
            'Data': '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/PNet_v1_2023_10_06/$SAMPLETAG/$SAMPLENAME/skims/Hto4b_0p8/PNet_*.root',
            'MC':   '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/$SAMPLENAME/skims/Hto4b_0p8/PNet_*.root' 
        },
    }
}

Luminosities_Inclusive = { # [<lumi>, <uncertainty in percent> ] in fb^-1
    Era_2016: [36.31, 1.2],
    Era_2017: [41.48, 2.3],
    Era_2018: [59.83, 2.5]
}

Luminosities_forGGFMode = { # [<lumi>, <uncertainty in percent> ] in fb^-1
    Era_2016: [36.31, 1.2],
    Era_2017: [41.48, 2.3],
    Era_2018: {
        'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4': [54.54, 2.5], # for HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_v* trigger. See ./data/luminosity/2018/output_brilcalc_314472-325175_UL18_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_Final.xlsx 
        'Trg_Combo_AK4AK8Jet_HT':                             [59.83, 2.5], # See ./data/luminosity/2018/Luminosity_HLTPaths.xlsx
        'HLT_IsoMu24':                                        [59.82, 2.5], # See ./data/luminosity/2018/output_brilcalc_314472-325175_UL18_HLT_IsoMu24_v.xlsx
        'HLT_IsoMu27':                                        [59.83, 2.5], # See ./data/luminosity/2018/output_brilcalc_314472-325175_UL18_HLT_IsoMu27_v.xlsx
        'Trg_Combo_Mu':                                       [59.83, 2.5], # See ./data/luminosity/2018/Luminosity_HLTPaths.xlsx
    }, 
}
Luminosities_forGGFMode_perEra = {
    Era_2018: {
        'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4': {
            'A':  8.736,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839
        },
        'Trg_Combo_AK4AK8Jet_HT': {
            'A': 14.027,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839,
        },
        'HLT_IsoMu24': {
            'A': 14.019,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839
        }, 
        'HLT_IsoMu27': {
            'A': 14.027,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839
        },  
        'Trg_Combo_Mu': {
            'A': 14.027,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839
        }, 
    }
}
Triggers_perEra = {
    Era_2018: {
        'Trg_Combo_AK4AK8Jet_HT': {
            'HLT_PFJet500':                                       ['L1_SingleJet180'], 
            'HLT_PFHT1050':                                       ['L1_HTT360er'],
            'HLT_AK8PFHT800_TrimMass50':                          ['L1_HTT360er'],
            'HLT_AK8PFJet500':                                    ['L1_SingleJet180'],
            'HLT_AK8PFJet400_TrimMass30':                         ['L1_SingleJet180'],
            'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4': ['L1_SingleJet180'],
        },
        'Trg_Combo_Mu': {
            'HLT_IsoMu24': ['L1_SingleMu22'],
            'HLT_IsoMu27': ['L1_SingleMu22', 'L1_SingleMu25'],
            'HLT_Mu50':    ['L1_SingleMu22', 'L1_SingleMu25'],
        }
    }
}


sFilesGoldenJSON = {
    Era_2016: '',
    Era_2017: '',    
    Era_2018: 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt',    
}

# List of recommended MET filters. https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
MET_Filters = {} 
MET_Filters[Era_2016] = {
    "Data": [
        "goodVertices",                       # primary vertex filter ("Flag_goodVertices")
        "globalSuperTightHalo2016Filter",     # beam halo filter ("Flag_globalSuperTightHalo2016Filter")
        "HBHENoiseFilter",                    # HBHE noise filter ("Flag_HBHENoiseFilter")
        "HBHENoiseIsoFilter",                 # HBHEiso noise filter ("Flag_HBHENoiseIsoFilter")
        "EcalDeadCellTriggerPrimitiveFilter", # ECAL TP filter ("EcalDeadCellTriggerPrimitiveFilter")
        "BadPFMuonFilter",                    # Bad PF Muon Filter ("Flag_BadPFMuonFilter")
        "BadPFMuonDzFilter",                  # Bad PF Muon Dz Filter ("Flag_BadPFMuonDzFilter")
        "eeBadScFilter",                      # ee badSC noise filter ("Flag_eeBadScFilter")
        "hfNoisyHitsFilter",                  # HF noisy hits filter ("Flag_hfNoisyHitsFilter")
    ]
}
MET_Filters[Era_2016]["MC"]  = MET_Filters[Era_2016]["Data"] 

MET_Filters[Era_2018] = {
    "Data": [
        "goodVertices",                       # primary vertex filter ("Flag_goodVertices")
        "globalSuperTightHalo2016Filter",     # beam halo filter ("Flag_globalSuperTightHalo2016Filter")
        "HBHENoiseFilter",                    # HBHE noise filter ("Flag_HBHENoiseFilter")
        "HBHENoiseIsoFilter",                 # HBHEiso noise filter ("Flag_HBHENoiseIsoFilter")
        "EcalDeadCellTriggerPrimitiveFilter", # ECAL TP filter ("Flag_EcalDeadCellTriggerPrimitiveFilter")
        "BadPFMuonFilter",                    # Bad PF Muon Filter ("Flag_BadPFMuonFilter")
        "BadPFMuonDzFilter",                  # Bad PF Muon Dz Filter ("Flag_BadPFMuonDzFilter")
        "hfNoisyHitsFilter",                  # HF noisy hits filter ("Flag_hfNoisyHitsFilter")
        "eeBadScFilter",                      # ee badSC noise filter ("Flag_eeBadScFilter")
        "ecalBadCalibFilter",                 # ECAL bad calibration filter update ("Flag_ecalBadCalibFilter")
    ],
}
MET_Filters[Era_2018]["MC"]  = MET_Filters[Era_2018]["Data"]
MET_Filters[Era_2017]        = MET_Filters[Era_2018]
# -----------------------------------------------------------------------------------


class JetIDs(enum.IntEnum):
    tightIDFailingLeptonVeto = 2
    tightIDPassingLeptonVeto = 6


kMCSamplesStitch_PhSpOverlapRemove = 'MCSamplesStitch_PhSpOverlapRemove'
kMCSamplesStitch_PhSpOverlapRewgt  = 'MCSamplesStitch_PhSpOverlapRewgt'

class MCSamplesStitchOptions(enum.Enum):
    PhSpOverlapRemove = 'PhSpOverlapRemove'
    PhSpOverlapRewgt  = 'PhSpOverlapRewgt' 


sFileLumiScalesPhSpOverlapRewgt = {
    Era_2018: {
        "inputFile":     "data/lumiScale/2018.root",
        "histogramName": "$SAMPLECATEGORY/$SAMPLECATEGORY_LumiScale_PhSpOverlapRewghted",
    }
}

HEM1516Issue2018_AffectedRunRange = [319077, 325175]
DataFractionAffectedBy2018HEM1516Issue = 0.7105 # factor = (luminosity for run >= 319077) / (2018 luminosity) = 38.7501 / 54.5365. Calculated for 2018 HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 trigger 

bTagWPs = { # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
    Era_2018: {
        'DeepCSV': { # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
            'L': 0.1208,
            'M': 0.4506,
        },
        'DDBvL': { # not provided for UL samples
            'M': 0.8, # taken from Si's code
        },
        'DDBvLV2': { # not provided for UL samples
            'M': 0.8, # taken from Si's code
        },
        'ParticleNetMD_XbbvsQCD': { # BTV-22-001 
            # https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=BTV-22-001&tp=an&id=2622&ancode=BTV-22-001
            # https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2021/005
            'L': 0.9172
        },
        'ParticleNetMD_Hto4b_Htoaa4bOverQCD': {
            # https://ssawant.web.cern.ch/ssawant/HToAA/DatavsMC/20231106_PNetSignificanceScan_Msd90to140/?match=ParticleNetMD_Hto4b_Htoaa4bOverQCD         
            'WP-80': 0.920,
            #'WP-60': 0.978,
            'WP-60': 0.975, # https://indico.cern.ch/event/1348321/?note=257291#31-saswati-nandan
            'WP-40': 0.992,
        }
    },
}


Corrections = {

    "PURewgt": {
        Era_2018: {
            "inputFile":     "data/correction/mc/PURewgt/PURewgts_2018.root",
            "histogramName": "MC_PURewgt"
        }

    },
        
    "HTRewgt" : { # ./data/correction/mc/HTSamplesStitch/HTSamplesStitchSF_2018.root
        "QCD_bGen": {
            "2018": {
                "FitFunctionFormat": "{p0} + ({p1} * (x - {HTBinMin}))",
                "HT100to200": "0.927235 + (0.001153 * (x - 100))",
                "HT200to300": "0.936276 + (0.000820 * (x - 200))",
                "HT300to500": "0.936325 + (0.000512 * (x - 300))",
                "HT500to700": "0.968147 + (0.000378 * (x - 500))",
                "HT700to1000": "0.931336 + (0.000229 * (x - 700))",
                "HT1000to1500": "0.956714 + (0.000122 * (x - 1000))",
                "HT1500to2000": "0.965210 + (0.000068 * (x - 1500))",
                "HT2000to3000": "1.006649 + (0.000030 * (x - 2000))"
            }
        }
    }, 

    "TopPtRewgt": { # Top pT reweights for ttbat sample. https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf#page=12
        "TuneCP5": {
            "FitFunctionFormat": "exp( {a} + ({b} * x) + ({c} * x * x) + ({d}/(x + {e})) )",
            "FitFunction": "exp( -2.02274e-01 + (1.09734e-04 * x) + (-1.30088e-07 * x * x) + (5.83494e+01/(x + 1.96252e+02)) )",
            "FitRange": [0, 3000],
        },
        "TuneCUETP": {
            "FitFunctionFormat": "{a} + ({b} * TanH({c} + ({d} * x) )",
            "FitFunction": "1.04554e+00 + (5.19012e-02 * TanH(-1.72927e+00 + (2.57113e-03 * x) )",
            "FitRange": [0, 3000],
        }
    },
    
    'ParticleNetMD_XbbvsQCD': { # Data-to-MC SFs for ParticleNetMD_XbbvsQCD:  BTV-22-001 
        # https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=BTV-22-001&tp=an&id=2622&ancode=BTV-22-001
        # https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2021/005
        Era_2018: {
            'L': {
                'pT_binEdges': [ 0,      450,      500,       600,    kPt_Max ],
                'SFs':         [      1,     0.921,     1.006,     1.001      ],
            },
        },

    },

}

massPseudoscalarA_windows_dict = OD([
    ('mA15Window', [14.2, 15.6]),
    ('mA30Window', [28.0, 32.0]),
    ('mA55Window', [50.0, 58.0]),
])



### Miscellaneous variables
HistogramNameExtensions_QCD = ['0bCat', '1bCat', '2bCat', '3bCat', '4bCat', '5bAndMoreCat'] 



### GEN-level variables 
PDGID_BottomQuark = 5
PDGID_TopQuark    = 6

MASS_BottomQuark = 4.18

GENPART_STATUSFLAGS_LIST = [
    "isPrompt",
    "isDecayedLeptonHadron",
    "isTauDecayProduct",
    "isPromptTauDecayProduct",
    "isDirectTauDecayProduct",
    "isDirectPromptTauDecayProduct",
    "isDirectHadronDecayProduct",
    "isHardProcess",
    "fromHardProcess",
    "isHardProcessTauDecayProduct",
    "isDirectHardProcessTauDecayProduct",
    "fromHardProcessBeforeFSR",
    "isFirstCopy",
    "isLastCopy",
    "isLastCopyBeforeFSR",
]

class GENPART_STATUSFLAGS(enum.IntEnum):
    isPrompt                           =  0
    isDecayedLeptonHadron              =  1
    isTauDecayProduct                  =  2
    isPromptTauDecayProduct            =  3
    isDirectTauDecayProduct            =  4
    isDirectPromptTauDecayProduct      =  5
    isDirectHadronDecayProduct         =  6
    isHardProcess                      =  7
    fromHardProcess                    =  8
    isHardProcessTauDecayProduct       =  9
    isDirectHardProcessTauDecayProduct = 10
    fromHardProcessBeforeFSR           = 11
    isFirstCopy                        = 12
    isLastCopy                         = 13
    isLastCopyBeforeFSR                = 14


config_Template = OD([
    ("nEventsToAnalyze", -1),
    ("era", ''),
    ("dataset", ''), 
    ("inputFiles", ''),
    ("outputFile", ''),
    ("sampleCategory", ''),
    ("isMC", False),
    #("Luminosity", 0),
    ("crossSection", 0),
    ("nEvents", -1),
    ("sumEvents", -1),
    ("downloadIpFiles", False),   
    ("server", ''),   
])
