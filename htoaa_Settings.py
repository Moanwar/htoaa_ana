
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

Era_2016 = '2016'
Era_2017 = '2017'
Era_2018 = '2018'

sFileSamplesInfo = {
    Era_2016: "Samples_2016UL.json",
    Era_2017: "Samples_2017UL.json",
    Era_2018: "Samples_2018UL.json"
}

Luminosities = { # [<lumi>, <uncertainty in percent> ] in fb^-1
    Era_2016: [36.31, 1.2],
    Era_2017: [41.48, 2.3],
    #Era_2018: [59.83, 2.5], # 2018 recommendation - inclusive 
    Era_2018: [54.54, 2.5], # for HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_v* trigger data/luminosity/2018/output_brilcalc_314472-325175_UL18_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_Final.xlsx 
}

sFilesGoldenJSON = {
    Era_2016: '',
    Era_2017: '',    
    Era_2018: 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt',    
}

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

Corrections = {
        
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

}

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
    },
}


### Miscellaneous
kLHE_HT_Max = 99999.0


GENPART_STATUSFLAGS = [
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
])
