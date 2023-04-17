
from collections import OrderedDict as OD

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

Luminosities = { # [<lumi>, <uncertainty in percent> ]
    Era_2016: [36.31, 1.2],
    Era_2017: [41.48, 2.3],
    Era_2018: [59.83, 2.5],
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
