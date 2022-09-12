
from collections import OrderedDict as OD

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



config_Template = OD([
    ("nEventsToAnalyze", -1),
    ("era", ''),
    ("inputFiles", ''),
    ("outputFile", ''),
    ("sampleCategory", ''),
    #("Luminosity", 0),
    ("crossSection", 0),
    ("nEvents", -1),
    ("sumEvents", -1),
])
