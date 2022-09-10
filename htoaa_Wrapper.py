

import json
import glob
import argparse
from collections import OrderedDict as OD

from htoaa_Settings import *
from htoaa_Samples import Samples2018


parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
parser.add_argument('-era', dest='era', type=str, default=Era_2018, choices=[Era_2016, Era_2017, Era_2018], required=False)
parser.add_argument('-run_mode', type=str, default='local', choices=['local', 'condor'])
parser.add_argument('-v', '--version', type=str, default=None, required=True)
parser.add_argument('-nFilesPerJob', type=int, default=1)
args=parser.parse_args()
print("args: {}".format(args))

era           = args.era
run_mode      = args.run_mode
nFilesPerJob  = args.nFilesPerJob
anaVersion    = args.version 






sAnalysis = "htoaa_Analysis.py"
sConfig   = "config_htoaa.json"

samplesList = None
samplesInfo = None
if era == Era_2018:
    samplesList = Samples2018
    with open(sFileSamplesInfo[era]) as fSamplesInfo:
        samplesInfo = json.load(fSamplesInfo)


print("samplesList: {}".format(samplesList))
print("\n\nsamplesInfo: {}".format(samplesInfo))


config = OD()



for sample_category, samples in samplesList.items():
    print("sample_category {}, samples {}".format(sample_category, samples))
    for sample in samples:
        print("\nsample: {}".format(sample))
        print("samplesInfo[sample]: {}".format(samplesInfo[sample]))

        
        
        print("files_0: {}".format(samplesInfo[sample][sampleFormat]))

        file1 = samplesInfo[sample][sampleFormat][0]
        print("files_1: {}".format(file1))
        print("files_2: {}".format(glob.glob(file1)))

        

