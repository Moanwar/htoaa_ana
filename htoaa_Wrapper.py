
import os
import json
import glob
import argparse
from collections import OrderedDict as OD
import numpy as np

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

pwd = os.getcwd()
DestinationDir = "./%s" % (anaVersion)
if   "/home/siddhesh/" in pwd: DestinationDir = "/home/siddhesh/Work/CMS/htoaa/analysis/%s" % (anaVersion)
elif "/afs/cern.ch/"   in pwd: DestinationDir = "/afs/cern.ch/work/s/ssawant/private/htoaa/analysis/%s" % (anaVersion)
if not os.path.exists(DestinationDir): os.mkdir( DestinationDir )



sAnalysis = "htoaa_Analysis.py"
sConfig   = "config_htoaa.json"

samplesList = None
samplesInfo = None
Luminosity  = None
if era == Era_2018:
    samplesList = Samples2018
    with open(sFileSamplesInfo[era]) as fSamplesInfo:
        samplesInfo = json.load(fSamplesInfo)
    Luminosity = Luminosities[era][0]


print("samplesList: {}".format(samplesList))
print("\n\nsamplesInfo: {}".format(samplesInfo))


config = config_Template

for sample_category, samples in samplesList.items():
    #print("sample_category {}, samples {}".format(sample_category, samples))
    for sample in samples:
        sampleInfo = samplesInfo[sample]
        fileList = sampleInfo[sampleFormat]
        files = []
        for iEntry in fileList:
            # file name with wildcard charecter *
            if "*" in iEntry:  files.extend( glob.glob( iEntry ) )
            else:              files.append( iEntry )
        sample_cossSection = sampleInfo["cross_section"]
        sample_nEvents     = sampleInfo["nEvents"]
        sample_sumEvents   = sampleInfo["sumEvents"]
            
        print("\nsample: {}".format(sample))
        print("samplesInfo[sample]: {}".format(samplesInfo[sample]))
        print("files ({}): {}".format(len(files), files))
        
        nSplits = int(len(files) / nFilesPerJob) + 1
        
        files_splitted = np.array_split(files, nSplits)
        print("files_splitted: {}".format(files_splitted))
        '''        
        files_splitted = []
        for iSplit in range(nSplits):
            idxStart = iSplit*nFilesPerJob
            idxEnd   = idxStart + nFilesPerJob if iSplit < (nSplits - 1) else None
            files_splitted.append( files[ idxStart : idxEnd ]  )
        print("files_splitted: {}".format(files_splitted))
        '''
        for iJob in range(len(files_splitted)):
            config["era"] = era
            config["inputFiles"] = list( files_splitted[iJob] )
            config["outputFile"] = '%s/analyze_htoaa_%s_0_%d' % (DestinationDir, sample, iJob)
            config["sampleCategory"] = sample_category
            #config["Luminosity"] = Luminosity
            config["crossSection"] = sample_cossSection
            config["nEvents"] = sample_nEvents
            config["sumEvents"] = sample_sumEvents

            print("config: {}".format(config))

            sConfig_to_use = "%s/%s" % (DestinationDir, sConfig.replace(".json", "_0_%d.json" % (iJob)))
            with open(sConfig_to_use, "w") as fConfig:
                json.dump( config,  fConfig, indent=4)
            
            

        

        

