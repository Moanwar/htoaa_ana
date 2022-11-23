
import os
import json
import glob
import argparse
from collections import OrderedDict as OD
import numpy as np

from htoaa_Settings import *
from htoaa_Samples import (
    Samples2018,
    kData
)


parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
parser.add_argument('-era', dest='era', type=str, default=Era_2018, choices=[Era_2016, Era_2017, Era_2018], required=False)
parser.add_argument('-run_mode',        type=str, default='local', choices=['local', 'condor'])
parser.add_argument('-v', '--version',  type=str, default=None, required=True)
parser.add_argument('-samples',         type=str, default=None, help='samples to run seperated by comma')
parser.add_argument('-nFilesPerJob',    type=int, default=10)
parser.add_argument('-forcefully',      action='store_true', default=False)
args=parser.parse_args()
print("args: {}".format(args))

era              = args.era
run_mode         = args.run_mode
nFilesPerJob     = args.nFilesPerJob
selSamplesToRun  = args.samples
anaVersion       = args.version
submitForcefully = args.forcefully

pwd = os.getcwd()
DestinationDir = "./%s" % (anaVersion)
if   "/home/siddhesh/" in pwd: DestinationDir = "/home/siddhesh/Work/CMS/htoaa/analysis/%s" % (anaVersion)
elif "/afs/cern.ch/"   in pwd: DestinationDir = "/afs/cern.ch/work/s/ssawant/private/htoaa/analysis/%s" % (anaVersion)
if not os.path.exists(DestinationDir): os.mkdir( DestinationDir )



sAnalysis = "htoaa_Analysis_wCoffea.py"  # "htoaa_Analysis.py"
sConfig   = "config_htoaa.json"

samplesList = None
samplesInfo = None
Luminosity  = None
if era == Era_2018:
    samplesList = Samples2018
    with open(sFileSamplesInfo[era]) as fSamplesInfo:
        samplesInfo = json.load(fSamplesInfo)
    Luminosity = Luminosities[era][0]

selSamplesToRun_list = []
if selSamplesToRun:
    selSamplesToRun_list = selSamplesToRun.split(',')

print("samplesList: {}".format(samplesList))
print("\n\nsamplesInfo: {}".format(samplesInfo))
print(f"\n\nselSamplesToRun_list: {selSamplesToRun_list}")


config = config_Template

for sample_category, samples in samplesList.items():
    #print("sample_category {}, samples {}".format(sample_category, samples))
    for sample in samples:
        if len(selSamplesToRun_list) > 0:
            skipThisSample = True
            for selSample in selSamplesToRun_list:
                if sample.startswith(selSample): skipThisSample = False
            if skipThisSample:
                continue

        print(f"sample_category: {sample_category}, sample: {sample}")
            
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

        
        nSplits = int(len(files) / nFilesPerJob) + 1 if nFilesPerJob > 0 else 1
        
        
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
            config["outputFile"] = '%s/analyze_htoaa_%s_0_%d.root' % (DestinationDir, sample, iJob)
            config["sampleCategory"] = sample_category
            config["isMC"] = (sample_category != kData)
            #config["Luminosity"] = Luminosity
            config["crossSection"] = sample_cossSection
            config["nEvents"] = sample_nEvents
            config["sumEvents"] = sample_sumEvents

            #outputFile_tmp = config["outputFile"].replace('.root', '_wCoffea.root')            
            if (not submitForcefully) and (os.path.exists(outputFile_tmp)):
                print(f"Skipping submission for {outputFile_tmp}")
                continue # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            configName = "%s" % sConfig.replace(".json", "_%s_0_%d.json" % (sample, iJob))
            sConfig_to_use = "%s/%s" % (DestinationDir, configName)
            print("config {}: {}".format(sConfig_to_use, config))
            with open(sConfig_to_use, "w") as fConfig:
                json.dump( config,  fConfig, indent=4)
            
            
            ## condor job scripts ----------------------------------------------------------------------------
            condorJobName = "%s_%s" % (sAnalysis.replace(".py", ""), configName.replace(".json", ""))
            condor_exec_file = '%s/condor_%s_exec.sh' % (DestinationDir, condorJobName) 
            if not os.path.isfile(condor_exec_file):
                with open(condor_exec_file, 'w') as f:
                    f.write("#!/bin/bash  \n\n")
                    f.write("cd %s \n" % pwd)
                    f.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n")
                    f.write("export SCRAM_ARCH=slc6_amd64_gcc700  \n")
                    f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")
                    #f.write("cd ")
                    f.write("export X509_USER_PROXY=/afs/cern.ch/user/s/ssawant/x509up_u108989  \n")
                    f.write("eval \n")
                    f.write("cd %s \n" % (pwd))
                    f.write("source /afs/cern.ch/user/s/ssawant/.bashrc \n")
                    f.write("which conda \n")
                    f.write("time conda env list \n")
                    f.write("conda activate ana_htoaa \n")
                    f.write("time conda env list \n")
                    
                    f.write("time conda list \n")
                    f.write("which python3 \n")
                    f.write("python3 -V \n")
                    #f.write(" \n")
                    f.write("conda activate ana_htoaa \n")
                    #f.write("time python3 %s/%s  %s \n" % (pwd,sAnalysis, sConfig_to_use))
                    f.write("time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 %s/%s  %s \n" % (pwd,sAnalysis, sConfig_to_use))
            
            condor_submit_file = '%s/condor_%s_submit.sh' % (DestinationDir, condorJobName)
            if not os.path.isfile(condor_submit_file):
                with open(condor_submit_file, 'w') as f:
                    f.write("universe = vanilla \n")
                    f.write("executable = %s \n" % condor_exec_file)
                    f.write("getenv = TRUE \n")
                    f.write("log = %s/condor_%s.log \n" % (DestinationDir, condorJobName))
                    f.write("output = %s/condor_%s.out \n" % (DestinationDir, condorJobName))
                    f.write("error = %s/condor_%s.error \n" % (DestinationDir, condorJobName))
                    f.write("notification = never \n")
                    f.write("should_transfer_files = YES \n")
                    f.write("when_to_transfer_output = ON_EXIT \n")
                    f.write("+JobFlavour = \"longlunch\" \n")
                    f.write("queue \n")
                    
            '''
            Job Flavours::
            espresso     = 20 minutes
            microcentury = 1 hour
            longlunch    = 2 hours
            workday      = 8 hours
            tomorrow     = 1 day
            testmatch    = 3 days
            nextweek     = 1 week
            '''
            
            os.system("chmod a+x %s" % condor_exec_file)
            os.system("chmod a+x %s" % condor_submit_file)
            if run_mode == 'condor':
                cmd1 = "condor_submit %s" % condor_submit_file
                print("Now:  %s " % cmd1)                
                os.system(cmd1)
            else:
                pass
                

        

