
'''

'''


import os
import sys
import subprocess
from pathlib import Path
import json
import glob
import argparse
from collections import OrderedDict as OD
import numpy as np
import time
from datetime import datetime
import copy
import enum

print(f"htoaa_Wraper:: here1 {datetime.now() = }")

from htoaa_Settings import *
print(f"htoaa_Wraper:: here2 {datetime.now() = }")
from htoaa_Samples import (
    Samples2018,
    kData, kQCDIncl, kQCD_bGen, kQCD_bEnrich
)
print(f"htoaa_Wraper:: here3 {datetime.now() = }")
from htoaa_CommonTools import (
    executeBashCommand
)
print(f"htoaa_Wraper:: here4 {datetime.now() = }")




#sAnalysis         = "htoaa_Analysis_wCoffea.py"  # "htoaa_Analysis.py"
sConfig           = "config_htoaa.json"
sRunCommandFile   = "1_RunCommand.txt"
sJobSubLogFile    = "1_JobSubmission.log"
sOpRootFile       = "analyze_htoaa_$SAMPLE_$STAGE_$IJOB.root"

printLevel = 2

#UserHomePath = os.path.expanduser("~")
UserHomePath = str(Path.home()) # Python 3.5+
UserName     = os.getlogin()
print(f"htoaa_Wraper:: here5 {datetime.now() = }")


class JobStatus(enum.Enum):
    NotSubmitted  = 'NotSubmitted'
    Finished      = 'Finished'
    Running       = 'Running'
    Failed_Misc   = 'Failed_Misc'
    Failed_Abort  = 'Failed_Abort'
    Failed_XRootD = 'Failed_XRootD'
    

def writeCondorExecFile(
        condor_exec_file,
        sConfig_to_use,
        sOpFile_to_use,
        EosDestinationDir_to_use,
        inpurFiles_to_use,
        server
):
    if not os.path.isfile(condor_exec_file):    
        with open(condor_exec_file, 'w') as f:
            f.write("#!/bin/bash  \n\n")
            #f.write("cd %s \n" % pwd)
            f.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n")
            #f.write("export SCRAM_ARCH=slc6_amd64_gcc700  \n")
            f.write("export SCRAM_ARCH=slc7_amd64_gcc10  \n")
            f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")
            #f.write("cd ")
            #f.write("export X509_USER_PROXY=/afs/cern.ch/user/s/ssawant/x509up_u108989  \n")
            f.write("export X509_USER_PROXY=%s/x509up_u108989  \n" % (UserHomePath)) # voms-proxy-init --rfc --voms cms -valid 192:00 --out ~/x509up_u108989 
            
            # Using x509 proxy without shipping it with the job  https://batchdocs.web.cern.ch/tutorial/exercise2e_proxy.html
            f.write("export X509_USER_PROXY=$1 \n")
            f.write("export EOS_MGM_URL=root://eoscms.cern.ch  \n")
            f.write("voms-proxy-info -all \n")
            f.write("voms-proxy-info -all -file $1 \n")
            
            #f.write("eval \n")
            #f.write("cd %s \n" % (pwd))
            f.write("\nsource %s/.bashrc \n" % (UserHomePath))
            f.write("which conda \n")
            f.write("time conda env list \n")
            myCondaEnv = 'myCondaEnv'
            if UserName.lower() == 'ssawant':
                myCondaEnv = 'ana_htoaa'
            f.write("conda activate %s \n" % (myCondaEnv))
            #f.write("time conda env list \n")

            #f.write("time conda list \n")
            #f.write("which python3 \n")
            #f.write("python3 -V \n")
            #f.write(" \n")
            #f.write("conda activate ana_htoaa \n")
            f.write("\ncp -r %s/* . \n" % (SourceCodeDir))

            f.write("printf \"pwd: \\n\" \n")
            f.write("pwd \n")
            f.write("printf \"ls: \\n\" \n")
            f.write("ls \n")
            f.write("echo \"$1 \" $1 \n")
            f.write("echo \"$2 \" $2 \n")
            
            f.write("\npythonPath=$(which python3) \n")
            #f.write("time python3 %s/%s  %s \n" % (pwd,sAnalysis, sConfig_to_use))
            #f.write("time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 %s/%s  %s \n" % (pwd,sAnalysis, sConfig_to_use))
            #f.write("time ${pythonPath} %s/%s  %s \n" % (pwd,sAnalysis, sConfig_to_use))
            f.write("time ${pythonPath} %s  $2 \n" % (sAnalysis))

            f.write("printf \"After execution pwd: \\n\" \n")
            f.write("pwd \n")
            f.write("printf \"ls: \\n\" \n")
            f.write("ls \n")
            
            cp_commandToUse = ''
            if server in ['lxplus']:
                cp_commandToUse = 'eos cp' # works on lxplus
            else:
                cp_commandToUse = 'cp'
            f.write("time %s %s %s   \n" % (cp_commandToUse, sOpFile_to_use, EosDestinationDir_to_use) )
            #f.write("rm -rf ./inputFiles \n")
            for sInputFile in inpurFiles_to_use:
                sFileLocal = './inputFiles/%s' %(os.path.basename(sInputFile))
                f.write("rm -rf %s \n" % (sFileLocal))
            #f.write(" \n")

        os.system("chmod a+x %s" % condor_exec_file)

    return;


def writeCondorSumitFile(
        condor_submit_file,
        condor_exec_file,
        sCondorLog_to_use,
        sCondorOutput_to_use,
        sCondorError_to_use,
        sIpConfig_to_use,
        increaseJobFlavour=False
):
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

    jobFlavours = OD([
        (0, 'espresso'),
        (1, 'microcentury'),
        (2, 'longlunch'),
        (3, 'workday'),
        (4, 'tomorrow'),
        (5, 'testmatch'),
        (6, 'nextweek'),
    ])
    iJobFlavour = 2 # 2, 'longlunch' 2 hours
    #iJobFlavour = 1 # 1, 'microcentury' 
    if increaseJobFlavour: iJobFlavour += 1
    
    
    #if not os.path.isfile(condor_submit_file):
    with open(condor_submit_file, 'w') as f:
        f.write("universe = vanilla \n")
        
        #f.write("x509userproxy = /afs/cern.ch/user/s/ssawant/x509up_u108989 \n")
        #f.write("use_x509userproxy = true \n")
        f.write("X509_USER_PROXY = %s/x509up_u108989  \n" % (UserHomePath))
        f.write("arguments = $(X509_USER_PROXY) %s \n"  % (sIpConfig_to_use) )       
        
        f.write("executable = %s \n" % condor_exec_file)
        f.write("getenv = TRUE \n")
        f.write("log = %s \n" % (sCondorLog_to_use))
        f.write("output = %s \n" % (sCondorOutput_to_use))
        f.write("error = %s \n" % (sCondorError_to_use))
        f.write("transfer_input_files = $(X509_USER_PROXY), %s \n" % (sIpConfig_to_use) )
        f.write("transfer_output_files = \"\" \n")
        f.write("notification = never \n")
        f.write("should_transfer_files = YES \n")
        f.write("when_to_transfer_output = ON_EXIT \n")
        
        #f.write("+JobFlavour = \"longlunch\" \n")
        f.write("+JobFlavour = \"%s\" \n" % (jobFlavours[iJobFlavour]))
        f.write("queue \n")


    os.system("chmod a+x %s" % condor_submit_file)
    return


def searchStringInFile(sFileName, searchString, nLinesToSearch, SearchFromEnd=True):
    lines = None
    
    with open(sFileName, 'r') as f:
    
        if SearchFromEnd:
            # https://stackoverflow.com/questions/260273/most-efficient-way-to-search-the-last-x-lines-of-a-file
            # f.seek(offset, whence)
            #        offset − This is the position of the read/write pointer within the file.
            #        whence − This is optional and defaults to 0 which means absolute file positioning, other values are 1 which means seek relative to the current position and 2 means seek relative to the file's end.
            f.seek(0, 2)  # Seek @ EOF
            fsize = f.tell() # Get size of file
            f.seek(max(fsize-1024, 0), 0) # Set pos @ last n chars
            lines = f.readlines()  # Read to end

    lines = lines[-1*nLinesToSearch : ] # Get last x lines
    searchStringFound = False
    for line in lines:
        if searchString in line:
            searchStringFound = True
            break

    return searchStringFound
            
        


if __name__ == '__main__':

    print(f"htoaa_Wraper:: here6 {datetime.now() = }"); sys.stdout.flush()
    print("htoaa_Wrapper:: main: {}".format(sys.argv)); sys.stdout.flush()
    
    parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
    parser.add_argument('-analyze',           type=str, default="htoaa_Analysis_GGFMode.py", choices=["htoaa_Analysis_GGFMode.py", "countSumEventsInSample.py"], required=True)
    parser.add_argument('-era', dest='era',   type=str, default=Era_2018,                    choices=[Era_2016, Era_2017, Era_2018], required=False)
    parser.add_argument('-run_mode',          type=str, default='condor',                    choices=['local', 'condor'])
    parser.add_argument('-v', '--version',    type=str, default=None,                        required=True)
    parser.add_argument('-samples',           type=str, default=None,                        help='samples to run seperated by comma')
    parser.add_argument('-excludeSamples',    type=str, default=None,                        help='samples to exclude seperated by comma')
    parser.add_argument('-nFilesPerJob',      type=int, default=5)
    parser.add_argument('-nResubMax',         type=int, default=80)
    parser.add_argument('-ResubWaitingTime',  type=int, default=15,                          help='Resubmit failed jobs after every xx minutes')
    parser.add_argument('-iJobSubmission',    type=int, default=0,                           help='Job submission iteration. Specify previous last job submittion iteration if script terminated for some reason.')
    parser.add_argument('-xrdcpIpAftNResub',  type=int, default=0,                           help='Download input files after n job failures')
    parser.add_argument('-server',            type=str, default='lxplus',                    choices=['lxplus', 'tifr'])
    parser.add_argument('-dryRun',            action='store_true', default=False)    
    args=parser.parse_args()
    print("args: {}".format(args))
    print(f"htoaa_Wraper:: here7 {datetime.now() = }"); sys.stdout.flush()

    sAnalysis               = args.analyze
    era                     = args.era
    run_mode                = args.run_mode
    nFilesPerJob            = args.nFilesPerJob
    selSamplesToRun         = args.samples
    selSamplesToExclude     = args.excludeSamples
    anaVersion              = args.version
    nResubmissionMax        = args.nResubMax
    ResubWaitingTime        = args.ResubWaitingTime
    iJobSubmission          = args.iJobSubmission
    xrdcpIpAftNResub        = args.xrdcpIpAftNResub
    server                  = args.server
    dryRun                  = args.dryRun

    SourceCodeDir     = os.getcwd()
    DestinationDir    = "../analysis/%s/%s" % (anaVersion, era)
    EosDestinationDir = "/eos/cms/store/user/%s/htoaa/analysis/%s/%s" % (UserName, anaVersion, era)

    os.chdir( SourceCodeDir )
    os.makedirs( DestinationDir, exist_ok=True )
    os.chdir( DestinationDir )
    DestinationDirAbsolute = os.getcwd() # save absolute path
    #os.makedirs( DestinationDirAbsolute, exist_ok=True )
    try:
        os.makedirs( EosDestinationDir, exist_ok=True )
    except:
        EosDestinationDir = DestinationDirAbsolute # if /eos area for user is not available then save histograms in DestinationDir

    os.chdir( SourceCodeDir )
    samplesList = None
    samplesInfo = None
    if era == Era_2018:
        samplesList = Samples2018 # htoaa_Samples.py
    with open(sFileSamplesInfo[era]) as fSamplesInfo:
        samplesInfo = json.load(fSamplesInfo) # Samples_Era.json
    selSamplesToRun_list = []
    if selSamplesToRun:
        selSamplesToRun_list = selSamplesToRun.split(',')
    selSamplesToExclude_list = []
    if selSamplesToExclude:
        selSamplesToExclude_list = selSamplesToExclude.split(',')

    ## Settings ---------------------------------------------------------------------------------

    ## MCSamplesStitchOptions.PhSpOverlapRewgt
    MCSamplesStitchOption                     = MCSamplesStitchOptions.PhSpOverlapRewgt 
    samples_wMCSamplesStitch_PhSpOverlapRewgt = [ kQCDIncl, kQCD_bGen, kQCD_bEnrich ]
    ## MCSamplesStitchOptions.PhSpOverlapRemove
    #MCSamplesStitchOption                     = MCSamplesStitchOptions.PhSpOverlapRemove 
    #samples_wMCSamplesStitch_PhSpOverlapRewgt = []

    #  Settings for GGF H->aa->4b analysis
    if sAnalysis in ["htoaa_Analysis_GGFMode.py"]:
        # exclude irrelevant samples from running
        selSamplesToExclude_list.extend( [
                "SUSY_VBFH_HToAATo4B", "SUSY_WH_WToAll_HToAATo4B", "SUSY_ZH_ZToAll_HToAATo4B", "SUSY_TTH_TTToAll_HToAATo4B", 
                "WJetsToLNu", "W1JetsToLNu", "W2JetsToLNu", "W3JetsToLNu", "W4JetsToLNu",
                "TTJets",  # TTJets LO, NLO Madgraph samples 
        ] )
    ## ------------------------------------------------------------------------------------------


    print("\nsamplesList: {}".format(json.dumps(samplesList, indent=4)))
    #print("\n\nsamplesInfo: {}".format(samplesInfo))
    print(f"\n\nselSamplesToRun_list: {selSamplesToRun_list}")
    print(f"selSamplesToExclude_list: {selSamplesToExclude_list}")
        
    sFileRunCommand = "%s/%s" % (DestinationDirAbsolute, sRunCommandFile)
    sFileJobSubLog  = "%s/%s" % (DestinationDirAbsolute, sJobSubLogFile)
    
    # save run command into a .txt tile
    with open(sFileRunCommand, 'a') as fRunCommand:
        datatime_now = datetime.now()
        sCommand = ' '.join(sys.argv)
        fRunCommand.write('%s    %s \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), sCommand))    

    
    jobSubmissionInfo_dict = {}

    allJobsSuccessful          = False
    OpRootFiles_Target         = None
    OpRootFilesAbsPath_Target  = None
    os.chdir( DestinationDirAbsolute )

    while iJobSubmission <= nResubmissionMax:

        print('\n\n%s \t Starting iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))

        OpRootFilesAbsPath_Target  = []
        OpRootFiles_Target         = []
        OpRootFiles_Exist          = []
        OpRootFiles_iJobSubmission = []
        jobStatus_dict             = {} # OD([])
        
        for sample_category, samples in samplesList.items():
            #print("sample_category {}, samples {}".format(sample_category, samples))
            sample_isMC = True
            for sampleSubString_toCheck in [kData, 'Run2016','Run2017', 'Run2018']:
                if sampleSubString_toCheck in sample_category:
                    sample_isMC = False
                    break

            for sample in samples:
                if len(selSamplesToRun_list) > 0:
                    skipThisSample = True
                    for selSample in selSamplesToRun_list:
                        if sample.startswith(selSample): skipThisSample = False
                    if skipThisSample:
                        continue

                if len(selSamplesToExclude_list) > 0:
                    skipThisSample = False
                    for selSample in selSamplesToExclude_list:
                        if sample.startswith(selSample): skipThisSample = True
                    if skipThisSample:
                        continue
                    
                #
                OpRootFileFinalDir = '%s/%s' % (EosDestinationDir, sample)
                JobLogsDir         = '%s/%s' % (DestinationDirAbsolute, sample)
                if not os.path.exists(OpRootFileFinalDir): os.makedirs( OpRootFileFinalDir, exist_ok=True )
                if not os.path.exists(JobLogsDir):         os.makedirs( JobLogsDir, exist_ok=True )
                os.chdir( JobLogsDir )
                    
                print(f"sample_category: {sample_category}, sample: {sample}")

                sampleInfo = samplesInfo[sample] # Samples_Era.json                
                fileList = sampleInfo[sampleFormat]
                files = []
                for iEntry in fileList:
                    # file name with wildcard charecter *
                    if "*" in iEntry:  files.extend( glob.glob( iEntry ) )
                    else:              files.append( iEntry )
                sample_dataset     = sampleInfo["dataset"]
                sample_cossSection = sampleInfo["cross_section"] if sample_isMC else None
                sample_nEvents     = sampleInfo["nEvents"]
                sample_sumEvents   = sampleInfo["sumEvents"] if sample_isMC else None

                if printLevel >= 6:
                    print("\nsample: {}".format(sample))
                    print("samplesInfo[{}]: {}".format(sample, samplesInfo[sample]))
                    print("files ({}): {}".format(len(files), files))


                nSplits = int(len(files) / nFilesPerJob) + 1 if (nFilesPerJob > 0) and (len(files) != nFilesPerJob) else 1


                files_splitted = np.array_split(files, nSplits)
                if printLevel >= 6:
                    print("files_splitted: {}".format(files_splitted))

                for iJob in range(len(files_splitted)):
                    JobStage = 0
                    
                    config = copy.deepcopy(config_Template)

                    # Job related files
                    #sOpRootFile_to_use      = '%s/%s' % (DestinationDir, sOpRootFile)
                    sOpRootFile_to_use      = '%s' % (sOpRootFile)
                    sOpRootFile_to_use      = sOpRootFile_to_use.replace('$SAMPLE', sample)
                    sOpRootFile_to_use      = sOpRootFile_to_use.replace('$STAGE', str(JobStage))
                    sOpRootFile_to_use      = sOpRootFile_to_use.replace('$IJOB', str(iJob))
                    sOpRootFileFinal_to_use = '%s/%s' % (OpRootFileFinalDir, sOpRootFile_to_use)
                    
                    sConfig_to_use          = sOpRootFile_to_use.replace('.root', '_config.json')
                    sCondorExec_to_use      = sOpRootFile_to_use.replace('.root', '_condor_exec.sh')
                    sCondorSubmit_to_use    = sOpRootFile_to_use.replace('.root', '_condor_submit.sh')
                    sCondorLog_to_use       = sOpRootFile_to_use.replace('.root', '_condor.log')
                    sCondorOutput_to_use    = sOpRootFile_to_use.replace('.root', '_condor.out')
                    sCondorError_to_use     = sOpRootFile_to_use.replace('.root', '_condor.error')

                    # Check if job related file exist or not
                    isConfigExist           = os.path.isfile(sConfig_to_use)
                    isOpRootFileExist       = os.path.isfile(sOpRootFileFinal_to_use)
                    isCondorExecExist       = os.path.isfile(sCondorExec_to_use)
                    isCondorSubmitExist     = os.path.isfile(sCondorSubmit_to_use)
                    isCondorLogExist        = os.path.isfile(sCondorLog_to_use)
                    isCondorOutputExist     = os.path.isfile(sCondorOutput_to_use)
                    isCondorErrorExist      = os.path.isfile(sCondorError_to_use)

                    if printLevel >= 3:
                        print(f"sOpRootFile_to_use: {sOpRootFile_to_use} ")
                        #print(f" {sConfig_to_use = }: {isConfigExist = } ")
                        #print(f" {sOpRootFileFinal_to_use = }: {isOpRootFileExist = } ")
                        #print(f" {sCondorExec_to_use = }: {isCondorExecExist = } ")
                        #print(f" {sCondorSubmit_to_use = }: {isCondorSubmitExist = } ")
                        #print(f" {sCondorLog_to_use = }: {isCondorLogExist = } ")
                        #print(f" {sCondorOutput_to_use = }: {isCondorOutputExist = } ")
                        #print(f" {sCondorError_to_use = }: {isCondorErrorExist = } ")
                    
                    # JobStatus
                    jobStatus = JobStatus.NotSubmitted # -1
                    #jobStatusForJobSubmission = [0, 3, 4, 5]
                    jobStatusForJobSubmission = [
                        JobStatus.NotSubmitted, #0
                        #JobStatus.Finished, #1
                        #JobStatus.Running, #2
                        JobStatus.Failed_Misc, #3
                        JobStatus.Failed_Abort, #4
                        JobStatus.Failed_XRootD, #5
                    ]

                    if not isConfigExist:
                        jobStatus = JobStatus.NotSubmitted #0 # job not yet submitted
                        if printLevel >= 3:
                            print(f"  jobStatus = 0")

                    elif isOpRootFileExist:
                        jobStatus = JobStatus.Finished #1 # job ran successfully
                        OpRootFiles_Exist.append(sOpRootFile_to_use)
                        if printLevel >= 3:
                            print(f"  jobStatus = 1")
                            
                    else:
                        if isCondorLogExist:
                            
                            if (searchStringInFile(                                    
                                    sFileName       = sCondorLog_to_use,
                                    searchString    = 'Job terminated',
                                    nLinesToSearch  = 3,
                                    SearchFromEnd   = True)):
                                # check wheter the job was terminated or not
                                jobStatus = JobStatus.Failed_Misc #3 # job failed due to some other error
                                if printLevel >= 3:
                                    print(f"  jobStatus = 3")
                                    
                                    
                                # check if job failed due to XRootD error
                                if (searchStringInFile(                                        
                                        sFileName       = sCondorError_to_use,
                                        searchString    = 'OSError: XRootD error: [ERROR]', 
                                        nLinesToSearch  = 150,
                                        SearchFromEnd   = True) or \
                                    searchStringInFile(
                                        sFileName       = sCondorError_to_use,
                                        searchString    = '[ERROR] Invalid redirect URL', 
                                        nLinesToSearch  = 150,
                                        SearchFromEnd   = True) ):
                                    jobStatus = JobStatus.Failed_XRootD #5 # job failed due to XRootD error
                                    if printLevel >= 3:
                                        print(f"  jobStatus = 5")
                                        

                            elif (searchStringInFile(
                                    sFileName       = sCondorLog_to_use,
                                    searchString    = 'Job was aborted',
                                    nLinesToSearch  = 3,
                                    SearchFromEnd   = True)):
                                # check wheter sCondorError does not exist due to Job was aborted
                                jobStatus = JobStatus.Failed_Abort #4 # job aborted
                                if printLevel >= 3:
                                    print(f"  jobStatus = 4")
                                    
                                    
                            else:
                                jobStatus = JobStatus.Running #2 # job is running
                                if printLevel >= 3:
                                    print(f"  jobStatus = 2")
                                
                                
                    OpRootFiles_Target.append(sOpRootFile_to_use)
                    OpRootFilesAbsPath_Target.append(sOpRootFileFinal_to_use)
                    if jobStatus in jobStatusForJobSubmission : # [0, 3, 4]:
                        OpRootFiles_iJobSubmission.append(sOpRootFile_to_use)

                    if jobStatus not in jobStatus_dict.keys():
                        jobStatus_dict[jobStatus] = [sOpRootFile_to_use]
                    else:
                        jobStatus_dict[jobStatus].append(sOpRootFile_to_use)
                    
                    if sOpRootFile_to_use not in jobSubmissionInfo_dict:
                        jobSubmissionInfo_dict[sOpRootFile_to_use] = {}
                        jobSubmissionInfo_dict[sOpRootFile_to_use]['nResubmissions'] = 0
                    else:
                        if jobStatus in jobStatusForJobSubmission :
                            jobSubmissionInfo_dict[sOpRootFile_to_use]['nResubmissions'] = jobSubmissionInfo_dict[sOpRootFile_to_use]['nResubmissions'] + 1
                    jobSubmissionInfo_dict[sOpRootFile_to_use]['JobStatusLast'] = jobStatus
                        
                    
                    
                        
                    if printLevel >= 0:
                        #print(f"\t {sOpRootFile_to_use}:: {jobStatus}, Config: {isConfigExist}, OpRootFile: {isOpRootFileExist}, CondorExec: {isCondorExecExist}, CondorSubmit: {isCondorSubmitExist}, CondorLog: {isCondorLogExist}, CondorOutput: {isCondorOutputExist}, CondorError: {isCondorErrorExist}"); sys.stdout.flush()
                        print(f"\t {sOpRootFile_to_use}:: {jobStatus}, Config: {isConfigExist}, OpRootFile: {isOpRootFileExist},  CondorLog: {isCondorLogExist}, CondorOutput: {isCondorOutputExist}, CondorError: {isCondorErrorExist}"); sys.stdout.flush()
                        

                    #if iJobSubmission == 0:
                    #if jobStatus == 0 or 1==1:
                    if jobStatus in jobStatusForJobSubmission:
                        config["era"] = era
                        config["dataset"]    = sample_dataset
                        #config["dataset"]    = list( sample_dataset )
                        config["inputFiles"] = list( files_splitted[iJob] )
                        config["outputFile"] = sOpRootFile_to_use 
                        config["sampleCategory"] = sample_category
                        config["isMC"] = sample_isMC 
                        config["nEvents"] = sample_nEvents
                        if sample_isMC:
                            config["crossSection"] = sample_cossSection
                            config["sumEvents"] = sample_sumEvents
                            
                            if MCSamplesStitchOption == MCSamplesStitchOptions.PhSpOverlapRewgt and \
                               sample_category in samples_wMCSamplesStitch_PhSpOverlapRewgt:
                                # MCSamplesStitch_PhSpOverlapRewgt: Read lumiScale from histogram saved in a ROOT file 
                                config["MCSamplesStitchOption"] = MCSamplesStitchOptions.PhSpOverlapRewgt.value 
                                config["MCSamplesStitchInputs"] = sFileLumiScalesPhSpOverlapRewgt[era]
                                
                        else:
                            del config["crossSection"]
                            del config["sumEvents"]
                        config["downloadIpFiles"] = True if jobSubmissionInfo_dict[sOpRootFile_to_use]['nResubmissions'] >= xrdcpIpAftNResub else False
                        config["server"] = server

                        if printLevel >= 4:
                            print("config {}: {}".format(sConfig_to_use, config))
                        with open(sConfig_to_use, "w") as fConfig:
                            json.dump( config,  fConfig, indent=4)


                        writeCondorExecFile(
                            sCondorExec_to_use,
                            sConfig_to_use,
                            sOpRootFile_to_use,
                            OpRootFileFinalDir,
                            config["inputFiles"],
                            server 
                        )


                    if jobStatus in jobStatusForJobSubmission: #[0, 3, 4]:
                        if jobStatus == [JobStatus.Failed_Misc, JobStatus.Failed_XRootD]: #[3, 5]:
                            # save previos .out and .error files with another names
                            sCondorOutput_vPrevious = sCondorOutput_to_use.replace('.out', '_v%d.out' % (iJobSubmission-1))
                            sCondorError_vPrevious  = sCondorError_to_use.replace('.error', '_v%d.error' % (iJobSubmission-1))
                            os.rename(sCondorOutput_to_use, sCondorOutput_vPrevious)
                            os.rename(sCondorError_to_use,  sCondorError_vPrevious)

                        increaseJobFlavour = False
                        #if jobStatus == 4 or jobSubmissionInfo_dict[sOpRootFile_to_use]['nResubmissions'] >= xrdcpIpAftNResub:
                        if jobStatus == JobStatus.Failed_Abort or jobSubmissionInfo_dict[sOpRootFile_to_use]['nResubmissions'] >= xrdcpIpAftNResub:
                            increaseJobFlavour = True
                            
                        writeCondorSumitFile(
                            sCondorSubmit_to_use,
                            sCondorExec_to_use,
                            sCondorLog_to_use,
                            sCondorOutput_to_use,
                            sCondorError_to_use,
                            sConfig_to_use,
                            increaseJobFlavour)



                    if jobStatus in [JobStatus.Finished, JobStatus.Running]: #[1, 2]:
                        # job is either running or succeeded
                        continue

                    '''
                    if jobStatus in [3]:
                        # job failed, but failure reason needs investigation
                        continue
                    '''
                    
                    if run_mode == 'condor':
                        cmd1 = "condor_submit %s" % sCondorSubmit_to_use 
                        
                        if not dryRun:
                            if printLevel >= 5:
                                print("Now:  %s " % cmd1)
                            os.system(cmd1)
                    else:
                        pass

                    


        # write JobSubmission status report in JobSubLog file
        with open(sFileJobSubLog, 'a') as fJobSubLog:
            fJobSubLog.write('%s \t iJobSubmission %d \t OpRootFiles_Target (%d):  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission, len(OpRootFiles_Target)))
            if iJobSubmission == 0:
                for f in OpRootFiles_Target:
                    fJobSubLog.write('\t %s \n' % (f))
            else:
                fJobSubLog.write('OpRootFiles_Exist %d out of %d. \n' % (len(OpRootFiles_Exist), len(OpRootFiles_Target)))
                fJobSubLog.write('OpRootFiles_iJobSubmission (%d): ' % (len(OpRootFiles_iJobSubmission)))
                for f in OpRootFiles_iJobSubmission:
                    fJobSubLog.write('\t %s \n' % (f))

                fJobSubLog.write('\n\nJob status wise output files: \n')
                for jobStatus in jobStatus_dict.keys():
                    fJobSubLog.write('\t jobStatus %s (%d) \n' % (str(jobStatus.value), len(jobStatus_dict[jobStatus])))
                    #if jobStatus in [0, 1]: continue
                    if jobStatus in [JobStatus.NotSubmitted, JobStatus.Finished]: continue
                    
                    for f in jobStatus_dict[jobStatus]:
                        fJobSubLog.write('\t\t %s \n' % (f))
                
            fJobSubLog.write('%s\n\n\n' % ('-'*10))
        
        
        jobStatus_list = [ (jobStatus.value, len(jobStatus_dict[jobStatus])) for jobStatus in jobStatus_dict.keys() ]
        print('\n\n\n%s \t iJobSubmission %d \t OpRootFiles_Exist %d out of %d. No. of jobs submitted in this resubmission: %d:  ' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission, len(OpRootFiles_Exist), len(OpRootFiles_Target), len(OpRootFiles_iJobSubmission)))
        print(f"jobStatus_list: {jobStatus_list} \n"); sys.stdout.flush()
        
            
        if dryRun:
            print('%s \t druRun with iJobSubmission: %d  \nTerminating...\n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
            exit(0)
            
        if len(OpRootFiles_Target) == len(OpRootFiles_Exist):
            allJobsSuccessful = True
            break
        else:
            time.sleep( ResubWaitingTime * 60 )
            iJobSubmission += 1


    fJobSubLog = open(sFileJobSubLog, 'a')
    fJobSubLog.write('%s \t Jobs are done. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
    print('%s \t Jobs are done. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))

    if allJobsSuccessful:
        print('%s \t All jobs run successfully. Now hadd root files.  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S")))
        os.chdir( EosDestinationDir )

        sOpRootFile_stage0 = sOpRootFile
        sOpRootFile_stage0 = sOpRootFile_stage0.replace('_$SAMPLE',  '')
        sOpRootFile_stage0 = sOpRootFile_stage0.replace('_$STAGE',   '')
        sOpRootFile_stage0 = sOpRootFile_stage0.replace('_$IJOB',    '')
        sOpRootFile_stage0 = sOpRootFile_stage0.replace('.root',     '*.root')
        
        sOpRootFile_stage1 = sOpRootFile
        sOpRootFile_stage1 = sOpRootFile_stage1.replace('_$SAMPLE',  '')
        sOpRootFile_stage1 = sOpRootFile_stage1.replace('_$STAGE',   '_stage1')
        sOpRootFile_stage1 = sOpRootFile_stage1.replace('_$IJOB',    '')

        nFilesPerBatchForHadd              = 100
        nBatchesForHadd                    = int(len(OpRootFilesAbsPath_Target) / nFilesPerBatchForHadd) + 1 if len(OpRootFilesAbsPath_Target) != nFilesPerBatchForHadd else 1
        OpRootFilesAbsPath_Target_splitted = np.array_split(OpRootFilesAbsPath_Target, nBatchesForHadd)
        sOpRootFile_stage1_batches         = []
        #print(f"\n\nNo. of {nFilesPerBatchForHadd} files splits in OpRootFilesAbsPath_Target_splitted: {len(OpRootFilesAbsPath_Target_splitted)}")
        print(f"\n\n{nFilesPerBatchForHadd = }. {len(OpRootFilesAbsPath_Target)} files split into {nBatchesForHadd} batches as {[len(iL) for iL in OpRootFilesAbsPath_Target_splitted]}")
        #print(f"{OpRootFilesAbsPath_Target_splitted = }")
        for iHadd in range(len(OpRootFilesAbsPath_Target_splitted)):
            sOpRootFile_stage1_toUse        = sOpRootFile_stage1.replace('.root', '_batch%d.root' % (iHadd))
            OpRootFilesAbsPath_Target_toUse = OpRootFilesAbsPath_Target_splitted[iHadd]
            sOpRootFile_stage1_batches.append(sOpRootFile_stage1_toUse)
            print(f"\n\n{iHadd = }, No. of files to hadd: {len(OpRootFilesAbsPath_Target_toUse)}")

            cmd_hadd = "time hadd -f %s" % (sOpRootFile_stage1_toUse)
            for opFileName in OpRootFilesAbsPath_Target_toUse:
                cmd_hadd += " %s" % (opFileName)
            
            cmd_hadd_stdout = executeBashCommand(cmd_hadd)
            fJobSubLog.write('\n\n{iHadd = } \t %s: \n%s \n' % (cmd_hadd, cmd_hadd_stdout))
            fJobSubLog.write('\n\n{iHadd = } \t %s: hadd %s is done.' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), sOpRootFile_stage1_toUse))

        # Now hadd sOpRootFile_stage1_batches to sOpRootFile_stage1
        print(f"\n\n\n Now hadd sOpRootFile_stage1_batches to sOpRootFile_stage1 ")

        cmd_hadd = "time hadd -f %s" % (sOpRootFile_stage1)
        for opFileName in sOpRootFile_stage1_batches:
            cmd_hadd += " %s" % (opFileName)
        
        cmd_hadd_stdout = executeBashCommand(cmd_hadd)
        fJobSubLog.write('\n\nHadd stage1_batches %s: \n%s \n' % (cmd_hadd, cmd_hadd_stdout))
        fJobSubLog.write('\n\n%s:add stage1_batches: hadd %s is done.' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), sOpRootFile_stage0))



        executeBashCommand("pwd")
        executeBashCommand("ls")

    fJobSubLog.close()
