
import os
import sys
from pathlib import Path
import json
import glob
import argparse
from collections import OrderedDict as OD
import numpy as np
import time
from datetime import datetime
import copy

from htoaa_Settings import *
from htoaa_Samples import (
    Samples2018,
    kData
)





sAnalysis         = "htoaa_Analysis_wCoffea.py"  # "htoaa_Analysis.py"
sConfig           = "config_htoaa.json"
sRunCommandFile   = "1_RunCommand.txt"
sJobSubLogFile    = "1_JobSubmission.log"
sOpRootFile       = "analyze_htoaa_$SAMPLE_$STAGE_$IJOB.root"

printLevel = 3

#UserHomePath = os.path.expanduser("~")
UserHomePath = str(Path.home()) # Python 3.5+

def writeCondorExecFile(condor_exec_file, sConfig_to_use):
    if not os.path.isfile(condor_exec_file):    
        with open(condor_exec_file, 'w') as f:
            f.write("#!/bin/bash  \n\n")
            f.write("cd %s \n" % pwd)
            f.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n")
            f.write("export SCRAM_ARCH=slc6_amd64_gcc700  \n")
            f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")
            #f.write("cd ")
            #f.write("export X509_USER_PROXY=/afs/cern.ch/user/s/ssawant/x509up_u108989  \n")
            f.write("export X509_USER_PROXY=%s/x509up_u108989  \n" % (UserHomePath))

            # Using x509 proxy without shipping it with the job  https://batchdocs.web.cern.ch/tutorial/exercise2e_proxy.html
            #f.write("export X509_USER_PROXY=$1 \n")
            #f.write("voms-proxy-info -all \n")
            #f.write("voms-proxy-info -all -file $1 \n")
            
            #f.write("eval \n")
            f.write("cd %s \n" % (pwd))
            f.write("source %s/.bashrc \n" % (UserHomePath))
            #f.write("which conda \n")
            #f.write("time conda env list \n")
            f.write("conda activate ana_htoaa \n")
            #f.write("time conda env list \n")

            #f.write("time conda list \n")
            #f.write("which python3 \n")
            #f.write("python3 -V \n")
            #f.write(" \n")
            #f.write("conda activate ana_htoaa \n")
            f.write("pythonPath=$(which python3) \n")
            #f.write("time python3 %s/%s  %s \n" % (pwd,sAnalysis, sConfig_to_use))
            #f.write("time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/ana_htoaa/bin/python3 %s/%s  %s \n" % (pwd,sAnalysis, sConfig_to_use))
            f.write("time ${pythonPath} %s/%s  %s \n" % (pwd,sAnalysis, sConfig_to_use))

        os.system("chmod a+x %s" % condor_exec_file)

    return;


def writeCondorSumitFile(condor_submit_file, condor_exec_file, sCondorLog_to_use, sCondorOutput_to_use, sCondorError_to_use, increaseJobFlavour=False):
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
    #iJobFlavour = 2 # 2, 'longlunch' 2 hours
    iJobFlavour = 1 # 1, 'microcentury' 
    if increaseJobFlavour: iJobFlavour += 1
    
    
    #if not os.path.isfile(condor_submit_file):
    with open(condor_submit_file, 'w') as f:
        f.write("universe = vanilla \n")
        
        #f.write("x509userproxy = /afs/cern.ch/user/s/ssawant/x509up_u108989 \n")
        #f.write("use_x509userproxy = true \n")
        f.write("X509_USER_PROXY = %s/x509up_u108989  \n" % (UserHomePath))
        f.write("arguments = $(X509_USER_PROXY)  \n")        
        
        f.write("executable = %s \n" % condor_exec_file)
        f.write("getenv = TRUE \n")
        f.write("log = %s \n" % (sCondorLog_to_use))
        f.write("output = %s \n" % (sCondorOutput_to_use))
        f.write("error = %s \n" % (sCondorError_to_use))
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

    print("htoaa_Wrapper:: main: {}".format(sys.argv)); sys.stdout.flush()
    
    parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
    parser.add_argument('-era', dest='era',   type=str, default=Era_2018, choices=[Era_2016, Era_2017, Era_2018], required=False)
    parser.add_argument('-run_mode',          type=str, default='condor', choices=['local', 'condor'])
    parser.add_argument('-v', '--version',    type=str, default=None, required=True)
    parser.add_argument('-samples',           type=str, default=None, help='samples to run seperated by comma')
    parser.add_argument('-nFilesPerJob',      type=int, default=10)
    parser.add_argument('-nResubMax',         type=int, default=80)
    parser.add_argument('-ResubWaitingTime',  type=int, default=15, help='Resubmit failed jobs after every xx minutes')
    parser.add_argument('-iJobSubmission',    type=int, default=0,  help='Job submission iteration. Specify previous last job submittion iteration if script terminated for some reason.')
    parser.add_argument('-dryRun',            action='store_true', default=False)
    args=parser.parse_args()
    print("args: {}".format(args))

    era              = args.era
    run_mode         = args.run_mode
    nFilesPerJob     = args.nFilesPerJob
    selSamplesToRun  = args.samples
    anaVersion       = args.version
    nResubmissionMax = args.nResubMax
    ResubWaitingTime = args.ResubWaitingTime
    iJobSubmission   = args.iJobSubmission
    dryRun           = args.dryRun

    pwd = os.getcwd()
    DestinationDir = "../analysis/%s" % (anaVersion)
    sFileRunCommand = "%s/%s" % (DestinationDir, sRunCommandFile)
    sFileJobSubLog  = "%s/%s" % (DestinationDir, sJobSubLogFile)
    if not os.path.exists(DestinationDir): os.mkdir( DestinationDir )

    # save run command into a .txt tile
    with open(sFileRunCommand, 'a') as fRunCommand:
        datatime_now = datetime.now()
        sCommand = ' '.join(sys.argv)
        fRunCommand.write('%s    %s \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), sCommand))



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
    #print("\n\nsamplesInfo: {}".format(samplesInfo))
    print(f"\n\nselSamplesToRun_list: {selSamplesToRun_list}")


    #config = config_Template


    #iJobSubmission = 0
    
    while iJobSubmission <= nResubmissionMax:

        print('\n\n%s \t Startiing iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))

        OpRootFiles_Target         = []
        OpRootFiles_Exist          = []
        OpRootFiles_iJobSubmission = []
        jobStatus_dict             = {} # OD([])
        
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
                sample_cossSection = sampleInfo["cross_section"] if (sample_category != kData) else None
                sample_nEvents     = sampleInfo["nEvents"]
                sample_sumEvents   = sampleInfo["sumEvents"] if (sample_category != kData) else None

                if printLevel >= 6:
                    print("\nsample: {}".format(sample))
                    print("samplesInfo[sample]: {}".format(samplesInfo[sample]))
                    print("files ({}): {}".format(len(files), files))


                nSplits = int(len(files) / nFilesPerJob) + 1 if nFilesPerJob > 0 else 1


                files_splitted = np.array_split(files, nSplits)
                if printLevel >= 6:
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
                    #config = config_Template.deepcopy()
                    config = copy.deepcopy(config_Template)

                    # Job related files
                    sOpRootFile_to_use    = '%s/%s' % (DestinationDir, sOpRootFile)
                    sOpRootFile_to_use    = sOpRootFile_to_use.replace('$SAMPLE', sample)
                    sOpRootFile_to_use    = sOpRootFile_to_use.replace('$STAGE', str(0))
                    sOpRootFile_to_use    = sOpRootFile_to_use.replace('$IJOB', str(iJob))

                    sConfig_to_use        = sOpRootFile_to_use.replace('.root', '_config.json')
                    sCondorExec_to_use    = sOpRootFile_to_use.replace('.root', '_condor_exec.sh')
                    sCondorSubmit_to_use  = sOpRootFile_to_use.replace('.root', '_condor_submit.sh')
                    sCondorLog_to_use     = sOpRootFile_to_use.replace('.root', '_condor.log')
                    sCondorOutput_to_use  = sOpRootFile_to_use.replace('.root', '_condor.out')
                    sCondorError_to_use   = sOpRootFile_to_use.replace('.root', '_condor.error')

                    # Check if job related file exist or not
                    isConfigExist        = os.path.isfile(sConfig_to_use)
                    isOpRootFileExist    = os.path.isfile(sOpRootFile_to_use)
                    isCondorExecExist    = os.path.isfile(sCondorExec_to_use)
                    isCondorSubmitExist  = os.path.isfile(sCondorSubmit_to_use)
                    isCondorLogExist     = os.path.isfile(sCondorLog_to_use)
                    isCondorOutputExist  = os.path.isfile(sCondorOutput_to_use)
                    isCondorErrorExist   = os.path.isfile(sCondorError_to_use)

                    if printLevel >= 3:
                        print(f"sOpRootFile_to_use: {sOpRootFile_to_use} ")
                    
                    # JobStatus
                    jobStatus = -1
                    jobStatusForJobSubmission = [0, 3, 4, 5]

                    if not isConfigExist:
                        jobStatus = 0 # job not yet submitted
                        if printLevel >= 3:
                            print(f"  jobStatus = 0")

                    elif isOpRootFileExist:
                        jobStatus = 1 # job ran successfully
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
                                jobStatus = 3 # job failed due to some other error
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
                                    jobStatus = 5 # job failed due to XRootD error
                                    if printLevel >= 3:
                                        print(f"  jobStatus = 5")
                                        

                            elif (searchStringInFile(
                                    sFileName       = sCondorLog_to_use,
                                    searchString    = 'Job was aborted',
                                    nLinesToSearch  = 3,
                                    SearchFromEnd   = True)):
                                # check wheter sCondorError does not exist due to Job was aborted
                                jobStatus = 4 # job aborted
                                if printLevel >= 3:
                                    print(f"  jobStatus = 4")
                                    
                                    
                            else:
                                jobStatus = 2 # job is running
                                if printLevel >= 3:
                                    print(f"  jobStatus = 2")
                                
                                

                    OpRootFiles_Target.append(sOpRootFile_to_use)
                    if jobStatus in jobStatusForJobSubmission : # [0, 3, 4]:
                        OpRootFiles_iJobSubmission.append(sOpRootFile_to_use)

                    if jobStatus not in jobStatus_dict.keys():
                        jobStatus_dict[jobStatus] = [sOpRootFile_to_use]
                    else:
                        jobStatus_dict[jobStatus].append(sOpRootFile_to_use)
                    

                    if printLevel >= 0:
                        print(f"\t {sOpRootFile_to_use}:: jobStatus: {jobStatus}, isConfigExist: {isConfigExist}, isOpRootFileExist: {isOpRootFileExist}, isCondorExecExist: {isCondorExecExist}, isCondorSubmitExist: {isCondorSubmitExist}, isCondorLogExist: {isCondorLogExist}, isCondorOutputExist: {isCondorOutputExist}, isCondorErrorExist: {isCondorErrorExist}"); sys.stdout.flush()
                        

                    #if iJobSubmission == 0:
                    if jobStatus == 0:
                        config["era"] = era
                        config["inputFiles"] = list( files_splitted[iJob] )
                        config["outputFile"] = sOpRootFile_to_use 
                        config["sampleCategory"] = sample_category
                        config["isMC"] = (sample_category != kData)
                        #config["Luminosity"] = Luminosity
                        config["nEvents"] = sample_nEvents
                        if (sample_category != kData):
                            config["crossSection"] = sample_cossSection
                            config["sumEvents"] = sample_sumEvents
                        else:
                            del config["crossSection"]
                            del config["sumEvents"]


                        print("config {}: {}".format(sConfig_to_use, config))
                        with open(sConfig_to_use, "w") as fConfig:
                            json.dump( config,  fConfig, indent=4)


                        writeCondorExecFile(sCondorExec_to_use, sConfig_to_use)


                    if jobStatus in jobStatusForJobSubmission: #[0, 3, 4]:
                        if jobStatus == [3, 5]:
                            # save previos .out and .error files with another names
                            sCondorOutput_vPrevious = sCondorOutput_to_use.replace('.out', '_v%d.out' % (iJobSubmission-1))
                            sCondorError_vPrevious  = sCondorError_to_use.replace('.error', '_v%d.error' % (iJobSubmission-1))
                            os.rename(sCondorOutput_to_use, sCondorOutput_vPrevious)
                            os.rename(sCondorError_to_use,  sCondorError_vPrevious)

                        increaseJobFlavour = False
                        if jobStatus == 4:
                            increaseJobFlavour = True
                            
                        writeCondorSumitFile(sCondorSubmit_to_use, sCondorExec_to_use, sCondorLog_to_use, sCondorOutput_to_use, sCondorError_to_use,
                                             increaseJobFlavour)



                    if jobStatus in [1, 2]:
                        # job is either running or succeeded
                        continue

                    '''
                    if jobStatus in [3]:
                        # job failed, but failure reason needs investigation
                        continue
                    '''
                    
                    if run_mode == 'condor':
                        cmd1 = "condor_submit %s" % sCondorSubmit_to_use 
                        print("Now:  %s " % cmd1)
                        if not dryRun:
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
                    fJobSubLog.write('\t jobStatus %d (%d) \n' % (jobStatus, len(jobStatus_dict[jobStatus])))
                    if jobStatus in [0, 1]: continue
                    
                    for f in jobStatus_dict[jobStatus]:
                        fJobSubLog.write('\t\t %s \n' % (f))
                
            fJobSubLog.write('%s\n\n\n' % ('-'*10))
        

        jobStatus_list = [ (jobStatus, len(jobStatus_dict[jobStatus])) for jobStatus in jobStatus_dict.keys() ]
        print('\n\n\n%s \t iJobSubmission %d \t OpRootFiles_Exist %d out of %d. No. of jobs submitted in this resubmission: %d:  ' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission, len(OpRootFiles_Exist), len(OpRootFiles_Target), len(OpRootFiles_iJobSubmission)))
        print(f"jobStatus_list: {jobStatus_list} \n"); sys.stdout.flush()
        
            
        if dryRun:
            print('%s \t druRun with iJobSubmission: %d  \nTerminating...\n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
            exit(0)
            
        if len(OpRootFiles_Target) == len(OpRootFiles_Exist):
            break
        else:
            time.sleep( ResubWaitingTime * 60 )
            iJobSubmission += 1


    with open(sFileJobSubLog, 'a') as fJobSubLog:
        fJobSubLog.write('%s \t Jobs are done. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))
    print('%s \t Jobs are done. iJobSubmission: %d  \n' % (datetime.now().strftime("%Y/%m/%d %H:%M:%S"), iJobSubmission))

