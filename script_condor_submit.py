#!/usr/bin/env python

'''

''' 



from collections import OrderedDict as OD
import json
import os
import subprocess
#import ROOT
#ROOT.PyConfig.IgnoreCommandLineOptions = True
import time
import datetime
import sys
import glob


#command0 = "source /afs/cern.ch/user/s/ssawant/.bashrc && time conda env list && time conda activate mlemv_fromSi && time conda env list && python --version && python -V && scl enable rh-python36 bash && python -V && time conda list && time python htoaa_BDT3.py "
#command0 = "source /afs/cern.ch/user/s/ssawant/.bashrc && time conda env list && time conda activate mlemv_fromSi && time conda env list && python --version && python3 -V && which python && which python3 && time conda list && time python3 htoaa_BDT3.py "
#command0 = "source /afs/cern.ch/user/s/ssawant/.bashrc && time conda env list && time conda activate mlemv_fromSi && time conda env list && python --version && python3 -V && time conda list && time python htoaa_BDT3.py "
#command0 = "source /afs/cern.ch/user/s/ssawant/.bashrc && time conda env list && time conda activate mlemv_fromSi && source activate /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/mlemv_fromSi && time conda env list && python --version && python3 -V && which python && which python3 && time conda list && echo $PATH && time python3 htoaa_BDT3.py "
#command0 = "source /afs/cern.ch/user/s/ssawant/.bashrc && time conda env list && time source /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/bin/activate mlemv_fromSi && time conda env list && python --version && python3 -V && which python && which python3 && time conda list && which conda && time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/mlemv_fromSi/bin/python3 htoaa_BDT3.py "
#command0 = "source /afs/cern.ch/user/s/ssawant/.bashrc && time conda env list && time source /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/bin/activate mlemv_fromSi && time conda env list && python --version && python3 -V && which python && which python3 && time conda list && which conda && time python3 htoaa_BDT3.py " # didn't work
#command0 = "source /afs/cern.ch/user/s/ssawant/.bashrc && time conda env list && time conda activate mlemv_fromSi && time conda env list && python --version && python3 -V && which python && which python3 && time conda list && which conda && time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/mlemv_fromSi/bin/python3 htoaa_BDT3.py "
commands = [
    "source /afs/cern.ch/user/s/ssawant/.bashrc",
    "time conda env list",
    "time conda activate mlemv_fromSi",
    "time conda env list",
    "python --version",
    "python3 -V",
    "which python",
    "which python3",
    "time conda list",
    "which conda",
    "time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/mlemv_fromSi/bin/python3 htoaa_BDT3.py"
]
scriptName0 = "htoaa_BDT3"



pwd = os.getcwd()

    
condor_exec_file = 'condor_exec_%s.sh' % (scriptName0)

with open(condor_exec_file, 'w') as f:
    f.write("#!/bin/bash  \n\n")

    if "t3storage3" in pwd:
        f.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n")
        f.write("export SCRAM_ARCH=slc6_amd64_gcc700  \n")
        f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")
        #f.write("cd ")
        f.write("export X509_USER_PROXY=/home/ssawant/x509up_u56558 \n")
    elif "cern" in pwd:
        #f.write("unset PYTHONPATH \n")
        #f.write("unset PYTHONHOME \n")

        #f.write("export PYTHONHOME=/usr/local  \n")

        '''
        # https://twiki.cern.ch/twiki/bin/view/Main/HomerWolfeCMSSWAndGDB : Do you see "ImportError: No module named site"
        Run scram tool info python, you'll see
        ...
        PYTHON_BASE=/afs/cern.ch/cms/slc5_amd64_gcc472/external/python/2.7.3-cms4
        ...
        or whatever. set PYTHONHOME equal to this value. 

        '''
        #f.write("export PYTHONHOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/python/2.7.14-omkpbe4  \n")

        f.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch \n")
        f.write("export SCRAM_ARCH=slc7_amd64_gcc700   \n")
        f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")
        f.write("export X509_USER_PROXY=/afs/cern.ch/user/s/ssawant/x509up_u108989  \n")

    #f.write("cd /afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/CMSSW_11_2_0/src/   \n")
    #f.write("cmsenv   \n")
    #f.write("cd /afs/cern.ch/work/s/ssawant/private/L1T_ServiceTasks/hcalPUsub_v4_20210510/myStudies/run_1/mimicMLJetRec/Run3_MC/tmp/  \n")
    #f.write("eval \n")
    f.write("cd %s \n" % (pwd))
    #f.write("%s \n" % (command0))
    for command_0 in commands:
        f.write("%s \n" % (command_0))


condor_submit_file = 'condor_submit_%s.sh' % (scriptName0)
with open(condor_submit_file, 'w') as f:
    f.write("universe = vanilla \n")
    f.write("executable = %s \n" % condor_exec_file)
    f.write("getenv = TRUE \n")
    f.write("log = %s.log \n" % (scriptName0))
    f.write("output = %s.out \n" % (scriptName0))
    f.write("error = %s.error \n" % (scriptName0))
    f.write("notification = never \n")
    f.write("should_transfer_files = YES \n")
    f.write("when_to_transfer_output = ON_EXIT \n")
    f.write("+JobFlavour = \"longlunch\" \n")
    #f.write("+JobFlavour = \"espresso\" \n")
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


cmd1 = "condor_submit %s" % condor_submit_file
print("Now:  %s " % cmd1)            
os.system(cmd1)
