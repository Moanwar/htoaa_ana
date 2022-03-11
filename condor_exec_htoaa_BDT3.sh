#!/bin/bash  

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch 
export SCRAM_ARCH=slc7_amd64_gcc700   
source /cvmfs/cms.cern.ch/cmsset_default.sh 

export X509_USER_PROXY=/afs/cern.ch/user/s/ssawant/x509up_u108989  
cd /afs/cern.ch/work/s/ssawant/private/htoaa/repo_Si/htoaa 
source /afs/cern.ch/user/s/ssawant/.bashrc 
time conda env list 
time conda activate mlemv_fromSi 
time conda env list 
python --version 
python3 -V 
which python 
which python3 
time conda list 
which conda 
time /afs/cern.ch/work/s/ssawant/private/softwares/anaconda3/envs/mlemv_fromSi/bin/python3 htoaa_BDT3.py 
