#!/bin/bash  

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch 
export SCRAM_ARCH=slc7_amd64_gcc10  
source /cvmfs/cms.cern.ch/cmsset_default.sh 

#export X509_USER_PROXY=/afs/cern.ch/user/s/ssawant/x509up_u108989  
#export X509_USER_PROXY=$1 
#export EOS_MGM_URL=root://eoscms.cern.ch  
#voms-proxy-info -all 
#voms-proxy-info -all -file $1 

source /afs/cern.ch/user/m/moanwar/.bashrc 
which conda 
time conda env list 
conda activate myCondaEnv
conda env list

#cp -r /afs/cern.ch/work/s/ssawant/private/htoaa/htoaa_b_ana_SS/* . 
printf "pwd: \n" 
pwd 
printf "ls: \n" 
ls 
echo "$1 " $1 
#echo "$2 " $2 

pythonPath=$(which python3) 
# --TopMDWP can be tight, med, lo, veto, inclusive
#time ${pythonPath} pnet_calib_trial.py $1 --TopMDWP lo  --jetsys jerDown 
#time ${pythonPath} pnet_calib_trial.py $1 --TopMDWP med
#time ${pythonPath} pnet_calib_trial.py $1 --TopMDWP veto
#time ${pythonPath} pnet_calib_trial.py $1 --TopMDWP tight
#time ${pythonPath} pnet_calib_trial.py $1 --TopMDWP lo
#time ${pythonPath} pnet_calib_trial.py $1 --TopMDWP lo  --jetsys jerUp
#time ${pythonPath} pnet_calib_trial.py $1 --TopMDWP lo  --jetsys jesDown
#time ${pythonPath} pnet_calib_trial.py $1 --TopMDWP lo  --jetsys jesUp
time ${pythonPath} htoaa_Analysis_dataOnly_VBF_v13.py $1 -j 2 -c 50000
printf "After execution pwd: \n" 
pwd 
printf "ls: \n" 
ls 
time cp *.root  /eos/user/m/moanwar/htoaa/run2_vbfdatav2/
rm -rf *.root
