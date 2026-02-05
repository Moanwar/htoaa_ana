#!/bin/bash

# Make Rucio request for dataset listed in './datasetNamesForRucioRequest.txt'
#  Run: Write dataset name line-by-line in ./datasetNamesForRucioRequest.txt'
#       . makeRucioRequest.sh
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookFileTransfer#Auto_Approval_of_Rules
# https://cms-rucio-webui.cern.ch/ 

fDatasetNamesList=datasetNamesForRucioRequest.txt

echo "voms-proxy-init --voms cms --valid 192:00"
voms-proxy-init --voms cms --valid 192:00
#bash
echo "source /cvmfs/cms.cern.ch/rucio/setup-py3.sh"
source /cvmfs/cms.cern.ch/rucio/setup-py3.sh
echo "export RUCIO_ACCOUNT=$USER"
export RUCIO_ACCOUNT=$USER



while IFS= read -r line; do
    echo "Dataset: $line"
    echo "rucio add-rule cms:$line 1 T2_CH_CERN --grouping 'ALL' --ask-approval --activity 'User AutoApprove' --lifetime 15550000 --comment 'Run 2 UL NanoAOD'"
    rucio add-rule cms:$line 1 T2_CH_CERN --grouping 'ALL' --ask-approval --activity 'User AutoApprove' --lifetime 15550000 --comment 'Run 2 UL NanoAOD'
done < "$fDatasetNamesList"

echo "rucio list-rules --account $USER"
rucio list-rules --account $USER


