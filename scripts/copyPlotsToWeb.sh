#!/bin/bash

####################################
## Script to copy plots to www
#  To run: for e.g.
#         time . copyPlotsToWeb.sh /eos/cms/store/user/ssawant/htoaa/analysis/20250805_DataMC_1 /eos/user/s/ssawant/www/HToAA/DatavsMC/20250805_DataMC
####################################

SourceDir=$1
DestinationDir=$2

echo "copy from ${SourceDir} to ${DestinationDir}"
Eras=("Run2" "2016preVFP" "2016postVFP" "2017" "2018" )
#Categories=("gg0l" "Vjj" "Zvv" "tt0l")
Categories=("Vjj" "Zvv" "tt0l")

PWD=`pwd`
for era in ${Eras[@]}; do 
    echo "era ${era}"
    for Cat in ${Categories[@]}; do        
        source="${SourceDir}/${era}/${Cat}/plots/*"
        destination="${DestinationDir}/${Cat}/${era}"

        echo "mkdir -p ${destination}"
        mkdir -p ${destination}
        echo "cp -r ${source} ${destination}"
        cp -r ${source} ${destination}

    done
done


echo "find ${DestinationDir}/ -type d | xargs -n 1 cp -v ${DestinationDir}/../index.php"
find ${DestinationDir}/ -type d | xargs -n 1 cp -v ${DestinationDir}/../index.php 



