#!/bin/bash

Year="UL18"
GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"


HLTPaths=(
    "HLT_IsoMu24"
    "HLT_IsoMu27"
    "HLT_Mu50"
    "HLT_OldMu100"
    "HLT_TkMu100"
    
)


source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env

# Example command:
# brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_PFJet500_v*" -o output_brilcalc_314472-325175_UL18_HLT_PFJet500_v.csv
for HLTPath in "${HLTPaths[@]}"
do 
    command="brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ${GoldenJSON} --hltpath \"${HLTPath}_v*\" -o output_brilcalc_${Year}_${HLTPath}_v.csv"
    printf "${command} \n"
    brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ${GoldenJSON} --hltpath "${HLTPath}_v*" -o output_brilcalc_${Year}_${HLTPath}_v.csv

    # 2018 HEM issues affected luminosity
    command="brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ${GoldenJSON} --hltpath \"${HLTPath}_v*\" --begin 319077 --end 325175 -o output_brilcalc_${Year}_${HLTPath}_v_run319077To325175.csv"
    printf "${command} \n"
    brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ${GoldenJSON} --hltpath "${HLTPath}_v*" --begin 319077 --end 325175 -o output_brilcalc_${Year}_${HLTPath}_v_run319077To325175.csv



done    


