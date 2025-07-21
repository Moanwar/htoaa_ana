#!/bin/bash

Year="UL17"
GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
: '
HLTPaths=(
    "HLT_PFJet500"
    "HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2"
    "HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2"
    "HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5"
    "HLT_PFHT1050"
    "HLT_AK8PFHT750_TrimMass50"
    "HLT_AK8PFHT800_TrimMass50"
    "HLT_AK8PFJet500"
    "HLT_AK8PFJet360_TrimMass30"
    "HLT_AK8PFJet380_TrimMass30"
    "HLT_AK8PFJet400_TrimMass30"

    "HLT_AK8PFJet330_PFAK8BTagCSV_p17"
    "HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33"
    "HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0"

    "HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1"
    "HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2"

    "HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1"
    "HLT_PFMET120_PFMHT120_IDTight_PFHT60"
    "HLT_PFMET120_PFMHT120_IDTight"
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60"
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"
    "HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60"
    "HLT_PFMETTypeOne120_PFMHT120_IDTight"
    "HLT_PFMET140_PFMHT140_IDTight"
    "HLT_PFMETTypeOne140_PFMHT140_IDTight"
    "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned"
    "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"

    "HLT_IsoMu24"
    "HLT_IsoMu27"
    "HLT_Mu50"
)
'

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

done    


