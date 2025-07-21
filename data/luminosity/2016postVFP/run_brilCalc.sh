#!/bin/bash

Year="UL16postVFP"
GoldenJSON="/afs/cern.ch/work/s/ssawant/private/htoaa/htoaa_b_ana_SS/data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016postVFP_Collisions16_JSON.txt"

: '
HLTPaths=(
    "HLT_PFJet450"
    "HLT_DiCentralPFJet430"
    "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5"
    "HLT_PFHT750_4JetPt50"
    "HLT_PFHT800"
    "HLT_PFHT900"

    "HLT_AK8PFJet360_TrimMass30"
    "HLT_AK8PFJet450"
    "HLT_AK8PFHT650_TrimR0p1PT0p03Mass50"
    "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50"

    "HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20"
    "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20"
    "HLT_PFHT400_SixJet30_DoubleBTagCSV_p056"
    "HLT_PFHT450_SixJet40_BTagCSV_p056"
    "HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20"

    "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6"
    "HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160"
    "HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6"
    "HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172"
    "HLT_DoubleJet90_Double30_TripleBTagCSV_p08"
    "HLT_QuadJet45_TripleBTagCSV_p087"

    "HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460"
    "HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500"
    "HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200"
    "HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240"

    "HLT_MET200"
    "HLT_PFMET110_PFMHT110_IDTight"
    "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight"
    "HLT_PFMET120_PFMHT120_IDTight"
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"
    "HLT_PFMET170_HBHECleaned"
    "HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight"
    
)
'

HLTPaths=(
    "HLT_IsoMu24"
    "HLT_IsoTkMu24"
    "HLT_IsoMu27"
    "HLT_Mu50"
    "HLT_TkMu50"
    
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


