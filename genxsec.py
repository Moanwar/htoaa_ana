#usage: cmsRun genxsec.py inputFiles=file:/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v2_2024_11_22/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-15_TuneCP5_13TeV_madgraph_pythia8/r1/PNet_v1_Skim_Lp_weight.root
#or
#cmsRun genxsec.py inputFiles=/store/mc/RunIISummer20UL18MiniAODv2/SUSY_GluGluH_01J_HToAATo4B_M-15_TuneCP5_13TeV_madgraph_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/033AE6B3-77E2-E749-A0F8-463D31167A46.root

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.parseArguments()

process = cms.Process('GENXSEC')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000000

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(options.inputFiles),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.genxsec = cms.EDAnalyzer('GenXSecAnalyzer')

process.path = cms.Path(process.genxsec)
