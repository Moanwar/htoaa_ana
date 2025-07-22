
from collections import OrderedDict as OD
import enum
import copy

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookXrootdService
# cms-xrd-global.cern.ch "global redirector"
# xrootd-cms.infn.it for Europe and Asia
xrootd_redirectorName  = "root://xrootd-cms.infn.it//"
# try / or // at the end
xrootd_redirectorNames = [
    "root://xrootd-cms.infn.it/",
    "root://cms-xrd-global.cern.ch/",
    "root://cmsxrootd.fnal.gov/",
    "root://cms-xrd-global.cern.ch:1094/",
]
sampleFormat = "nanoAOD" 

### Miscellaneous constants
kPt_Max     = 99999.0
kLHE_HT_Max = 99999.0
NanoAODFileSize_Min = 0.3 # in MB
#------------------------------------

Era_2016        = '2016'
Era_2016preVFP  = '2016preVFP'
Era_2016postVFP = '2016postVFP'
Era_2017        = '2017'
Era_2018        = '2018'
Era_Run2        = 'Run2'

class DatasetToAnalyze(enum.Enum):
    FullRun2 = 'FullRun2'
    SingleYear = 'SingleYear'

### Set DatasetToAnalyze.SingleYear: to analyze a single year, DatasetToAnalyze.FullRun2: to analyse full Run2 data
# This is important for co-related/de-correlated systematic uncertainties.
kDatasetToAnalyze = DatasetToAnalyze.FullRun2 # DatasetToAnalyze.SingleYear, DatasetToAnalyze.FullRun2

sFileSamplesInfo = {
    #Era_2016: "Samples_2016UL.json",
    Era_2016preVFP:  "Samples_2016preVFPUL.json",
    Era_2016postVFP: "Samples_2016postVFPUL.json",
    Era_2017:        "Samples_2017UL.json",
    Era_2018:        "Samples_2018UL.json"
}

# Refer https://docs.google.com/spreadsheets/d/1xDLsr3ikLJxuMPNiSRs79YjTzbN64RetXL3A-tL6-hY/edit?usp=sharing
# /eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/QCD*/r1/PNet_*.root
sPathSkimmedNanoAODs = {
    Era_2016preVFP: {
        'skim_v2': {
            'Data': '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2016/data/PNet_v2_2024_11_22/$SAMPLENAME/r1_$ERATAG/PNet_*.root',
            'MC':   '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2016APV/MC/PNet_v2_2024_11_22/$SAMPLENAME/r*/PNet_*.root' 
        },
    },
    Era_2016postVFP: {
        'skim_v2': {
            'Data': '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2016/data/PNet_v2_2024_11_22/$SAMPLENAME/r1_$ERATAG/PNet_*.root',
            'MC':   '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2016/MC/PNet_v2_2024_11_22/$SAMPLENAME/r*/PNet_*.root' 
        },
    },
    Era_2017: {
        'skim_v2': {
            'Data': '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2017/data/PNet_v2_2024_11_22/$SAMPLENAME/r1_$ERATAG/PNet_*.root',
            'MC':   '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2017/MC/PNet_v2_2024_11_22/$SAMPLENAME/r*/PNet_*.root' 
        },
    },
    Era_2018: {
        'skim_v1': {
            'Data': '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/PNet_v1_2023_10_06/$SAMPLETAG/$SAMPLENAME/r*/PNet_*.root',
            'MC':   '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/$SAMPLENAME/r1/PNet_*.root' 
        },
        #'skim_Hto4b_0p8': {
        #    'Data': '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/PNet_v1_2023_10_06/$SAMPLETAG/$SAMPLENAME/skims/Hto4b_0p8/PNet_*.root',
        #    'MC':   '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/$SAMPLENAME/skims/Hto4b_0p8/PNet_*.root' 
        #},
        'skim_v2': {
            'Data': '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/PNet_v2_2024_11_22/$SAMPLENAME/r1_$ERATAG/PNet_*.root',
            'MC':   '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v2_2024_11_22/$SAMPLENAME/r*/PNet_*.root' 
        },
    }
}


sFilesGoldenJSON = {
    Era_2016:        'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt',   # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
    Era_2016preVFP:  'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt',
    Era_2016postVFP: 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt',
    Era_2017:        'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt', # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt   
    Era_2018:        'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt',   # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
}

YearsAndEras_dict = {
    Era_2016:        ['B-ver2_HIPM', 'C-HIPM', 'D-HIPM', 'E-HIPM', 'F-HIPM', 'F', 'G', 'H'],
    Era_2016preVFP:  ['B-ver2_HIPM', 'C-HIPM', 'D-HIPM', 'E-HIPM', 'F-HIPM'],
    Era_2016postVFP: ['F', 'G', 'H'],
    Era_2017:        ['B', 'C', 'D', 'E', 'F'],
    Era_2018:        ['A', 'B', 'C', 'D'],
}
Luminosities_Inclusive = { # [<lumi>, <uncertainty in percent> ] in fb^-1
    Era_2016:        [36.31, 1.2],
    Era_2016preVFP:  [19.5 , 1.2],
    Era_2016postVFP: [16.8 , 1.2],
    Era_2017:        [41.48, 2.3],
    Era_2018:        [59.83, 2.5]
}
Luminosities_perTrigger = {
    Era_2016preVFP: {
        'HLT_PFJet450':                                                     [19.498, 1.2],
        'HLT_DiCentralPFJet430':                                            [16.710, 1.2],
        'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5':                               [19.498, 1.2],
        'HLT_PFHT750_4JetPt50':                                             [19.498, 1.2],
        'HLT_PFHT800':                                                      [19.498, 1.2],
        'HLT_PFHT900':                                                      [19.498, 1.2],

        'HLT_AK8PFJet360_TrimMass30':                                       [19.498, 1.2],
        'HLT_AK8PFJet450':                                                  [16.710, 1.2],
        'HLT_AK8PFHT650_TrimR0p1PT0p03Mass50':                              [13.941, 1.2],
        'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50':                              [19.498, 1.2],

        'HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20':                     [13.941, 1.2],
        'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20':                     [19.498, 1.2],
        'HLT_PFHT400_SixJet30_DoubleBTagCSV_p056':                          [19.498, 1.2],
        'HLT_PFHT450_SixJet40_BTagCSV_p056':                                [19.498, 1.2],
        'HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20':                  [13.941, 1.2],

        'HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6': [19.498, 1.2],
        'HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160':           [19.498, 1.2],
        'HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6': [19.498, 1.2],
        'HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172':           [19.498, 1.2],
        'HLT_DoubleJet90_Double30_TripleBTagCSV_p08':                       [19.498, 1.2], # 2016:36.47. Got error with brilcalc. Hence guess
        'HLT_QuadJet45_TripleBTagCSV_p087':                                 [19.498, 1.2],

        'HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460':                            [16.721, 1.2],
        'HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500':                            [19.498, 1.2],
        'HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200':                        [16.721, 1.2],
        'HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240':                        [19.498, 1.2],

        'HLT_MET200':                                                       [19.498, 1.2],
        'HLT_PFMET110_PFMHT110_IDTight':                                    [19.498, 1.2],
        'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight':                            [19.498, 1.2],
        'HLT_PFMET120_PFMHT120_IDTight':                                    [19.498, 1.2],
        'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight':                            [19.498, 1.2],
        'HLT_PFMET170_HBHECleaned':                                         [19.498, 1.2],
        'HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight':         [19.498, 1.2],

        'HLT_IsoMu24':                                                      [19.498, 1.2],
        'HLT_IsoTkMu24':                                                    [19.498, 1.2],
        'HLT_IsoMu27':                                                      [19.498, 1.2],
        'HLT_Mu50':                                                         [19.498, 1.2],
        'HLT_TkMu50':                                                       [16.710, 1.2],
        
    },
    Era_2016postVFP: {
        'HLT_PFJet450':                                                     [16.812, 1.2],
        'HLT_DiCentralPFJet430':                                            [16.812, 1.2],
        'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5':                               [16.812, 1.2],
        'HLT_PFHT750_4JetPt50':                                             [ 8.072, 1.2],
        'HLT_PFHT800':                                                      [ 8.072, 1.2],
        'HLT_PFHT900':                                                      [16.812, 1.2],

        'HLT_AK8PFJet360_TrimMass30':                                       [16.812, 1.2],
        'HLT_AK8PFJet450':                                                  [16.812, 1.2],
        'HLT_AK8PFHT650_TrimR0p1PT0p03Mass50':                              [ 6.159, 1.2],
        'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50':                              [16.812, 1.2],

        'HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20':                     [ 6.159, 1.2],
        'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20':                     [16.812, 1.2],
        'HLT_PFHT400_SixJet30_DoubleBTagCSV_p056':                          [16.812, 1.2],
        'HLT_PFHT450_SixJet40_BTagCSV_p056':                                [16.812, 1.2],
        'HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20':                  [ 6.159, 1.2],

        'HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6': [16.617, 1.2],
        'HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160':           [16.617, 1.2],
        'HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6': [16.812, 1.2],
        'HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172':           [16.812, 1.2],
        'HLT_DoubleJet90_Double30_TripleBTagCSV_p08':                       [16.812, 1.2], # 2016:36.47. Got error with brilcalc. Hence guess
        'HLT_QuadJet45_TripleBTagCSV_p087':                                 [16.812, 1.2],

        'HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460':                            [ 8.527, 1.2],
        'HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500':                            [16.812, 1.2],
        'HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200':                        [ 8.527, 1.2],
        'HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240':                        [16.812, 1.2],

        'HLT_MET200':                                                       [16.812, 1.2],
        'HLT_PFMET110_PFMHT110_IDTight':                                    [16.190, 1.2],
        'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight':                            [16.190, 1.2],
        'HLT_PFMET120_PFMHT120_IDTight':                                    [16.812, 1.2],
        'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight':                            [16.812, 1.2],
        'HLT_PFMET170_HBHECleaned':                                         [16.812, 1.2],
        'HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight':         [16.812, 1.2],

        'HLT_IsoMu24':                                                      [16.812, 1.2],
        'HLT_IsoTkMu24':                                                    [16.812, 1.2],
        'HLT_IsoMu27':                                                      [16.812, 1.2],
        'HLT_Mu50':                                                         [16.812, 1.2],
        'HLT_TkMu50':                                                       [16.812, 1.2],

    },
    Era_2017: {
        'HLT_PFJet500':                                                     [41.478, 2.3], #[41.54, 2.3],
        'HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2':                       [36.675, 2.3], #[36.75, 2.3],
        'HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2':                   [27.122, 2.3], #[27.13, 2.3],
        'HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5':                             [32.112, 2.3], #[32.13, 2.3],
        'HLT_PFHT1050':                                                     [41.478, 2.3], #[41.54, 2.3],
        'HLT_AK8PFHT750_TrimMass50':                                        [30.897, 2.3], #[30.96, 2.3],
        'HLT_AK8PFHT800_TrimMass50':                                        [36.421, 2.3], #[36.49, 2.3],
        'HLT_AK8PFJet500':                                                  [41.478, 2.3], #[41.54, 2.3],
        'HLT_AK8PFJet360_TrimMass30':                                       [28.230, 2.3], #[28.3 , 2.3],
        'HLT_AK8PFJet380_TrimMass30':                                       [31.150, 2.3], #[31.22, 2.3],
        'HLT_AK8PFJet400_TrimMass30':                                       [36.675, 2.3], #[36.75, 2.3],
        'HLT_AK8PFJet330_PFAK8BTagCSV_p17':                                 [ 7.728, 2.3], #[ 7.73, 2.3],

        'HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33':              [36.264, 2.3], #[36.34, 2.3],
        'HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0':        [36.675, 2.3], #[36.75, 2.3],

        'HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1':             [ 7.728, 2.3], #[7.73, 2.3],
        'HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2':                       [ 7.728, 2.3], #[7.73, 2.3],

        'HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1':                    [36.675, 2.3], #[36.75, 2.3],
        'HLT_PFMET120_PFMHT120_IDTight_PFHT60':                             [36.675, 2.3], #[36.75, 2.3],
        'HLT_PFMET120_PFMHT120_IDTight':                                    [40.610, 2.3], #[40.67, 2.3],
        'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60':                     [36.675, 2.3], #[36.75, 2.3],
        'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight':                            [40.610, 2.3], #[40.67, 2.3],
        'HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60':                      [36.675, 2.3], #[36.75, 2.3],
        'HLT_PFMETTypeOne120_PFMHT120_IDTight':                             [40.610, 2.3], #[40.67, 2.3],
        'HLT_PFMET140_PFMHT140_IDTight':                                    [41.478, 2.3], 
        'HLT_PFMETTypeOne140_PFMHT140_IDTight':                             [41.478, 2.3], #[40.67, 2.3],
        'HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned':                         [36.675, 2.3], #[36.75, 2.3],
        'HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight':         [41.478, 2.3], 
        
        'HLT_IsoMu24':                                                      [37.997, 2.3], #[38.06, 2.3],
        'HLT_IsoMu27':                                                      [41.478, 2.3], #[41.54, 2.3],
        'HLT_Mu50':                                                         [41.478, 2.3], #[41.54, 2.3],
        'HLT_OldMu100':                                                     [36.675, 2.3],
        'HLT_TkMu100':                                                      [36.675, 2.3],
    },
    Era_2018: {
        'HLT_PFHT1050':                                                     [59.827, 2.5],
        'HLT_PFJet500':                                                     [59.827, 2.5],
        'HLT_AK8PFHT800_TrimMass50':                                        [59.827, 2.5],
        'HLT_AK8PFJet500':                                                  [59.827, 2.5],
        'HLT_AK8PFJet400_TrimMass30':                                       [59.827, 2.5],
        'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4':               [54.536, 2.5],

        'HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71':          [54.537, 2.5],
        'HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5':    [59.828, 2.5],

        'HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1':       [54.537, 2.5],
        'HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2':                 [54.537, 2.5],

        'HLT_PFMET120_PFMHT120_IDTight_PFHT60':                             [59.820, 2.5],
        'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight':                            [59.828, 2.5],
        'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60':                     [59.82, 2.5],
        'HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1':                [54.537, 2.5],
        'HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned':                         [59.828, 2.5],
        'HLT_PFMETTypeOne140_PFMHT140_IDTight':                             [59.828, 2.5],

        'HLT_IsoMu24':                                                      [59.820, 2.5],
        'HLT_IsoMu27':                                                      [59.827, 2.5],
        'HLT_Mu50':                                                         [59.827, 2.5],
        'HLT_OldMu100':                                                     [59.827, 2.5],
        'HLT_TkMu100':                                                      [59.827, 2.5],

        'HLT_Ele32_WPTight_Gsf':                                            [59.828, 2.5],
        'HLT_Ele35_WPTight_Gsf_L1EGMT':                                     [59.828, 2.5],
        'HLT_Ele115_CaloIdVT_GsfTrkIdT':                                    [59.828, 2.5],
        'HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165':                            [59.828, 2.5],
        
    },
}
Luminosities_TotalPerYear = { # [<lumi>, <uncertainty in percent> ] in fb^-1
    Era_2016: { # [36.31, 1.2]
        'Trg_Combo_AK4AK8Jet_HT' :     [36.31, 1.2], # https://twiki.cern.ch/twiki/bin/view/CMS/HLTPathsRunIIList#2017
        'Trg_Combo_AK4AK8Jet_HT_VBF' : [36.31, 1.2],
        'Trg_Combo_MET' :              [36.31, 1.2],
        'Trg_Combo_Mu' :               [36.31, 1.2],
        
    },
    Era_2016preVFP: { # [19.5 , 1.2]
        'Trg_Combo_AK4AK8Jet_HT' :     [19.5 , 1.2], # https://twiki.cern.ch/twiki/bin/view/CMS/HLTPathsRunIIList#2017
        'Trg_Combo_AK4AK8Jet_HT_VBF' : [19.5 , 1.2],
        'Trg_Combo_MET' :              [19.5 , 1.2],
        'Trg_Combo_Mu' :               [19.5 , 1.2],
        
    },
    Era_2016postVFP: { # [16.8 , 1.2]
        'Trg_Combo_AK4AK8Jet_HT' :     [16.8 , 1.2], # https://twiki.cern.ch/twiki/bin/view/CMS/HLTPathsRunIIList#2017
        'Trg_Combo_AK4AK8Jet_HT_VBF' : [16.8 , 1.2],
        'Trg_Combo_MET' :              [16.8 , 1.2],
        'Trg_Combo_Mu' :               [16.8 , 1.2],
        
    },
    Era_2017: { # [41.48, 2.3]
        'Trg_Combo_AK4AK8Jet_HT' :     [41.54, 2.3], # https://twiki.cern.ch/twiki/bin/view/CMS/HLTPathsRunIIList#2017
        'Trg_Combo_AK4AK8Jet_HT_VBF' : [41.54, 2.3],
        'Trg_Combo_MET' :              [40.67, 2.3],
        'Trg_Combo_Mu' :               [41.54, 2.3],
        
    },
    Era_2018: {
        'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4': [54.54, 2.5], # for HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_v* trigger. See ./data/luminosity/2018/output_brilcalc_314472-325175_UL18_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_Final.xlsx 
        'Trg_Combo_AK4AK8Jet_HT':                             [59.83, 2.5], # https://docs.google.com/spreadsheets/d/19ot4nFlhiJoD6v81qhgyKSYjE5PhgqgT2dz98PNIWg0/edit?usp=sharing
        'Trg_Combo_AK4AK8Jet_HT_VBF':                         [59.83, 2.5], # https://docs.google.com/spreadsheets/d/19ot4nFlhiJoD6v81qhgyKSYjE5PhgqgT2dz98PNIWg0/edit?usp=sharing
        'Trg_Combo_AK4AK8Jet_HT_MET':                         [59.83, 2.5], # Copied from Trg_Combo_AK4AK8Jet_HT. Needs to be checked.
        'Trg_Combo_MET':                                      [59.83, 2.5], # https://docs.google.com/spreadsheets/d/19ot4nFlhiJoD6v81qhgyKSYjE5PhgqgT2dz98PNIWg0/edit?usp=sharing
        'HLT_IsoMu24':                                        [59.82, 2.5], # See ./data/luminosity/2018/output_brilcalc_314472-325175_UL18_HLT_IsoMu24_v.xlsx
        'HLT_IsoMu27':                                        [59.83, 2.5], # See ./data/luminosity/2018/output_brilcalc_314472-325175_UL18_HLT_IsoMu27_v.xlsx
        'Trg_Combo_Mu':                                       [59.83, 2.5], # See ./data/luminosity/2018/Luminosity_HLTPaths.xlsx
    }, 
}
Luminosities_TotalPerYear_perEra = {
    Era_2018: {
        'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4': {
            'A':  8.736,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839
        },
        'Trg_Combo_AK4AK8Jet_HT': {
            'A': 14.027,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839,
        },
        'Trg_Combo_AK4AK8Jet_HT_VBF': {
            'A': 14.027,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839,
        },
        'Trg_Combo_AK4AK8Jet_HT_MET': { # Copied from Trg_Combo_AK4AK8Jet_HT. Needs to be checked.
            'A': 14.027,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839,
        },
        'Trg_Combo_MET': { 
            'A': 14.027,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839,
        },
        'HLT_IsoMu24': {
            'A': 14.019,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839
        }, 
        'HLT_IsoMu27': {
            'A': 14.027,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839
        },  
        'Trg_Combo_Mu': {
            'A': 14.027,
            'B':  7.067,
            'C':  6.895,
            'D': 31.839
        }, 
    }
}
Triggers_perEra = {
    Era_2016: {
        'Trg_Combo_AK4AK8Jet_HT_VBF': {
            'JetHT': { # JetHT primary dataset
                'HLT_PFJet450':                                    [], # 36.47 / 36.47
                'HLT_DiCentralPFJet430':                           [], # 33.64 / 36.47
                'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5':              [], # 36.47 / 36.47   # GGH->aa->4b efficiency: 28%
                'HLT_PFHT750_4JetPt50':                            [], # 27.71 / 36.47
                'HLT_PFHT800':                                     [], # 27.71 / 36.47
                'HLT_PFHT900':                                     [], # 36.47 / 36.47
                #
                'HLT_AK8PFJet360_TrimMass30':                      [], # 36.47 / 36.47
                'HLT_AK8PFJet450':                                 [], # 33.64 / 36.47
                'HLT_AK8PFHT650_TrimR0p1PT0p03Mass50':             [], # 20.2 / 36.47
                'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50':             [], # 36.47 / 36.47   # GGH->aa->4b efficiency: 61%
                #
                'HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20':    [], # 20.2 / 36.47
                'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20':    [], # 36.47 / 36.47
                'HLT_PFHT400_SixJet30_DoubleBTagCSV_p056':         [], # 36.47 / 36.47
                'HLT_PFHT450_SixJet40_BTagCSV_p056':               [], # 36.47 / 36.47
                'HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20': [], # 20.2 / 36.47
            },
            'BTagCSV': { # BTagCSV primary dataset
                'HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6': [], # 36.27 / 36.47
                'HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160': [], # 36.27 / 36.47
                'HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6': [], # 36.47 / 36.47
                'HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172': [], # 36.47 / 36.47
                'HLT_DoubleJet90_Double30_TripleBTagCSV_p08': [], # 36.47 / 36.47
                'HLT_QuadJet45_TripleBTagCSV_p087': [], # 36.47 / 36.47
                # VBF triggers
                'HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460': [], # 25.36 / 36.47
                'HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500': [], # 36.47 / 36.47
                'HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200': [], # 25.36 / 36.47
                'HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240': [], # 36.47 / 36.47
            },
        },
        'Trg_Combo_MET': {
            'MET': { # MET primary dataset
                'HLT_MET200':                                               [], # 36.47 / 36.47
                'HLT_PFMET110_PFMHT110_IDTight':                            [], # 35.83 / 36.47
                'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight':                    [], # 35.83 / 36.47
                'HLT_PFMET120_PFMHT120_IDTight':                            [], # 36.47 / 36.47
                'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight':                    [], # 36.47 / 36.47
                'HLT_PFMET170_HBHECleaned':                                 [], # 36.47 / 36.47
                'HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight': [], # 36.47 / 36.47
            },
        },
        'Trg_Combo_Mu': {
            'SingleMuon': { # SingleMuon primary dataset
                'HLT_IsoMu24':   [], # ['L1_SingleMu22'], # 36.47
                'HLT_IsoTkMu24': [], 
                'HLT_IsoMu27':   [], # ['L1_SingleMu22', 'L1_SingleMu25'], # 36.47
                'HLT_Mu50':      [], # ['L1_SingleMu22', 'L1_SingleMu25'], # 36.47
                'HLT_TkMu50':    [],
            },
        },
    },
    Era_2017: {
        'Trg_Combo_AK4AK8Jet_HT_VBF': {
            'JetHT': { # JetHT primary dataset
                'HLT_PFJet500':                                                  [], # 41.54 
                'HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2':                    [], # 36.75  
                'HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2':                [], # 
                'HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5':                          [], # 
                'HLT_PFHT1050':                                                  [], # 41.54           
                'HLT_AK8PFHT750_TrimMass50':                                     [], # 30.96 / 41.54   # GGH->aa->4b efficiency: 43%
                'HLT_AK8PFHT800_TrimMass50':                                     [], # 36.49 / 41.54   # GGH->aa->4b efficiency: 37%
                'HLT_AK8PFJet500':                                               [], # 41.54            
                'HLT_AK8PFJet360_TrimMass30':                                    [], # 28.30 / 41.54  # GGH->aa->4b efficiency: 59%
                'HLT_AK8PFJet380_TrimMass30':                                    [], # 31.22 / 41.54  # GGH->aa->4b efficiency: 52%
                'HLT_AK8PFJet400_TrimMass30':                                    [], # 36.75 / 41.54  # GGH->aa->4b efficiency: 45%     
            },     
            'BTagCSV': { # BTagCSV primary dataset
                #'HLT_AK8PFJet330_PFAK8BTagCSV_p17':                             [], # 7.73 / 41.54  BTagCSV dataset
                #
                'HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33':          [], # 36.34 / 41.54  BTagCSV dataset
                'HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0':    [], # 36.75 / 41.54   BTagCSV dataset           
                # VBF triggers:
                #'HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1':          [], # 7.73 / 41.54  BTagCSV dataset
                #'HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2':                    [], # 7.73 / 41.54  BTagCSV dataset
            },
        },       
        'Trg_Combo_MET': { 
            'MET': { # MET primary dataset
                'HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1':            [], # 36.75 / 41.54     
                'HLT_PFMET120_PFMHT120_IDTight_PFHT60':                     [], # 36.75 / 41.54   
                'HLT_PFMET120_PFMHT120_IDTight':                            [], # 40.67 / 41.54 
                'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60':             [], # 36.75 / 41.54 
                'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight':                    [], # 40.67 / 41.54 
                'HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60':              [], # 36.75 / 41.54    
                'HLT_PFMETTypeOne120_PFMHT120_IDTight':                     [], # 40.67 / 41.54 
                'HLT_PFMET140_PFMHT140_IDTight': [], # 41.54 / 41.54 
                'HLT_PFMETTypeOne140_PFMHT140_IDTight': [], # 41.54 / 41.54 
                'HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned':                 [], # 36.75 / 41.54 
                'HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight': [], # 41.54 / 41.54 
            },            
        },        
        'Trg_Combo_Mu': {
            'SingleMuon': { # SingleMuon primary dataset
                'HLT_IsoMu24': [], # ['L1_SingleMu22'], # 38.06 / 41.54
                'HLT_IsoMu27': [], # ['L1_SingleMu22', 'L1_SingleMu25'], # 41.54
                'HLT_Mu50':    [], # ['L1_SingleMu22', 'L1_SingleMu25'], # 41.54
                'HLT_OldMu100': [],
                'HLT_TkMu100': [],  
            },
        },        
    },
    Era_2018: {
        'Trg_Combo_AK4AK8Jet_HT_VBF': {
            'JetHT': { # JetHT primary dataset
                'HLT_PFJet500':                                                  [], # ['L1_SingleJet180'], 
                'HLT_PFHT1050':                                                  [], # ['L1_SingleJet180', 'L1_HTT360er'],
                'HLT_AK8PFHT800_TrimMass50':                                     [], # ['L1_SingleJet180', 'L1_HTT360er'],
                'HLT_AK8PFJet500':                                               [], # ['L1_SingleJet180'],
                'HLT_AK8PFJet400_TrimMass30':                                    [], # ['L1_SingleJet180'],
                'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4':            [], # ['L1_SingleJet180'],
                #
                'HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71':       [], # ['L1_DoubleJet112er2p3_dEta_Max1p6', 'L1_DoubleJet150er2p5'],
                'HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5': [], # ['L1_HTT320er', 'L1_HTT360er', 'L1_HTT400er', 'L1_ETT2000', 'L1_HTT320er_QuadJet_70_55_40_40_er2p4', 'L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3' ],
                #
                'HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1':    [], # ['L1_SingleJet180', 'L1_HTT320er', 'L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5'],
                'HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2':              [], # ['L1_SingleJet180', 'L1_HTT320er', 'L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5'],      
            },      
        },
        'Trg_Combo_MET': { # https://indico.cern.ch/event/1424480/#17-andrew-brinkerhoff
            'MET': { # MET primary dataset
                'HLT_PFMET120_PFMHT120_IDTight_PFHT60':               [], # ['L1_ETMHF90_HTT60er', 'L1_ETMHF100_HTT60er', 'L1_ETMHF110_HTT60er'], 
                #'HLT_PFMET120_PFMHT120_IDTight': [],
                'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60':       [], # ['L1_ETMHF90_HTT60er', 'L1_ETMHF100_HTT60er', 'L1_ETMHF110_HTT60er'],
                #'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight': [],
                'HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1':  [], # ['L1_ETMHF100', 'L1_ETMHF110', 'L1_ETMHF120', 'L1_ETMHF130'], 
                'HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned':           [], # ['L1_ETMHF100', 'L1_ETMHF110', 'L1_ETMHF120', 'L1_ETMHF130'], 
                #'HLT_PFMET200_HBHE_BeamHaloCleaned': [],
                'HLT_PFMETTypeOne140_PFMHT140_IDTight':               [], # ['L1_ETMHF100', 'L1_ETMHF110', 'L1_ETMHF120', 'L1_ETMHF130'], 
                #'HLT_PFMET140_PFMHT140_IDTight': [],
                #'HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight': [], # ZH->aa->4b+MET efficiency 99%
            },
        },        
        'Trg_Combo_Mu': {
            'SingleMuon': { # SingleMuon primary dataset
                'HLT_IsoMu24': [], # ['L1_SingleMu22'],
                'HLT_IsoMu27': [], # ['L1_SingleMu22', 'L1_SingleMu25'],
                'HLT_Mu50':    [], # ['L1_SingleMu22', 'L1_SingleMu25'],
                'HLT_OldMu100': [],
                'HLT_TkMu100': [],                
            },
        },
    }
}
Triggers_perEra[Era_2016preVFP ] = Triggers_perEra[Era_2016]
Triggers_perEra[Era_2016postVFP] = Triggers_perEra[Era_2016]
'''
## Andrew's slides on triggers https://indico.cern.ch/event/1479951/contributions/6234638/attachments/2968060/5241665/2024_11_15_HToAATo4B_selection_catgories_NanoAODTools.pdf#page=11
trigFat :
((HLT_PFJet500 || HLT_AK8PFJet500 || HLT_AK8PFJet400_TrimMass30 || HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4) && (L1_SingleJet180)) 
|| 
((HLT_AK8PFHT800_TrimMass50 || HLT_PFHT1050) && (L1_SingleJet180 || L1_HTT360er))

trigBtag : 
(HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 && (L1_DoubleJet112er2p3_dEta_Max1p6 || L1_DoubleJet150er2p5)) 
|| 
(HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5 && (L1_HTT320er || L1_HTT360er || L1_HTT400er || L1_ETT2000 || L1_HTT320er_QuadJet_70_55_40_40_er2p4 || L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3))

trigVBF :
(HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2 || HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1) && (L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5 || L1_HTT320er || L1_SingleJet180)

trigMET :
 ((HLT_PFMET120_PFMHT120_IDTight_PFHT60 || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60) && (L1_ETMHF90_HTT60er || L1_ETMHF100_HTT60er || L1_ETMHF110_HTT60er)) 
|| 
((HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1 || HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned || HLT_PFMETTypeOne140_PFMHT140_IDTight) && (L1_ETMHF100 || L1_ETMHF110 || L1_ETMHF120 || L1_ETMHF130))
'''
'''
## add individual HLT+L1 paths into Triggers_perEra to check Data-MC for individual HLT+L1 triggers. ^^^^ Its broken now
for era_ in Triggers_perEra:
    triggerCombs_0 = list(Triggers_perEra[era_].keys())
    for trg_ in triggerCombs_0:
        for hlt_, l1t_list in Triggers_perEra[era_][trg_].items():
            if hlt_ not in Triggers_perEra[era_]:
                Triggers_perEra[era_][hlt_] = { hlt_: l1t_list }
'''

## 2018 HEM1516 issue
HEM1516Issue2018_AffectedRunRange = [319077, 325175]
DataFractionAffectedBy2018HEM1516Issue = 0.7105 # factor = (luminosity for run >= 319077) / (2018 luminosity) = 38.7501 / 54.5365. Calculated for 2018 HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 trigger 
Weight_HEM1516Issue2018_perTrigger = {
    # Google-sheet: https://docs.google.com/spreadsheets/d/19ot4nFlhiJoD6v81qhgyKSYjE5PhgqgT2dz98PNIWg0/edit?usp=sharing
    'HLT_PFHT1050':	                                                 0.3523,
    'HLT_PFJet500':                	                                 0.3523,
    'HLT_AK8PFHT800_TrimMass50':	                                 0.3523,
    'HLT_AK8PFJet500':	                                             0.3523,
    'HLT_AK8PFJet400_TrimMass30':	                                 0.3523,
    'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4':	         0.2895,
                                    
    'HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71':	     0.2895,
    'HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5': 0.3523,
                                    
    'HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1': 	 0.2895,
    'HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2':	             0.2895,
                                    
    'HLT_PFMET120_PFMHT120_IDTight_PFHT60': 	                     0.3522,
    'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight':	                     0.3523,
    'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60':	                 0.3522,
    'HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1':	         0.2895,
    'HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned':	                     0.3523,
    'HLT_PFMETTypeOne140_PFMHT140_IDTight':	                         0.3523,    

    'HLT_IsoMu24':                                                   0.3522, # 1 - (38.750 / 59.820)
    'HLT_IsoMu27':                                                   0.3523, # 1 - (38.750 / 59.827)
    'HLT_Mu50':                                                      0.3523, # 1 - (38.750 / 59.827)
    'HLT_OldMu100':                                                  0.3523, # 1 - ( 38.750 / 59.827)
    'HLT_TkMu100':                                                   0.3523, # 1 - ( 38.750 / 59.827)    
    
}


# List of recommended MET filters. https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
MET_Filters = {} 
MET_Filters[Era_2016] = {
    "Data": [
        "goodVertices",                       # primary vertex filter ("Flag_goodVertices")
        "globalSuperTightHalo2016Filter",     # beam halo filter ("Flag_globalSuperTightHalo2016Filter")
        "HBHENoiseFilter",                    # HBHE noise filter ("Flag_HBHENoiseFilter")
        "HBHENoiseIsoFilter",                 # HBHEiso noise filter ("Flag_HBHENoiseIsoFilter")
        "EcalDeadCellTriggerPrimitiveFilter", # ECAL TP filter ("EcalDeadCellTriggerPrimitiveFilter")
        "BadPFMuonFilter",                    # Bad PF Muon Filter ("Flag_BadPFMuonFilter")
        "BadPFMuonDzFilter",                  # Bad PF Muon Dz Filter ("Flag_BadPFMuonDzFilter")
        "eeBadScFilter",                      # ee badSC noise filter ("Flag_eeBadScFilter")
        "hfNoisyHitsFilter",                  # HF noisy hits filter ("Flag_hfNoisyHitsFilter")
    ]
}
MET_Filters[Era_2016]["MC"]  = MET_Filters[Era_2016]["Data"] 

MET_Filters[Era_2018] = {
    "Data": [
        "goodVertices",                       # primary vertex filter ("Flag_goodVertices")
        "globalSuperTightHalo2016Filter",     # beam halo filter ("Flag_globalSuperTightHalo2016Filter")
        "HBHENoiseFilter",                    # HBHE noise filter ("Flag_HBHENoiseFilter")
        "HBHENoiseIsoFilter",                 # HBHEiso noise filter ("Flag_HBHENoiseIsoFilter")
        "EcalDeadCellTriggerPrimitiveFilter", # ECAL TP filter ("Flag_EcalDeadCellTriggerPrimitiveFilter")
        "BadPFMuonFilter",                    # Bad PF Muon Filter ("Flag_BadPFMuonFilter")
        "BadPFMuonDzFilter",                  # Bad PF Muon Dz Filter ("Flag_BadPFMuonDzFilter")
        #"hfNoisyHitsFilter", ##                 # HF noisy hits filter ("Flag_hfNoisyHitsFilter")
        "eeBadScFilter",                      # ee badSC noise filter ("Flag_eeBadScFilter")
        "ecalBadCalibFilter",                 # ECAL bad calibration filter update ("Flag_ecalBadCalibFilter")
    ],
}
MET_Filters[Era_2018]["MC"]  = MET_Filters[Era_2018]["Data"]
MET_Filters[Era_2017]        = MET_Filters[Era_2018]

MET_Filters[Era_2016preVFP]  = MET_Filters[Era_2016]
MET_Filters[Era_2016postVFP] = MET_Filters[Era_2016]
# -----------------------------------------------------------------------------------


class JetIDs(enum.IntEnum):
    tightIDFailingLeptonVeto = 2
    tightIDPassingLeptonVeto = 6


kMCSamplesStitch_PhSpOverlapRemove = 'MCSamplesStitch_PhSpOverlapRemove'
kMCSamplesStitch_PhSpOverlapRewgt  = 'MCSamplesStitch_PhSpOverlapRewgt'

class MCSamplesStitchOptions(enum.Enum):
    PhSpOverlapRemove = 'PhSpOverlapRemove'
    PhSpOverlapRewgt  = 'PhSpOverlapRewgt' 


sFileLumiScalesPhSpOverlapRewgt = {
    Era_2018: {
        "inputFile":     "data/lumiScale/2018.root",
        "histogramName": "$SAMPLECATEGORY/$SAMPLECATEGORY_LumiScale_PhSpOverlapRewghted",
    }
}


bTagWPs = {}
bTagWPs[Era_2018] = {
    'AK4DeepJet': { # https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
        'L': 0.0490,
        'M': 0.2783,
        'T': 0.7100
    },
    'DeepCSV': { # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        'L': 0.1208,
        'M': 0.4506,
    },
    'DDBvL': { # not provided for UL samples
        'M': 0.8, # taken from Si's code
    },
    'DDBvLV2': { # not provided for UL samples
        'M': 0.8, # taken from Si's code
    },
    'ParticleNetMD_XbbvsQCD': { 
        'VL': 0.75,
        # BTV-22-001 
        # https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=BTV-22-001&tp=an&id=2622&ancode=BTV-22-001
        # https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2021/005
        'L': 0.9172
    },
    'ParticleNetMD_Hto4b_Htoaa4bOverQCD': { # earlier ParticleNetMD_Hto4b_Htoaa4bOverQCD
        # https://ssawant.web.cern.ch/ssawant/HToAA/DatavsMC/20231106_PNetSignificanceScan_Msd90to140/?match=ParticleNetMD_Hto4b_Htoaa4bOverQCD         
        'WP-80': 0.920,
        #'WP-60': 0.978,
        'WP-60': 0.975, # https://indico.cern.ch/event/1348321/?note=257291#31-saswati-nandan
        'WP-40': 0.992,
        'WP-95': 0.80, # sideband minimum threshold for WP60 <--> Assumption
        'WP-99': 0.50, # sideband minimum threshold for WP80 <--> Assumption
    },
    'PNet_Xto4bv1_Htoaa4bOverQCD': { # earlier ParticleNetMD_Hto4b_Htoaa4bOverQCD
        # https://ssawant.web.cern.ch/ssawant/HToAA/DatavsMC/20231106_PNetSignificanceScan_Msd90to140/?match=ParticleNetMD_Hto4b_Htoaa4bOverQCD         
        'WP-80': 0.920,
        #'WP-60': 0.978,
        'WP-60': 0.975, # https://indico.cern.ch/event/1348321/?note=257291#31-saswati-nandan
        'WP-40': 0.992,
        'WP-95': 0.80, # sideband minimum threshold for WP60 <--> Assumption
        'WP-99': 0.50, # sideband minimum threshold for WP80 <--> Assumption
    },
    'PNet_Xto4bv2_Htoaa4b': {
        # Htoaato4b channel: https://mattermost.web.cern.ch/cms-exp/pl/icc97qchspnfprkiozptnqkksy  
        # X4b_v2 = (FatJet_PNet_X4b_v2a_Haa4b_score + FatJet_PNet_X4b_v2b_Haa4b_score) / 2.0 
        'SRWP-40':  0.96,   'SBWP-40':  0.84,    # fake rate 0.1%
        'SRWP-45a': 0.955,  'SBWP-45a': 0.80,    # fake rate ??%
        'SRWP-45b': 0.950,  'SBWP-45b': 0.77,    # fake rate ??%
        'SRWP-50':  0.945,  'SBWP-50':  0.74,    # fake rate 0.2%
        'SRWP-60':  0.93,   'SBWP-60':  0.66,    # fake rate 0.3%
        'SRWP-65':  0.92,   'SBWP-65':  0.60,    # fake rate 0.4%
        'SRWP-70':  0.90,   'SBWP-70':  0.50,    # fake rate 0.5%
        'SRWP-80':  0.84,   'SBWP-80':  0.40,    # fake rate 0.1%
        'SRWP-95':  0.66, # 'SBWP-60'
    },       
    'PNet_Xto4bv2a_Htoaa4b': {
        # Andrew, Hichem, Siddhesh chat: https://mattermost.web.cern.ch/cms-exp/pl/j6dnq8aid7nadnhb1w4sypqt3r 
        'SRWP-40': 0.968,  'SBWP-40': 0.84,    # fake rate 0.1%
        'SRWP-60': 0.944,  'SBWP-60': 0.66,    # fake rate 0.3%
        'SRWP-80': 0.868,  'SBWP-80': 0.40,    # fake rate 0.1%
    }, 
    'PNet_Xto4bv2b_Htoaa4b': {
        # Andrew, Hichem, Siddhesh chat: https://mattermost.web.cern.ch/cms-exp/pl/j6dnq8aid7nadnhb1w4sypqt3r 
        'SRWP-40': 0.952,  'SBWP-40': 0.84,    # fake rate 0.1%
        'SRWP-60': 0.916,  'SBWP-60': 0.66,    # fake rate 0.3%
        'SRWP-80': 0.814,  'SBWP-80': 0.40,    # fake rate 0.1%
    }, 
    'PNet_Xto34bv2_Htoaa4b': {
        # Andrew, Hichem, Siddhesh chat: https://mattermost.web.cern.ch/cms-exp/pl/48j5369mdbnyf8o9anb96ebjra 
        'SRWP-40':  0.89,   'SBWP-40':  0.55,    # fake rate 0.6%
        'SRWP-50':  0.85,   'SBWP-50':  0.41,    # fake rate 1.1%
        'SRWP-60':  0.78,   'SBWP-60':  0.27,    # fake rate 1.9%
        'SRWP-65':  0.73,   'SBWP-65':  0.21,    # fake rate 2.6%
        'SRWP-70':  0.67,   'SBWP-70':  0.14,    # fake rate 3.5%
        'SRWP-80':  0.53,   'SBWP-80':  0.03,    # fake rate 6.3%
    },                 
}

bTagWPs[Era_2017] = copy.deepcopy(bTagWPs[Era_2018])
bTagWPs[Era_2017]['AK4DeepJet'].update({ # https://btv-wiki.docs.cern.ch/ScaleFactors/Run2UL2017/
    'L': 0.0532,
    'M': 0.3040,
    'T': 0.7476
})
bTagWPs[Era_2017]['DeepCSV'].update({ # https://btv-wiki.docs.cern.ch/ScaleFactors/Run2UL2017/
    'L': 0.1355,
    'M': 0.4506,
})

bTagWPs[Era_2016preVFP] = copy.deepcopy(bTagWPs[Era_2018])
bTagWPs[Era_2016preVFP]['AK4DeepJet'].update({ # https://btv-wiki.docs.cern.ch/ScaleFactors/Run2UL2016preVFP/
    'L': 0.0508,
    'M': 0.2598,
    'T': 0.6502
})
bTagWPs[Era_2016preVFP]['DeepCSV'].update({ # https://btv-wiki.docs.cern.ch/ScaleFactors/Run2UL2016preVFP/
    'L': 0.2027,
    'M': 0.6001,
})

bTagWPs[Era_2016postVFP] = copy.deepcopy(bTagWPs[Era_2018])
bTagWPs[Era_2016postVFP]['AK4DeepJet'].update({ # https://btv-wiki.docs.cern.ch/ScaleFactors/Run2UL2016postVFP/
    'L': 0.0480,
    'M': 0.2489,
    'T': 0.6377
})
bTagWPs[Era_2016postVFP]['DeepCSV'].update({ # https://btv-wiki.docs.cern.ch/ScaleFactors/Run2UL2016postVFP/
    'L': 0.1918,
    'M': 0.5847,
})

topTagWPs = { # https://twiki.cern.ch/twiki/bin/viewauth/CMS/ParticleNetSFs
    Era_2016preVFP: {
        'PNetTvsQCD': {
            'T': 0.957,
        },
    },
    Era_2016postVFP: {
        'PNetTvsQCD': {
            'T': 0.958,
        },
    },
    Era_2017: {
        'PNetTvsQCD': {
            'T': 0.970,
        },
    },
    Era_2018: {
        'PNetTvsQCD': {
            'T': 0.97,
        },
    },    
}


## b-tag efficiency in MC for b-tag SF application
bTagSFEfficiencyDict = {}
bTagSFEfficiencyDict[Era_2018] = { # 'AK4DeepJet' WP-M
    'inputFile':    'data/correction/mc/BtagSF/2018/jetBtagEfficiency.root',
    'histogramName': {
        'b-flavour':     'hJetBtagEffi_b_TT_Presel', #'hJetBtagEffi_b_QCD_TT_Presel',
        'c-flavour':     'hJetBtagEffi_c_TT_Presel', #'hJetBtagEffi_c_QCD_TT_Presel',
        'light-flavour': 'hJetBtagEffi_l_TT_Presel', #'hJetBtagEffi_l_QCD_TT_Presel',           
    },
    'pTAxisRange': [20, 1000],
}
#
bTagSFEfficiencyDict[Era_2017       ] = copy.deepcopy( bTagSFEfficiencyDict[Era_2018] )
bTagSFEfficiencyDict[Era_2016preVFP ] = copy.deepcopy( bTagSFEfficiencyDict[Era_2018] )
bTagSFEfficiencyDict[Era_2016postVFP] = copy.deepcopy( bTagSFEfficiencyDict[Era_2018] )
bTagSFEfficiencyDict[Era_2017       ]['inputFile'] = 'data/correction/mc/BtagSF/2017/jetBtagEfficiency.root'
bTagSFEfficiencyDict[Era_2016preVFP ]['inputFile'] = 'data/correction/mc/BtagSF/2016preVFP/jetBtagEfficiency.root'
bTagSFEfficiencyDict[Era_2016postVFP]['inputFile'] = 'data/correction/mc/BtagSF/2016postVFP/jetBtagEfficiency.root'


Corrections = {

    "PURewgt": {
        Era_2018: {
            "inputFile":     "data/correction/mc/PURewgt/PURewgts_2018.root",
            "histogramName": "MC_PURewgt"
        }

    },
        
    "HTRewgt" : { # ./data/correction/mc/HTSamplesStitch/HTSamplesStitchSF_2018.root
        "QCD_bGen": {
            Era_2016preVFP: {
                "FitFunctionFormat": "1",
                "HT100to200":   "1",
                "HT200to300":   "1",
                "HT300to500":   "1",
                "HT500to700":   "1",
                "HT700to1000":  "1",
                "HT1000to1500": "1",
                "HT1500to2000": "1",
                "HT2000to3000": "1"
            },
            Era_2016postVFP: {
                "FitFunctionFormat": "1",
                "HT100to200":   "1",
                "HT200to300":   "1",
                "HT300to500":   "1",
                "HT500to700":   "1",
                "HT700to1000":  "1",
                "HT1000to1500": "1",
                "HT1500to2000": "1",
                "HT2000to3000": "1"
            },
            Era_2017: {
                "FitFunctionFormat": "1",
                "HT100to200":   "1",
                "HT200to300":   "1",
                "HT300to500":   "1",
                "HT500to700":   "1",
                "HT700to1000":  "1",
                "HT1000to1500": "1",
                "HT1500to2000": "1",
                "HT2000to3000": "1"
            },
            Era_2018: {
                "FitFunctionFormat": "{p0} + ({p1} * (x - {HTBinMin}))",
                "HT100to200":   "0.927235 + (0.001153 * (x - 100))",
                "HT200to300":   "0.936276 + (0.000820 * (x - 200))",
                "HT300to500":   "0.936325 + (0.000512 * (x - 300))",
                "HT500to700":   "0.968147 + (0.000378 * (x - 500))",
                "HT700to1000":  "0.931336 + (0.000229 * (x - 700))",
                "HT1000to1500": "0.956714 + (0.000122 * (x - 1000))",
                "HT1500to2000": "0.965210 + (0.000068 * (x - 1500))",
                "HT2000to3000": "1.006649 + (0.000030 * (x - 2000))"
            },
        }
    }, 

    "TopPtRewgt": { # Top pT reweights for ttbat sample. https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf#page=12
        "TuneCP5": {
            "FitFunctionFormat": "exp( {a} + ({b} * x) + ({c} * x * x) + ({d}/(x + {e})) )",
            "FitFunction": "exp( -2.02274e-01 + (1.09734e-04 * x) + (-1.30088e-07 * x * x) + (5.83494e+01/(x + 1.96252e+02)) )",
            "FitRange": [0, 3000],
        },
        "TuneCUETP": {
            "FitFunctionFormat": "{a} + ({b} * TanH({c} + ({d} * x) )",
            "FitFunction": "1.04554e+00 + (5.19012e-02 * TanH(-1.72927e+00 + (2.57113e-03 * x) )",
            "FitRange": [0, 3000],
        }
    },

    "HiggsPtRewgt": {
        "GGH_HToAATo4B": {
            'inputFile':     'data/correction/mc/HiggsPtRewgt/2018/ggHHiggsPtRewgt_HToAATo4B_TH1D.root', #'/eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts/2018/ggHHiggsPtRewgt_HToAATo4B.root', 
            'histogramName': 'hGenHiggsPt_Nom_Wgt_Hqt',
            'xAxisRange': [20, 800], # pT(GenHiggs)
        },
        "VBFH_HToAATo4B": {
            'inputFile':     'data/correction/mc/HiggsPtRewgt/2018/VBFHHiggsPtRewgt_HToAATo4B_TH1D.root', # /eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts/2018/VBFHHiggsPtRewgt_HToAATo4B.root
            'histogramName': 'hGenHiggsPt_Nom_Wgt_NLO',
            'xAxisRange': [20, 650], # pT(GenHiggs)
        },
        "WH_HToAATo4B": {
            'inputFile':     'data/correction/mc/HiggsPtRewgt/WH_ZH_wgts_Hichem/WH_2D_weight_18.root', # Hichem's file share on 01/07/2025: /afs/cern.ch/user/h/hboucham/public/commonFiles/VH_pt_SF/
            'histogramName': 'WH_weights_histo',
            'xAxisRange': [-1.10, 1.10], # log2( (2*pT_H) / (pT_H + pT_W) )
            'yAxisRange': [ 7.0 , 9.6 ], # log2( pT_H )
        },
        "ZH_HToAATo4B": {
            'inputFile':     'data/correction/mc/HiggsPtRewgt/WH_ZH_wgts_Hichem/ZH_2D_weight_18.root', # Hichem's file share on 01/07/2025: /afs/cern.ch/user/h/hboucham/public/commonFiles/VH_pt_SF/
            'histogramName': 'ZH_weights_histo',
            'xAxisRange': [-1.10, 1.10], # log2( (2*pT_H) / (pT_H + pT_W) )
            'yAxisRange': [ 7.0 , 9.6 ], # log2( pT_H )
        },
        "TTH_HToAATo4B": {
            'inputFile':     'data/correction/mc/HiggsPtRewgt/2018/ttHHiggsPtRewgt_HToAATo4B_TH1D.root', # /eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts/2018/ttHHiggsPtRewgt_HToAATo4B.root.
            'histogramName': 'hGenHiggsPt_Nom_Wgt_NLO',
            'xAxisRange': [20, 650], # pT(GenHiggs)
        },

    },

    "TrigEffi": {
        'Hadronic': { # SF as a fuction of (pT, msoft-drop)
            Era_2016preVFP: {
                'inputFile': {
                    'pTIncl': 'data/correction/mc/TrgEffSF/Hadronic/GGF_triggerSF_soup_1D_Inc_pT_Pre2016.json'
                },  
                'corrSetName': 'gg0l_triggerSF2016',
            },
            Era_2016postVFP: {
                'inputFile':   {
                    'pTIncl': 'data/correction/mc/TrgEffSF/Hadronic/GGF_triggerSF_soup_1D_Inc_pT_Post2016.json'
                },
                'corrSetName': 'gg0l_triggerSF2016',
            },
            Era_2017: {
                'inputFile':   {
                    'pTIncl': 'data/correction/mc/TrgEffSF/Hadronic/GGF_triggerSF_soup_1D_Inc_pT_2017.json'
                },
                'corrSetName': 'gg0l_triggerSF2017',
            },
            Era_2018: {
                'inputFile':   {
                    'pTIncl': 'data/correction/mc/TrgEffSF/Hadronic/GGF_triggerSF_soup_1D_Inc_pT_2018.json'
                },
                'corrSetName': 'gg0l_triggerSF2018',
            },             
        },
        'MET': { # SF as a fuction of (pT, msoft-drop)
            Era_2016preVFP: {
                'inputFile':   'data/correction/mc/TrgEffSF/MET/met_triggerSF_Pre2016.json',
                'corrSetName': 'met_triggerSF2016',
            },
            Era_2016postVFP: {
                'inputFile':   'data/correction/mc/TrgEffSF/MET/met_triggerSF_Post2016.json',
                'corrSetName': 'met_triggerSF2016',
            },
            Era_2017: {
                'inputFile':   'data/correction/mc/TrgEffSF/MET/met_triggerSF_2017.json',
                'corrSetName': 'met_triggerSF2017',
            },
            Era_2018: {
                'inputFile':   'data/correction/mc/TrgEffSF/MET/met_triggerSF_2018.json',
                'corrSetName': 'met_triggerSF2018',
            },             
        },
        
    },
    
    'ParticleNetMD_XbbvsQCD': { # Data-to-MC SFs for ParticleNetMD_XbbvsQCD:  BTV-22-001 
        # https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=BTV-22-001&tp=an&id=2622&ancode=BTV-22-001
        # https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2021/005
        Era_2018: {
            'L': {
                'pT_binEdges': [ 0,      450,      500,       600,    kPt_Max ],
                'SFs':         [      1,     0.921,     1.006,     1.001      ],
            },
        },

    },

}

SystNameConvUp   = 'Up'
SystNameConvDown = 'Down'
# 'PU':                 'CMS_pileup_$YEAR',
SystNameConvs = {
    ## Experimental 
    'PU':                 'CMS_pileup_$YEAR',
    'JetTrigEffi':        'CMS_eff_j_trigger_$YEAR',
    'MetTrigEffi':        'CMS_eff_met_trigger_$YEAR',
    'EleTrigEffi':        'CMS_eff_e_trigger_$YEAR',
    'MuTrigEffi':         'CMS_eff_m_trigger_$YEAR',
    'L1Prefire':          'CMS_l1_ecal_prefiring_$YEAR',

    '2018HEM1516Issue':   'CMS_HEM_2018',

    'Btag':               'CMS_btag_fixedWP_comb_bc_$YEAR',    
    'BtagCorr':           'CMS_btag_fixedWP_comb_bc_correlated',
    'BtagUncorr':         'CMS_btag_fixedWP_comb_bc_uncorrelated_$YEAR',

    'AK8JetJES':          'CMS_scale_fj_$YEAR', 
    'AK8JetJER':          'CMS_res_fj_$YEAR',
    'AK4JetJES':          'CMS_scale_j_$YEAR', 
    'AK4JetJER':          'CMS_res_j_$YEAR', 
    'METJES':             'CMS_scale_met_$YEAR', 
    'METJER':             'CMS_res_met_$YEAR', 
    'METUnclE':           'CMS_scale_met_unclustered_energy_$YEAR', 

    ## Theoretical 
    'ggHPtRewgt':         'higgs_pt_reweighting_ggH',
    'VBFHPtRewgt':        'higgs_pt_reweighting_qqH',
    'WHPtRewgt':          'higgs_pt_reweighting_WH',
    'ZHPtRewgt':          'higgs_pt_reweighting_ZH',
    'ttHPtRewgt':         'higgs_pt_reweighting_ttH',
    'LPRewgt':            'CMS_eff_fj_LundPlan_reweighting',

    'TopPtReWeight':      'top_pt_reweighting',
    
    'ISR':                'ps_isr',
    'FSR':                'ps_fsr',
    'QCDFactr':           'QCDscale_fac',
    'QCDRenorm':          'QCDscale_ren',
    'PDF':                'pdf_99',
    
}

massPseudoscalarA_windows_dict = OD([
    ('mA15Window', [14.2, 15.6]),
    ('mA30Window', [28.0, 32.0]),
    ('mA55Window', [50.0, 58.0]),
])
massHiggs_windows_dict = {
    'MsoftdropHiggsWindow': [110, 140],
    'PNet_massH_Hto4b_HiggsWindow': [110, 140],
}



### Miscellaneous variables
SplitQCDInGENCats = False
HistogramNameExtensions_QCD = ['0bCat', '1bCat', '2bCat', '3bCat', '4bAndMoreCat'] # ['0bCat', '1bCat', '2bCat', '3bCat', '4bCat', '5bAndMoreCat']
 



### GEN-level variables 
PDGID_BottomQuark = 5
PDGID_TopQuark    = 6

MASS_BottomQuark = 4.18

GENPART_STATUSFLAGS_LIST = [
    "isPrompt",
    "isDecayedLeptonHadron",
    "isTauDecayProduct",
    "isPromptTauDecayProduct",
    "isDirectTauDecayProduct",
    "isDirectPromptTauDecayProduct",
    "isDirectHadronDecayProduct",
    "isHardProcess",
    "fromHardProcess",
    "isHardProcessTauDecayProduct",
    "isDirectHardProcessTauDecayProduct",
    "fromHardProcessBeforeFSR",
    "isFirstCopy",
    "isLastCopy",
    "isLastCopyBeforeFSR",
]

class GENPART_STATUSFLAGS(enum.IntEnum):
    isPrompt                           =  0
    isDecayedLeptonHadron              =  1
    isTauDecayProduct                  =  2
    isPromptTauDecayProduct            =  3
    isDirectTauDecayProduct            =  4
    isDirectPromptTauDecayProduct      =  5
    isDirectHadronDecayProduct         =  6
    isHardProcess                      =  7
    fromHardProcess                    =  8
    isHardProcessTauDecayProduct       =  9
    isDirectHardProcessTauDecayProduct = 10
    fromHardProcessBeforeFSR           = 11
    isFirstCopy                        = 12
    isLastCopy                         = 13
    isLastCopyBeforeFSR                = 14


config_Template = OD([
    ("nEventsToAnalyze", -1),
    ("era", ''),
    ("dataset", ''), 
    ("inputFiles", ''),
    ("outputFile", ''),
    ("sampleCategory", ''),
    ("isMC", False),
    #("Luminosity", 0),
    ("crossSection", 0),
    ("nEvents", -1),
    ("sumEvents", -1),
    ("downloadIpFiles", False),   
    ("server", ''),   
])
