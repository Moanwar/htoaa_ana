from collections import OrderedDict as OD
import enum
Era_2018 = '2018'

MET_Filters = {}
MET_Filters[Era_2018] = {
    "Data": [
        "goodVertices",                       # primary vertex filter ("Flag_goodVertices")
        "globalSuperTightHalo2016Filter",     # beam halo filter ("Flag_globalSuperTightHalo2016Filter")
        "HBHENoiseFilter",                    # HBHE noise filter ("Flag_HBHENoiseFilter")
        "HBHENoiseIsoFilter",                 # HBHEiso noise filter ("Flag_HBHENoiseIsoFilter")
        "EcalDeadCellTriggerPrimitiveFilter", # ECAL TP filter ("Flag_EcalDeadCellTriggerPrimitiveFilter")
        "BadPFMuonFilter",                    # Bad PF Muon Filter ("Flag_BadPFMuonFilter")
        "BadPFMuonDzFilter",                  # Bad PF Muon Dz Filter ("Flag_BadPFMuonDzFilter")
        "hfNoisyHitsFilter",                  # HF noisy hits filter ("Flag_hfNoisyHitsFilter")
        "eeBadScFilter",                      # ee badSC noise filter ("Flag_eeBadScFilter")
        "ecalBadCalibFilter",                 # ECAL bad calibration filter update ("Flag_ecalBadCalibFilter")
    ],
}

MET_Filters[Era_2018]["MC"]  = MET_Filters[Era_2018]["Data"]


class JetIDs(enum.IntEnum):
    tightIDFailingLeptonVeto = 2
    tightIDPassingLeptonVeto = 6

Luminosities_forGGFMode = { Era_2018: {'Trg_Combo_AK4AK8Jet_HT': [59.83, 2.5]}, }

Triggers_perEra = {
    Era_2018: {
        'Trg_Combo_AK4AK8Jet_HT': {
            'HLT_PFJet500':                                       ['L1_SingleJet180'],
            'HLT_PFHT1050':                                       ['L1_HTT360er'],
            'HLT_AK8PFHT800_TrimMass50':                          ['L1_HTT360er'],
            'HLT_AK8PFJet500':                                    ['L1_SingleJet180'],
            'HLT_AK8PFJet400_TrimMass30':                         ['L1_SingleJet180'],
            'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4': ['L1_SingleJet180'],
            'HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5': ['L1_HTT320er', 'L1_HTT360er', 'L1_HTT400er', 'L1_ETT2000', 'L1_HTT320er_QuadJet_70_55_40_40_er2p4', 'L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3' ],
            'HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71': ['L1_DoubleJet112er2p3_dEta_Max1p6','L1_DoubleJet150er2p5'],
            'HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2': ['L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5','L1_HTT320er','L1_SingleJet180'],
            'HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1': ['L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5','L1_HTT320er','L1_SingleJet180'],
        }
    }
}


