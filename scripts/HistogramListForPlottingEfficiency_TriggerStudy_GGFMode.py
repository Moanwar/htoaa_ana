import os
import numpy as np
from collections import OrderedDict as OD
import json
import copy

sXRange = "xAxisRange"; sYRange = "yAxisRange";
sXLabel = 'xAxisLabel'; sYLabel = 'yAxisLabel';
sXScale = 'xAxisScale';
sNRebinX = 'nRebinX';  sNRebinY = 'nRebinY';
sHistosToOverlay = 'histosToOverlay'
sHistosToHadd = 'histosToHadd'
sIpFileNameNice = 'ipFileNameNice'
sHistName   = 'histogramName'
sHistNumerator = 'Numerator' #'histogramNumerator'
sHistDenominator = 'Denominator' #'histogramDenominator'
sSelTagNameNice = 'sSelTagNameNice'
sPlotTag = 'plotTag'



sIpFiles = OD([
    # (<file name to refer>, <file path+name>)
    ("fIp1", '/eos/cms/store/user/ssawant/htoaa/analysis/20231123_TrgEffiSF_PNetMD_Hto4b_Htoaa4bOverQCDWP80/2018/analyze_htoaa_stage1.root')
])
sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20231123_TrgEffiSF_PNetMD_Hto4b_Htoaa4bOverQCDWP80/2018/plots_effi'

sampleCategory_dict = OD([
    ('SingleMuon_Run2018ABCD', ['SingleMuon_Run2018A', 'SingleMuon_Run2018B', 'SingleMuon_Run2018C', 'SingleMuon_Run2018D', ]),
    ('TTToSemiLeptonic_powheg', ['TTToSemiLeptonic_powheg']),
    ('MC', [
        'QCD_0bCat', 'QCD_1bCat', 'QCD_2bCat', 'QCD_3bCat', 'QCD_4bCat', 'QCD_5bAndMoreCat',  
        'TTToHadronic_powheg', 'TTToSemiLeptonic_powheg', 'TTTo2L2Nu_powheg', "SingleTop", 
        'ZJetsToQQ_HT', 
        #"DYJets_M-10to50_Incl_LO", "DYJets_M-50_Incl_LO", #"DYJets_M-10to50_Incl_NLO", "DYJets_M-50_Incl_NLO", 
        "DYJets_M-50_Incl_NLO", #"DYJets_M-10to50_Incl_LO", "DYJets_HT_LO", 
        'WJetsToQQ_HT', 
        'WJetsToLNu_HT_LO', #'WJetsToLNu_Incl_NLO', 
    ]),
    #('TT_W', ['TTToSemiLeptonic_powheg', 'WJetsToLNu_HT_LO', ])
])
#sampleCategory_dict = OD([
#    ('TTToSemiLeptonic_powheg', ['TTToSemiLeptonic_powheg']),
#])
sampleCategory_Data_forEffiSF = 'SingleMuon_Run2018ABCD'
#sampleCategory_MC_forEffiSF   = 'TTToSemiLeptonic_powheg'
sampleCategory_MC_forEffiSF_list   = ['TTToSemiLeptonic_powheg', 'MC'] #['TTToSemiLeptonic_powheg', 'MC', 'TT_W']

efficiencyHistogramNameNice_dict = OD([
    ('hLeadingFatJetPt', {sXLabel: 'hLeadingFatJetPt', sXRange: [180, 1000],  sNRebinX: 4,}),
    ('hLeadingFatJetPt_msoftdropGt60_btagHbbGtnp1', {sXLabel: 'LeadingFatJetPt_msoftdropGt60_btagHbbGtnp1 [GeV]', sXRange: [180, 1000],  sNRebinX: 4,}),
    ('hLeadingFatJetPt_msoftdropGt60_PNetMD_Hto4b_Htoaa4bOverQCDWP80', {sXLabel: 'LeadingFatJetPt_msoftdropGt60_PNetMD_Hto4b_Htoaa4bOverQCDWP80 [GeV]', sXRange: [180, 1000],  sNRebinX: 4,}),
    ('hLeadingFatJetEta', {sXLabel: 'hLeadingFatJetEta', sXRange: [-3.5, 3.5],  sNRebinX: 2,}),
    ('hLeadingFatJetPhi', {sXLabel: 'hLeadingFatJetPhi', sXRange: [-3.14, 3.14],  sNRebinX: 2,}),
    ('hLeadingFatJetMass', {sXLabel: 'hLeadingFatJetMass', sXRange: [0, 300],  sNRebinX: 5,}),
    ('hLeadingFatJetMSoftDrop', {sXLabel: 'hLeadingFatJetMSoftDrop', sXRange: [0, 300],  sNRebinX: 5,}),
    ('hLeadingFatJetMass_pTGt400_btagHbbGtnp1', {sXLabel: 'hLeadingFatJetMass_pTGt400_btagHbbGtnp1', sXRange: [0, 300],  sNRebinX: 5,}),
    ('hLeadingFatJetMSoftDrop_pTGt400_btagHbbGtnp1', {sXLabel: 'hLeadingFatJetMSoftDrop_pTGt400_btagHbbGtnp1', sXRange: [0, 300],  sNRebinX: 5,}),
    ('hLeadingFatJetMSoftDrop_pTGt400_PNetMD_Hto4b_Htoaa4bOverQCDWP80', {sXLabel: 'hLeadingFatJetMSoftDrop_pTGt400_PNetMD_Hto4b_Htoaa4bOverQCDWP80', sXRange: [0, 300],  sNRebinX: 5,}),
    ('hLeadingFatJetBtagCSVV2', {sXLabel: 'hLeadingFatJetBtagCSVV2', sXRange: [0, 1],  sNRebinX: 2,}),
    ('hLeadingFatJetBtagDDBvLV2', {sXLabel: 'hLeadingFatJetBtagDDBvLV2', sXRange: [0, 1],  sNRebinX: 2,}),
    ('hLeadingFatJetBtagDeepB', {sXLabel: 'hLeadingFatJetBtagDeepB', sXRange: [0, 1],  sNRebinX: 2,}),
    ('hLeadingFatJetBtagHbb', {sXLabel: 'hLeadingFatJetBtagHbb', sXRange: [-1, 1],  sNRebinX: 2,}),
    ('hLeadingFatJetParticleNetMD_XbbOverQCD', {sXLabel: 'hLeadingFatJetParticleNetMD_XbbOverQCD', sXRange: [0, 1],  sNRebinX: 2,}),
    ('hLeadingFatJetBtagHbb_pTGt400_msoftdropGt60', {sXLabel: 'hLeadingFatJetBtagHbb_pTGt400_msoftdropGt60', sXRange: [-1, 1],  sNRebinX: 2,}),
    ('hLeadingFatJetParticleNetMD_XbbOverQCD_pTGt400_msoftdropGt60', {sXLabel: 'hLeadingFatJetParticleNetMD_XbbOverQCD_pTGt400_msoftdropGt60', sXRange: [0, 1],  sNRebinX: 2,}),
    ('hLeadingFatJetPNetMD_Hto4b_Htoaa4bOverQCD_pTGt400_msoftdropGt60', {sXLabel: 'hLeadingFatJetPNetMD_Hto4b_Htoaa4bOverQCD_pTGt400_msoftdropGt60', sXRange: [0, 1],  sNRebinX: 4,}),
    #('', {sXLabel: '',  sNRebinX: ,}),

    
])

efficiencyHistogramNameNice_dict = OD([
#    ('hLeadingFatJetPt_msoftdropGt60_btagHbbGtnp1', {sXLabel: 'LeadingFatJetPt_msoftdropGt60_btagHbbGtnp1 [GeV]', sXRange: [180, 1000],  sNRebinX: 4,}),
#    ('hdR_leadingMuon_leadingFatJet', {sXLabel: 'hdR_leadingMuon_leadingFatJet',  sNRebinX: 2,}),
    ('hLeadingFatJetPt_msoftdropGt60_PNetMD_Hto4b_Htoaa4bOverQCDWP80', {sXLabel: 'LeadingFatJetPt_msoftdropGt60_PNetMD_Hto4b_Htoaa4bOverQCDWP80 [GeV]', sXRange: [180, 1000],  sNRebinX: 4,}),
    
])

#selectionConditionTag_forNumerator_list   = ['SR_TrgAK8330_M30_BDBnp4', 'SR_Trg2AK4116_DCSVp71', 'SR_TrgAK8400_M30', 'SR_TrgAK8500', 'SR_TrgComb2', 'SR_TrgComb4' ]
#selectionConditionTag_forDenominator = 'SR'
selectionConditionTags_forEffiRatio_dict = OD([
    #('SR_TrgAK8330_M30_BDBnp4', {sHistDenominator: 'SR', sHistNumerator: 'SR_TrgAK8330_M30_BDBnp4', sSelTagNameNice: 'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 (A)'}),
    #('SR_Trg2AK4116_DCSVp71', {sHistDenominator: 'SR', sHistNumerator: 'SR_Trg2AK4116_DCSVp71', sSelTagNameNice: 'HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 (B)'}),
    #('SR_TrgAK8400_M30', {sHistDenominator: 'SR', sHistNumerator: 'SR_TrgAK8400_M30', sSelTagNameNice: 'HLT_AK8PFJet400_TrimMass30 (C)'}),
    #('SR_TrgAK8500', {sHistDenominator: 'SR', sHistNumerator: 'SR_TrgAK8500', sSelTagNameNice: 'HLT_AK8PFJet500 (D)'}),
    #('SR_TrgComb2', {sHistDenominator: 'SR', sHistNumerator: 'SR_TrgComb2', sSelTagNameNice: 'A + B'}),
    #('SR_TrgComb4', {sHistDenominator: 'SR', sHistNumerator: 'SR_TrgComb4', sSelTagNameNice: 'A + B + C + D'}),
    ('SR_Trg_PFJet500', {sHistDenominator: 'SR', sHistNumerator: 'SR_Trg_PFJet500', sSelTagNameNice: 'HLT_PFJet500'}),
    ('SR_Trg_PFHT1050', {sHistDenominator: 'SR', sHistNumerator: 'SR_Trg_PFHT1050', sSelTagNameNice: 'HLT_PFHT1050'}),
    ('SR_Trg_AK8PFHT800_TrimMass50', {sHistDenominator: 'SR', sHistNumerator: 'SR_Trg_AK8PFHT800_TrimMass50', sSelTagNameNice: 'HLT_AK8PFHT800_TrimMass50'}),
    ('SR_Trg_AK8PFJet500', {sHistDenominator: 'SR', sHistNumerator: 'SR_Trg_AK8PFJet500', sSelTagNameNice: 'HLT_AK8PFJet500'}),
    ('SR_Trg_AK8PFJet400_TrimMass30', {sHistDenominator: 'SR', sHistNumerator: 'SR_Trg_AK8PFJet400_TrimMass30', sSelTagNameNice: 'HLT_AK8PFJet400_TrimMass30'}),
    ('SR_Trg_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4', {sHistDenominator: 'SR', sHistNumerator: 'SR_Trg_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4', sSelTagNameNice: 'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4'}),
    ('SR_Trg_Combo_AK4AK8Jet_HT', {sHistDenominator: 'SR', sHistNumerator: 'SR_Trg_Combo_AK4AK8Jet_HT', sSelTagNameNice: 'All'}),
])


print(f"sampleCategory_dict: {json.dumps(sampleCategory_dict, indent=4)}")
print(f"efficiencyHistogramNameNice_dict: {json.dumps(efficiencyHistogramNameNice_dict, indent=4)}")
print(f"selectionConditionTags_forEffiRatio_dict: {json.dumps(selectionConditionTags_forEffiRatio_dict, indent=4)}")


efficiencyHistograms_dict = OD() # dictionary holds information for efficiency plot. For example, see below commented part
'''
efficiencyHistograms_dict['hEfficiency_LeadingFatJetPt_msoftdropGt60_btagHbbGtnp1'] = {
    sHistNumerator: [
        {sIpFileNameNice: "fIp1", sHistName: 'evt/SingleMuon_Run2018A/hLeadingFatJetPt_msoftdropGt60_btagHbbGtnp1_SR_noweight'}, # add multiple histograms if you want to hadd
    ],
    sHistDenominator: [
        {sIpFileNameNice: "fIp1", sHistName: 'evt/SingleMuon_Run2018A/hLeadingFatJetPt_msoftdropGt60_btagHbbGtnp1_sel_JetID_noweight'}, # add multiple histograms if you want to hadd
    ],
    sXLabel: 'LeadingFatJetPt_msoftdropGt60_btagHbbGtnp1 [GeV]',
    sNRebinX: 5,
}
'''
for efficiencyHistogramNameNice, histoDetails_dict in efficiencyHistogramNameNice_dict.items():
    efficiencyHistograms_dict[efficiencyHistogramNameNice] = copy.deepcopy(histoDetails_dict)

    for sampleCategoryNameNice in sampleCategory_dict:
        efficiencyHistograms_dict[efficiencyHistogramNameNice][sampleCategoryNameNice] = OD()
        
        for sSelTagName, sSelTag_details_dict in selectionConditionTags_forEffiRatio_dict.items():
            efficiencyHistograms_dict[efficiencyHistogramNameNice][sampleCategoryNameNice][sSelTagName] = OD()
            sSelTagForNumerator   = sSelTag_details_dict[sHistNumerator]
            sSelTagForDenominator = sSelTag_details_dict[sHistDenominator]
            efficiencyHistograms_dict[efficiencyHistogramNameNice][sampleCategoryNameNice][sSelTagName][sSelTagNameNice] = sSelTag_details_dict[sSelTagNameNice]

            #efficiencyHistName_toUse = 'hEfficiency_%s_%s' % (efficiencyHistogramNameNice, sampleCategoryNameNice)
            #efficiencyHistograms_dict[efficiencyHistName_toUse] = copy.deepcopy(efficiencyHistogramNameNice_dict[efficiencyHistogramNameNice])
            #efficiencyHistograms_dict[efficiencyHistName_toUse][sPlotTag] = sampleCategoryNameNice
            #print(f"{efficiencyHistogramNameNice = }, {sampleCategoryNameNice = }, {efficiencyHistName_toUse = }")
            isMC = True
            for sDataStr in ['SingleMuon', 'JetHT']:
                if sDataStr in sampleCategoryNameNice: isMC = False
            sSystematics = 'central' if isMC else 'noweight'

            # append histogroms from the same category, which needs to be hadded, into a list
            efficiencyHistograms_dict[efficiencyHistogramNameNice][sampleCategoryNameNice][sSelTagName][sHistNumerator]   = []
            efficiencyHistograms_dict[efficiencyHistogramNameNice][sampleCategoryNameNice][sSelTagName][sHistDenominator] = []            
            for sampleCategory in sampleCategory_dict[sampleCategoryNameNice]: 
                efficiencyHistograms_dict[efficiencyHistogramNameNice][sampleCategoryNameNice][sSelTagName][sHistNumerator].append({
                    sIpFileNameNice: "fIp1",
                    sHistName:       'evt/%s/%s_%s_%s' % (sampleCategory, efficiencyHistogramNameNice, sSelTagForNumerator, sSystematics)
                })
                efficiencyHistograms_dict[efficiencyHistogramNameNice][sampleCategoryNameNice][sSelTagName][sHistDenominator].append({
                    sIpFileNameNice: "fIp1",
                    sHistName:       'evt/%s/%s_%s_%s' % (sampleCategory, efficiencyHistogramNameNice, sSelTagForDenominator, sSystematics)
                })
            






