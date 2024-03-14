import os
import numpy as np
from collections import OrderedDict as OD

sXRange = "xAxisRange"; sYRange = "yAxisRange";
sXLabel = 'xAxisLabel'; sYLabel = 'yAxisLabel';
sXScale = 'xAxisScale';
sNRebin = "nRebin"
sHistosToOverlay = 'histosToOverlay'
sHistosToHadd = 'histosToHadd'
sIpFileNameNice = 'ipFileNameNice'
sHistName   = 'histogramName'

mA_list_full = [8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.5, 13.0, 13.5, 14.0, 16.0, 17.0, 18.5, 21.5, 23.0, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5]
#mA_list = [ 11.0, 11.5, 12.5, 13.0, 13.5, 14.0, 16.0,]
mA_dict = OD([
    ('mA_8To12', [ 8.5, 9.0, 10.0, 11.0, 12.5,]),
    ('mA_11To16', [ 11.0, 11.5, 12.5, 13.0, 13.5, 14.0, 16.0,]),
    ('mA_16To30', [ 16.0, 17.0, 18.5, 21.5, 23.0, 27.5]),
    ('mA_32To62', [ 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5]),     
    ('mA_selective', [ 8.5, 10.5, 17.0, 37.5, 62.5]),

])

sIpFiles = OD()

sAnaVersion = 'CheckNewSignal'

for mA_ in mA_list_full:
    sMA_ = '%.1f' % (mA_)
    sIpFiles['mA_%s'%(sMA_)] = "/eos/cms/store/user/ssawant/htoaa/analysis/20240222_CheckNewSignalMC/2018/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-%s_TuneCP5_13TeV_madgraph_pythia8/analyze_htoaa_SUSY_GluGluH_01J_HToAATo4B_Pt150_M-%s_TuneCP5_13TeV_madgraph_pythia8_0_0.root" % (sMA_, sMA_)

#sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240222_CheckNewSignalMC/2018/plots' 
sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240228_CheckNewSignalMC/2018/plots' 

sSelection = 'sel_leadingFatJetPt'

if not os.path.exists(sOpDir):
    os.makedirs(sOpDir)

histograms_dict = OD()

'''
for mA_ in mA_list:
    sMA_ = '%.1f' % (mA_)
    histograms_dict["hmassA_Signal_mA%s"%(sMA_)] = {
        sXLabel: 'm(A) [GeV]', sYLabel: 'Events',
        #sXRange: [0, 40],
        #sNRebin: 2,
        sHistosToOverlay: OD([# ('h1', [{histoToOverlay1}]),  ('h2', [{histoToOverlay2}]), ('h3', [{histoTohadd3p1}, {histoTohadd3p2}, ...]), ('hi',[{}]), ...]
            ("hmA_%s"%(sMA_), [
                {sIpFileNameNice: 'mA_%s'%(sMA_), sHistName: 'evt/SUSY_GluGluH_01J_HToAATo4B_M-%s_HPtAbv150/hLeadingFatJetParticleNet_massA_Hto4b_avg_v013_Presel_central'%(sMA_)}, # SRWP40
            ]),
        ])
    }
'''

for s_mA_name, mA_list in mA_dict.items():
    sHistName_0 = "hmassA_Signal_%s"%(s_mA_name)
    mA_min = min(mA_list)
    mA_max = max(mA_list)
    mA_span = mA_max - mA_min
    xRangeMin = max(mA_min - mA_span*0.2, 5)
    xRangeMax = min(mA_max + mA_span*0.2, 70)
    histograms_dict[sHistName_0] = {
        sXLabel: 'm(A) [GeV]', sYLabel: 'Events',
        sXRange: [xRangeMin, xRangeMax],
        #sNRebin: 2,
        sHistosToOverlay: OD(),
    }
    for mA_ in mA_list:
        sMA_ = '%.1f' % (mA_)
        histograms_dict[sHistName_0][sHistosToOverlay]["hmA_%s"%(sMA_)] = [
            {sIpFileNameNice: 'mA_%s'%(sMA_), sHistName: 'evt/SUSY_GluGluH_01J_HToAATo4B_M-%s_HPtAbv150/hLeadingFatJetParticleNet_massA_Hto4b_avg_v013_%s_central'%(sMA_, sSelection)}, # SRWP40
        ]


for s_mA_name, mA_list in mA_dict.items():
    sHistName_0 = "hmassH_Signal_%s"%(s_mA_name)
    mA_min = min(mA_list)
    mA_max = max(mA_list)
    mA_span = mA_max - mA_min
    histograms_dict[sHistName_0] = {
        sXLabel: 'm(H) [GeV]', sYLabel: 'Events',
        sXRange: [50, 200],
        sNRebin: 4,
        sHistosToOverlay: OD(),
    }
    for mA_ in mA_list:
        sMA_ = '%.1f' % (mA_)
        histograms_dict[sHistName_0][sHistosToOverlay]["hmA_%s"%(sMA_)] = [
            {sIpFileNameNice: 'mA_%s'%(sMA_), sHistName: 'evt/SUSY_GluGluH_01J_HToAATo4B_M-%s_HPtAbv150/hLeadingFatJetParticleNet_massH_Hto4b_avg_v0123_%s_central'%(sMA_, sSelection)}, # SRWP40
        ]


for s_mA_name, mA_list in mA_dict.items():
    sHistName_0 = "hLeadingFJPNetMD_Signal_Hto4b_Htoaa4bOverQCD_%s"%(s_mA_name)
    mA_min = min(mA_list)
    mA_max = max(mA_list)
    mA_span = mA_max - mA_min
    histograms_dict[sHistName_0] = {
        sXLabel: 'LeadingFJPNetMD_Signal__Hto4b_Htoaa4bOverQCD', sYLabel: 'Events',
        sXRange: [0, 1], #[0.8, 1],
        sNRebin: 10,
        sHistosToOverlay: OD(),
    }
    for mA_ in mA_list:
        sMA_ = '%.1f' % (mA_)
        histograms_dict[sHistName_0][sHistosToOverlay]["hmA_%s"%(sMA_)] = [
            {sIpFileNameNice: 'mA_%s'%(sMA_), sHistName: 'evt/SUSY_GluGluH_01J_HToAATo4B_M-%s_HPtAbv150/hLeadingFatJetParticleNetMD_Hto4b_Htoaa4bOverQCD_%s_central'%(sMA_, sSelection)}, # SRWP40
        ]