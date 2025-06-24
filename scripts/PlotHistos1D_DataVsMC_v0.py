# %%
#from IPython.display import display, HTML
#display(HTML("<style>.container { width:100% !important; }</style>"))

# to run .ipynb from terminal, run: ipython -c "%run <notebook>.ipynb"

# %%
import os, sys
from collections import OrderedDict as OD
import enum
from parse import *
import math
import numpy as np
#import uproot3
import uproot as uproot
import hist

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import mplhep as hep

from hist.intervals import ratio_uncertainty

#sys.path.insert(1, '../') # to import file from other directory (../ in this case)
sys.path.append( os.path.abspath('../') )
print(f"{os.path.abspath('../') = }")

from htoaa_Settings import *
from htoaa_CommonTools import (
    rebinTH1, rebinTH2, variableRebinTH1,
)

class DataBlindingOptions(enum.Enum):
    BlindPartially = '(partially blind)'
    BlindFully     = '(blind)'
    Unblind        = ' '


#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20240627_gg0l_1/2018/analyze_htoaa_stage1.root'
#sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240627_gg0l_1/2018/plots_tmp'
#sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250122_gg0l_DataMC_1/2018/analyze_htoaa_stage1.root'
#sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250122_gg0l_DataMC_1/2018/plots'
#sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250122_gg0l_DataMC_1/2018/plots_1'
sIpFile = '/eos/cms/store/user/ssawant/htoaa/analysis/20250213_gg0l_H34bCat/2018/analyze_htoaa_stage1.root'
sOpDir  = '/eos/cms/store/user/ssawant/htoaa/analysis/20250213_gg0l_H34bCat/2018/plots_1'


CAT = 'gg0lIncl_H34b' # 'gg0l', 'VBFjj', 'Wlv', 'Zll', 'Zvv',  'Vjj'. 'ZvvIncl','ZvvLo', 'ZvvHi', 'gg0lIncl', 'gg0lLo', 'gg0lHi', 'tt0l', 
# 'tt0l_ge1NonHFatJet_0BExtra', 'tt0l_ge1NonHFatJet_1BExtra', 'tt0l_ge1NonHFatJet_ge2BExtra', 'tt0l_0NonHFatJet_ge2B'

sSubcategory = '%s_' % (CAT) if CAT in ['ZvvIncl', 'ZvvLo', 'ZvvHi', 'gg0lIncl', 'gg0lLo', 'gg0lHi', 'tt0l_ge1NonHFatJet_0BExtra'] else '' # 'ZvvIncl' # 'gg0l', 'VBFjj', 'Wlv', 'Zll', 'Zvv',  'Vjj'
#selectionTags = ['%sSBWP80to40' % (sSubcategory)]
selectionTags = ['gg0lIncl_H34b_Xto34bv2_SRWP80'] # ['Presel', 'gg0lIncl_SRWP60', 'gg0lIncl_SBWP60'] #['tt0l_ge1NonHFatJet_1BExtra_Hi_SBWP95to60', 'tt0l_ge1NonHFatJet_1BExtra_Med_SBWP95to60', 'tt0l_ge1NonHFatJet_1BExtra_Lo_SBWP95to60']


#from HistogramListForPlottingDataVsMC_TriggerStudy_GGFMode import *
#from HistogramListForPlottingDataVsMC_Analysis_GGFMode import *
#from HistogramListForPlottingDataVsMC_Analysis_VHHadronicMode import *
#from HistogramListForPlottingDataVsMC_Analysis_ZH_4b2nu import *
#from HistogramListForPlottingDataVsMC_Analysis_Example import *

if 'gg0l'   in CAT: from HistogramListForPlottingDataVsMC_Analysis_GGFMode               import *
if 'Vjj'    in CAT: from HistogramListForPlottingDataVsMC_Analysis_VHHadronicMode        import *
if 'Zvv'    in CAT: from HistogramListForPlottingDataVsMC_Analysis_ZH_4b2nu              import *
if 'tt0l'   in CAT: from HistogramListForPlottingDataVsMC_Analysis_ttHHadronicMode       import *

cmsWorkStatus                  = 'Work in Progress'
era                            = '2018'
luminosity_total               = Luminosities_TotalPerYear[era][HLT_toUse][0] # 54.54  #59.83
dataBlindOption                = DataBlindingOptions.BlindPartially # DataBlindingOptions.BlindPartially , DataBlindingOptions.BlindFully , DataBlindingOptions.Unblind
#significantThshForDataBlinding = 4 # 0.125 # blind data in bins with S/sqrt(B) > significantThshForDataBlinding while running with dataBlindOption = DataBlindingOptions.BlindPartially
significantThshForDataBlinding = 10 # for significance Z

sOpDir = '%s/%s' % (sOpDir, CAT)
if not os.path.exists(sOpDir):
    os.makedirs(sOpDir)
    
fIpFile = uproot.open(sIpFile)

# %% [markdown]
# 

# %%
def getNonZeroMin(arr):
    min_ = 1e20
    #a_   = arr[np.nonzero(arr)]
    a_ = arr[ np.argwhere(arr > 0) ]
    if len(a_) > 0:
        min_ = np.min( a_ )
    return min_


# %%
# Function to draw box error bars
# https://matplotlib.org/stable/gallery/statistics/errorbars_and_boxes.html#sphx-glr-gallery-statistics-errorbars-and-boxes-py
def make_error_boxes(ax, xdata, ydata, xerror, yerror, 
                     facecolor='lightgrey',
                     edgecolor='none', alpha=0.5, hatch='////', linewidth=0
                     #kwagrs_
                     ):

    # Loop over data points; create box from errors at each point
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Rectangle.html
    # matplotlib.patches.Rectangle(xy, width, height, *, angle=0.0, rotation_point='xy', **kwargs)
    #errorboxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
    #              for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T)]
    errorboxes = [Rectangle((x - xe, y - ye), 2*xe, 2*ye)
                  for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T)]

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor, hatch=hatch, linewidth=linewidth)

    # Add collection to axes
    ax.add_collection(pc)

    artists = None
    # Plot errorbars
    #artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
    #                      fmt='none', ecolor=facecolor)

    return artists


## Calculate significance
def calSignificance1(S, B):
    significance = np.where(
        B > 1e-10,
        np.sqrt( 2 * ((S+B)*np.log(1 + (S/B)) - S) ),
        np.full_like(S, 1e-6)
    )
    return significance

def calSignificance2(S, B, Bvariance):
    denom = np.sqrt(B + Bvariance)
    significance = np.where(
        denom > 0,
        S / denom,
        np.full_like(S, 1e-6)
    )
    return significance

# %%
#colors_bkg_list = ['blue', 'orange', 'brown'] # ["#9b59b6", "#e74c3c", "#34495e", "#2ecc71"] #['lightcoral', 'burlywood', 'cyan', 'saddlebrown', 'slateblue', 'lightpink', 'darkkhaki', 'antiquewhite', 'limegreen', 'violet', 'firebrick', 'darkorchid', 'tan', 'olive', 'purple']

colors_bkg_list_NonCMS = [ 
    # ['color', <transperent>, '<fill pattern>']
    ["#3f90da",    0.7,  ''],
    ["#ffa90e",    0.7,  ''],

    ['lightcoral',    0.7,  ''],
    ['cyan',          0.7,  '' ],
    ['burlywood',     0.7,  '' ],     
    ['slateblue',     0.7,  '' ],
    ['saddlebrown',   0.7,  '' ],
    ['lightpink',     0.7,  'xx' ],
    ['darkkhaki',     0.9,  '' ],
    ['antiquewhite',  0.9,  '//' ],
    ['limegreen',     0.6,  '' ],
    ['violet',        0.4,  'oo' ],
    ['lightskyblue',  0.4,  '||' ],    
    ['firebrick',     0.4,  '--' ],
    ['rosybrown',     0.4,  '..' ],
    ['darkorchid',    0.4,  '' ],
    ['tan',           0.9,  '' ],
    ['olive',         0.9,  '' ],
    ['purple',        0.9,  ''],

    ['gainsboro',      0.7,  '.'],
    ['rosybrown',      0.7,  ''],
    ['cadetblue',      0.7,  'o'],
    ['oldlace',        0.7,  ''],
    ['palevioletred',  0.7,  ''],
    ['sandybrown',     0.7,  ''],

    ['limegreen',     0.6,  'xx' ],
    ['violet',        0.4,  '--' ],
    ['lightskyblue',  0.4,  '.' ],    
    ['firebrick',     0.4,  '//' ],
    ['rosybrown',     0.4,  '||' ],    
]

colors_sig_list_NonCMS = [
    # ['color', <transperent>, '<fill pattern>', ]
    ['blue',          0.9,  ''],
    ['red',           0.9,  ''],
    ['green',         0.9,  ''],
    ['magenta',       0.9,  ''],
    ['orange',        0.9,  ''],
]

## CMS color schemes: https://gitlab.cern.ch/cms-analysis/analysisexamples/plotting-demo/-/blob/master/1-tutorial_CAT_recommendations.ipynb?ref_type=heads
# 6-color scheme: ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]
# 10-color scheme: "#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]
colors_bkg_list = [ 
    # ['color', <transperent>, '<fill pattern>']
    ["#3f90da",    1,  ''],
    ["#ffa90e",    1,  ''],
    ["#94a4a2",    1,  ''],
    ["#a96b59",    1,  ''],
    ["#b9ac70",    1,  ''],
    ["#717581",    1,  ''],
    ["#92dadd",    1,  ''],   

]

colors_sig_list = [
    # ['color', <transperent>, '<fill pattern>', ]
    ["#bd1f01",    1,  ''],
    ["#832db6",    1,  ''],
    ["#e76300",    1,  ''],
    
]

#errps = {'hatch':'////', 'facecolor':'none', 'lw': 0, 'edgecolor': 'k', 'alpha': 0.5}
errps = {'hatch':'////', 'facecolor':'none', 'linewidth': 0, 'edgecolor': 'k', 'alpha': 0.5}

PlotRatioPlot = True
PlotSignificancePlot = False #True


hep.style.use("CMS")

for sData, ExpData_list in ExpData_dict.items():
    luminosity_toUse = 0
    for ExpData_component in ExpData_list:
        DatasetEra_         = ExpData_component.split(era)[1][0] # 'JetHT_Run2018A'.split('2018')[1][0]
        luminosity_forEra_  = Luminosities_TotalPerYear_perEra[era][HLT_toUse][DatasetEra_]
        luminosity_toUse   += luminosity_forEra_
        print(f"{ExpData_list = }, {DatasetEra_ = }, {luminosity_forEra_ = } ")
    luminosity_Scaling_toUse = round(luminosity_toUse, 2) / round(luminosity_total, 2)
    luminosity_toUse = round(luminosity_toUse, 1)
    print(f"{sData}: {ExpData_list}, {luminosity_toUse = }, {luminosity_total = },  {luminosity_Scaling_toUse = }")

    for selectionTag in selectionTags:    
        #dataBlindOption_toUse = dataBlindOption if selectionTag != 'SR' else DataBlindingOptions.BlindPartially

        for histo_name in histograms_dict.keys():
            dataBlindOption_toUse = dataBlindOption
            if 'ParticleNet_massA_Hto4b' in histo_name:
                dataBlindOption_toUse = DataBlindingOptions.BlindFully

            histo_name_toUse = '%s_%s' % (histo_name, selectionTag)
            for systematic in systematics_list:
                for yAxisScale in ['linearY', 'logY']: # ['linearY', 'logY']
                    xAxisRange = histograms_dict[histo_name][sXRange] if sXRange in histograms_dict[histo_name].keys() else None
                    yAxisRange = histograms_dict[histo_name][sYRange] if sYRange in histograms_dict[histo_name].keys() else None
                    xAxisLabel = histograms_dict[histo_name][sXLabel] if sXLabel in histograms_dict[histo_name].keys() else None
                    yAxisLabel = histograms_dict[histo_name][sYLabel] if sYLabel in histograms_dict[histo_name].keys() else None
                    nRebinX    = histograms_dict[histo_name][sNRebinX] if sNRebinX in histograms_dict[histo_name].keys() else 1
                    nRebinY    = histograms_dict[histo_name][sNRebinY] if sNRebinY in histograms_dict[histo_name].keys() else 1
                    XRebinning = histograms_dict[histo_name][sXRebinning] if sXRebinning in histograms_dict[histo_name].keys() else None
                    YRebinning = histograms_dict[histo_name][sYRebinning] if sYRebinning in histograms_dict[histo_name].keys() else None
                    if yAxisRange and yAxisRange[0] > yAxisRange[1]:
                        yAxisRange = None                        

                    nHistoDimemsions = None
                    yAxisRange_cal      = [1e20, -1e10]
                    yRatioAxisRange_cal = [1e20, -1e10]
                    ySignfAxisRange_cal = [1e20, -1e10]                    
                    xError = np.array([])
                    hData = None
                    hBkgTot_values = None
                    hBkgTot_variance = None
                    hStack_values_list = np.array([]) 
                    hStack_edges = np.array([])
                    hStack_centers = np.array([])
                    sStack_list = []
                    nBkgTot = 0
                    hBkgTot = None
                    significance_list = [] #np.array([])

                    sEventYieldTable = ''

                    print(f"\n\n {histo_name_toUse = }, {systematic = }, {yAxisScale = }, ")
                    #fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(8,10), sharex='col', gridspec_kw={'height_ratios': [3, 1]}, subplot_kw={'ymargin': 0.4})
                    ###fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(8,10), sharex='col', gridspec_kw={'height_ratios': [4, 1], 'hspace': 0})
                    #fig, ax = plt.subplots(ncols=1, nrows=3, figsize=(8,10), sharex='col')
                    #print(f"fig: {fig}, axs: {axs}")

                    if PlotRatioPlot and (not PlotSignificancePlot):
                        fig, (axTop, axRatio) = plt.subplots(2, 1, gridspec_kw=dict(height_ratios=[3, 1], hspace=0.1), sharex=True)
                    if PlotRatioPlot and PlotSignificancePlot:
                        fig, (axTop, axRatio, axSignf) = plt.subplots(3, 1, gridspec_kw=dict(height_ratios=[3.5, 0.5, 0.5], hspace=0.1), sharex=True)

                    #fig1, ax1 = plt.subplots()
                    
                    histos_dict = OD()
                    mask_DataBlindedBins = None

                    
                    #if len(MCBkg_list) > 0:
                    if len(list(MCBkg_dict.keys())) > 0:
                        hBkg_list = []
                        sBkg_list = []
                        hBkg_integral_list = []
                        for i_, (MCBkgNameShort, MCBkg_list) in enumerate(MCBkg_dict.items()):
                            h = None
                            for dataset in MCBkg_list:
                                histo_name_toUse_full = 'evt/%s/%s_%s' % (dataset, histo_name_toUse, systematic)
                                #print(f"{histo_name_toUse_full = }")
                                h_i = fIpFile[histo_name_toUse_full].to_hist()
                                nHistoDimemsions = len(h_i.axes)
                                if nHistoDimemsions == 2 and yAxisScale == 'logY': break  # No need to plot 2-D hist with logY
                                if isinstance(XRebinning, list) or isinstance(XRebinning, (np.ndarray, np.generic)):
                                    h_i = variableRebinTH1(h_i, XRebinning)  if nHistoDimemsions == 1 else h_i
                                else:
                                    h_i = rebinTH1(h_i, nRebinX) if nHistoDimemsions == 1 else rebinTH2(h_i, nRebinX, nRebinY)
                                    #h_i = h_i.rebin(nRebinX) if nHistoDimemsions == 1 else rebinTH2(h_i, nRebinX, nRebinY)

                                if dataset == MCBkg_list[0]:  h = h_i
                                else:                         h = h + h_i
                            
                                

                            h = h * luminosity_Scaling_toUse

                            hBkgTot = h 
                            if i_ == 0: 
                                hBkgTot = h
                            else:
                                hBkgTot = hBkgTot + h 

                            nTot_ = h.values().sum()
                            hBkg_list.append(h)
                            sBkg_list.append(MCBkgNameShort)
                            hBkg_integral_list.append(nTot_)

                            histos_dict[MCBkgNameShort] = h 
                            if not isinstance(mask_DataBlindedBins, np.ndarray):
                                mask_DataBlindedBins = np.full_like(h.values(), False)

                            if abs(nTot_ - 0) < 1e-10: continue

                            if nHistoDimemsions == 1:
                                mask_XRange = ((h.axes.centers[0] >= xAxisRange[0]) & (h.axes.centers[0] <= xAxisRange[1])) if xAxisRange else np.full_like(h.values(), True)
                                yMin_ = getNonZeroMin(h.values()[mask_XRange])
                                yMax_ = np.max(h.values()[mask_XRange])
                                if yMin_ < yAxisRange_cal[0]:
                                    yAxisRange_cal[0] = yMin_
                                if yMax_ > yAxisRange_cal[1]:
                                    yAxisRange_cal[1] = yMax_                        

                        # No need to plot 2-D hist with logY
                        if nHistoDimemsions == 2 and yAxisScale == 'logY': 
                            plt.close(fig)
                            continue 


                        # sort histograms in decreasing yield
                        isReverseSortForStack = True
                        idx_hBkg_sortedByIntegral = sorted(range(len(hBkg_integral_list)), key=lambda i: hBkg_integral_list[i], reverse=isReverseSortForStack)            

                        hStack_list = [ hBkg_list[idx] for idx in idx_hBkg_sortedByIntegral ]  
                        sStack_list = [ sBkg_list[idx] for idx in idx_hBkg_sortedByIntegral ]  

                        hStack_values_list    = np.array( [ h.values() for h in hStack_list ] )
                        hStack_variance_list  = np.array( [ h.variances() for h in hStack_list ] )
                        hStack_error_list     = np.array( [ np.sqrt(h.variances()) for h in hStack_list ] )
                        hStack_edges          = hStack_list[0].axes[0].edges
                        hStack_centers        = hStack_list[0].axes[0].centers
                        xError                = (hStack_list[0].axes[0].edges[1:] - hStack_list[0].axes[0].edges[0:-1]) / 2 if len(xError) == 0 else xError

                        hBkgTot_values        = np.sum(hStack_values_list, axis=0)
                        hBkgTot_variance      = np.sum(hStack_variance_list, axis=0)

                        # Update yRange for hStackBkg -------
                        if nHistoDimemsions == 1:
                            #mask_XRange = ((h.axes.centers[0] >= xAxisRange[0]) & (h.axes.centers[0] <= xAxisRange[1])) if xAxisRange else np.full_like(h.values(), True)
                            yMin_ = getNonZeroMin(hBkgTot_values[mask_XRange])
                            yMax_ = np.max(hBkgTot_values[mask_XRange])
                            if yMin_ < yAxisRange_cal[0]:
                                yAxisRange_cal[0] = yMin_
                            if yMax_ > yAxisRange_cal[1]:
                                yAxisRange_cal[1] = yMax_    

                        nHists = len(MCBkg_list)
                        colors_toUse = [ colors_bkg_list[i][0] for i in range(nHists) ]
                        alpha_toUse  = [ colors_bkg_list[i][1] for i in range(nHists) ]
                        hatch_toUse  = [ colors_bkg_list[i][2] for i in range(nHists) ]

                       
                        if nHistoDimemsions == 1: # 1-D histogram
                            hep.histplot(
                                hStack_values_list, 
                                bins=hStack_edges, 
                                ax=axTop, 
                                histtype='fill', 
                                stack=True, 
                                label=sStack_list, 
                                color=colors_toUse,
                                #alpha=alpha_toUse,
                                #hatch=hatch_toUse,
                                sort='yield'
                                )

                            #hep.histplot(hBkgTot_values, histtype='band', ax=axTop, **errps)   
                            make_error_boxes(
                                ax=axTop, 
                                xdata=hStack_centers, 
                                ydata=hBkgTot_values, 
                                xerror=xError, 
                                yerror=np.sqrt(hBkgTot_variance), 
                                **errps
                                )   
                            
                        elif nHistoDimemsions == 2 and 1==0: # 2-D histogram  
                            hep.hist2dplot(
                                hBkgTot_values,
                                xbins=hStack_list[0].axes[0].edges,
                                ybins=hStack_list[0].axes[1].edges,
                                #labels='Bkg_total',
                                cmin=getNonZeroMin(hStack_list[0].values()),
                                ax=axTop
                            )   
                        
                        # No. of events in total background
                        nBkgTot = np.sum(hBkgTot_values)



                    if len(MCSig_list) > 0:
                        hSig_list = []
                        sSig_list = []
                        hSig_integral_list = []
                        for iSig, dataset in enumerate(MCSig_list):
                            histo_name_toUse_full = 'evt/%s/%s_%s' % (dataset, histo_name_toUse, systematic)
                            h = fIpFile[histo_name_toUse_full].to_hist()
                            h = rebinTH1(h, nRebinX) if nHistoDimemsions == 1 else rebinTH2(h, nRebinX, nRebinY)

                            h = h * luminosity_Scaling_toUse

                            nTot_ = h.values().sum()
                            hSig_list.append(h)
                            sSig_list.append(dataset)
                            hSig_integral_list.append(h.values().sum())
                            #print(f"{histo_name_toUse_full} integral: {h.values().sum()}")

                            histo_edges = h.axes[0].edges
                            xError      = (h.axes[0].edges[1:] - h.axes[0].edges[0:-1]) / 2 if len(xError) == 0 else xError


                            label_MCSig = dataset
                            label_MCSig = sLableSig[iSig]
                            if abs(scale_MCSig - 1) > 1e-6:
                                if scale_MCSig >= 1:
                                    label_MCSig = '%s x %d' % (label_MCSig, scale_MCSig)
                                else:
                                    label_MCSig = '%s x %g' % (label_MCSig, scale_MCSig)
                                
                            if nHistoDimemsions == 1:
                                mask_XRange = ((h.axes.centers[0] >= xAxisRange[0]) & (h.axes.centers[0] <= xAxisRange[1])) if xAxisRange else np.full_like(h.values(), True)
                                yMin_ = getNonZeroMin(h.values()[mask_XRange])
                                yMax_ = np.max(h.values()[mask_XRange])
                                if yMin_ < yAxisRange_cal[0]:
                                    yAxisRange_cal[0] = yMin_
                                if yMax_ > yAxisRange_cal[1]:
                                    yAxisRange_cal[1] = yMax_                        

                            # plot signal
                            if nHistoDimemsions == 1:
                                hep.histplot(
                                    h.values() * scale_MCSig, 
                                    bins=histo_edges, 
                                    ax=axTop, 
                                    yerr=np.sqrt(h.variances()) * scale_MCSig, 
                                    histtype='step', #'errorbar', 
                                    label=label_MCSig,
                                    color=colors_sig_list[iSig][0],                             
                                    #marker='o',
                                    #markerfacecolor=colors_sig_list[iSig][0],
                                    #markersize=3
                                    )

                            nSig = np.sum(h.values())
                            #print(f"{iSig = }, {dataset = } {nSig = }, {nBkgTot = }")
                            # S/sqrt(B) or S/sqrt(S+B)
                            if nSig > 0 and nBkgTot > 0:
                                #S_ = h.values() / nSig
                                #B_ = np.sqrt(hBkgTot_values / nBkgTot)
                                #significance_i = np.divide(S_, B_, where=B_!=0, out=np.zeros(B_.shape))
                                #B_ = hBkgTot_values / nBkgTot
                                #SB_ = np.sqrt( S_ + B_ )
                                #significance_i = np.divide(S_, SB_, where=SB_!=0, out=np.zeros(SB_.shape))

                                S_ = h.values()
                                B_ = np.sqrt(hBkgTot_values)
                                #significance_i = np.divide(S_, B_, where=B_!=0, out=np.zeros(B_.shape))
                                significance_i = calSignificance1(S_, hBkgTot.values())
                                #significance_i = calSignificance2(S_, hBkgTot.values(), hBkgTot.variances())

                                # set high significant when S_ > 0 and B_ = 0
                                significance_i = np.where(
                                    np.logical_and(S_ > 0, hBkgTot_values < 1e-6),
                                    np.full(B_.shape, 10000),
                                    significance_i)
                                significance_list.append(significance_i)

                        
                        significanceMax = np.array(significance_list)
                        #print(f"{significanceMax = }")
                        #significanceMax = np.sum(significanceMax, axis=0)
                        #significanceMax = np.divide(significanceMax, len(MCSig_list) )
                        significanceMax = np.max(significanceMax, axis=0)
                        #print(f"{significanceMax = }")

                        #print(f"significanceMax (max: {np.max(significanceMax)}): {significanceMax}")                        




                    #print(f"\nAfter MCSig {yAxisRange_cal = }")
                    
                    if dataBlindOption_toUse in [DataBlindingOptions.Unblind, DataBlindingOptions.BlindPartially]: #sData:
                        hData = None
                        for ExpData_component in ExpData_list:
                            histo_name_toUse_full = 'evt/%s/%s_%s' % (ExpData_component, histo_name_toUse, systematics_forData)
                            h = fIpFile[histo_name_toUse_full].to_hist()
                            if hData == None: hData = h
                            else:             hData = hData + h

                        hData = rebinTH1(hData, nRebinX) if nHistoDimemsions == 1 else rebinTH2(hData, nRebinX, nRebinY)
                        xError = (hData.axes[0].edges[1:] - hData.axes[0].edges[0:-1]) / 2

                        if nHistoDimemsions == 1:
                            mask_XRange = ((hData.axes.centers[0] >= xAxisRange[0]) & (hData.axes.centers[0] <= xAxisRange[1])) if xAxisRange else np.full_like(hData.values(), True)
                            yMin_ = getNonZeroMin(hData.values()[mask_XRange])
                            yMax_ = np.max(hData.values()[mask_XRange])
                            if yMin_ < yAxisRange_cal[0]:
                                yAxisRange_cal[0] = yMin_
                            if yMax_ > yAxisRange_cal[1]:
                                yAxisRange_cal[1] = yMax_
                            #print(f"Data: {yMin_ = }, {yMin_}")

                        hData_values_toUse = hData.values()
                        hData_errors_toUse = np.sqrt(hData.variances())
                        histos_dict['Data'] = hData

                        
                        hData_values_toUse = np.where(
                            hData_values_toUse >= 1,
                            hData_values_toUse,
                            np.full(len(hData_values_toUse), -1),
                        )
                        hData_errors_toUse = np.where(
                            hData_values_toUse >= 1,
                            hData_errors_toUse,
                            np.full(len(hData_values_toUse), 0),
                        )
                        

                        # blind data with high S/sqrt(B) bins
                        #print(f"{len(significanceMax) = }")
                        if dataBlindOption_toUse in [DataBlindingOptions.BlindPartially] and \
                            len(significanceMax):
                            # inflate significantThshForDataBlinding for higher S/sqrt(B) when histogram is rebinned, 
                            # so that blinding of data is independent of rebinning
                            #significantThshForDataBlinding_toUse = significantThshForDataBlinding * math.sqrt(nRebinX)
                            significantThshForDataBlinding_toUse = significantThshForDataBlinding
                            mask_DataBlindedBins = (significanceMax > significantThshForDataBlinding_toUse)
                            hData_values_toUse = np.where(
                                (significanceMax > significantThshForDataBlinding_toUse),
                                np.full(len(hData_values_toUse), 0),
                                hData_values_toUse
                            )
                            hData_errors_toUse = np.where(
                                (significanceMax > significantThshForDataBlinding_toUse),
                                np.full(len(hData_values_toUse), 0),
                                hData_errors_toUse
                            )
                            #print(f"{(significanceMax > significantThshForDataBlinding_toUse) =}")
                            #print(f"Data blinding x-values: { hData.axes[0].centers[(significanceMax > significantThshForDataBlinding_toUse)] }")
                            #print(f"hData_values_toUse ({len(hData_values_toUse)}): {hData_values_toUse}")

                            hData_values_toUse = np.where(
                                mask_DataBlindedBins,
                                np.full(len(hData_values_toUse), -1),
                                hData_values_toUse
                            )

                        #print(f"{hData_values_toUse = }")
                        if nHistoDimemsions == 1:
                            #hep.histplot(hData.values(), bins=hData.axes[0].edges, ax=axTop, yerr=np.sqrt(hData.variances()), histtype='errorbar', color='black', label='Data')
                            hep.histplot(
                                hData_values_toUse, 
                                bins=hData.axes[0].edges, 
                                ax=axTop, 
                                yerr=hData_errors_toUse, 
                                histtype='errorbar', 
                                color='black', 
                                label='%s %s' % (sData, dataBlindOption_toUse.value),
                                capsize=2,
                                )

                            # highlight blinded bins
                            axTop.plot(
                                hData.axes[0].centers[mask_DataBlindedBins],
                                np.zeros_like(hData.axes[0].centers)[mask_DataBlindedBins],
                                label='Data blinded bins',
                                color='red', 
                                marker='x',
                                markerfacecolor='red',
                                markersize=8
                            )    

                        elif nHistoDimemsions == 2 and 1==0: # 2-D histogram  
                            hep.hist2dplot(
                                hData_values_toUse,
                                xbins=hData.axes[0].edges,
                                ybins=hData.axes[1].edges,
                                #labels='Bkg_total',
                                cmin=getNonZeroMin(hData_values_toUse),
                                ax=axRatio
                            )                                              

                        #print(f"hData integral: {hData.values().sum()}")


                        # Ratio plot ---------------------------------------------------------       
                        ratio_values = np.divide(hData_values_toUse, hBkgTot_values, where=hBkgTot_values!=0, out=np.full(hData.shape[0], -1, dtype=float))
                        ratio_values_toUse = np.divide(hData_values_toUse, hBkgTot_values, where=hBkgTot_values!=0, out=np.full(hData.shape[0], -9999, dtype=float))
                        ratio_error  = hData_errors_toUse            
                        ratio_error  = np.divide(ratio_error, hBkgTot_values, where=hBkgTot_values!=0, out=np.zeros(hData.shape))
                        ratio_syst   = np.sqrt(hBkgTot_variance)
                        ratio_syst   = np.divide(ratio_syst, hBkgTot_values, where=hBkgTot_values!=0, out=np.zeros(hData.shape))
                        ratio_syst_CMS = ratio_uncertainty(hData_values_toUse, hBkgTot_values, 'poisson-ratio')

                        #print(f"{ratio_syst      = }")
                        #print(f"{ratio_syst_CMS = }")

                        #print(f"{list(zip(ratio_syst, ratio_syst_CMS[0], ratio_syst_CMS[1])) = }")

                        #print(f"ratio_values ({ratio_values.shape}): {ratio_values}")
                        if nHistoDimemsions == 1:
                            yMin_ = getNonZeroMin( ratio_values[mask_XRange] - ratio_error[mask_XRange])
                            yMax_ = np.max( ratio_values[mask_XRange] + ratio_error[mask_XRange])
                            if yMin_ < yRatioAxisRange_cal[0]:
                                yRatioAxisRange_cal[0] = yMin_
                            if yMax_ > yRatioAxisRange_cal[1]:
                                yRatioAxisRange_cal[1] = yMax_                          
                        
                        if nHistoDimemsions == 1:
                            hep.histplot(
                                ratio_values_toUse, 
                                bins=hData.axes[0].edges, 
                                ax=axRatio, 
                                yerr=ratio_error, 
                                histtype='errorbar', 
                                color='black', 
                                label='Data',
                                capsize=2,
                                )
                            #if xAxisRange: axRatio.set_xlim(xAxisRange[0], xAxisRange[1])

                            # plot totoal background error bars only for ratio plot
                            '''
                            make_error_boxes(
                                ax=axRatio, 
                                xdata=hData.axes[0].centers, 
                                ydata=np.full(len(hData.axes[0].centers), 1), 
                                xerror=xError, 
                                yerror=ratio_syst, 
                                facecolor='grey',
                                edgecolor='none', 
                                alpha=0.5
                                )
                            '''
                            '''
                            make_error_boxes(
                                ax=axRatio, 
                                xdata=hData.axes[0].centers, 
                                ydata=np.full(len(hData.axes[0].centers), 1), 
                                xerror=xError, 
                                yerror=ratio_syst, 
                                **errps
                                )
                            '''
                            axRatio.stairs(1+ratio_syst_CMS[1], edges=hData.axes[0].edges, baseline=1-ratio_syst_CMS[0], **errps)
                            
                            # highlight blinded bins
                            axRatio.plot(
                                hData.axes[0].centers[mask_DataBlindedBins],
                                np.ones_like(hData.axes[0].centers)[mask_DataBlindedBins],
                                label='Data blinded',
                                color='red', 
                                marker='x',
                                markerfacecolor='red',
                                markersize=8
                            )
                            
                        elif nHistoDimemsions == 2: # 2-D histogram  
                            hep.hist2dplot(
                                ratio_values,
                                xbins=hData.axes[0].edges,
                                ybins=hData.axes[1].edges,
                                #labels='Bkg_total',
                                cmin=yRatioLimit[0], cmax=yRatioLimit[1],
                                ax=axTop
                            )    

                    if yAxisScale == 'linearY' and dataBlindOption_toUse != DataBlindingOptions.BlindFully and 1==0:
                        sEventYieldTable = ''
                        dataName_tmp_ = ''
                        for dataName, histo_ in histos_dict.items():
                            #print(f"{dataName = }, {histo_dict['values'   ].shape = }, {mask_DataBlindedBins.shape = }")
                            nEvents_  = histo_.values()[~ mask_DataBlindedBins].sum()
                            variance_ = histo_.variances()[~ mask_DataBlindedBins].sum()
                            sEventYieldTable += '%s \t %g \t %g \t %g \n' % (dataName, nEvents_, math.sqrt(variance_), variance_)
                            dataName_tmp_ = dataName
                        print(f"Blinded x points: {histos_dict[dataName_tmp_].axes[0].centers[mask_DataBlindedBins] = }")
                        print(f"\n\n\n Event yield table {histo_name_toUse}: \n{sEventYieldTable}\n\n")

                    
                    if PlotSignificancePlot and len(significance_list) > 0:
                        for i_, significance_i in enumerate(significance_list):
                            hep.histplot(
                                significance_i, 
                                bins=hBkgTot.axes[0].edges, 
                                ax=axSignf, 
                                histtype='step', #'errorbar', 
                                #label=label_MCSig,
                                color=colors_sig_list[i_][0],                             
                                #marker='o',
                                #markerfacecolor=colors_sig_list[iSig][0],
                                #markersize=3
                            ) 
                            if nHistoDimemsions == 1:
                                mask_XRange = ((hBkgTot.axes.centers[0] >= xAxisRange[0]) & (hBkgTot.axes.centers[0] <= xAxisRange[1])) if xAxisRange else np.full_like(hBkgTot.values(), True)
                                yMin_ = getNonZeroMin(significance_i[mask_XRange])
                                yMax_ = np.max(significance_i[mask_XRange])
                                if yMin_ < ySignfAxisRange_cal[0]:
                                    ySignfAxisRange_cal[0] = yMin_
                                if yMax_ > ySignfAxisRange_cal[1]:
                                    ySignfAxisRange_cal[1] = yMax_                        



                    
                    # Upper plot cosmetics ---------
                    if xAxisRange: axTop.set_xlim(xAxisRange[0], xAxisRange[1])
                    print(f"\nAt the end {yAxisRange_cal = }")
                    if yAxisRange: axTop.set_ylim(yAxisRange[0], yAxisRange[1])
                    elif nHistoDimemsions == 1:          
                        #yMaxOffset = 10**(math.log10(yAxisRange_cal[1] / abs(yAxisRange_cal[0])) * 0.4) if yAxisScale == 'logY' else 1.6
                        yMaxOffset = 10**(math.log10(yAxisRange_cal[1] / abs(yAxisRange_cal[0])) * 0.55) if yAxisScale == 'logY' else 2.0
                        print(f"{yMaxOffset = }, {yAxisRange_cal[1] * yMaxOffset = }, \t\t {abs(yAxisRange_cal[0]) * logYMinScaleFactor = }")
                        if yAxisScale == 'logY':
                            yAxisRange_cal[0] = abs(yAxisRange_cal[0]) * logYMinScaleFactor
                            yAxisRange_cal[1] = yAxisRange_cal[1] * yMaxOffset
                        else:
                            yAxisRange_cal[0] = yAxisRange_cal[0]
                            yAxisRange_cal[1] = yAxisRange_cal[1] * yMaxOffset
                        print(f"\nAt the end updated {yAxisRange_cal = } \t {yAxisScale = }")
                        axTop.set_ylim(yAxisRange_cal[0], yAxisRange_cal[1])
                    if xAxisLabel:                              axTop.set_xlabel(xAxisLabel)
                    if (PlotRatioPlot or PlotSignificancePlot): axTop.set_xlabel("")
                    if yAxisLabel:                              axTop.set_ylabel(yAxisLabel)       
                    if yAxisScale == 'logY': axTop.set_yscale('log', base=10)
                    handles_, labels_ = axTop.get_legend_handles_labels()         
                    #axTop.legend(reversed(handles_), reversed(labels_), fontsize=14, loc='best', ncol=2, bbox_to_anchor=(-0.1, 0.65, 1.1, 0.36))
                    #axTop.legend(reversed(handles_), reversed(labels_), title='Category: %s'%(CAT), loc='best', ncol=2)
                    axTop.legend(reversed(handles_), reversed(labels_), loc='best', ncol=2)
                    #axTop.legend(reversed(handles_), reversed(labels_), )

                     #axTop.set_ymargin(1.)
                    #axTop.grid()

                    # Ratio plot cosmetics ---------
                    if yRatioAxisRange_cal[0] < yRatioLimit[0]: yRatioAxisRange_cal[0] = yRatioLimit[0]
                    if yRatioAxisRange_cal[1] > yRatioLimit[1]: yRatioAxisRange_cal[1] = yRatioLimit[1]                    
                    yRatioAxisRange_cal_maxDeviation = max(abs(yRatioAxisRange_cal[0] - 1), abs(yRatioAxisRange_cal[1] - 1))
                    yRatioAxisRange_cal[0] = 1 - yRatioAxisRange_cal_maxDeviation
                    yRatioAxisRange_cal[1] = 1 + yRatioAxisRange_cal_maxDeviation
                    yRatioAxisRange_cal[0] = max(yRatioAxisRange_cal[0], 0)
                    if xAxisRange: axRatio.set_xlim(xAxisRange[0], xAxisRange[1]) 
                    axRatio.set_ylim(yRatioAxisRange_cal[0], yRatioAxisRange_cal[1])
                    print(f"{yRatioAxisRange_cal = }") 

                    if xAxisLabel: axRatio.set_xlabel(xAxisLabel)
                    if PlotSignificancePlot: axRatio.set_xlabel("")
                    #axRatio.set_ylabel('Data/MC')
                    axRatio.set_ylabel(r'$\frac{Data}{MC}$')
                    
                    axRatio.axhline(y=1, ls='--', color='k')
                    #axRatio.grid()

                    # Significance plot cosmetics ---------
                    if PlotSignificancePlot:
                        if ySignfAxisRange_cal[0] < ySignfLimit[0]: ySignfAxisRange_cal[0] = ySignfLimit[0]
                        if ySignfAxisRange_cal[1] > ySignfLimit[1]: ySignfAxisRange_cal[1] = ySignfLimit[1] 
                        axSignf.set_ylim(ySignfAxisRange_cal[0], ySignfAxisRange_cal[1])
                        if xAxisRange: axSignf.set_xlim(xAxisRange[0], xAxisRange[1]) 
                        if xAxisLabel: axSignf.set_xlabel(xAxisLabel)
                        axSignf.set_ylabel('Sign.')
                        #axSignf.set_yscale('log')


                    

                    isData = True if dataBlindOption_toUse != DataBlindingOptions.BlindFully else False
                    fontsize_toUse = 18 if isData else 15
                    #hep.cms.label(ax=axTop, data=isData, year=era, lumi=luminosity_toUse, label=cmsWorkStatus, fontsize=fontsize_toUse)
                    hep.cms.label(ax=axTop, data=isData, year=era, lumi=luminosity_toUse, label=cmsWorkStatus)
                    #hep.cms.label("Work in Progress", ax=axTop, data=isData, year=era, lumi=luminosity_toUse, )

                    
                    axTop.text(0.8, 0.51, 'Cat. %s'%(CAT), #selectionTag,
                            fontsize=18, fontstyle='italic',
                            horizontalalignment='center',
                            verticalalignment='center',
                            transform=axTop.transAxes
                            )
                    
                    
                    sOpDir_toUse = '%s/%s' % (sOpDir, selectionTag)
                    if not os.path.exists(sOpDir_toUse):
                        os.makedirs(sOpDir_toUse)

                    #fig.savefig('%s/%s_%s_%s_%s.png' % (sOpDir_toUse,histo_name_toUse.replace('_%s'%selectionTag, ''),systematic,sData, yAxisScale), transparent=False, dpi=80, bbox_inches="tight")
                    fig.savefig('%s/%s_%s_%s.png' % (sOpDir_toUse,histo_name_toUse.replace('_%s'%selectionTag, ''),systematic, yAxisScale), transparent=False, dpi=80, bbox_inches="tight")
    

                    plt.close(fig)

                    


