
'''
Calculate lumiScale weights for samples with Overlap Phase Space.  
'''

import os
import sys
from pathlib import Path
import json
import glob
import argparse
from collections import OrderedDict as OD
import numpy as np
import time
from datetime import datetime
import copy

import ROOT as R

from htoaa_Settings import *
from htoaa_Samples import (
    Samples2018,
    kQCDIncl, kQCD_bGen, kQCD_bEnrich
)
from htoaa_CommonTools import (
    getSampleHTRange,
    calculate_lumiScale
)


printLevel = 10


def setHistogramSampleComponenetsVsHT(sHistoName, sample_HTBinEdges, sample_components):
    h_ = R.TH2D(
        sHistoName, "",
        len(sample_HTBinEdges)-1, np.array(sample_HTBinEdges),
        len(sample_components), -0.5, len(sample_components)-0.5
    )
    h_.GetXaxis().SetTitle('HT [GeV]')
    h_.GetYaxis().SetTitle('Samples')
    for yBin_, sample_component in enumerate(sample_components):
        print(f"setHistogramSampleComponenetsVsHT():: {sHistoName = }, {yBin_+1 = }, {sample_component = }, ")
        h_.GetYaxis().SetBinLabel(yBin_+1, sample_component)
    return h_
    



if __name__ == "__main__":
     
    print("calculateLumiscale_wPhSpOverlapRewght.py:: main: {}".format(sys.argv)); sys.stdout.flush()

    parser = argparse.ArgumentParser(description='calculateLumiscale_wPhSpOverlapRewght')
    parser.add_argument('-era', dest='era',   type=str, default=Era_2018,                    choices=[Era_2016, Era_2017, Era_2018], required=False)
    args=parser.parse_args()
    print("args: {}".format(args))

    era              = args.era


    samples_wPhSpOverlap = OD([
        ('QCD', {
            'InclusiveSample':   kQCDIncl,
            'ExclusiveSamples': [kQCD_bGen, kQCD_bEnrich],
            'HTBinEdges': [50, 100, 200, 300, 500, 700, 1000, 1500, 2000, kLHE_HT_Max]
        }),
    ])

    sFOut_dict = {
        Era_2016: "",
        Era_2017: "",        
        Era_2018: "data/lumiScale/2018.root"
    }

    
    
    samplesList = None
    samplesInfo = None
    Luminosity  = None
    if era == Era_2018:
        samplesList = Samples2018 # htoaa_Samples.py
        
    with open(sFileSamplesInfo[era]) as fSamplesInfo:
        samplesInfo = json.load(fSamplesInfo) # Samples_Era.json
    Luminosity = Luminosities[era][0]



    if printLevel >= 100:
        print(f"samplesInfo ({type(samplesInfo)}) ({len(samplesInfo)}): {json.dumps(samplesInfo, indent=4)}")


    sFOut = sFOut_dict[era]
    os.makedirs(os.path.dirname(sFOut), exist_ok=True)
    fOut = R.TFile(sFOut, 'RECREATE')
    
    for sample_wPhSpOverlap in samples_wPhSpOverlap:
        sample_Inclusive  = samples_wPhSpOverlap[sample_wPhSpOverlap]['InclusiveSample']
        samples_Exclusive = samples_wPhSpOverlap[sample_wPhSpOverlap]['ExclusiveSamples']
        sample_HTBinEdges = samples_wPhSpOverlap[sample_wPhSpOverlap]['HTBinEdges']

        sample_components = [sample_Inclusive] + samples_Exclusive
        

        hSampleSumEvents    = setHistogramSampleComponenetsVsHT("%s_SumEvents"      % (sample_wPhSpOverlap), sample_HTBinEdges, sample_components)
        hSampleCrossSection = setHistogramSampleComponenetsVsHT("%s_CrossSection" % (sample_wPhSpOverlap), sample_HTBinEdges, sample_components)
        hSampleLumiScale    = setHistogramSampleComponenetsVsHT("%s_LumiScale_PhSpOverlapRewghted" % (sample_wPhSpOverlap), sample_HTBinEdges, sample_components)

        # Set SumEvents and Cross-section
        for sample_component in sample_components:
            samples_name_ = samplesList[sample_component]
            for sample_name_ in samples_name_:
                sample_HT_Min_, sample_HT_Max_ = getSampleHTRange(sample_name_)

                sumEvents_    = samplesInfo[sample_name_]["sumEvents"]
                crossSection_ = samplesInfo[sample_name_]["cross_section"]
                
                HTbin_     = hSampleSumEvents.GetXaxis().FindBin(sample_HT_Min_)
                sampleBin_ = hSampleSumEvents.GetYaxis().FindBin(sample_component)
                hSampleSumEvents   .SetBinContent(HTbin_, sampleBin_, sumEvents_)
                hSampleCrossSection.SetBinContent(HTbin_, sampleBin_, crossSection_)

        # Calculate lumiScale
        for iHTBin in range(len(sample_HTBinEdges) - 1):
            HTBin_min_                               = sample_HTBinEdges[iHTBin]
            HTbin_                                   = hSampleSumEvents.GetXaxis().FindBin(HTBin_min_)

            # Inclusive sample
            sampleInclusiveBin_                      = hSampleSumEvents.GetYaxis().FindBin(sample_Inclusive)
            sumEvents_SampleInclusive                = hSampleSumEvents   .GetBinContent(HTbin_, sampleInclusiveBin_)
            crossSection_SampleInclusive             = hSampleCrossSection.GetBinContent(HTbin_, sampleInclusiveBin_)
            crossSection_SampleInclusiveRemnant      = crossSection_SampleInclusive

            # Exclusive samples
            for sample_Exclusive in samples_Exclusive:
                sampleBin_                           = hSampleSumEvents.GetYaxis().FindBin(sample_Exclusive)
                sumEvents_SampleExclusive            = hSampleSumEvents   .GetBinContent(HTbin_, sampleBin_)
                crossSection_SampleExclusive         = hSampleCrossSection.GetBinContent(HTbin_, sampleBin_)
                crossSection_SampleInclusiveRemnant -= crossSection_SampleExclusive
                if sumEvents_SampleExclusive <= 0: continue

                PhSpFraction     = crossSection_SampleExclusive / crossSection_SampleInclusive
                sumEvents_total  = sumEvents_SampleExclusive + (PhSpFraction * sumEvents_SampleInclusive)
                lumiScale_       = calculate_lumiScale(Luminosity, crossSection_SampleExclusive, sumEvents_total)

                hSampleLumiScale.SetBinContent(HTbin_, sampleBin_, lumiScale_)

            # PhSp: Remnant of Inclusive sample
            PhSpFraction     = crossSection_SampleInclusiveRemnant / crossSection_SampleInclusive
            sumEvents_total  = (PhSpFraction * sumEvents_SampleInclusive)
            lumiScale_       = calculate_lumiScale(Luminosity, crossSection_SampleInclusiveRemnant, sumEvents_total)

            hSampleLumiScale.SetBinContent(HTbin_, sampleInclusiveBin_, lumiScale_)

        
        # Rename hSampleLumiScale bin label 'Inclusive' to 'Inclusive_Remnant'
        sampleInclusiveBin_                      = hSampleLumiScale.GetYaxis().FindBin(sample_Inclusive)
        hSampleLumiScale.GetYaxis().SetBinLabel(sampleInclusiveBin_, sample_Inclusive + '_Remnant')
            
            
            
        fOut.cd()
        fOut.mkdir(sample_wPhSpOverlap)
        fOut.cd(sample_wPhSpOverlap)
        hSampleSumEvents.Write()
        hSampleCrossSection.Write()
        hSampleLumiScale.Write()
       

        #for sample_PhSpCombination in sample_PhSpCombinations:
        #    for sample_PhSp in sample_PhSpCombination:
                

        
    
    fOut.Close()
    print(f"Wrote lumiScale in {sFOut}")
