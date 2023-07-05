#htoaa analysis main code

import os
import sys
import subprocess
import json
import glob
from collections import OrderedDict as OD
import time
import tracemalloc
import math
import numpy as np
from copy import copy, deepcopy
#import uproot
import uproot3 as uproot

'''
Count sum of positive and negative gen-weight events in sample

References:
  * Coffea framework used for TTGamma analysis: https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py
  * Coffea installation: /home/siddhesh/anaconda3/envs/ana_htoaa/lib/python3.10/site-packages/coffea
'''

#import coffea.processor as processor
from coffea import processor, util
from coffea.nanoevents import schemas
from coffea.nanoevents.methods import nanoaod, vector
from coffea.analysis_tools import PackedSelection, Weights
#import hist
from coffea import hist
import awkward as ak
import uproot
#from dask.distributed import Client
from particle import Particle # For PDG particle listing https://github.com/scikit-hep/particle



from htoaa_Settings import *
from htoaa_CommonTools import (
    GetDictFromJsonFile, calculate_lumiScale, setXRootDRedirector, getNanoAODFile,
    xrdcpFile
)
from htoaa_Samples import (
    kData, kQCD_bEnrich, kQCD_bGen, kQCDIncl
)




# use GOldenJSON

 
printLevel = 0
nEventToReadInBatch =  0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze = -1 # 1000 # 100000 # -1


sWeighted = "Wtd: "


    

# -----------------------------------------------------------------------------------
def printVariable(sName, var):
    printInDetail=True
    #if nEventsToAnalyze == -1: printInDetail = False
    if str(type(var)) in ['numpy.ndarray', "<class 'numpy.ndarray'>"]: printInDetail = False # as gave error
    #print(f"printInDetail: {printInDetail} {sName} ({type(var)}) ({len(var)}): {var}")
    if not printInDetail:
        print(f"{sName} ({type(var)}) ({len(var)}): {var}")
    else:
        print(f"{sName} ({type(var)}) ({len(var)}): {var.to_list()}")

# -----------------------------------------------------------------------------------

    
    
class HToAATo4bProcessor(processor.ProcessorABC):
    def __init__(self, datasetInfo={}):

        ak.behavior.update(nanoaod.behavior)

        self.datasetInfo = datasetInfo
        
        
        #dataset_axis = hist.axis.StrCategory(name="dataset", label="", categories=[], growth=True)
        #muon_axis = hist.axis.Regular(name="massT", label="Transverse Mass [GeV]", bins=50, start=15, stop=250)
        dataset_axis    = hist.Cat("dataset", "Dataset")
        systematic_axis = hist.Cat("systematic", "Systematic Uncertatinty")

        cutFlow_axis  = hist.Bin("CutFlow",   r"Cuts",            21, -0.5, 20.5)

        sXaxis      = 'xAxis'
        sXaxisLabel = 'xAxisLabel'
        sYaxis      = 'yAxis'
        sYaxisLabel = 'yAxisLabel'
        histos = OD([
            # ('histogram_name',  {sXaxis: hist.Bin() axis,  sXaxisLabel: "histogram axis label"})
            ('hCutFlow',                                  {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),
            ('hCutFlowWeighted',                          {sXaxis: cutFlow_axis,    sXaxisLabel: 'Cuts'}),
            
            
        ])        

        self._accumulator = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int)
        })

        for histName, histAttributes in histos.items():
            #hXaxis = histAttributes[sXaxis].copy()
            hXaxis = deepcopy(histAttributes[sXaxis])
            hXaxis.label = histAttributes[sXaxisLabel]

            if sYaxis not in histAttributes.keys():
                # TH1
                self._accumulator.add({
                    histName: hist.Hist(
                        "Counts",
                        dataset_axis,
                        hXaxis, #nObject_axis,
                        systematic_axis,
                    )
                })
            else:
                # TH2
                hYaxis = deepcopy(histAttributes[sYaxis])
                hYaxis.label = histAttributes[sYaxisLabel]
                
                self._accumulator.add({
                    histName: hist.Hist(
                        "Counts",
                        dataset_axis,
                        hXaxis, #nObject_axis,
                        hYaxis,
                        systematic_axis,
                    )
                })
                



    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata["dataset"] # dataset label
        self.datasetInfo[dataset]['isSignal'] = False
        self.datasetInfo[dataset]['isQCD'] = False

        if printLevel >= 20:
            print(f"nEvents: {len(events)}")
        if printLevel >= 20:
            print(f"\n events.fields: {events.fields}")
            print(f"\n events.GenPart.fields: {events.GenPart.fields}")
            
        if nEventsToAnalyze != -1:
            print(f"\n (run:ls:event): {ak.zip([events.run, events.luminosityBlock, events.event])}")
              
            
            output = self.accumulator.identity()
            systematics_shift = [None]
            for _syst in systematics_shift:
                output += self.process_shift(events, _syst)
        else:
            output = self.process_shift(events, None)

        return output


    
    def process_shift(self, events, shift_syst=None):
        
        output = self.accumulator.identity()
        dataset = events.metadata["dataset"] # dataset label
        #print(f"dataset: {dataset}")

        

        ones_list  = np.ones(len(events))
        trues_list = np.ones(len(events), dtype=bool)


          
        ##################
        # OBJECT SELECTION
        ##################

        
        ##################
        # EVENT VARIABLES
        ##################
        
            
        #####################
        # EVENT SELECTION
        #####################
        selection = PackedSelection()

        sel_posGenWgt = None
        sel_negGenWgt = None        
        
        if self.datasetInfo[dataset]["isMC"]:
            selection.add("posGenWgt", (events.genWeight > 0))
        
            selection.add("negGenWgt", (events.genWeight < 0))

            sel_posGenWgt           = selection.all("posGenWgt")
            sel_negGenWgt           = selection.all("negGenWgt")
        
        ################
        # EVENT WEIGHTS
        ################

        # create a processor Weights object, with the same length as the number of events in the chunk
        weights     = Weights(len(events))
        weights_gen = Weights(len(events))



        if self.datasetInfo[dataset]["isMC"]:
            weights_gen.add(
                "genWeight",
                weight=np.copysign(np.ones(len(events)), events.genWeight)
            )
            
            
        
        ###################
        # FILL HISTOGRAMS
        ###################

        systList = []
        if self.datasetInfo[dataset]['isMC']:
            if shift_syst is None:
                systList = [
                    "central"
                ]
            else:
                systList = [shift_syst]
        else:
            systList = ["noweight"]

            
        output['cutflow']['all events'] += len(events)
        output['cutflow'][sWeighted+'all events'] += weights.weight().sum() 
        for n in selection.names:
            sel_i = selection.all(n)
            output['cutflow'][n] += sel_i.sum()
            output['cutflow'][sWeighted+n] += weights.weight()[sel_i].sum()
           

        for syst in systList:

            # find the event weight to be used when filling the histograms
            weightSyst = syst
            
            # in the case of 'central', or the jet energy systematics, no weight systematic variation is used (weightSyst=None)
            if syst in ["central", "JERUp", "JERDown", "JESUp", "JESDown"]:
                weightSyst = None

            
            
            if syst == "noweight":
                evtWeight = np.ones(len(events))
            else:
                evtWeight     = weights.weight(weightSyst)
                evtWeight_gen = weights_gen.weight(weightSyst)

               

            # all events
            iBin = 0
            output['hCutFlow'].fill(
                dataset=dataset,
                CutFlow=(ones_list * iBin),
                systematic=syst
            )
            output['hCutFlowWeighted'].fill(
                dataset=dataset,
                CutFlow=(ones_list * iBin),
                systematic=syst,
                weight=evtWeight
            )

            # MC ----------------------------------------------
            if self.datasetInfo[dataset]['isMC'] :
                # proper GenWeights
                iBin = 1
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list * iBin),
                    systematic=syst
                )
                output['hCutFlowWeighted'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list * iBin),
                    systematic=syst,
                    weight=evtWeight_gen
                )

                # events w/ positive GenWeights
                iBin = 2
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[sel_posGenWgt] * iBin),
                    systematic=syst
                )
                output['hCutFlowWeighted'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[sel_posGenWgt] * iBin),
                    systematic=syst,
                    weight=evtWeight_gen[sel_posGenWgt]
                )

                # events w/ negative GenWeights
                iBin = 3
                output['hCutFlow'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[sel_negGenWgt] * iBin),
                    systematic=syst
                )
                output['hCutFlowWeighted'].fill(
                    dataset=dataset,
                    CutFlow=(ones_list[sel_negGenWgt] * iBin),
                    systematic=syst,
                    weight=evtWeight_gen[sel_negGenWgt]
                )
       
        return output


    def postprocess(self, accumulator):
        #pass
        return accumulator



    
                       

def printWithType(sX, X):
    #print(f"{sX} ({type(X)}): {X}")
    print(f"{sX} : {X}")


        
    
    
if __name__ == '__main__':
    print("htoaa_Analysis:: main: {}".format(sys.argv)); sys.stdout.flush()


    if len(sys.argv) != 2:
        print("htoaa_Analysis:: Command-line config file missing.. \t **** ERROR **** \n")

    sConfig = sys.argv[1]

    
    config = GetDictFromJsonFile(sConfig)
    print("Config {}: \n{}".format(sConfig, json.dumps(config, indent=4)))

    lumiScale = 1
    sInputFiles         = config["inputFiles"]
    sOutputFile         = config["outputFile"]
    sample_dataset      = config["dataset"] 
    sample_category     = config['sampleCategory']
    isMC                = config["isMC"]
    era                 = config['era']
    downloadIpFiles     = config['downloadIpFiles'] if 'downloadIpFiles' in config else False
    '''
    # lumiScale is not applied for countSumEventsInSample.py
    if isMC:
        luminosity          = Luminosities[era][0]
        sample_crossSection = config["crossSection"]
        sample_nEvents      = config["nEvents"]
        sample_sumEvents    = config["sumEvents"] if config["sumEvents"] != -1 else sample_nEvents
        if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
        lumiScale = calculate_lumiScale(luminosity=luminosity, crossSection=sample_crossSection, sumEvents=sample_sumEvents)    
    '''    
    #branchesToRead = htoaa_nanoAODBranchesToRead
    #print("branchesToRead: {}".format(branchesToRead))
    sample_dataset = sample_dataset[0] if isinstance(sample_dataset, list) else sample_dataset # dataset is list of datasets w/ same sample name, as they are considered together recently. Those set of datasets are extension of the same samples.
    
    print(f"isMC: {isMC}, lumiScale: {lumiScale}")
    sInputFiles_toUse = []
    for sInputFile in sInputFiles:
        if "*" in sInputFile:  sInputFiles_toUse.extend( glob.glob( sInputFile ) )
        else:                  sInputFiles_toUse.append( sInputFile )
    sInputFiles = sInputFiles_toUse
    print(f"Initial sInputFiles ({len(sInputFiles)}) (type {type(sInputFiles)}):");
    for sInputFile in sInputFiles:
        print(f"\t{sInputFile}");  sys.stdout.flush()

    for iFile in range(len(sInputFiles)):     
        sInputFile = sInputFiles[iFile]
        sFileLocal = './inputFiles/%s' %(os.path.basename(sInputFile))   
        sInputFiles[iFile] = getNanoAODFile(
            fileName = sInputFile, 
            useLocalFileIfExists = True, 
            downloadFile = True, 
            fileNameLocal = './inputFiles/%s' %(os.path.basename(sInputFile)), 
            nTriesToDownload = 3)
        


    startTime = time.time()
    
    tracemalloc.start()
    



    #client = Client("tls://localhost:8786")
    #executor = processor.DaskExecutor(client=client)
    chunksize = nEventToReadInBatch
    maxchunks = None if nEventsToAnalyze == -1 else int(nEventsToAnalyze/nEventToReadInBatch)
    nWorkers  = 4 if nEventsToAnalyze == -1 else 1
    print(f"nEventsToAnalyze: {nEventsToAnalyze},  nEventToReadInBatch: {nEventToReadInBatch}, chunksize: {chunksize},  maxchunks: {maxchunks},  nWorkers: {nWorkers}")
    run = processor.Runner(
        #executor=executor,
        executor=processor.FuturesExecutor(workers=nWorkers),
        schema=schemas.NanoAODSchema,
        savemetrics=True,
        chunksize=chunksize,  #3 ** 20,  ## Governs the number of times LeptonJetProcessor "process" is called
        maxchunks=maxchunks
    )

    output, metrics = run(
        fileset={sample_category: sInputFiles},
        #fileset={"QCD": ["/home/siddhesh/Work/CMS/htoaa/analysis/tmp/20BE2B12-EFF6-8645-AB7F-AFF6A624F816.root"]},
        treename="Events",
        processor_instance=HToAATo4bProcessor(
            datasetInfo={
                "era": era, 
                sample_category: {"isMC": isMC, "lumiScale": lumiScale}
            }
        )
    )

    print(f"metrics: {metrics}")


    if 'cutflow' in output.keys():
        print("Cutflow::")
        #for key, value in output['cutflow'].items():
        for key in output['cutflow'].keys():
            #print(key, value)
            if key.startswith(sWeighted): continue

            print("%10f\t%10d\t%s" % (output['cutflow'][sWeighted+key], output['cutflow'][key], key))


    
    if sOutputFile is not None:
        if not sOutputFile.endswith('.root'): sOutputFile += '.root'
        #sDir1 = 'evt/%s' % (sample_category)
        #sample_dataset_toUse = sample_dataset[1:] if sample_dataset.startswith('/') else sample_dataset
        #sDir1 = 'evt/%s' % (sample_dataset_toUse) #
        sample_dataset_toUse = sample_dataset[0] if isinstance(sample_dataset, list) else sample_dataset
        sample_name_toUse = sample_dataset_toUse.split('/')[1]
        sDir1 = 'evt/%s' % (sample_name_toUse)
        print(f"Output file: {sOutputFile}, directory: {sDir1} ")
        
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                #print(f"key: {key},  value ({type(value)}): {value}")

                if not isinstance(value, hist.Hist): continue
                #print(f"1: key {key}, value ({type(value)})     Hist: {type(hist.Hist)},    isinstance(value, hist.Hist): {isinstance(value, hist.Hist)}") # value: {value}")



                for _dataset in value.axis('dataset').identifiers():
                    #print(f"_dataset ({type(_dataset)}): {_dataset}")

                    for _syst in value.axis('systematic').identifiers():
                        #print(f"_syst ({type(_syst)}): {_syst}")

                        h1 = value.integrate('dataset',_dataset).integrate('systematic',_syst).to_hist()
                        #print(f"h1 ({type(h1)}): {h1}")
                        #print("histoNameToSave: %s/%s_%s" % (sDir1, key, _syst))

                        #fOut['%s/%s_%s_%s' % (sDir1, key, _dataset, _syst)] = h1
                        fOut['%s/%s_%s' % (sDir1, key, _syst)] = h1
                
                #fOut['%s%s' % (sDir1, key)] = value
                #fOut['%s%s' % (sDir1, key)] = hist.export1d(value)
                #fOut['%s%s' % (sDir1, key)] = h1

                
        

        #util.save(output, sOutputFile)
            
        
        print("Wrote to sOutputFile {}".format(sOutputFile))
   















    
    current_memory, peak_memory = tracemalloc.get_traced_memory() # https://medium.com/survata-engineering-blog/monitoring-memory-usage-of-a-running-python-program-49f027e3d1ba
    print(f"\n\nMemory usage:: current {current_memory / 10**6}MB;  peak {peak_memory / 10**6}MB")

    endTime = time.time()
    totalTime = endTime - startTime
    totalTime_hr  = int(totalTime/60/60)
    totalTime_min = totalTime - float(totalTime_hr * 60)
    totalTime_min = int(totalTime_min/60)
    totalTime_sec = totalTime - float(totalTime_hr * 60*60) - float(totalTime_min * 60)
    print(f"Total run time: {totalTime_hr}h {totalTime_min}m {totalTime_sec}s = {totalTime}sec ")
    
