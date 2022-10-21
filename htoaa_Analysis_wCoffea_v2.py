#htoaa analysis main code

import os
import sys
import json
from collections import OrderedDict as OD
import time
import tracemalloc
import math

#import coffea.processor as processor
from coffea import processor, util
from coffea.nanoevents import schemas
import awkward as ak
import hist
import uproot
#from dask.distributed import Client


from htoaa_Settings import *
from htoaa_CommonTools import GetDictFromJsonFile, calculate_lumiScale
'''
from htoaa_Settings import *
from htoaa_NanoAODBranches import htoaa_nanoAODBranchesToRead
from htoaa_CommonTools import GetDictFromJsonFile, DfColLabel_convert_bytes_to_string, calculate_lumiScale
from htoaa_CommonTools import cut_ObjectMultiplicity, cut_ObjectPt, cut_ObjectEta, cut_ObjectPt_1
'''




printLevel = 5
nEventToReadInBatch = 200 # 0.5*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze =  200 # -1 # 1000 # 100000 # -1
#pd.set_option('display.max_columns', None)

#print("".format())


class Processor(processor.ProcessorABC):
    def __init__(self):
        dataset_axis = hist.axis.StrCategory(name="dataset", label="", categories=[], growth=True)
        #muon_axis = hist.axis.Regular(name="massT", label="Transverse Mass [GeV]", bins=50, start=15, stop=250)
        
        self.output = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int),
            'hnFatJet_level0': hist.Hist.new.Reg(20, -0.5, 19.5, name="No. of FatJet").Weight(),
            'hnFatJet_level1': hist.Hist.new.Reg(20, -0.5, 19.5, name="No. of FatJet").Weight(),
        })



    def process(self, events):
        nFatJet = 1
        FatJetPtThsh = 220
        FatJetEtaThsh = 2.4
        
        #dataset = events.metadata["dataset"]
    
        # Keep track of muons and electrons by tagging them 0/1.
        #electrons = ak.with_field(events.Electron, 11, 'flavor')
        #muons     = ak.with_field(events.Muon, 13, 'flavor')
        
        if printLevel >= 5:
            print(f"events.fields: {events.fields}")
            print(f"events.FatJet.fields: {events.FatJet.fields}")
            print(f"events.FatJet.pt: {events.FatJet.pt}")


        self.output['cutflow']['All events'] += len(events)

        '''
        jets.cut(jets['FatJet_pt'] > 170)
        jets.cut(jets['FatJet_eta'].abs() < 2.4)
        jets.cut(jets['FatJet_btagDDBvL'] > 0.8)
        jets.cut(jets['FatJet_btagDeepB'] > 0.4184)
        jets.cut(jets['FatJet_msoftdrop'] > 90)
        jets.cut(jets['FatJet_msoftdrop'] < 200)#<= 200)
        jets.cut(jets['FatJet_mass'] > 90)
        #jets.cut(jets['FatJet_mass'] <= 200)
        other.cut(other['PV_npvsGood'] >= 1)
        '''

        self.output['hnFatJet_level0'].fill( ak.to_numpy(ak.num(events.FatJet)) )
        
        #events = events[(
        #    (ak.num(events.FatJet) >= 1) )]
        #self.output['cutflow']['nFatJet >= 1'] += len(events)

        if printLevel >= 5:
            print(f"events.FatJet.pt _0 ({len(events)}): {events.FatJet.pt.to_list()}")

        #events = events[(
        #    (events.FatJet.pt[:, 0] > FatJetPtThsh) )]
        #self.output['cutflow']['LeadingFatJetPt > %s' % (FatJetPtThsh)] += len(events)

        if printLevel >= 5:
            print(f"\n\nevents.FatJet.pt _1 ({len(events)}): {events.FatJet.pt.to_list()}")
            cut1  = (events.FatJet.pt > FatJetPtThsh)
            print(f"\n\ncut1 ({len(cut1)}): {cut1.to_list()}")
            #events1 = events.FatJet[cut1]
            events1 = events[ak.any(cut1, axis=1)]
            print(f"\n\nevents1 ({len(events1)}): {events1.to_list()}")
            print(f"\n\nevents1.FatJet.pt _1 ({len(events1)}): {events1.FatJet.pt.to_list()}")

        events = events[(
            (events.FatJet.pt > FatJetPtThsh) )]
        self.output['cutflow']['LeadingFatJetPt _1 > %s' % (FatJetPtThsh)] += len(events)
        
        if printLevel >= 5:
            print(f"events.FatJet.pt _2 ({len(events)}): {events.FatJet.pt.to_list()}")

        
        events = events[(
            (abs(events.FatJet.eta) < 2.4) )]
        self.output['cutflow']['FatJetAbsEta < 2.4'] += len(events)

        if printLevel >= 5:
            print(f"nEvents FatJetEta cut ({len(events)}) ")
            print(f"events: {events.to_list()}")
            
        self.output['hnFatJet_level1'].fill( ak.to_numpy(ak.num(events.FatJet)) )

        
        return self.output


    def postprocess(self, accumulator):
        pass



    
def getLorentzVectorFromAwkArray(awkArray, ptObjectName, etaObjectName, phiObjectName, massObjectName):
    '''
    return ak.Array(
        awkArray,
        with_name='Momentum4D'
    )
    '''
    v1 = ak.zip( {
        'pt':   awkArray[ptObjectName],
        'eta':  awkArray[etaObjectName],
        'phi':  awkArray[phiObjectName],
        'mass': awkArray[massObjectName],
    })
    if printLevel >= 12:
        print(f"getLorentzVectorFromAwkArray(): v1 ({type(v1)}): {v1.to_list()}")

    return ak.Array(v1, with_name='Momentum4D')
                       




if __name__ == '__main__':
    print("htoaa_Analysis:: main: {}".format(sys.argv)); sys.stdout.flush()


    if len(sys.argv) != 2:
        print("htoaa_Analysis:: Command-line config file missing.. \t **** ERROR **** \n")

    sConfig = sys.argv[1]

    
    config = GetDictFromJsonFile(sConfig)
    print("Config {}: \n{}".format(sConfig, json.dumps(config, indent=4)))
    
    sInputFiles  = config["inputFiles"]
    sOutputFile  = config["outputFile"]
    sample_category = config['sampleCategory']
    era = config['era']
    luminosity = Luminosities[era][0]
    sample_crossSection = config["crossSection"]
    sample_nEvents = config["nEvents"]
    sample_sumEvents = config["sumEvents"] if config["sumEvents"] != -1 else sample_nEvents
    if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
    lumiScale = calculate_lumiScale(luminosity=luminosity, crossSection=sample_crossSection, sumEvents=sample_sumEvents)
    #branchesToRead = htoaa_nanoAODBranchesToRead
    #print("branchesToRead: {}".format(branchesToRead))

    for iFile in range(len(sInputFiles)):
        if sInputFiles[iFile].startswith("/store/"): # LFN: Logical File Name
            sInputFiles[iFile] = xrootd_redirectorName + sInputFiles[iFile]
    print(f"sInputFiles ({len(sInputFiles)}): {sInputFiles}"); sys.stdout.flush()

    startTime = time.time()
    
    tracemalloc.start()
    



    #client = Client("tls://localhost:8786")
    #executor = processor.DaskExecutor(client=client)
    chunksize = nEventToReadInBatch
    maxchunks = None if nEventsToAnalyze == -1 else int(nEventsToAnalyze/nEventToReadInBatch)
    print(f"nEventsToAnalyze: {nEventsToAnalyze},  nEventToReadInBatch: {nEventToReadInBatch}, chunksize: {chunksize},  maxchunks: {maxchunks}")
    run = processor.Runner(
        #executor=executor,
        executor=processor.FuturesExecutor(workers=4),
        schema=schemas.NanoAODSchema,
        savemetrics=True,
        chunksize=chunksize,  #3 ** 20,  ## Governs the number of times LeptonJetProcessor "process" is called
        maxchunks=maxchunks
    )

    output, metrics = run(
        fileset={sample_category: sInputFiles},
        #fileset={"QCD": ["/home/siddhesh/Work/CMS/htoaa/analysis/tmp/20BE2B12-EFF6-8645-AB7F-AFF6A624F816.root"]},
        treename="Events",
        processor_instance=Processor()
    )

    print(f"metrics: {metrics}")


    if 'cutflow' in output.keys():
        print("Cutflow::")
        for key, value in output['cutflow'].items():
            print(key, value)    



    if sOutputFile is not None:
        if not sOutputFile.endswith('.root'): sOutputFile += '.root'
        sOutputFile = sOutputFile.replace('.root', '_wCoffea.root')
        sDir1 = 'evt/%s' % (sample_category)
        
        with uproot.recreate(sOutputFile) as fOut:
            for key, value in output.items():
                #if not (key.startswith('h') or key != 'cutflow'): continue
                if not isinstance(value, hist.Hist): continue
                #print(f"1: key {key}, value ({type(value)})     Hist: {type(hist.Hist)},    isinstance(value, hist.Hist): {isinstance(value, hist.Hist)}") # value: {value}")

                
                fOut['%s%s' % (sDir1, key)] = value

            
        
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
    
