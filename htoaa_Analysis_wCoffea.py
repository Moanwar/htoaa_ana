#htoaa analysis main code

import os
import sys
import json
from collections import OrderedDict as OD
import time
import tracemalloc
import math

import coffea.processor as processor
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




printLevel = 0
nEventToReadInBatch =  1.0*10**6 # 2500000 #  1000 # 2500000
nEventsToAnalyze =  -1 # 1000 # 100000 # -1
#pd.set_option('display.max_columns', None)

#print("".format())


class Processor(processor.ProcessorABC):
    def __init__(self):
        dataset_axis = hist.axis.StrCategory(name="dataset", label="", categories=[], growth=True)
        #muon_axis = hist.axis.Regular(name="massT", label="Transverse Mass [GeV]", bins=50, start=15, stop=250)
        
        self.output = processor.dict_accumulator({
            #'massT': hist.Hist(dataset_axis, muon_axis),
            #'hMuonPt': hist.Hist(dataset_axis,
            #                     hist.axis.Regular(name="hMuonPt", label="Muon pT [GeV]", bins=50, start=0, stop=100)
            #                     ),
            'cutflow': processor.defaultdict_accumulator(int),
            'hnMuon': hist.Hist.new.Reg(20, -0.5, 19.5, name="No. of muons").Weight(),
            'hMuonPt': hist.Hist.new.Reg(100, 0, 100, name="Muon pT").Weight(),
            'hMuonPt_lead': hist.Hist.new.Reg(100, 0, 100, name="Leading muon pT").Weight(),
            'hMuonPt_sublead': hist.Hist.new.Reg(100, 0, 100, name="Sub-leading muon pT").Weight(),
            'hMuonEta': hist.Hist.new.Reg(40, -3, 3, name="Muon eta").Weight(),
            'hMuonPhi': hist.Hist.new.Reg(40, -3.14, 3.14, name="Muon phi").Weight(),
            ##
            'hnElectron': hist.Hist.new.Reg(20, -0.5, 19.5, name="No. of electron").Weight(),
            'hElectronPt': hist.Hist.new.Reg(100, 0, 100, name="Electron pT").Weight(),
            'hElectronPt_lead': hist.Hist.new.Reg(100, 0, 100, name="Leading Electron pT").Weight(),
            'hElectronPt_sublead':  hist.Hist.new.Reg(100, 0, 100, name="Sub-leading Electron pT").Weight(),
            'hElectronEta':         hist.Hist.new.Reg(40, -3, 3, name="Electron eta").Weight(),
            'hElectronPhi': hist.Hist.new.Reg(40, -3.14, 3.14, name="Electron phi").Weight(),
            #
            'hmass_2Mu': hist.Hist.new.Reg(100, 0, 200, name="mass(2 mu) [GeV]").Weight(),
            'hmass_2Ele': hist.Hist.new.Reg(100, 0, 200,  name="mass(2 ele) [GeV]").Weight(),
            'hdR_2Mu': hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 mu)").Weight(),
            'hdR_2Ele': hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 ele)").Weight(),
            'hdR_2Mu_2Ele': hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 mu, 2 ele)").Weight(),
            'h2ddR_2Mu_vs_dR_2Ele': (
                hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 mu)")
                .Reg(100, 0, 3.14, name="deltaR(2 ele)")
                .Weight() ),            
        })



    def process(self, events):
        nMuonsToSelect           = 2
        nElectronsToSelect       = 2
        muonPtThrshsToSelect     = [4.0, 3.0]
        electronPtThrshsToSelect = [4.0, 3.0]
                    
        
        #dataset = events.metadata["dataset"]
    
        # Keep track of muons and electrons by tagging them 0/1.
        #electrons = ak.with_field(events.Electron, 11, 'flavor')
        #muons     = ak.with_field(events.Muon, 13, 'flavor')


        self.output['cutflow']['All events'] += len(events)
        
        if printLevel >= 5:
            print(f"events.fields: {events.fields}")
            print(f"events.Muon.fields: {events.Muon.fields}")

        if printLevel >= 13:
            #print(f"events ({len(events)}) ({type(events)}): {events}")
            print(f"nEvents ({len(events)}) ")
            #print(f"events.Muon ({type(events.Muon)}) {events.Muon}")
        if printLevel >= 10:  print(f"nEvents 0 ({len(events)}) ")
        #if printLevel >= 10:  print(f"events ({len(events)}) {events.to_list()}")
        

        events = events[ (
            (ak.num(events.Muon)     >= nMuonsToSelect) &
            (ak.num(events.Electron) >= nElectronsToSelect) )]

        self.output['cutflow']['Muon, electron multiplicity cut'] += len(events)

        if printLevel >= 12:  print(f"nEvents nMu, nEle cut ({len(events)}) ")
        #if printLevel >= 10:  print(f"events ({len(events)}) {events.to_list()}")
        if printLevel >= 13:
            #print(f"events.nMuon: {ak.num(events.Muon).to_list()}")
            #print(f"events.nElectron: {ak.num(events.Electron).to_list()}")
            print(f"[nMuon, nElectron]: {list(zip(ak.num(events.Muon), ak.num(events.Electron)))}")
            

        if printLevel >= 12:
            #print(f"events.Muon.pt ({type(events.Muon.pt)}): {events.Muon.pt.to_list()}")
            #print(f"events.Electron.pt ({type(events.Electron.pt)}): {events.Electron.pt.to_list()}")
            #print(f"events.Muon ({type(events.Muon)}): {events.Muon.to_list()}")
            #print(f"events.Electron ({type(events.Electron)}): {events.Electron.to_list()}")
            print(f"[nMuon, nElectron]: {list(zip(ak.num(events.Muon), ak.num(events.Electron)))}")


        if printLevel >= 10:
            print(f"events.Muon.pt before cut({type(events.Muon.pt)}): {events.Muon.pt.to_list()}")


        if printLevel >= 10:
            print(f"events.Muon.pt[:, 0] {events.Muon.pt[:, 0].to_list()}")
            #print(f"(events.Muon.pt[0] > muonPtThrshsToSelect[0]): {(events.Muon.pt[0] > muonPtThrshsToSelect[0])}")


    
        events = events[ (
            (events.Muon.pt[:, 0] > muonPtThrshsToSelect[0]) &
            (events.Muon.pt[:, 1] > muonPtThrshsToSelect[1]) &
            (abs(events.Muon.eta[:, 0]) < 2.4) &
            (abs(events.Muon.eta[:, 1]) < 2.4) ) ]

        self.output['cutflow']['Muon pT, eta cut'] += len(events)

        if printLevel >= 10:  print(f"nEvents mu pT eta cut ({len(events)}) ")
        #if printLevel >= 10:  print(f"events mu pT eta cut: ({len(events)}) {events.to_list()}")

        if printLevel >= 12:
            print(f"events.Muon.pt after cut ({type(events.Muon.pt)}): {events.Muon.pt.to_list()}")

        
        if printLevel >= 12:
            print(f"events.Electron.pt ({type(events.Electron.pt)}): {events.Electron.pt.to_list()}")

            


    
        events = events[ (
            (events.Electron.pt[:, 0] > electronPtThrshsToSelect[0]) &
            (events.Electron.pt[:, 1] > electronPtThrshsToSelect[1]) &
            (abs(events.Electron.eta[:, 0]) < 2.5) &
            (abs(events.Electron.eta[:, 1]) < 2.5) ) ]

        self.output['cutflow']['Electron pT, eta cut'] += len(events)

        if printLevel >= 10:  print(f"nEvents ele pT eta cut ({len(events)}) ")
        if printLevel >= 12:
            print(f"events.Muon.charge ({type(events.Muon.charge)}): {events.Muon.charge}")
            print(f"events.Muon.charge[:, :nMuonsToSelect] ({type(events.Muon.charge)}): {events.Muon.charge[:, :nMuonsToSelect]}")
            print(f"ak.sum(events.Muon.charge, axis=-1): {ak.sum(events.Muon.charge[:, :nMuonsToSelect], axis=-1)}")
            

        events = events[ (
            (ak.sum(events.Muon.charge[:, :nMuonsToSelect], axis=-1) == 0) & 
            (ak.sum(events.Electron.charge[:, :nElectronsToSelect], axis=-1) == 0) ) ]

        self.output['cutflow']['Muon, electron charge sum cut'] += len(events)

        if printLevel >= 10:  print(f"nEvents mu, ele chargeSum cut ({len(events)}) ")
        if printLevel >= 12:
            print(f"events.Muon.charge ({type(events.Muon.charge)}): {events.Muon.charge}")
            print(f"events.Muon.charge[:, :nMuonsToSelect] ({type(events.Muon.charge)}): {events.Muon.charge[:, :nMuonsToSelect]}")
            print(f"ak.sum(events.Muon.charge, axis=-1): {ak.sum(events.Muon.charge[:, :nMuonsToSelect], axis=-1)}")

        if printLevel >= 13:
            print(f"events.Muon.pt:  ({type(events.Muon.pt)}): {events.Muon.pt.to_list()}")
            print(f"events.Muon.energy: {events.Muon.energy}")
            print(f"events.Muon[:, 0].energy: {events.Muon[:, 0].energy}")


        if printLevel >= 10:
            print(f"(events.Muon[:, 0] + events.Muon[:, 1]).mass: {(events.Muon[:, 0] + events.Muon[:, 1]).mass}")
            
        
        self.output['hnMuon'].fill( ak.to_numpy(ak.num(events.Muon)) )
        self.output['hMuonPt'].fill( ak.flatten(events.Muon.pt) )
        self.output['hMuonPt_lead'].fill( ak.to_numpy(events.Muon.pt[:, 0]) )
        self.output['hMuonPt_sublead'].fill( ak.to_numpy(events.Muon.pt[:, 1]) )
        self.output['hMuonEta'].fill( ak.flatten(events.Muon.eta) )
        self.output['hMuonPhi'].fill( ak.flatten(events.Muon.phi) )
        #
        self.output['hnElectron'].fill( ak.to_numpy(ak.num(events.Electron)) )
        self.output['hElectronPt'].fill( ak.flatten(events.Electron.pt) )
        self.output['hElectronPt_lead'].fill( ak.to_numpy(events.Electron.pt[:, 0]) )
        self.output['hElectronPt_sublead'].fill( ak.to_numpy(events.Electron.pt[:, 1]) )
        self.output['hElectronEta'].fill( ak.flatten(events.Electron.eta) )
        self.output['hElectronPhi'].fill( ak.flatten(events.Electron.phi) )
        #
        self.output['hmass_2Mu'].fill( (events.Muon[:, 0] + events.Muon[:, 1]).mass )
        self.output['hmass_2Ele'].fill( (events.Electron[:, 0] + events.Electron[:, 1]).mass )
        self.output['hdR_2Mu'].fill( events.Muon[:, 0].delta_r(events.Muon[:, 1])  )
        self.output['hdR_2Ele'].fill( events.Electron[:, 0].delta_r(events.Electron[:, 1])  )
        self.output['hdR_2Mu_2Ele'].fill( (events.Muon[:, 0] + events.Muon[:, 1]).delta_r(events.Electron[:, 0] + events.Electron[:, 1]) )
        self.output['h2ddR_2Mu_vs_dR_2Ele'].fill( events.Muon[:, 0].delta_r(events.Muon[:, 1]), events.Electron[:, 0].delta_r(events.Electron[:, 1]) )

        
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
    print("htoaa_Analysis:: main: {}".format(sys.argv))


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
        sDir1 = 'evt/'
        
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
    
