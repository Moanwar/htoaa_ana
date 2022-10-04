#htoaa analysis main code

import os
import sys
import json
from collections import OrderedDict as OD
import time

import pandas as pd
import uproot
#import uproot3 as uproot
import numpy as np
import math
import awkward as ak
import matplotlib.pyplot as plt
#import mplhep as hep # for histograms like CMS template
import hist # for histogram
import vector

from collections import Counter
import linecache
import os
import tracemalloc

from htoaa_Settings import *
from htoaa_NanoAODBranches import htoaa_nanoAODBranchesToRead
from htoaa_CommonTools import GetDictFromJsonFile, DfColLabel_convert_bytes_to_string, calculate_lumiScale
from htoaa_CommonTools import cut_ObjectMultiplicity, cut_ObjectPt, cut_ObjectEta, cut_ObjectPt_1


printLevel = 0
nEventToReadInBatch =  2500000 # 10000 # 2500000
nEventsToAnalyze = -1 # 100000
pd.set_option('display.max_columns', None)

#print("".format())



def display_top(snapshot, key_type='lineno', limit=3):
    '''
    Memory consumption monitoring function.
    https://stackoverflow.com/questions/552744/how-do-i-profile-memory-usage-in-python

    How to use it:
        tracemalloc.start()
        counts = Counter()

        .... # your commands to check for memory

        snapshot = tracemalloc.take_snapshot()
        display_top(snapshot, key_type='lineno', limit=10)
    
    '''
    
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        print("#%s: %s:%s: %.1f KiB"
              % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))

    return;



def data_read_and_select(sInputFiles, branchesToRead):
    data_all = pd.DataFrame()
    for sInputFile in sInputFiles:
        if printLevel >= 10: print("sInputFile: {}".format(sInputFile))
        with uproot.open(sInputFile) as fInputFile:            
            tree = fInputFile['Events']
            #tree = fInputFile.get('Events')
            #if printLevel >= 1: print("file: {}, nEvents: {}".format(sInputFile, uproot.numentries(sInputFile, 'Events')))
            if printLevel >= 1: print("file: {}, nEvents: {}".format(sInputFile, uproot.num_entries(sInputFile, 'Events')))

            for data_i in tree.iterate(branchesToRead,
                                       flatten=False, outputtype=pd.DataFrame,
                                       entrysteps=nEventToReadInBatch):

                ## select data... Put data selection cuts here --------------------------------------------------

                #fail = data_i[ ~ np.vectorize(cut_ObjectMultiplicity)(data_i.nMuon, nObjects_min=2, nObjects_max=99) ].index
                fail = data_i[ ~ np.vectorize(cut_ObjectMultiplicity)(data_i.nMuon, nObjects_min=3) ].index
                #fail = data_i[ ~ np.vectorize(cut_ObjectMultiplicity)(data_i.nMuon,  nObjects_max=1) ].index
                #data_i.drop(index=fail, inplace=True)
                data_i.drop(fail, inplace=True)
                if printLevel >= 10: print("\n\ndata_i selected 0: {}".format(data_i))

                #PtThrshs=[[25, 15, 10]]
                #print("PtThrshs_0 ({}): {}".format(type(PtThrshs), PtThrshs))
                #fail = data_i[ ~ np.vectorize(cut_ObjectPt)(data_i.Muon_pt, PtThrshs=PtThrshs)].index
                fail = data_i[ ~ np.vectorize(cut_ObjectPt)(data_i.Muon_pt, PtThrsh_Lead=25, PtThrsh_Sublead=15, PtThrsh_Third=10)].index
                data_i.drop(fail, inplace=True)
                print("\n\ndata_i selected 1: {}".format(data_i))

                fail = data_i[ ~ np.vectorize(cut_ObjectEta)(data_i.Muon_eta, EtaThrsh=2.3, nObjects=3) ].index
                data_i.drop(index=fail, inplace=True)
                
 
                data_i['mll_tmp'] = np.vectorize(calc_mll)(data_i.Muon_pt, data_i.Muon_eta, data_i.Muon_phi)
                print("\n\ndata_i selected 2: {}".format(data_i))
                
                if data_i.empty: continue
                data_all = data_all.append( data_i )


    return data_all;

                
def data_read_and_select_wNumpy(sInputFiles, branchesToRead, sOutFile=None):
    data_all = {}    
    for sInputFile in sInputFiles:
        if printLevel >= 10: print("sInputFile: {}".format(sInputFile))
        with uproot.open(sInputFile) as fInputFile:            
            tree = fInputFile['Events']
            #tree = fInputFile.get('Events')
            #if printLevel >= 1: print("file: {}, nEvents: {}".format(sInputFile, uproot.numentries(sInputFile, 'Events')))
            if printLevel >= 1: print("file: {}, nEvents: {}".format(sInputFile, tree.num_entries))
            if printLevel >= 100:
                print("tree.show(): {}".format(tree.show()))
                print("tree.keys(): {}, tree.typenames(): {}".format(tree.keys(), tree.typenames()))



            #tracemalloc.start()
            #counts = Counter()

            for data_i in tree.iterate(branchesToRead,
                                       library="np",
                                       step_size=100000):

                if printLevel >= 12:
                    print("data_i.keys() ({}): {}".format(type(data_i.keys()), data_i.keys()))
                    print("type(data_i): {}".format(type(data_i)))
                    print("type(data_i[data_i.keys()[0]]): {}".format(type(data_i[list(data_i)[0]])))
                if printLevel >= 12: 
                    #print("\ndata_i 0 ({}): \nkeys: {}".format(type(data_i), data_i.keys())); sys.stdout.flush();
                    print("\ndata_i 0 ({}):  {}".format(type(data_i), data_i)); sys.stdout.flush();
                
                    #print("data_i['nMuon'] ({})  ({}): {}".format(type(data_i['nMuon']), len(data_i['nMuon']), data_i['nMuon']))

                mask_nMu = data_i['nMuon'] >= 3
                if printLevel >= 20: 
                    print("\n\nmask: {}, \n\ndata_i['nMuon'][mask]: {}".format(mask_nMu, data_i['nMuon'][mask_nMu] ))
                    print("\n\ndata_i['Muon_pt'][mask: {}".format(data_i['Muon_pt'][mask_nMu]))

                
                apply_numpy_cut(data_i, mask_nMu);
                if printLevel >= 12: 
                    print("\ndata_i selected 1 ({}):  {}".format(type(data_i), data_i)); sys.stdout.flush();

                    
                mask_nEle = data_i['nElectron'] >= 2
                apply_numpy_cut(data_i, mask_nEle);
                if printLevel >= 10: 
                    print("\ndata_i selected 2 ({}):  {}".format(len(data_i[list(data_i)[0]]), data_i)); sys.stdout.flush();
                

                data_all = append_to_dictOfNpArrays(data_all, data_i)
                if printLevel >= 10:
                    #print("\ndata_all :  {}".format( data_all)); sys.stdout.flush();
                    print("\ndata_all ({}):  {}".format(len(data_all[list(data_all)[0]]), data_all)); sys.stdout.flush();


            if printLevel >= 9:                
                print("\ndata_all ({}):  {}".format(len(data_all[list(data_all)[0]]), data_all)); sys.stdout.flush();
                #print("\n\nnp.stack(data_all['Muon_pt'], axis=0): {}".format(np.stack(data_all['Muon_pt'], axis=0)))
                #print("Leading muon pt: {}".format(np.stack(data_all['Muon_pt'], axis=0)[:, 1]))

                for treeBranch in data_all.keys():
                    print("Br {}: type: {}, \t br[0] type: {}".format(treeBranch, type(data_all[treeBranch]), type(data_all[treeBranch][0]) ))
                print("ak.flatten(data_all['Muon_pt']: {}".format(ak.flatten(data_all['Muon_pt'])))




    ## Make histograms with selected dataset
    


    return




def apply_numpy_cut(data, mask):    
    for treeBranch in data.keys():
        #data[treeBranch] = np.delete(data[treeBranch], ~mask, axis=0) # this command used slightly more memory than the next command
        data[treeBranch] = data[treeBranch][mask]
                
    return 

def append_to_dictOfNpArrays(data_sum, data_i):
    if not data_sum:
        data_sum = data_i
    else:
        for treeBranch in data_i.keys():
            data_sum[treeBranch] = np.concatenate([data_sum[treeBranch], data_i[treeBranch]], axis=0)
    #print("data_sum: {}".format(data_sum))
    return data_sum
    
### -------------------------------------------------------------------------------------------------------------------------------------------


                
def data_read_and_select_wAkward(sInputFiles, branchesToRead, sOutFile=None):
    nMuonsToSelect           = 2
    nElectronsToSelect       = 2
    muonPtThrshsToSelect     = [4.0, 3.0]
    electronPtThrshsToSelect = [4.0, 3.0]
                    
    vector.register_awkward()
    
    data_all = None    
    for sInputFile in sInputFiles:
        if printLevel >= 10: print("sInputFile: {}".format(sInputFile))
        with uproot.open(sInputFile) as fInputFile:            
            tree = fInputFile['Events']
            #tree = fInputFile.get('Events')
            #if printLevel >= 1: print("file: {}, nEvents: {}".format(sInputFile, uproot.numentries(sInputFile, 'Events')))
            if printLevel >= 0: print("file: {}, nEvents: {}".format(sInputFile, tree.num_entries)); sys.stdout.flush();
            if printLevel >= 100:
                print("tree.show(): {}".format(tree.show()))
                print("tree.keys(): {}, tree.typenames(): {}".format(tree.keys(), tree.typenames()))


            htoaa_nanoAODBranchesToRead_info = [
                'run',
                'luminosityBlock',
                'event',
                #'nMuon',
                #'nElectron',
            ]

            htoaa_nanoAODBranchesToRead_muon = [
                'nMuon',
                'Muon_pt',
                'Muon_eta',
                'Muon_phi',
                'Muon_mass',
                'Muon_charge',
            ]

            htoaa_nanoAODBranchesToRead_electron = [
                'nElectron',
                'Electron_pt',
                'Electron_eta',
                'Electron_phi',
                'Electron_eCorr',
                'Electron_mass',
                'Electron_charge',
            ]
            


            branchesToRead_1 = htoaa_nanoAODBranchesToRead_info + htoaa_nanoAODBranchesToRead_muon + htoaa_nanoAODBranchesToRead_electron
            nEventsToAnalyze_1 = tree.num_entries if nEventsToAnalyze == -1 else nEventsToAnalyze
            
            #tracemalloc.start()
            #counts = Counter()
            iIteration = 0
            for data_i in tree.iterate(branchesToRead_1,
                                       library="ak",
                                       entry_start=0, entry_stop=nEventsToAnalyze_1,
                                       step_size=nEventToReadInBatch,
                                       ):
                '''
                uproot tree --> akward array:  list of event-level dict
                For e.g.
                [
                {“nMuon”: <>,  “Muon_pt”: []}, # event_0 
                {“nMuon”: <>,  “Muon_pt”: []}, # event_1
                ...
                ]
                '''

                if printLevel >= 0: print("nAnalyzedEvents: {}".format(iIteration * nEventToReadInBatch )); sys.stdout.flush();
                iIteration += 1
                
                # num_events_data_i = len(data_i)


                if printLevel >= 12:
                    #print(f"iIteration {iIteration}: data_i ({ak.count(data_i)}): {data_i.to_list()}")
                    print(f"iIteration {iIteration}: data_i ({len(data_i)}): {data_i[:]['Muon_pt'].to_list()}")


                    
                #mask_nMu_1 = data_i[:]['nMuon'] >= 3
                mask_nMu = np.vectorize(cut_ObjectMultiplicity)(data_i[:]['nMuon'], nObjects_min=nMuonsToSelect) 
                if ak.sum(mask_nMu) == 0: continue 
                data_i = data_i[mask_nMu]

                
                mask_nEle = np.vectorize(cut_ObjectMultiplicity)(data_i[:]['nElectron'], nObjects_min=nElectronsToSelect)
                if ak.sum(mask_nEle) == 0: continue 
                data_i = data_i[mask_nEle]

                if printLevel >= 12: 
                    print("\ndata_i selected 1 ({}) ({}):  {}".format(type(data_i), len(data_i), ak.to_list(data_i))); sys.stdout.flush();
                    
                
                #mask_MuPt = ((data_i[:, 'Muon_pt', 0] > 5) ) # & 
                #             (data_i[:, 'Muon_pt', 1] > 3) )
                mask_MuPt = cutAwkArray_ObjectPt(awkArray=data_i, objectName='Muon_pt', PtThrshs=muonPtThrshsToSelect)
                if ak.sum(mask_MuPt) == 0: continue   
                data_i = data_i[mask_MuPt]


                #mask_ElePt = ((data_i[:, 'Electron_pt', 0] > 20) & 
                #              (data_i[:, 'Electron_pt', 1] > 15)  )
                mask_ElePt = cutAwkArray_ObjectPt(awkArray=data_i, objectName='Electron_pt', PtThrshs=electronPtThrshsToSelect)
                if ak.sum(mask_ElePt) == 0: continue   
                data_i = data_i[mask_ElePt]

                
                if printLevel >= 12:
                    #print("\ndata_i selected 3 ({}) ({}):  {}".format(type(data_i), len(data_i), ak.to_list(data_i))); sys.stdout.flush();
                    print("data_i[:, 'Muon_eta', :) ({}): {}".format(type(data_i[:, 'Muon_eta', :]), data_i[:, 'Muon_eta', :]))
                    #print("abs(data_i[:, 'Muon_eta', :) ): {}".format(abs(data_i[:, 'Muon_eta', :])))

                    
                #mask_MuEta = abs(data_i[:, 'Muon_eta', :]) < 2.5
                #mask_MuEta = ((abs(data_i[:, 'Muon_eta', 0]) < 2.4) & 
                #              (abs(data_i[:, 'Muon_eta', 1]) < 2.4) )
                mask_MuEta = cutAwkArray_ObjectEta(awkArray=data_i, objectName='Muon_eta', EtaThrsh=2.4, nObjects=nMuonsToSelect)
                if ak.sum(mask_MuEta) == 0: continue   
                data_i = data_i[mask_MuEta]
                
                
                #mask_EleEta = ((abs(data_i[:, 'Electron_eta', 0]) < 2.4) & 
                #               (abs(data_i[:, 'Electron_eta', 1]) < 2.4) )
                mask_EleEta = cutAwkArray_ObjectEta(awkArray=data_i, objectName='Electron_eta', EtaThrsh=2.5, nObjects=nElectronsToSelect)
                if ak.sum(mask_EleEta) == 0: continue   
                data_i = data_i[mask_EleEta]


                mask_MuChargeSum = cutAwkArray_chargeSum(awkArray=data_i, objectName='Muon_charge', nObjects=nMuonsToSelect, chargeSumCondition=0)
                if ak.sum(mask_MuChargeSum) == 0: continue
                data_i = data_i[mask_MuChargeSum]

                mask_EleChargeSum = cutAwkArray_chargeSum(awkArray=data_i, objectName='Electron_charge', nObjects=nElectronsToSelect, chargeSumCondition=0)
                if ak.sum(mask_EleChargeSum) == 0: continue
                data_i = data_i[mask_EleChargeSum]



                
                
                if len(data_i) > 0:
                    if data_all is None:
                        data_all = data_i
                    else:
                        if printLevel >= 14: 
                            print("\nconcatenatedata_i selected 2p1 ({}) ({}):  {}".format(type(data_i), len(data_i), data_i)); sys.stdout.flush();
                            print("\ndata_all  ({}) ({}):  {}".format(type(data_all), len(data_all), ak.to_list(data_all))); sys.stdout.flush();
                        #ak.concatenate(data_all, data_i, axis=0)
                        #ak.concatenate([data_all, data_i], axis=0)
                        data_all = ak.concatenate([data_all, data_i], axis=0)
                        #print("here1")

                if printLevel >= 10:
                    #print("\ndata_all :  {}".format( data_all)); sys.stdout.flush();
                    print("\ndata_all  ({}) ({}):  {}".format(type(data_all), len(data_all), data_all)); sys.stdout.flush();
                    #print("\ndata_all_1  ({}) ({}):  {}".format(type(data_all_1), len(data_all_1), ak.to_list(data_all_1))); sys.stdout.flush();


            if printLevel >= 10:                
                print("\ndata_all final  ({}) ({}):  {}".format(type(data_all), len(data_all), ak.to_list(data_all))); sys.stdout.flush();
                #print("\n\nnp.stack(data_all['Muon_pt'], axis=0): {}".format(np.stack(data_all['Muon_pt'], axis=0)))
                #print("Leading muon pt: {}".format(np.stack(data_all['Muon_pt'], axis=0)[:, 1]))
                #print("\ndata_all_1 final ({}) ({}):  {}".format(type(data_all_1), len(data_all_1), data_all_1)); sys.stdout.flush();



    if printLevel >= 11:
        mu_pt = ak.flatten(data_all[:]['Muon_pt'])
        print("data_all[:]['Muon_pt']  ({})   ({}) : {}".format(type(mu_pt), len(mu_pt), mu_pt))
    #hnMuon, _ = np.histogram(aNMuon1, bins=50, range=(0, 100))

    
    if printLevel >= 1  and len(data_all) > 0:                
        print("\ndata_all final  ({}) ({}):  {}".format(type(data_all), len(data_all), ak.to_list(data_all))); sys.stdout.flush();

        print(f"\n\nSelected events ({len(data_all)}): {list(zip(data_all['run'], data_all['luminosityBlock'],data_all['event'] ))}")





    



        
                
    if sOutFile is not None  and len(data_all) > 0:
        if not sOutFile.endswith('.root'): sOutFile += '.root'
        sOutFile = sOutFile.replace('.root', '_wAwkwardArray.root')



        if printLevel >=12:
            data_tmp1 = data_all[:, ['Muon_pt', 'Muon_eta', 'Muon_phi', 'Muon_mass'], :]
            print(f"\n\ndata_tmp1 ({len(data_tmp1)}): {data_tmp1.to_list()}")
            print(f"data_all['Muon_pt']: {data_all['Muon_pt']}")
            print(f"aw.Array(data_all, with_name='Momentum4D'): {ak.Array(data_all, with_name='Momentum4D')}")

            #vMu = ak.Array(data_all, with_name='Momentum4D')
            vMu = getLorentzVectorFromAwkArray(data_tmp1, ptObjectName='Muon_pt', etaObjectName='Muon_eta', phiObjectName='Muon_phi', massObjectName='Muon_mass')
            print(f"vMu ({type(vMu)}) :{vMu}  list: {vMu.to_list()}")

            vMuCombinations = ak.combinations(vMu, 2, axis=1)
            print(f"ak.combinations(vMu, 2, axis=1) ({type(vMuCombinations)}): {vMuCombinations.to_list()}")

            vMus = ak.unzip(vMuCombinations)
            print(f"ak.unzip(vMuCombinations) {type(vMus)}: {vMus}")

            print(f"vMus[0]: {type(vMus[0])} {vMus[0].to_list()}")
            print(f"vMus[1]: {type(vMus[1])} {vMus[1].to_list()}")

            print(f"(vMus[0] + vMus[1]).mass: {(vMus[0] + vMus[1]).mass.to_list()}")
            print(f"(vMus[0] + vMus[1]).mass _1: {(vMus[0] + vMus[1]).mass[:, 0].to_list()}")

            print(f"(vMu[:, 0] + vMu[:, 1])  : {(vMu[:, 0] + vMu[:, 1]).to_list()}")
            print(f"(vMu[:, 0] + vMu[:, 1]).mass  : {(vMu[:, 0] + vMu[:, 1]).mass.to_list()}")

            print(f"vMu[:, 0].deltaR(vMu[:, 1].to_list(): {vMu[:, 0].deltaR(vMu[:, 1]).to_list()}")



        vMuons     = getLorentzVectorFromAwkArray(data_all, ptObjectName='Muon_pt',     etaObjectName='Muon_eta',     phiObjectName='Muon_phi',     massObjectName='Muon_mass')
        vElectrons = getLorentzVectorFromAwkArray(data_all, ptObjectName='Electron_pt', etaObjectName='Electron_eta', phiObjectName='Electron_phi', massObjectName='Electron_mass')


        
        with uproot.recreate(sOutFile) as fOut:
            #fOut['evt/hnMuon'] = np.histogram(ak.to_numpy(data_all[:]['nMuon']), bins=20, range=(-0.5, 19.5))

            hnMuon = hist.Hist.new.Reg(20, -0.5, 19.5, name="No. of muons").Weight()
            hnMuon.fill(ak.to_numpy(data_all[:, 'nMuon']))
            fOut['evt/hnMuon'] = hnMuon


            hMuonPt = hist.Hist.new.Reg(100, 0, 100, name="Muon pT").Weight()
            hMuonPt.fill(ak.to_numpy( ak.flatten(data_all[:, 'Muon_pt']) ))
            fOut['evt/hMuonPt'] = hMuonPt

            hMuonPt_lead = hist.Hist.new.Reg(100, 0, 100, name="Leading muon pT").Weight()
            hMuonPt_lead.fill(ak.to_numpy( data_all[:, 'Muon_pt', 0] ))
            fOut['evt/hLeadingMuonPt'] = hMuonPt_lead

            hMuonPt_sublead = hist.Hist.new.Reg(100, 0, 100, name="Sub-leading muon pT").Weight()
            hMuonPt_sublead.fill(ak.to_numpy( data_all[:, 'Muon_pt', 1] ))
            fOut['evt/hSubleadingMuonPt'] = hMuonPt_sublead

            '''
            hMuonPt_thridlead = hist.Hist.new.Reg(100, 0, 100, name="Third-leading muon pT").Weight()
            hMuonPt_thridlead.fill(ak.to_numpy( data_all[:, 'Muon_pt', 2] ))
            fOut['evt/hThirdleadingMuonPt'] = hMuonPt_thridlead
            '''

            ## Eta --
            hMuonEta = hist.Hist.new.Reg(40, -3, 3, name="Muon eta").Weight()
            hMuonEta.fill(ak.to_numpy( ak.flatten(data_all[:, 'Muon_eta']) ))
            fOut['evt/hMuonEta'] = hMuonEta

            
            ## Phi --
            hMuonPhi = hist.Hist.new.Reg(40, -3.14, 3.14, name="Muon phi").Weight()
            hMuonPhi.fill(ak.to_numpy( ak.flatten(data_all[:, 'Muon_phi']) ))
            fOut['evt/hMuonPhi'] = hMuonPhi


            ### Electron
            hnElectron = hist.Hist.new.Reg(20, -0.5, 19.5, name="No. of electron").Weight()
            hnElectron.fill(ak.to_numpy(data_all[:, 'nElectron']))
            fOut['evt/hnElectron'] = hnElectron
            
            hElectronPt = hist.Hist.new.Reg(100, 0, 100, name="Electron pT").Weight()
            hElectronPt.fill(ak.to_numpy( ak.flatten(data_all[:, 'Electron_pt']) ))
            fOut['evt/hElectronPt'] = hElectronPt

            hElectronPt_lead = hist.Hist.new.Reg(100, 0, 100, name="Leading Electron pT").Weight()
            hElectronPt_lead.fill(ak.to_numpy( data_all[:, 'Electron_pt', 0] ))
            fOut['evt/hLeadingElectronPt'] = hElectronPt_lead

            hElectronPt_sublead = hist.Hist.new.Reg(100, 0, 100, name="Sub-leading Electron pT").Weight()
            hElectronPt_sublead.fill(ak.to_numpy( data_all[:, 'Electron_pt', 1] ))


            ## Eta --
            hElectronEta = hist.Hist.new.Reg(40, -3, 3, name="Electron eta").Weight()
            hElectronEta.fill(ak.to_numpy( ak.flatten(data_all[:, 'Electron_eta']) ))
            fOut['evt/hElectronEta'] = hElectronEta

            
            ## Phi --
            hElectronPhi = hist.Hist.new.Reg(40, -3.14, 3.14, name="Electron phi").Weight()
            hElectronPhi.fill(ak.to_numpy( ak.flatten(data_all[:, 'Electron_phi']) ))
            fOut['evt/hElectronPhi'] = hElectronPhi




            hmass_2Mu = hist.Hist.new.Reg(100, 0, 200, name="mass(2 mu) [GeV]").Weight()
            hmass_2Mu.fill( (vMuons[:, 0] + vMuons[:, 1]).mass )
            fOut['evt/hmass_2Mu'] = hmass_2Mu

            hmass_2Ele = hist.Hist.new.Reg(100, 0, 200,  name="mass(2 ele) [GeV]").Weight()
            hmass_2Ele.fill( (vElectrons[:, 0] + vElectrons[:, 1]).mass )
            fOut['evt/hmass_2Ele'] = hmass_2Ele

            hdR_2Mu = hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 mu)").Weight()
            hdR_2Mu.fill( vMuons[:, 0].deltaR(vMuons[:, 1]) )
            fOut['evt/hdR_2Mu'] = hdR_2Mu

            hdR_2Ele = hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 ele)").Weight()
            hdR_2Ele.fill( vElectrons[:, 0].deltaR(vElectrons[:, 1]) )
            fOut['evt/hdR_2Ele'] = hdR_2Ele            

            hdR_2Mu_2Ele = hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 mu, 2 ele)").Weight()
            hdR_2Mu_2Ele.fill( (vMuons[:, 0] + vMuons[:, 1]).deltaR(vElectrons[:, 0] + vElectrons[:, 1]) )
            fOut['evt/hdR_2Mu_2Ele'] = hdR_2Mu_2Ele


            h2ddR_2Mu_vs_dR_2Ele = (
                hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 mu)")
                .Reg(100, 0, 3.14, name="deltaR(2 ele)")
                .Weight() )
            h2ddR_2Mu_vs_dR_2Ele.fill( vMuons[:, 0].deltaR(vMuons[:, 1]), vElectrons[:, 0].deltaR(vElectrons[:, 1]) )
            fOut['evt/h2ddR_2Mu_vs_dR_2Ele'] = h2ddR_2Mu_vs_dR_2Ele
            

            print("Wrote to sOutFile {}".format(sOutFile))

    

    return


def cutAwkArray_ObjectPt(awkArray, objectName, PtThrshs=[]):
    mask = None
    for iParticle, PtThrsh in enumerate(PtThrshs):
        mask_i = awkArray[:, objectName, iParticle] > PtThrsh
        if mask is None:
            mask = mask_i
        else:
            mask = (mask & mask_i)

        if printLevel >= 12:
            print(f"mask_i: {mask_i}")
            print(f"mask: {mask}")

    return mask;
    

def cutAwkArray_ObjectEta(awkArray, objectName, EtaThrsh, nObjects):
    '''
    Check Objects abs(Eta) is greater than thresholds set in EtaThrshs list.

    Return:
        True: All objects' Eta is about respective threshold
        False: Else false
    '''   

    mask = None
    for iParticle in range(nObjects):
        mask_i = abs(awkArray[:, objectName, iParticle]) < EtaThrsh
        if mask is None:
            mask = mask_i
        else:
            mask = (mask & mask_i)

        if printLevel >= 12:
            print(f"mask_i: {mask_i}")
            print(f"mask: {mask}")

    return mask;




def cutAwkArray_chargeSum(awkArray, objectName, nObjects, chargeSumCondition):
    if printLevel >= 12:
        print(f"awkArray[:, {objectName}, :{nObjects}]:  {awkArray[:, objectName, :nObjects] } ")
        print(f"awk.sum(awkArray[:, {objectName}, :{nObjects}], axis=-1):  {ak.sum(awkArray[:, objectName, :nObjects], axis=-1)} ")
        print(f"awk.sum(awkArray[:, {objectName}, :{nObjects}], axis=-1) == chargeSumCondition:  {(ak.sum(awkArray[:, objectName, :nObjects], axis=-1) == chargeSumCondition)} ")


    return (ak.sum(awkArray[:, objectName, :nObjects], axis=-1) == chargeSumCondition)


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
                       
### -------------------------------------------------------------------------------------------------------------------------------------------


















def data_read_and_select_wPandas(sInputFiles, branchesToRead, sOutFile=None):
    nMuonsToSelect           = 2
    nElectronsToSelect       = 2
    muonPtThrshsToSelect     = [4.0, 3.0]
    electronPtThrshsToSelect = [4.0, 3.0]
    
    data_all = None

    collection_info  = Collection('info')
    collection_mu    = Collection('Muon')
    collection_ele   = Collection('Electron')
    events           = Events(collections=[
        collection_info,
        collection_mu,
        collection_ele
    ])
    
    for sInputFile in sInputFiles:
        if printLevel >= 10: print("sInputFile: {}".format(sInputFile))
        with uproot.open(sInputFile) as fInputFile:            
            tree = fInputFile['Events']
            if printLevel >= 0: print("file: {}, nEvents: {}".format(sInputFile, tree.num_entries)); sys.stdout.flush();

            htoaa_nanoAODBranchesToRead_info = [
                'run',
                'luminosityBlock',
                'event',
                #'nMuon',
                #'nElectron',
            ]

            htoaa_nanoAODBranchesToRead_muon = [
                'nMuon',
                'Muon_pt',
                'Muon_eta',
                'Muon_phi',
                'Muon_mass',
                'Muon_charge',
            ]

            htoaa_nanoAODBranchesToRead_electron = [
                'nElectron',
                'Electron_pt',
                'Electron_eta',
                'Electron_phi',
                'Electron_eCorr',
                'Electron_mass',
                'Electron_charge',
            ]
            
            #tracemalloc.start()
            #counts = Counter()
            nEventsToAnalyze_1 = tree.num_entries if nEventsToAnalyze == -1 else nEventsToAnalyze
            for iEvtStart in range(0, nEventsToAnalyze_1, nEventToReadInBatch):
                if printLevel >= 1: print("\n\nRead entries from {} to {}".format(iEvtStart, (iEvtStart+nEventToReadInBatch-1)))
                data_subset = None

                data_subset_info = tree.arrays(
                    htoaa_nanoAODBranchesToRead_info,
                    entry_start=iEvtStart, entry_stop=(iEvtStart+nEventToReadInBatch),
                    library='pd'
                )
                collection_info.appendCollection(data_subset_info)

                data_subset_mu = tree.arrays(
                    htoaa_nanoAODBranchesToRead_muon,
                    entry_start=iEvtStart, entry_stop=(iEvtStart+nEventToReadInBatch),
                    library='pd'
                )
                collection_mu.appendCollection(data_subset_mu)

                data_subset_ele = tree.arrays(
                    htoaa_nanoAODBranchesToRead_electron,
                    entry_start=iEvtStart, entry_stop=(iEvtStart+nEventToReadInBatch),
                    library='pd'
                )
                collection_ele.appendCollection(data_subset_ele)

                if printLevel >= 13:
                    print("\n\n data_subset_info ({}): \n{}".format(len(data_subset_info.index), data_subset_info))
                    print("\n data_subset_mu ({}): \n{}".format(len(data_subset_mu.loc[pd.IndexSlice[:, 0], :].index), data_subset_mu.to_string()))
                    print("\n data_subset_ele ({}): \n{}".format(len(data_subset_ele.loc[pd.IndexSlice[:, 0], :].index), data_subset_ele))
                if printLevel >= 13:
                    events.print()



                mask_pass_all_cut = None

                #mask_pass_nMu_cut = collection_mu.df['nMuon'] >= 2
                mask_pass_nMu_cut = cutDf_ObjectMultiplicity(df=collection_mu.df, objectName='nMuon', nObjects_min=nMuonsToSelect)
                mask_pass_all_cut = mask_pass_nMu_cut


                #mask_pass_nEle_cut = collection_ele.df['nElectron'] >= 2
                mask_pass_nEle_cut = cutDf_ObjectMultiplicity(df=collection_ele.df, objectName='nElectron', nObjects_min=nElectronsToSelect)
                mask_pass_all_cut = (mask_pass_all_cut & mask_pass_nEle_cut)                
                if printLevel >= 12:
                    print("mask_pass_nMu_cut  ({}): {}".format(len(mask_pass_nMu_cut.index), mask_pass_nMu_cut))
                    print("mask_pass_nEle_cut ({}): {}".format(len(mask_pass_nEle_cut.index), mask_pass_nEle_cut))
                    print("mask_pass_all_cut  ({}): {}".format(len(mask_pass_all_cut.index), mask_pass_all_cut))

                    
                #events.dropEventsFailingMask(mask_pass_cut = mask_pass_all_cut)
                #print("\n Dayaset selected after nMuon, nElectron cut: :: Events: \n"); events.print()

                if printLevel >= 12:
                    #print("\n data_subset_mu _2 ({}): \n{}".format(len(data_subset_mu.loc[pd.IndexSlice[:, 0], :].index), data_subset_mu.to_string()))
                    print("\n data_subset_mu  _2 ({}): \n{}".format(getDataframeNEvents(data_subset_mu), data_subset_mu.to_string()))
                    print("\n data_subset_ele _2 ({}): \n{}".format(getDataframeNEvents(data_subset_ele), data_subset_ele.to_string()))
                
                
                '''
                mask_pass_MuPt_cut = (
                    (collection_mu.df.loc[pd.IndexSlice[:, 0], 'Muon_pt'] > 15) &
                    (collection_mu.df.loc[pd.IndexSlice[:, 1], 'Muon_pt'] > 10)
                )
                mask_pass_all_cut = mask_pass_all_cut & mask_pass_MuPt_cut


                mask_pass_MuEta_cut = (
                    (abs(collection_mu.df.loc[pd.IndexSlice[:, 0], 'Muon_eta']) < 2.5) &
                    (abs(collection_mu.df.loc[pd.IndexSlice[:, 1], 'Muon_eta']) < 2.5)
                )
                '''


                mask_pass_MuPt_cut = cutDf_ObjectPt(df=collection_mu.df, objectName='Muon_pt', PtThrshs=muonPtThrshsToSelect)
                mask_pass_all_cut = (mask_pass_all_cut & mask_pass_MuPt_cut)
                 
                if printLevel >= 12:
                    print("mask_pass_MuPt_cut  ({}): {}".format(len(mask_pass_MuPt_cut.index), mask_pass_MuPt_cut))
                    
                mask_pass_ElePt_cut = cutDf_ObjectPt(df=collection_ele.df, objectName='Electron_pt', PtThrshs=electronPtThrshsToSelect)
                mask_pass_all_cut = (mask_pass_all_cut & mask_pass_ElePt_cut)

                
                mask_pass_MuEta_cut = cutDf_ObjectEta(df=collection_mu.df, objectName='Muon_eta', EtaThrsh=2.4, nObjects=nMuonsToSelect)
                mask_pass_all_cut = (mask_pass_all_cut & mask_pass_MuEta_cut)               
                
                mask_pass_EleEta_cut = cutDf_ObjectEta(df=collection_ele.df, objectName='Electron_eta', EtaThrsh=2.5, nObjects=nElectronsToSelect)
                mask_pass_all_cut = (mask_pass_all_cut & mask_pass_EleEta_cut)                


                
                #events.dropEventsFailingMask(mask_pass_cut = mask_pass_all_cut); events.print()

                
                
                mask_pass_MuChargeSum_cut = cutDf_chargeSum(df=collection_mu.df, objectName='Muon_charge', nObjects=nMuonsToSelect, chargeSumCondition=0)
                mask_pass_all_cut = (mask_pass_all_cut & mask_pass_MuChargeSum_cut)
                
                mask_pass_EleChargeSum_cut = cutDf_chargeSum(df=collection_ele.df, objectName='Electron_charge', nObjects=nElectronsToSelect, chargeSumCondition=0)
                mask_pass_all_cut = (mask_pass_all_cut & mask_pass_EleChargeSum_cut)
                   
                if printLevel >= 12:
                    print("mask_pass_all_cut ({}): {}".format(len(mask_pass_all_cut.index), mask_pass_all_cut.to_string()))

                    
                if getDataframeNEvents(mask_pass_all_cut) == 0: continue




                events.dropEventsFailingMask(mask_pass_cut = mask_pass_all_cut)

                
                if printLevel >= 6:
                    print("\n Dayaset selected untilnow: :: Events: \n")
                    events.print()

                

                                

                

    if printLevel >= 5:
        print("\n\n\n Final dataset:: Events: \n")
        events.print()



    nMuons = nMuonsToSelect
    nElectrons = nElectronsToSelect
    vMuons = []
    for iParticle in range(nMuons):
        vMu = getLorentzVectorFromDf(df=collection_mu.df, indexLevel1=iParticle,
                                     ptObjectName='Muon_pt', etaObjectName='Muon_eta', phiObjectName='Muon_phi', massObjectName='Muon_mass')
        vMuons.append( vMu )


    vElectrons = []
    for iParticle in range(nElectrons):
        vEle = getLorentzVectorFromDf(df=collection_ele.df, indexLevel1=iParticle,
                                     ptObjectName='Electron_pt', etaObjectName='Electron_eta', phiObjectName='Electron_phi', massObjectName='Electron_mass')
        vElectrons.append( vEle )



        
    if printLevel >= 5:
        print("nEvents after selection: {}".format(events.nEvents()))


    if sOutFile is not None  and events.nEvents() > 0:
        if not sOutFile.endswith('.root'): sOutFile += '.root'
        sOutFile = sOutFile.replace('.root', '_wPandas.root')
        
        with uproot.recreate(sOutFile) as fOut:
            #fOut['evt/hnMuon'] = np.histogram(ak.to_numpy(data_all[:]['nMuon']), bins=20, range=(-0.5, 19.5))

            hnMuon = hist.Hist.new.Reg(20, -0.5, 19.5, name="No. of muons").Weight()
            hnMuon.fill( collection_mu.df['nMuon'].to_numpy() )
            fOut['evt/hnMuon'] = hnMuon


            hMuonPt = hist.Hist.new.Reg(100, 0, 100, name="Muon pT").Weight()
            hMuonPt.fill( collection_mu.df['Muon_pt'].to_numpy() )
            fOut['evt/hMuonPt'] = hMuonPt

            hMuonPt_lead = hist.Hist.new.Reg(100, 0, 100, name="Leading muon pT").Weight()
            hMuonPt_lead.fill( collection_mu.df.loc[pd.IndexSlice[:, 0], 'Muon_pt'].to_numpy() ) 
            fOut['evt/hLeadingMuonPt'] = hMuonPt_lead

            hMuonPt_sublead = hist.Hist.new.Reg(100, 0, 100, name="Sub-leading muon pT").Weight()
            hMuonPt_sublead.fill( collection_mu.df.loc[pd.IndexSlice[:, 1], 'Muon_pt'].to_numpy() )
            fOut['evt/hSubleadingMuonPt'] = hMuonPt_sublead

            '''
            hMuonPt_thridlead = hist.Hist.new.Reg(100, 0, 100, name="Third-leading muon pT").Weight()
            hMuonPt_thridlead.fill(ak.to_numpy( data_all[:, 'Muon_pt', 2] ))
            fOut['evt/hThirdleadingMuonPt'] = hMuonPt_thridlead
            '''

            ## Eta --
            hMuonEta = hist.Hist.new.Reg(40, -3, 3, name="Muon eta").Weight()
            hMuonEta.fill( collection_mu.df['Muon_eta'].to_numpy() )
            fOut['evt/hMuonEta'] = hMuonEta

            
            ## Phi --
            hMuonPhi = hist.Hist.new.Reg(40, -3.14, 3.14, name="Muon phi").Weight()
            hMuonPhi.fill( collection_mu.df['Muon_phi'].to_numpy() )
            fOut['evt/hMuonPhi'] = hMuonPhi


            ### Electron

            hElectronPt = hist.Hist.new.Reg(100, 0, 100, name="Electron pT").Weight()
            hElectronPt.fill( collection_ele.df['Electron_pt'].to_numpy() )
            fOut['evt/hElectronPt'] = hElectronPt

            hElectronPt_lead = hist.Hist.new.Reg(100, 0, 100, name="Leading Electron pT").Weight()
            hElectronPt_lead.fill( collection_ele.df.loc[pd.IndexSlice[:, 0], 'Electron_pt'].to_numpy() ) 
            fOut['evt/hLeadingElectronPt'] = hElectronPt_lead

            hElectronPt_sublead = hist.Hist.new.Reg(100, 0, 100, name="Sub-leading Electron pT").Weight()
            hElectronPt_sublead.fill( collection_ele.df.loc[pd.IndexSlice[:, 1], 'Electron_pt'].to_numpy() )
            fOut['evt/hSubleadingElectronPt'] = hElectronPt_sublead

            '''
            hElectronPt_thridlead = hist.Hist.new.Reg(100, 0, 100, name="Third-leading Electron pT").Weight()
            hElectronPt_thridlead.fill(ak.to_numpy( data_all[:, 'Electron_pt', 2] ))
            fOut['evt/hThirdleadingElectronPt'] = hElectronPt_thridlead
            '''

            ## Eta --
            hElectronEta = hist.Hist.new.Reg(40, -3, 3, name="Electron eta").Weight()
            hElectronEta.fill( collection_ele.df['Electron_eta'].to_numpy() )
            fOut['evt/hElectronEta'] = hElectronEta

            
            ## Phi --
            hElectronPhi = hist.Hist.new.Reg(40, -3.14, 3.14, name="Electron phi").Weight()
            hElectronPhi.fill( collection_ele.df['Electron_phi'].to_numpy() )
            fOut['evt/hElectronPhi'] = hElectronPhi



            hmass_2Mu = hist.Hist.new.Reg(100, 0, 200, name="mass(2 mu) [GeV]").Weight()
            hmass_2Mu.fill( (vMuons[0] + vMuons[1]).mass )
            fOut['evt/hmass_2Mu'] = hmass_2Mu

            hmass_2Ele = hist.Hist.new.Reg(100, 0, 200,  name="mass(2 ele) [GeV]").Weight()
            hmass_2Ele.fill( (vElectrons[0] + vElectrons[1]).mass )
            fOut['evt/hmass_2Ele'] = hmass_2Ele

            hdR_2Mu = hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 mu)").Weight()
            hdR_2Mu.fill( vMuons[0].deltaR(vMuons[1]) )
            fOut['evt/hdR_2Mu'] = hdR_2Mu

            hdR_2Ele = hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 ele)").Weight()
            hdR_2Ele.fill( vElectrons[0].deltaR(vElectrons[1]) )
            fOut['evt/hdR_2Ele'] = hdR_2Ele            

            hdR_2Mu_2Ele = hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 mu, 2 ele)").Weight()
            hdR_2Mu_2Ele.fill( (vMuons[0] + vMuons[1]).deltaR(vElectrons[0] + vElectrons[1]) )
            fOut['evt/hdR_2Mu_2Ele'] = hdR_2Mu_2Ele


            h2ddR_2Mu_vs_dR_2Ele = (
                hist.Hist.new.Reg(100, 0, 3.14, name="deltaR(2 mu)")
                .Reg(100, 0, 3.14, name="deltaR(2 ele)")
                .Weight() )
            h2ddR_2Mu_vs_dR_2Ele.fill( vMuons[0].deltaR(vMuons[1]), vElectrons[0].deltaR(vElectrons[1]) )
            fOut['evt/h2ddR_2Mu_vs_dR_2Ele'] = h2ddR_2Mu_vs_dR_2Ele
            

            print("Wrote to sOutFile {}".format(sOutFile))

    


class Collection:
    def __init__(self, name):
        self.name = name
        self.df = None #pd.DataFrame()

    def appendCollection(self, df_new):
        if len(df_new.index) == 0: return; # check if df_new is empty

        if self.df is None:
            self.df = df_new
        else:
            self.df = pd.concat([self.df, df_new])

    def print(self):
        #if len(df_new.index) == 0:
            
        '''
        nEvents = None
        index_df = self.df
        if not isinstance(self.df.index, pd.MultiIndex):
            nEvents = len(self.df.index)
        else:
            nEvents = len(self.df.loc[pd.IndexSlice[:, 0], :].index)
        '''
        print("Collection {}: nEvents {}.  \n{}".format(self.name, getDataframeNEvents(self.df), self.df))

        

class Events:
    def __init__(self, collections = []):
        self.collections = OD()
        for collection in collections:
            self.collections[collection.name] = collection

            
    def print(self):
        print("Events::print():")
        for collection in self.collections.values():
            collection.print()

            
    def dropEventsFailingMask(self, mask_pass_cut):
        
        if printLevel >= 14:
            print("Events::dropEventsFailingMask():: mask_pass_cut ({}): {}".format(getDataframeNEvents(mask_pass_cut), mask_pass_cut.to_string()))
            print("Events::dropEventsFailingMask():: mask_pass_cut[mask_pass_cut] ({}): {}".format(getDataframeNEvents(mask_pass_cut[mask_pass_cut]), mask_pass_cut[mask_pass_cut].to_string()))
            #print("Events::dropEventsFailingMask():: ".format())
            mask_pass_events = mask_pass_cut.groupby(axis=0, level=0).all()
            mask_pass_events_1 = mask_pass_events[mask_pass_events]
            selected_events = mask_pass_events_1.index
            print("Events::dropEventsFailingMask():: mask_pass_cut.groupby(axis=0, level=0).all() ({}): {}".format(type(mask_pass_events), mask_pass_events.to_string()))
            print("Events::dropEventsFailingMask():: mask_pass_cut.groupby(axis=0, level=0).all() selected events ({}): {}".format(type(mask_pass_events_1), mask_pass_events_1.to_string()))
            print("Events::dropEventsFailingMask():: mask_pass_cut.groupby(axis=0, level=0).all() selected events ({}): {} \n{}".format(type(selected_events), selected_events, set(selected_events)))
            #print("Events::dropEventsFailingMask()::  mask_pass_cut[ mask_pass_cut.groupby(axis=0, level=0).all() ] ({}): {}".format(getDataframeNEvents(mask_pass_cut[ mask_pass_cut.groupby(axis=0, level=0).all() ]), mask_pass_cut[ mask_pass_cut.groupby(axis=0, level=0).all() ].to_string() ) )
            
            

        indexLevel0_allCollections = set() # list of 'index level 0' ('entries') in all  collections combined. It represents list of events.
        indexLevel0_collectionwise = OD()
        for name, collection in self.collections.items():
            #indexLevel0_allCollections += set( getDataframeIndexLevel0(collection.df) )
            indexLevel0_collectionwise[name] = set( getDataframeIndexLevel0(collection.df) )
            indexLevel0_allCollections.update( indexLevel0_collectionwise[name] )

        indexLevel0_eventsPassingCut = mask_pass_cut.groupby(axis=0, level=0).all() # Boolean pd.Series with 'index=index level0 of mask_pass_cut' seting True when all subentries in mask_pass_cut are true
        indexLevel0_eventsPassingCut = set( indexLevel0_eventsPassingCut[indexLevel0_eventsPassingCut].index ) # get entry when all subentries in 'mask_pass_cut' are true
        indexLevel0_eventsToDrop =  indexLevel0_allCollections - indexLevel0_eventsPassingCut 
        
        if printLevel >= 11:
            print("Events::dropEventsFailingMask():: indexLevel0_allCollections ({}): {}".format(len(indexLevel0_allCollections), indexLevel0_allCollections))
            print("Events::dropEventsFailingMask():: indexLevel0_eventsPassingCut ({}): {}".format(len(indexLevel0_eventsPassingCut), indexLevel0_eventsPassingCut))
            print("Events::dropEventsFailingMask():: indexLevel0_eventsToDrop ({}): {}".format(len(indexLevel0_eventsToDrop), indexLevel0_eventsToDrop))
            

            
        for name, collection in self.collections.items():
            if printLevel >= 12:
                print("Events::dropEventsFailingMask():: Print collection {} _0".format(name))
                collection.print()

            indexLevel0_eventsToDropFromCollection = indexLevel0_eventsToDrop.intersection( indexLevel0_collectionwise[name] )
            if printLevel >= 12:
                print("Events::dropEventsFailingMask():: indexLevel0_eventsToDropFromCollection ({}): {}".format(len(indexLevel0_eventsToDropFromCollection), indexLevel0_eventsToDropFromCollection))
            if len(indexLevel0_eventsToDropFromCollection) == 0: continue

                                                    
            if not isinstance(collection.df, pd.MultiIndex):
                collection.df.drop(list(indexLevel0_eventsToDropFromCollection), axis=0,          inplace=True)
            else:
                collection.df.drop(list(indexLevel0_eventsToDropFromCollection), axis=0, level=0, inplace=True)

            if printLevel >= 12:
                print("Events::dropEventsFailingMask():: Print collection {} _1".format(name))
                collection.print()
        

        return

    
    def nEvents(self):
        #return getDataframeNEvents( self.collections.items()[0][1].df )
        return getDataframeNEvents( next(iter(self.collections.items()))[1].df ) # next(iter(self.collections.items()))[1] : get first collection



def getDataframeNEvents(df):
    nEvents = 0
    if len(df.index) == 0: return 0
    
    #index_df = df
    #print("getDataframeNEvents():: df ({}): {}".format(type(df), df))
    #print("getDataframeNEvents():: df.index ({}): {}".format(type(df.index), df.index))
    if not isinstance(df.index, pd.MultiIndex):
        nEvents = len(df.index)
    else:
        if isinstance(df, pd.Series):
            nEvents = len(df.loc[pd.IndexSlice[:, 0]].index)
        else:
            nEvents = len(df.loc[pd.IndexSlice[:, 0], :].index)
    return nEvents


def getDataframeIndexLevel0(df):
    indexLevel0 = []
    if len(df.index) == 0: return indexLevel0
    
    if not isinstance(df.index, pd.MultiIndex):
        indexLevel0 = df.index
    else:
        if isinstance(df, pd.Series):
            indexLevel0 = df.loc[pd.IndexSlice[:, 0]   ].index.get_level_values(level=0)
        else:
            indexLevel0 = df.loc[pd.IndexSlice[:, 0], :].index.get_level_values(level=0)
    return indexLevel0


def cutDf_ObjectMultiplicity(df, objectName, nObjects_min=None, nObjects_max=None):
    '''
    Check if "nObjects are with nObjects_min and nObjects_max", if they are specified.
    If either of nObjects_min and nObjects_max are not specified, then corresponding condition is not checked.

    Return:
        pd.Series of boolean, with 
            True: if nObjects passes the condition.
            False: if nObjects fails the condition
    '''
    mask = mask_low = mask_up  = None
    if nObjects_min is not None:
        mask_low = (df[objectName] >= nObjects_min)
    if nObjects_max is not None:
        mask_up  = (df[objectName] <= nObjects_max)

    if (nObjects_min is not None) and (nObjects_max is not None):
        mask = (mask_low and mask_up)
    elif (nObjects_min is not None):
        mask = mask_low
    else:
        mask = mask_up

    return mask





def cutDf_ObjectPt(df, objectName, PtThrshs=[]): 
    '''
    Check Objects Pt is above their resepctive thresholds.
    If PtThrsh_<rank> is not set, then their Pt condition is not checked.

    Return:
      pd.Series of boolean, with 
        True: All objects' Pt is about respective threshold
        False: Else false
    '''
    if len(PtThrshs) == 0:
        raise Exception("cutDf_ObjectPt():: PtThrshs not provided \t\t **** ERROR ****")
    
    mask_failingCut = df[objectName] != df[objectName] # set mask for all entry-subentry rows to False
    #print("cutDf_ObjectPt():: mask before  ({}): {} \n{}".format(getDataframeNEvents(mask),   mask.tolist(), mask.to_string()))
    for iParticle, PtTrsh in enumerate(PtThrshs):
        mask_failingCut_i = df.loc[pd.IndexSlice[:, iParticle], objectName] < PtTrsh # set mask to True for failing particles
        # mask = (mask & mask_i)        
        mask_failingCut = (mask_failingCut |
                (df.loc[pd.IndexSlice[:, iParticle], objectName] < PtTrsh) ) # set entry-subentry combination True for failing particles **** IMPORTANT ****
        
        
        if printLevel >= 12:
            print("cutDf_ObjectPt():: iParticle {}, PtTrsh {}".format(iParticle, PtTrsh))
            print("cutDf_ObjectPt():: mask   ({}): {} \n{}".format(getDataframeNEvents((~mask_failingCut)),   (~mask_failingCut).tolist(), (~mask_failingCut).to_string()))
            print("cutDf_ObjectPt():: mask_i ({}): {} \n{}".format(len((~mask_failingCut_i).index), (~mask_failingCut_i).tolist(), (~mask_failingCut_i).to_string()))

    
    return ( ~ mask_failingCut);



def cutDf_ObjectEta(df, objectName, EtaThrsh, nObjects):
    '''
    Check Objects abs(Eta) is greater than thresholds set in EtaThrshs list.

    Return:
        True: All objects' Eta is about respective threshold
        False: Else false
    '''   
    
    mask_failingCut = df[objectName] != df[objectName] # set mask for all entry-subentry rows to False 
    for iObject in range(nObjects):
        mask_failingCut = (mask_failingCut  |
                (abs(df.loc[pd.IndexSlice[:, iObject], objectName]) > EtaThrsh) ) # set entry-subentry combination True for failing particles **** IMPORTANT ****
        if printLevel >= 12:
            print("iObject, mask ({}): {}".format(iObject, getDataframeNEvents((~mask_failingCut)), (~mask_failingCut).to_string()))
            
    return ( ~ mask_failingCut);



def cutDf_chargeSum(df, objectName, nObjects, chargeSumCondition):
    if printLevel >= 12:
        df_tmp = df.groupby(axis=0, level=0)
        print("cutDf_chargeSum(): df.group(axis=0, level=0) ({}): \n{}".format(type(df_tmp), df_tmp))
        print("cutDf_chargeSum(): df.group(axis=0, level=0) _1 ({}): \n{}".format(type(df.groupby(axis=0, level=0)), df.groupby(axis=0, level=0)[objectName].sum()))

        print(f"df.loc[pd.IndexSlice[:, range(nObjects)], objectName]: {df.loc[pd.IndexSlice[:, range(nObjects)], objectName]}")
        #print(f"df.loc[pd.IndexSlice[:, range(nObjects)], objectName].sum(axis=0, level=0): {df.loc[pd.IndexSlice[:, range(nObjects)], objectName].sum(axis=0, level=0)}")
        #print(f"df.loc[pd.IndexSlice[:, range(1)], objectName]: {df.loc[pd.IndexSlice[:, range(1)], objectName].sum(axis=0, level=0)}")

        #print(f"(df.loc[pd.IndexSlice[:, range(nObjects)], objectName].sum(axis=0, level=0)== chargeSumCondition): {(df.loc[pd.IndexSlice[:, range(nObjects)], objectName].sum(axis=0, level=0) == chargeSumCondition)}")

        print(f"(df.loc[pd.IndexSlice[:, range(nObjects)], objectName]).groupby(axis=0, level=0).sum():  { (df.loc[pd.IndexSlice[:, range(nObjects)], objectName]).groupby(axis=0, level=0).sum()}")
        print(f"( (df.loc[pd.IndexSlice[:, range(nObjects)], objectName]).groupby(axis=0, level=0).sum() == chargeSumCondition):  {( (df.loc[pd.IndexSlice[:, range(nObjects)], objectName]).groupby(axis=0, level=0).sum() == chargeSumCondition)}")
        
    #return (df.groupby(axis=0, level=0)[objectName].sum() == chargeSumCondition)
    return ( (df.loc[pd.IndexSlice[:, range(nObjects)], objectName]).groupby(axis=0, level=0).sum() == chargeSumCondition) 
                       
    

def getLorentzVectorFromDf(df, indexLevel1, ptObjectName, etaObjectName, phiObjectName, massObjectName):
    return vector.array({'pt':  df.loc[pd.IndexSlice[:, indexLevel1], ptObjectName],
                         'eta': df.loc[pd.IndexSlice[:, indexLevel1], etaObjectName],
                         'phi': df.loc[pd.IndexSlice[:, indexLevel1], phiObjectName],
                         'mass':df.loc[pd.IndexSlice[:, indexLevel1], massObjectName]})
### -------------------------------------------------------------------------------------------------------------------------------------------



def calc_mll(lep_pts,lep_etas,lep_phis):
    theta_0 = 2*math.atan(math.exp(-lep_etas[0]))
    theta_1 = 2*math.atan(math.exp(-lep_etas[1]))
    #theta_2 = 2*math.atan(math.exp(-lep_etas[2]))
    #theta_3 = 2*math.atan(math.exp(-lep_etas[3]))
    p_0 = lep_pts[0]/math.sin(theta_0)
    p_1 = lep_pts[1]/math.sin(theta_1)
    #p_2 = lep_pts[2]/math.sin(theta_2)
    #p_3 = lep_pts[3]/math.sin(theta_3)
    pz_0 = p_0*math.cos(theta_0)
    pz_1 = p_1*math.cos(theta_1)
    #pz_2 = p_2*math.cos(theta_2)
    #pz_3 = p_3*math.cos(theta_3)
    px_0 = p_0*math.sin(theta_0)*math.cos(lep_phis[0])
    px_1 = p_1*math.sin(theta_1)*math.cos(lep_phis[1])
    #px_2 = p_2*math.sin(theta_2)*math.cos(lep_phis[2])
    #px_3 = p_3*math.sin(theta_3)*math.cos(lep_phis[3])
    py_0 = p_0*math.sin(theta_0)*math.sin(lep_phis[0])
    py_1 = p_1*math.sin(theta_1)*math.sin(lep_phis[1])
    #py_2 = p_2*math.sin(theta_2)*math.sin(lep_phis[2])
    #py_3 = p_3*math.sin(theta_3)*math.sin(lep_phis[3])
    sumpz = pz_0 + pz_1 #+ pz_2 + pz_3
    sumpx = px_0 + px_1 #+ px_2 + px_3
    sumpy = py_0 + py_1 #+ py_2 + py_3
    sumE = p_0 + p_1 #+ p_2 + p_3
    mllll = sumE**2 - sumpz**2 - sumpx**2 - sumpy**2
    if mllll < 0: mllll = 0
    return math.sqrt(mllll)# /1000 #/1000 to go from MeV to GeV


#def plot_data(data):
    
    


if __name__ == '__main__':
    print("htoaa_Analysis:: main: {}".format(sys.argv))

    if len(sys.argv) != 2:
        print("htoaa_Analysis:: Command-line config file missing.. \t **** ERROR **** \n")

    sConfig = sys.argv[1]

    ''' Json file does not support comments. Hence a simple tweak using GetDictFromJsonFile(sConfig)
    config = OD()
    with open(sConfig) as fConfig:
        config = json.load(fConfig)
    '''
    
    config = GetDictFromJsonFile(sConfig)
    print("Config {}: \n{}".format(sConfig, json.dumps(config, indent=4)))
    
    sInputFiles  = config["inputFiles"]
    sOutputFile  = config["outputFile"]
    era = config['era']
    luminosity = Luminosities[era][0]
    sample_crossSection = config["crossSection"]
    sample_nEvents = config["nEvents"]
    sample_sumEvents = config["sumEvents"] if config["sumEvents"] != -1 else sample_nEvents
    if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
    lumiScale = calculate_lumiScale(luminosity=luminosity, crossSection=sample_crossSection, sumEvents=sample_sumEvents)
    branchesToRead = htoaa_nanoAODBranchesToRead
    print("branchesToRead: {}".format(branchesToRead))


    startTime = time.time()
    
    tracemalloc.start()
    #counts = Counter()

    #data_read_and_select_wNumpy(sInputFiles, branchesToRead)

    data_read_and_select_wAkward(sInputFiles, branchesToRead, sOutputFile)
    #data_read_and_select_wPandas(sInputFiles, branchesToRead, sOutputFile)

    #snapshot = tracemalloc.take_snapshot()
    #display_top(snapshot, key_type='lineno', limit=10)

    
    current_memory, peak_memory = tracemalloc.get_traced_memory() # https://medium.com/survata-engineering-blog/monitoring-memory-usage-of-a-running-python-program-49f027e3d1ba
    print(f"\n\nMemory usage:: current {current_memory / 10**6}MB;  peak {peak_memory / 10**6}MB")

    endTime = time.time()
    totalTime = endTime - startTime
    totalTime_hr  = int(totalTime/60/60)
    totalTime_min = totalTime - float(totalTime_hr * 60)
    totalTime_min = int(totalTime_min/60)
    totalTime_sec = totalTime - float(totalTime_hr * 60*60) - float(totalTime_min * 60)
    print(f"Total run time: {totalTime_hr}h {totalTime_min}m {totalTime_sec}s = {totalTime}sec ")
