#htoaa analysis main code

import os
import sys
import json
from collections import OrderedDict as OD
import pandas as pd
import uproot
#import uproot3 as uproot
import numpy as np
import math
import awkward as ak
import matplotlib.pyplot as plt
#import mplhep as hep # for histograms like CMS template
import hist # for histogram

from collections import Counter
import linecache
import os
import tracemalloc

from htoaa_Settings import *
from htoaa_NanoAODBranches import htoaa_nanoAODBranchesToRead
from htoaa_CommonTools import GetDictFromJsonFile, DfColLabel_convert_bytes_to_string, calculate_lumiScale
from htoaa_CommonTools import cut_ObjectMultiplicity, cut_ObjectPt, cut_ObjectEta, cut_ObjectPt_1


printLevel = 10
nEventToReadInBatch = 1000 # 2500000
nEventsToAnalyze = 5000
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
    



                
def data_read_and_select_wAkward(sInputFiles, branchesToRead, sOutFile=None):
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



            #tracemalloc.start()
            #counts = Counter()
            iIteration = 0
            for data_i in tree.iterate(branchesToRead,
                                       library="ak",
                                       step_size=nEventToReadInBatch):
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
                    #print("\ndata_i 0 ({}): \nkeys: {}".format(type(data_i), data_i.keys())); sys.stdout.flush();
                    print("\ndata_i 0 ({}):  {}".format(type(data_i), data_i)); sys.stdout.flush();
                    print("data_i[0]: {}".format(data_i[0]))

                    print("ak.count(data_i): {}".format(ak.count(data_i)))
                    print("ak.num(data_i, axis=0): {}".format(ak.num(data_i, axis=0)))
                    print("len(data_i): {}".format(len(data_i)))
                    print("\ndata_i[:]: {}".format(data_i[:]))
                    print("\ndata_i[:]['nMuon'] ({}): {}".format(type(data_i[:]['nMuon']), data_i[:]['nMuon']))
                    print("\ndata_i[:]['Muon_pt'] ({}): {}".format(type(data_i[:]['Muon_pt']), data_i[:]['Muon_pt']))
                    #print("data_i['nMuon'] ({})  ({}): {}".format(type(data_i['nMuon']), len(data_i['nMuon']), data_i['nMuon']))

                    
                mask_nMu_1 = data_i[:]['nMuon'] >= 3
                mask_nMu = np.vectorize(cut_ObjectMultiplicity)(data_i[:]['nMuon'], nObjects_min=3)
                if printLevel >= 12:
                    print("mask_nMu ({}): {}".format(len(mask_nMu), list(mask_nMu)))
                    print("np.sum(mask_nMu): {}".format(np.sum(mask_nMu)))
                    print("mask_nMu_1 ({}): {}".format(len(mask_nMu_1), list(mask_nMu_1)))
                    print("np.sum(mask_nMu_1): {}".format(np.sum(mask_nMu_1)))
                    print("\n\nmask: {}, \n\ndata_i['nMuon'][mask]: {}".format(mask_nMu, data_i[mask_nMu]['nMuon'] ))
                    print("\n\ndata_i[mask]['Muon_pt']: {}".format(data_i[mask_nMu]['Muon_pt']))


                if ak.sum(mask_nMu) == 0: continue 
                data_i = data_i[mask_nMu]

                
                mask_nEle = np.vectorize(cut_ObjectMultiplicity)(data_i[:]['nElectron'], nObjects_min=2)

                if printLevel >= 12: 
                    print("\ndata_i selected 1 ({}) ({}):  {}".format(type(data_i), len(data_i), ak.to_list(data_i))); sys.stdout.flush();
                    
                
                    
                    
                #mask_nEle = data_i[:]['nElectron'] >= 2
                if printLevel >= 12:
                    print("mask_nEle: {}".format(list(mask_nEle)))
                    print("np.sum(mask_nEle): {}".format(np.sum(mask_nEle)))
                    print("\n\nmask: {}, \n\ndata_i['nElectron'][mask]: {}".format(mask_nEle, data_i[mask_nEle]['nElectron'] ))
                    print("\n\ndata_i[mask]['Electron_pt']: {}".format(data_i[mask_nEle]['Electron_pt']))

                if ak.sum(mask_nEle) == 0: continue   
                data_i = data_i[mask_nEle]
                
                if printLevel >= 12: 
                    #print("\ndata_i selected 2 ({}) ({}):  {}".format(type(data_i), len(data_i), ak.to_list(data_i))); sys.stdout.flush();
                    print("\ndata_i selected 2 ({}) ({}):  {}".format(type(data_i), len(data_i), ak.to_list(data_i))); sys.stdout.flush();
                    #print("len(data_i): {}, len(data_all): {}".format(len(data_i), len(data_all)))

                if printLevel >= 15:
                    print("ak.to_list(data_i[:]['Muon_pt']): {}".format(ak.to_list(data_i[:]['Muon_pt'])))
                    print("data_i[:]['Muon_pt'][1]): {}".format(data_i[:]['Muon_pt'][1]))
                    print("data_i[:]['Muon_pt'][1]): {}".format((data_i[:]['Muon_pt'])[:][1]))
                    print("data_i[0]['Muon_pt']: {}, \t  data_i[0]['Muon_pt'][1]): {}".format(data_i[0]['Muon_pt'], data_i[0]['Muon_pt'][1]))
                    print("data_i[1]['Muon_pt']: {}, \t  data_i[1]['Muon_pt'][1]): {}".format(data_i[1]['Muon_pt'], data_i[1]['Muon_pt'][1]))

                    print("data_i[:, 'Muon_pt']: {}".format(data_i[:, 'Muon_pt']))
                    print("data_i[:, 'Muon_pt', 1]: {}".format(data_i[:, 'Muon_pt', 1]))

                
                    
                #mask_MuPt = np.vectorize(cut_ObjectPt)(data_i[:, 'Muon_pt'], PtThrsh_Lead=25, PtThrsh_Sublead=20, PtThrsh_Third=15, PtThrsh_Fourth=10)
                #mask_MuPt = np.vectorize(cut_ObjectPt_1)(data_i[:, 'Muon_pt'], [25, 20, 15, 10])
                mask_MuPt = ((data_i[:, 'Muon_pt', 0] > 25) & 
                             (data_i[:, 'Muon_pt', 1] > 15) & 
                             (data_i[:, 'Muon_pt', 2] > 10) )
                #print("mask_MuPt ({}) ({}): {}  \t\t sum: {}".format(type(mask_MuPt), len(mask_MuPt), ak.to_list(mask_MuPt), ak.sum(mask_MuPt)))
                if ak.sum(mask_MuPt) == 0: continue   
                data_i = data_i[mask_MuPt]
                if printLevel >= 12:
                    print("\ndata_i selected 3 ({}) ({}):  {}".format(type(data_i), len(data_i), ak.to_list(data_i))); sys.stdout.flush();
                    print("data_i[:, 'Muon_eta', :) ({}): {}".format(type(data_i[:, 'Muon_eta', :]), data_i[:, 'Muon_eta', :]))
                    print("abs(data_i[:, 'Muon_eta', :) ): {}".format(abs(data_i[:, 'Muon_eta', :])))

                    
                #mask_MuEta = abs(data_i[:, 'Muon_eta', :]) < 2.5
                mask_MuEta = ((abs(data_i[:, 'Muon_eta', 0]) < 2.5) & 
                              (abs(data_i[:, 'Muon_eta', 1]) < 2.5) & 
                              (abs(data_i[:, 'Muon_eta', 2]) < 2.5) )
                if ak.sum(mask_MuEta) == 0: continue   
                data_i = data_i[mask_MuEta]
                if printLevel >= 12:
                    print("\ndata_i selected 4 ({}) ({}):  {}".format(type(data_i), len(data_i), ak.to_list(data_i))); sys.stdout.flush();


                mask_ElePt = ((data_i[:, 'Electron_pt', 0] > 20) & 
                              (data_i[:, 'Electron_pt', 1] > 15)  )
                if ak.sum(mask_ElePt) == 0: continue   
                data_i = data_i[mask_ElePt]

                
                mask_EleEta = ((abs(data_i[:, 'Electron_eta', 0]) < 2.4) & 
                               (abs(data_i[:, 'Electron_eta', 1]) < 2.4) )
                if ak.sum(mask_EleEta) == 0: continue   
                data_i = data_i[mask_EleEta]

                
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

    
    if printLevel >= 0  and len(data_all) > 0:                
        print("\ndata_all final  ({}) ({}):  {}".format(type(data_all), len(data_all), ak.to_list(data_all))); sys.stdout.flush();

                
    if sOutFile is not None  and len(data_all) > 0:
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
            
            hMuonPt_thridlead = hist.Hist.new.Reg(100, 0, 100, name="Third-leading muon pT").Weight()
            hMuonPt_thridlead.fill(ak.to_numpy( data_all[:, 'Muon_pt', 2] ))
            fOut['evt/hThirdleadingMuonPt'] = hMuonPt_thridlead


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


            print("Wrote to sOutFile {}".format(sOutFile))

    return


                
def data_read_and_select_wPandas_v0(sInputFiles, branchesToRead, sOutFile=None):
    data_all = None
    data_all_info = None
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
            ]

            htoaa_nanoAODBranchesToRead_electron = [
                'nElectron',
                'Electron_pt',
                'Electron_eta',
                'Electron_phi'
            ]
            
            #tracemalloc.start()
            #counts = Counter()
            for iEvtStart in range(0, 500, nEventToReadInBatch):
                data_subset = None

                data_subset_info = tree.arrays(
                    htoaa_nanoAODBranchesToRead_info,
                    entry_start=iEvtStart, entry_stop=(iEvtStart+nEventToReadInBatch),
                    library='pd'
                )
                if (data_subset_info is not None) and (len(data_subset_info.index) > 0):
                    if data_all_info is None:
                        data_all_info = data_subset_info
                    else:
                        data_all_info = pd.concat([data_all_info, data_subset_info])
                        

                for iTreeBranchSet, treeBranches_list in enumerate([
                        htoaa_nanoAODBranchesToRead_info + htoaa_nanoAODBranchesToRead_muon,
                        htoaa_nanoAODBranchesToRead_info + htoaa_nanoAODBranchesToRead_electron
                        #htoaa_nanoAODBranchesToRead_info + htoaa_nanoAODBranchesToRead_muon + htoaa_nanoAODBranchesToRead_electron
                        #htoaa_nanoAODBranchesToRead_info
                        #htoaa_nanoAODBranchesToRead_info, htoaa_nanoAODBranchesToRead_muon, htoaa_nanoAODBranchesToRead_electron                        
                ]):
                    #treeBranches_list = set(treeBranches_list)
                    if printLevel >= 10:
                        print("\n\ntreeBranches_list: {}".format(treeBranches_list))
                        
                    data_subset_byTrBrch = tree.arrays(
                        treeBranches_list,
                        entry_start=iEvtStart, entry_stop=(iEvtStart+nEventToReadInBatch),
                        library='pd'
                    )
                    '''
                    print("data_subset_byTrBrch 0 {}: \n{}".format(type(data_subset_byTrBrch), data_subset_byTrBrch)); sys.stdout.flush();
                    print("data_subset_byTrBrch.index ({}): {}".format(type(data_subset_byTrBrch.index), data_subset_byTrBrch.index)); sys.stdout.flush();
                    index_tmp = list(data_subset_byTrBrch.index)
                    print("index_tmp ({}): {}".format(type(index_tmp), index_tmp))
                    index_multi = pd.MultiIndex.from_product([list(data_subset_byTrBrch.index), [0]], names=["entry", "subentry"])
                    print("index_multi ({}): {}".format(type(index_multi), index_multi))
                    #data_subset_byTrBrch = pd.DataFrame( data_subset_byTrBrch, index=pd.MultiIndex(levels=[list(data_subset_byTrBrch.index), []], names=['entry', 'subentry']) )
                    #data_subset_byTrBrch = pd.DataFrame( data_subset_byTrBrch, index=[() for   ], names=['entry', 'subentry']) )
                    data_subset_byTrBrch = pd.DataFrame( data_subset_byTrBrch[treeBranches_list], index=index_multi )
                    print("data_subset_byTrBrch 1 {}: \n{}".format(type(data_subset_byTrBrch), data_subset_byTrBrch)); sys.stdout.flush();
                    '''
                    if printLevel >= 10:
                        #print("\ndata_subset_byTrBrch.index ({}): {}".format(type(data_subset_byTrBrch.index), data_subset_byTrBrch.index))
                        #print("data_subset_byTrBrch.columns ({}): {}".format(type(data_subset_byTrBrch.columns), data_subset_byTrBrch.columns))
                        print("data_subset_byTrBrch ({}): \n{}".format(type(data_subset_byTrBrch), data_subset_byTrBrch))

                    if len(data_subset_byTrBrch.index) > 0:
                        if data_subset is None:
                            data_subset = data_subset_byTrBrch 
                        else:
                            #data = data.join(data_part)
                            #data_subset = data_subset.join(data_subset_byTrBrch, how='outer')
                            data_subset = data_subset.join(data_subset_byTrBrch, how='outer', rsuffix='_%d' % (iTreeBranchSet))                            
                        if printLevel >= 10:
                            print("\n\ndata_subset_i: \n{}".format(data_subset))

                if printLevel >= 10:
                    print("\n\ndata_subset: \n{}".format(data_subset))

                            
                if (data_subset is not None) and (len(data_subset.index) > 0):
                    if data_all is None:
                        data_all = data_subset
                    else:
                        data_all = pd.concat([data_all, data_subset])
                        
                if printLevel >= 10:
                    print("\n\ndata_all_i: \n{}".format(data_all))

                        
                        
            if printLevel >= 9:
                print("\n\ndata_all: \n{}".format(data_all))

            mask_nMu = data_all.nMuon >= 2
            data_all = data_all.loc[ mask_nMu ]
            if printLevel >= 9:
                print("mask_nMu: {}".format(mask_nMu))
                print("\n\n data_all.loc[ data_all.nMuon > 2 ]**: \n{}".format(data_all))
                print("data_all[(:, 0)].index: {}.  len(data_all[(:, 0)].index): {} ".format(data_all.loc[(slice(None), 0), :].index, len(data_all.loc[(slice(None), 0), :].index)))

            #mask_MuPt = data_all.loc[(slice(None), 0), 'Muon_pt'] > 20
            mask_MuPt = data_all.loc[ pd.IndexSlice[:, 0], 'Muon_pt'] > 20
            mask_MuPt_1 = mask_MuPt[mask_MuPt].index
            mask_MuPt_2 = mask_MuPt_1.get_level_values(0) # mask_MuPt_1.droplevel(level=1)
            if printLevel >= 9:
                print("\nmask_MuPt ({}): {}".format(type(mask_MuPt), mask_MuPt))
                print("\n mask_MuPt.index: {}".format(mask_MuPt.index))
                print("\n mask_MuPt.index_1: {}".format(mask_MuPt[mask_MuPt].index))
                
                print("mask_MuPt_1 ({}): {}".format(type(mask_MuPt_1), mask_MuPt_1))
                
                print("mask_MuPt_2 ({}): {}".format(type(mask_MuPt_2), mask_MuPt_2))
                #mask_MuPt_3 = list(mask_MuPt_2)
                #print("mask_MuPt_3 ({}): {}".format(type(mask_MuPt_3), mask_MuPt_3))
            #data_all = data_all.loc[(slice(mask_MuPt_3), slice(None)), :]
            #data_all = data_all.loc[ pd.IndexSlice[mask_MuPt, :], pd.IndexSlice[:] ]
            #data_all = data_all.loc[ pd.IndexSlice[mask_MuPt_1, :], pd.IndexSlice[:] ]
            #data_all = data_all.loc[ pd.IndexSlice[mask_MuPt_1, :], : ]
            data_all = data_all.loc[ pd.IndexSlice[mask_MuPt_2, :], : ]
            if printLevel >= 9:
                print("data_all ({}): \n{}".format(len(data_all.loc[pd.IndexSlice[:, 0], :].index), data_all))
                #print("\n\n data_allmask_MuPt.loc[ data_all.nMuon > 2 ]**: \n{}".format(data_all))


            mask_1 = data_all.loc[pd.IndexSlice[:, 0], :].index
            mask_2 = data_all.loc[pd.IndexSlice[:, 0], :].index.get_level_values(level=0)
            if printLevel >= 9:
                print("mask_1 ({}): {}".format(type(mask_1), mask_1))
                print("mask_2 ({}): {}".format(type(mask_2), mask_2))


            mask_info_1 = data_all_info.index
            if printLevel >= 9:
                print("mask_info_1 ({}): {}".format(type(mask_info_1), mask_info_1))
                print("data_all_info ({}): {}".format(len(data_all_info.index), data_all_info))

            data_all_info = data_all_info.loc[pd.IndexSlice[mask_2], :]
            mask_info_1 = data_all_info.index
            if printLevel >= 9:
                print("mask_info_1 after selection ({}): {}".format(type(mask_info_1), mask_info_1))
                print("data_all_info after selection ({}): {}".format(len(data_all_info.index), data_all_info))
           
            
            #data_all.loc[]




def data_read_and_select_wPandas(sInputFiles, branchesToRead, sOutFile=None):
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
            ]

            htoaa_nanoAODBranchesToRead_electron = [
                'nElectron',
                'Electron_pt',
                'Electron_eta',
                'Electron_phi'
            ]
            
            #tracemalloc.start()
            #counts = Counter()
            for iEvtStart in range(0, nEventsToAnalyze, nEventToReadInBatch):
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
                mask_pass_nMu_cut = cutDf_ObjectMultiplicity(df=collection_mu.df, objectName='nMuon', nObjects_min=2)
                mask_pass_all_cut = mask_pass_nMu_cut


                #mask_pass_nEle_cut = collection_ele.df['nElectron'] >= 2
                mask_pass_nEle_cut = cutDf_ObjectMultiplicity(df=collection_ele.df, objectName='nElectron', nObjects_min=2)
                mask_pass_all_cut = (mask_pass_all_cut & mask_pass_nEle_cut)                
                if printLevel >= 12:
                    print("mask_pass_nMu_cut  ({}): {}".format(len(mask_pass_nMu_cut.index), mask_pass_nMu_cut))
                    print("mask_pass_nEle_cut ({}): {}".format(len(mask_pass_nEle_cut.index), mask_pass_nEle_cut))
                    print("mask_pass_all_cut  ({}): {}".format(len(mask_pass_all_cut.index), mask_pass_all_cut))

                    
                events.dropEventsFailingMask(mask_pass_cut = mask_pass_all_cut)


                if printLevel >= 10:
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


                #mask_pass_MuPt_cut = cutDf_ObjectPt(df=collection_mu.df, objectName='Muon_pt', PtThrshs=[20.0, 3.5])
                #mask_pass_MuPt_cut = cutDf_ObjectPt(df=collection_mu.df, objectName='Muon_pt', PtThrshs=[20.0])
                mask_pass_MuPt_cut = cutDf_ObjectPt(df=collection_mu.df, objectName='Muon_pt', PtThrshs=[5.0, 3.0])
                mask_pass_all_cut = (mask_pass_all_cut & mask_pass_MuPt_cut)
                 
                if printLevel >= 10:
                    print("mask_pass_MuPt_cut  ({}): {}".format(len(mask_pass_MuPt_cut.index), mask_pass_MuPt_cut))
                    


                    
                    
                if printLevel >= 10:
                    print("mask_pass_all_cut ({}): {}".format(len(mask_pass_all_cut.index), mask_pass_all_cut.to_string()))

                    
                if getDataframeNEvents(mask_pass_all_cut) == 0: continue
                

                events.dropEventsFailingMask(mask_pass_cut = mask_pass_all_cut)

                
                
                print("\n Dayaset selected untilnow: :: Events: \n")
                events.print()

                


    print("\n\n\n Final dataset:: Events: \n")
    events.print()


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
        
        if printLevel >= 10:
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
        #indexLevel0_allCollections = list(set(indexLevel0_allCollections)) # remove duplicate

        #indexLevel0_eventsPassingCut = set( getDataframeIndexLevel0(mask_pass_cut[mask_pass_cut]) )
        indexLevel0_eventsPassingCut = mask_pass_cut.groupby(axis=0, level=0).all() # Boolean pd.Series with 'index=index level0 of mask_pass_cut' seting True when all subentries in mask_pass_cut are true
        indexLevel0_eventsPassingCut = set( indexLevel0_eventsPassingCut[indexLevel0_eventsPassingCut].index ) # get entry when all subentries in 'mask_pass_cut' are true
        indexLevel0_eventsToDrop =  indexLevel0_allCollections - indexLevel0_eventsPassingCut 
        
        if printLevel >= 10:
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



def cutDf_ObjectPt_1(df, objectName, PtThrshs=[]): 
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

    mask = None
    #mask_tmp = df['nMuon'] != df['nMuon']
    mask_tmp = df[objectName] != df[objectName]
    #mask_tmp = df[objectName].astype(int) != df[objectName].astype(int)
    if printLevel >= 10:
        print("cutDf_ObjectPt():: mask_tmp   ({}): {} \n{}\n\n".format(type(mask_tmp), getDataframeNEvents(mask_tmp),   mask_tmp.tolist(), mask_tmp.to_string()))
        print("cutDf_ObjectPt():: mask_tmp.index: {}\n".format(mask_tmp.index))

    mask = mask_tmp
    for iParticle, PtTrsh in enumerate(PtThrshs):
        mask_i = df.loc[pd.IndexSlice[:, iParticle], objectName] > PtTrsh
        #if mask is None: mask = mask_i
        #else:            mask = (mask & mask_i)
        #else:            mask = mask.join(mask_i, how='outer')

        mask = (mask | mask_i)
        
        
        if printLevel >= 10:
            print("cutDf_ObjectPt():: iParticle {}, PtTrsh {}".format(iParticle, PtTrsh))
            print("cutDf_ObjectPt():: mask   ({}): {} \n{}".format(getDataframeNEvents(mask),   mask.tolist(), mask.to_string()))
            print("cutDf_ObjectPt():: mask_i ({}): {} \n{}".format(len(mask_i.tolist()), mask_i.tolist(), mask_i.to_string()))

    return mask;



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

    mask = None
    #mask_tmp = df['nMuon'] != df['nMuon']
    mask_tmp = df[objectName] == df[objectName]
    #mask_tmp = df[objectName].astype(int) != df[objectName].astype(int)
    if printLevel >= 10:
        print("cutDf_ObjectPt():: mask_tmp   ({}): {} \n{}\n\n".format(type(mask_tmp), getDataframeNEvents(mask_tmp),   mask_tmp.tolist(), mask_tmp.to_string()))
        print("cutDf_ObjectPt():: mask_tmp.index: {}\n".format(mask_tmp.index))

    mask = mask_tmp
    for iParticle, PtTrsh in enumerate(PtThrshs):
        mask_i = df.loc[pd.IndexSlice[:, iParticle], objectName] > PtTrsh
        #if mask is None: mask = mask_i
        #else:            mask = (mask & mask_i)
        #else:            mask = mask.join(mask_i, how='outer')

        #mask = ~ (mask &  mask_i)
        mask = (mask &  mask_i)
        
        
        if printLevel >= 10:
            print("cutDf_ObjectPt():: iParticle {}, PtTrsh {}".format(iParticle, PtTrsh))
            print("cutDf_ObjectPt():: mask   ({}): {} \n{}".format(getDataframeNEvents(mask),   mask.tolist(), mask.to_string()))
            print("cutDf_ObjectPt():: mask_i ({}): {} \n{}".format(len(mask_i.tolist()), mask_i.tolist(), mask_i.to_string()))

    return mask;



def cutDf_ObjectEta(df, objectName, EtaThrsh, nObjects):
    '''
    Check Objects abs(Eta) is greater than thresholds set in EtaThrshs list.

    Return:
        True: All objects' Eta is about respective threshold
        False: Else false
    '''    
    condition = True
    for iObject in range(nObjects):
        if abs(objects_Eta[iObject]) > EtaThrsh:
            condition = False
            break
    return condition



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


    #tracemalloc.start()
    #counts = Counter()

    #data_read_and_select_wNumpy(sInputFiles, branchesToRead)

    #data_read_and_select_wAkward(sInputFiles, branchesToRead, sOutputFile)
    data_read_and_select_wPandas(sInputFiles, branchesToRead, sOutputFile)

    #snapshot = tracemalloc.take_snapshot()
    #display_top(snapshot, key_type='lineno', limit=10)

    print("Done")
    exit(0)

    
    ### Event loop starts here --------------------------------------------------------------------------------------------------
    data_all = pd.DataFrame()
    for sInputFile in sInputFiles:
        if printLevel >= 10: print("sInputFile: {}".format(sInputFile))
        with uproot.open(sInputFile) as fInputFile:            
            tree = fInputFile['Events']
            #tree = fInputFile.get('Events')
            #if printLevel >= 1: print("file: {}, nEvents: {}".format(sInputFile, uproot.numentries(sInputFile, 'Events')))
            if printLevel >= 1: print("file: {}, nEvents: {}".format(sInputFile, uproot.num_entries('%s:Events' % sInputFile)))
            print("tree.num_entries(): {}".format(tree.num_entries))

            #for data_i in tree.iterate(branchesToRead,
            #                           flatten=False, outputtype=pd.DataFrame,
            #                           entrysteps=nEventToReadInBatch):
            for data_i in tree.iterate(branchesToRead,
                                       library="pd", how=dict,
                                       #library="np",
                                       #library="ak",
                                       #max_num_elements=4000, 
                                       step_size=1000):

                
                print("\ndata_i 0 ({}): {}".format(type(data_i), data_i)); sys.stdout.flush();
                #print("\ndata_i 0 ({}): {}".format(type(data_i), data_i.to_list())); sys.stdout.flush();
                #print("\n\ndata_i 0 {},  \n{}".format(data_i.shape, data_i.columns)); sys.stdout.flush();
                #print("\n\ndta_i.index ({}): {}".format(type(data_i.index), data_i.index))


                data_i = data_i.loc[data_i.nMuon >= 3]
                print("\n\ndata_i selected 0: {}".format(data_i)); sys.stdout.flush();

                data_i = data_i[data_i.nElectron >= 2]
                print("\n\ndata_i selected 1: {}".format(data_i)); sys.stdout.flush();

                continue
            
                
                #fail = data_i[ ~ np.vectorize(cut_ObjectMultiplicity)(data_i.nMuon, nObjects_min=2, nObjects_max=99) ].index
                fail = data_i[ ~ np.vectorize(cut_ObjectMultiplicity)(data_i.nMuon, nObjects_min=3) ].index
                #fail = data_i[ ~ np.vectorize(cut_ObjectMultiplicity)(data_i.nMuon,  nObjects_max=1) ].index
                #data_i.drop(index=fail, inplace=True)
                data_i.drop(fail, inplace=True)
                print("\n\ndata_i selected 0 {},  \n{}".format(data_i.shape, data_i.columns)); sys.stdout.flush();
                print("\n\ndata_i selected 0: {}".format(data_i))


                fail = data_i[ ~ np.vectorize(cut_ObjectMultiplicity)(data_i.nElectron, nObjects_min=2) ].index
                data_i.drop(fail, inplace=True)
                print("\n\ndata_i selected 0p1 {},  \n{}".format(data_i.shape, data_i.columns)); sys.stdout.flush();
                print("\n\ndata_i selected 0p1: {}".format(data_i))
                
                

                print("data_i.Muon_pt: {}".format(data_i.Muon_pt))
                #print("".format())

                continue
                
                #PtThrshs=[[25, 15, 10]]
                #print("PtThrshs_0 ({}): {}".format(type(PtThrshs), PtThrshs))
                #fail = data_i[ ~ np.vectorize(cut_ObjectPt)(data_i.Muon_pt, PtThrshs=PtThrshs)].index
                fail = data_i[ ~ np.vectorize(cut_ObjectPt)(data_i.Muon_pt, PtThrsh_Lead=25, PtThrsh_Sublead=15, PtThrsh_Third=10)].index
                data_i.drop(fail, inplace=True)
                print("\n\ndata_i selected 1 {},  \n{}".format(data_i.shape, data_i.columns)); sys.stdout.flush();
                print("\n\ndata_i selected 1: {}".format(data_i))

                fail = data_i[ ~ np.vectorize(cut_ObjectEta)(data_i.Muon_eta, EtaThrsh=2.3, nObjects=3) ].index
                data_i.drop(index=fail, inplace=True)
                
 
                data_i['mll_tmp'] = np.vectorize(calc_mll)(data_i.Muon_pt, data_i.Muon_eta, data_i.Muon_phi)
                print("\n\ndata_i selected 2 {},  \n{}".format(data_i.shape, data_i.columns)); sys.stdout.flush();
                print("\n\ndata_i selected 2: {}".format(data_i))
                
                if data_i.empty: continue
                data_all = data_all.append( data_i )


    exit(0)
    fIn = uproot.open(sInputFiles[0])
    tree = fInputFile['Events']
    
    print("data_all.Muon_pt: {}".format(data_all.Muon_pt))
    print("ak.flatten(tree['Muon_pt']]: {}".format(ak.flatten(tree['Muon_pt'])))
    #print("ak.flatten(data_all.Muon_pt): {}".format(ak.flatten(data_all.Muon_pt)))
    #h_Muon_pt = np.histogram(ak.flatten(data_all.Muon_pt), bins=40, range=(0, 100))
    h_Muon_pt = np.histogram(ak.flatten(tree["Muon_pt"]), bins=40, range=(0, 100))
    #h_Muon_pt = np.histogram(data_all.nMuon, bins=20, range=(0, 20))


    with uproot.recreate(sOutputFile) as fOutputFile:
        fOutputFile["muon_pt"] = h_Muon_pt
    
