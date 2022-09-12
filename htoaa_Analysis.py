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

from htoaa_Settings import *
from htoaa_CommonTools import GetDictFromJsonFile, DfColLabel_convert_bytes_to_string, calculate_lumiScale
from htoaa_NanoAODBranches import htoaa_nanoAODBranchesToRead

printLevel = 10
pd.set_option('display.max_columns', None)

def read_data(sInputFiles, branchesToRead):
    data_all = pd.DataFrame()
    for sInputFile in sInputFiles:
        if printLevel >= 10: print("sInputFile: {}".format(sInputFile))
        with uproot.open(sInputFile) as fInputFile:
            #tree = fInputFile['Events']
            tree = fInputFile.get('Events')

            print("nEvents: {}".format(uproot.numentries(sInputFile, 'Events')))
            
            if printLevel >= 10:
                '''
                print("fInputFile.keys(): {}, \nfInputFile.classnames(): {}".format(fInputFile.keys(), fInputFile.classnames()))
                print("type(fInputFile['Events']: {}".format(type(fInputFile['Events'])))
                print("type(fInputFile['Events']: {}".format(type(fInputFile[b'Events;1'])))
                print("type(fInputFile.get('Events'): {}".format(type(fInputFile.get('Events'))))
                print("type(tree): {}".format(type(tree)))

                print("tree.keys(): {}".format(tree.keys()))
                '''
                #data_i = tree.arrays(branchesToRead, library="pd")
                #data_i = tree.arrays(branchesToRead)
                #print("type(data_i): {}".format(type(data_i)))
                #data_i = pd.DataFrame( tree.arrays(branchesToRead, entry_stop=1000) )                
                data_i = pd.DataFrame( tree.arrays(branchesToRead) )
                #data_i.set_index([b'run', b'luminosityBlock', b'event'])
                #data_i = pd.DataFrame( tree.arrays(branchesToRead), columns=branchesToRead )

                
                '''
                for br_name in branchesToRead:
                    print("br_name: {}".format(br_name))
                    #data_i[br_name] = pd.DataFrame( tree.arrays( br_name ) )
                    data_i[br_name] = pd.Series( tree.arrays( br_name ) )
                '''
                print("type(data_i): {}".format(type(data_i)))
                print("data_i.columns: {}".format(data_i.columns))
                print("data_i.shape: {}".format(data_i.shape))
                #print("data_i.describe: {}".format(data_i.describe()))
                print("data_i: {}".format(data_i))

                data_i = DfColLabel_convert_bytes_to_string(data_i)

                print("data_i.columns: {}".format(data_i.columns))
                for col in data_i.columns:
                    print("col {}, type(col) {}".format(col, type(col)))

                '''
                rename_cols = { bytes(br_name, 'utf-8') : br_name for br_name in branchesToRead }
                print("rename_cols: {}".format(rename_cols))
                data_i.rename(columns=rename_cols, inplace=True)
                '''
                '''
                for br_name in branchesToRead:
                    #data_i.rename(columns={b'%s' % (br_name): '%s' % (br_name)}, inplace=True)
                    data_i.rename(columns={bytes(br_name, 'utf-8'): '%s' % (br_name)}, inplace=True)
                '''
                #print("data_i after renaming columns: {}".format(data_i))

                '''
                #print("mask: {}".format(data_i.nMuon > 0))
                mask = data_i['nMuon'] > 0
                print("mask{} : {}".format(type(mask), mask))
                print("mask.sum(): {}".format(mask.sum()))
                print("mask.sum(): {}".format((data_i['nMuon'] > 0).sum()))
                
                #data_i.drop(data_i[b'nMuon'] <= 1, inplace=True)
                print("data_i.loc[mask]): {}".format(data_i.loc[mask]))
                print("data_i.loc[data_i[b'nMuon'] > 0]): {}".format(data_i.loc[data_i['nMuon'] > 0]))
                print("data_i[mask]: {}".format(data_i[mask]))
                print("data_i[mask].index: {}".format(data_i[mask].index))
                #data_i.drop(index=data_i[~mask].index, inplace=True)
                #data_i.drop(index=data_i[data_i[b'nMuon'] == 0].index, inplace=True)
                '''
                data_i.drop(index=data_i.loc[data_i['nMuon'] <= 1].index, inplace=True)
                #print("data_i: {}".format(data_i))

                data_all = data_all.append( data_i )
                #print("data_all: {}".format(data_all))

                #print("\n\n\ndata_i.index ({}): {}".format(len(data_i), data_i.index.to_list()))
                #print("\n\ndata_all.index ({}): {}\n\n".format(len(data_all), data_all.index.to_list()))


    '''
    print("data_all['Muon_pt'][:]: {}".format(data_all['Muon_pt'][:]))
    print("\ndata_all['Muon_pt'][:][0]: {}".format(data_all['Muon_pt'][data_all.index][0]))
    print("\ndata_all['Muon_pt'][:][1]: {}".format(data_all['Muon_pt'][:][1]))
    '''

    print("data_all.head(): {}".format(data_all.head()))

    data_all['mll'] = np.vectorize(calc_mll)(data_all.Muon_pt, data_all.Muon_eta, data_all.Muon_phi)

    print("data_all.head(): {}".format(data_all.head()))
    
    #for i in range(2):
    #    print("data_all['Muon_pt'][{}]: {}".format(i, data_all[:, 'Muon_pt'][i]))
    
            
    return;


def read_data1(sInputFiles, branchesToRead):
    data_all = pd.DataFrame()
    for sInputFile in sInputFiles:
        if printLevel >= 10: print("sInputFile: {}".format(sInputFile))
        with uproot.open(sInputFile) as fInputFile:            
            #tree = fInputFile['Events']
            tree = fInputFile.get('Events')

            print("nEvents: {}".format(uproot.numentries(sInputFile, 'Events')))

            iteration = 0
            for data_i in tree.iterate(branchesToRead,
                                       flatten=False, outputtype=pd.DataFrame,
                                       entrysteps=2500000):
                print("data_i {} ({}): {}".format(iteration, type(data_i), data_i))
                iteration += 1


    

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
    
    sInputFiles = config["inputFiles"]
    era = config['era']
    luminosity = Luminosities[era][0]
    sample_crossSection = config["crossSection"]
    sample_nEvents = config["nEvents"]
    sample_sumEvents = config["sumEvents"] if config["sumEvents"] != -1 else sample_nEvents
    if sample_sumEvents == -1: sample_sumEvents = 1 # Case when sumEvents is not calculated
    lumiScale = calculate_lumiScale(luminosity=luminosity, crossSection=sample_crossSection, sumEvents=sample_sumEvents)

    # read_data(sInputFiles=sInputFiles, branchesToRead=htoaa_nanoAODBranchesToRead)
    read_data1(sInputFiles=sInputFiles, branchesToRead=htoaa_nanoAODBranchesToRead)
