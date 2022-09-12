#htoaa analysis main code

import os
import sys
import json
from collections import OrderedDict as OD
import pandas as pd
import uproot
#import uproot3 as uproot

from htoaa_Settings import *
from htoaa_CommonTools import GetDictFromJsonFile, calculate_lumiScale
from htoaa_NanoAODBranches import htoaa_nanoAODBranchesToRead

printLevel = 10

def read_data(sInputFiles, branchesToRead):
    data_all = pd.DataFrame()
    for sInputFile in sInputFiles:
        if printLevel >= 10: print("sInputFile: {}".format(sInputFile))
        with uproot.open(sInputFile) as fInputFile:
            #tree = fInputFile['Events']
            tree = fInputFile.get('Events')

            print("nEvents: {}".format(uproot.numentries(sInputFile, 'Events')))
            
            if printLevel >= 10:
                print("fInputFile.keys(): {}, \nfInputFile.classnames(): {}".format(fInputFile.keys(), fInputFile.classnames()))
                print("type(fInputFile['Events']: {}".format(type(fInputFile['Events'])))
                print("type(fInputFile['Events']: {}".format(type(fInputFile[b'Events;1'])))
                print("type(fInputFile.get('Events'): {}".format(type(fInputFile.get('Events'))))
                print("type(tree): {}".format(type(tree)))

                print("tree.keys(): {}".format(tree.keys()))

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

                #print("mask: {}".format(data_i.nMuon > 0))
                mask = data_i[b'nMuon'] > 0
                print("mask{} : {}".format(type(mask), mask))
                print("mask.sum(): {}".format(mask.sum()))
                print("mask.sum(): {}".format((data_i[b'nMuon'] > 0).sum()))
                
                #data_i.drop(data_i[b'nMuon'] <= 1, inplace=True)
                print("data_i.loc[mask]): {}".format(data_i.loc[mask]))
                print("data_i.loc[data_i[b'nMuon'] > 0]): {}".format(data_i.loc[data_i[b'nMuon'] > 0]))
                print("data_i[mask]: {}".format(data_i[mask]))
                print("data_i[mask].index: {}".format(data_i[mask].index))
                #data_i.drop(index=data_i[~mask].index, inplace=True)
                #data_i.drop(index=data_i[data_i[b'nMuon'] == 0].index, inplace=True)
                data_i.drop(index=data_i.loc[data_i[b'nMuon'] == 0].index, inplace=True)
                print("data_i: {}".format(data_i))

                data_all = data_all.append( data_i )
                print("data_all: {}".format(data_all))

                
            
            
    return;
    



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

    read_data(sInputFiles=sInputFiles, branchesToRead=htoaa_nanoAODBranchesToRead)
