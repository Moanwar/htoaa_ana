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
from htoaa_NanoAODBranches import htoaa_nanoAODBranchesToRead
from htoaa_CommonTools import GetDictFromJsonFile, DfColLabel_convert_bytes_to_string, calculate_lumiScale
from htoaa_CommonTools import cut_ObjectMultiplicity, cut_ObjectPt, cut_ObjectEta


printLevel = 10
nEventToReadInBatch = 2500000
pd.set_option('display.max_columns', None)



def data_read_and_select(sInputFiles, branchesToRead):
    data_all = pd.DataFrame()
    for sInputFile in sInputFiles:
        if printLevel >= 10: print("sInputFile: {}".format(sInputFile))
        with uproot.open(sInputFile) as fInputFile:            
            tree = fInputFile['Events']
            #tree = fInputFile.get('Events')
            if printLevel >= 1: print("file: {}, nEvents: {}".format(sInputFile, uproot.numentries(sInputFile, 'Events')))

            for data_i in tree.iterate(branchesToRead,
                                       flatten=False, outputtype=pd.DataFrame,
                                       entrysteps=nEventToReadInBatch):

                ## select data... Put data selection cuts here --------------------------------------------------

                #fail = data_i[ ~ np.vectorize(cut_ObjectMultiplicity)(data_i.nMuon, nObjects_min=2, nObjects_max=99) ].index
                fail = data_i[ ~ np.vectorize(cut_ObjectMultiplicity)(data_i.nMuon, nObjects_min=3) ].index
                #fail = data_i[ ~ np.vectorize(cut_ObjectMultiplicity)(data_i.nMuon,  nObjects_max=1) ].index
                #data_i.drop(index=fail, inplace=True)
                data_i.drop(fail, inplace=True)
                print("\n\ndata_i selected 0: {}".format(data_i))

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


def plot_data(data):
    
    


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
    data_selected = data_read_and_select(sInputFiles=sInputFiles, branchesToRead=htoaa_nanoAODBranchesToRead)
    plot_data( data_selected )
