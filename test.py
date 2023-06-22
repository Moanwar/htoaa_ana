import os
import sys
from datetime import datetime
import time

print(f"htoaa_Analysis_GGFMode:: here1 {datetime.now() = }"); sys.stdout.flush();
import subprocess
import json
from urllib.request import urlopen
import numpy as np 

from coffea.lumi_tools import LumiMask

'''
with urlopen("https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt") as fDataGoldenJSON:
    dataLSSelGoldenJSON = json.load(fDataGoldenJSON)    
    #print(f"0 : {dataLSSelGoldenJSON = }")
    dataLSSelGoldenJSON = {int(k): v for k, v in dataLSSelGoldenJSON.items()}
    #print(f"1 : {dataLSSelGoldenJSON = }")
'''
from htoaa_CommonTools import (
    selectRunLuminosityBlock, selectRunLuminosityBlock_ApprochEventBase
)


#runNumber_list       = np.array([315322, 315322, 315323, 315324, 315357, 315357, 315357, 315357, 315357])
#luminosityBlock_list = np.array([    25,   1339,     45,   1909,    700,    750,    800,    775,     40])
runNumber_list       = np.random.randint(315257, 325172, 1000000)
luminosityBlock_list = np.random.randint(1, 2400, 1000000)

#print(f" { = } ")
print(f" {runNumber_list = } ")
print(f" {luminosityBlock_list = } "); sys.stdout.flush();

golden_json_path = "data/goldenJsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"

print(f"htoaa_Analysis_GGFMode:: here2 {datetime.now() = }"); sys.stdout.flush();
start_time = time.time()
mask_run_ls = selectRunLuminosityBlock_ApprochEventBase(
    #dataLSSelGoldenJSON = dataLSSelGoldenJSON, 
    golden_json_path = golden_json_path, 
    runNumber_list = runNumber_list, 
    luminosityBlock_list = luminosityBlock_list)

print("selectRunLuminosityBlock_ApprochEventBase: --- %s seconds ---" % (time.time() - start_time)); sys.stdout.flush();
print(f"htoaa_Analysis_GGFMode:: here3 {datetime.now() = }")


start_time = time.time()
lumi_mask = LumiMask(golden_json_path)(runNumber_list, luminosityBlock_list)
#print(f"lumi_mask ({type(lumi_mask)}) ({len(lumi_mask)}): {lumi_mask}")
print("lumiMask: --- %s seconds ---" % (time.time() - start_time)); sys.stdout.flush();
print(f"htoaa_Analysis_GGFMode:: here4 {datetime.now() = }")

compare_ = (mask_run_ls ==lumi_mask)
print(f" {mask_run_ls = },  {mask_run_ls.shape = },  {np.sum(mask_run_ls) = }  ")
print(f" {lumi_mask = },  {lumi_mask.shape = },  {np.sum(lumi_mask) = }  ")
print(f" {compare_ = },  {compare_.shape = },  {np.sum(compare_) = }  ")
print(f"test ended"); sys.stdout.flush();