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
import uproot as uproot
#import uproot3 as uproot
#import parse
from parse import *
import enum

import hist

from coffea import hist as histC


class TestStrEnum(enum.StrEnum):
    Test1 = "test1"
    Test2 = "test2"

class TestEnum(enum.Enum):
    Test1 = 1
    Test2 = 2
    

if __name__ == '__main__':
    sFile = '../data/lumiScale/2018.root'
    sHisto = 'QCD/QCD_LumiScale_PhSpOverlapRewghted'

    fIn = uproot.open(sFile)
    h = fIn[sHisto]

    print(f"{h = } ")
    print(f"{type(h) = } ")
    print(f"{h.classname = } ")
    print(f"{h._members = } ")
    
    print(f"{h.axes = } ")
    print(f"{h.axes[0] = } ")
    print(f"{h.axes[0].edges = } ")

    #print(f"{h.values() = }")

    print(f"{h.to_hist() = }")
    h1 = h.to_hist()
    print(f"{h1 = }")
    print(f"{h1.axes = }")
    print(f"{h1.axes[0] = }")
    print(f"{h1.axes[0][0] = }")
    print(f"{h1.axes[1][0] = }")
    print(f" {len(h1.axes[0]) = }")

    print(f"\n {h1[1, 0] = } ")
    print(f" {h1[1, 1] = } ")
    print(f" {h1[1, 'QCD_bGen'] = } ")
    print(f" {h1[120.5j, 'QCD_bGen'] = } ")

    HT = [121.8, 345.4, 3010.6]

    ls = np.ones(len(HT))

    for iHTBin in range(len(h1.axes[0])):
        print(f"{iHTBin = }, {h1.axes[0][iHTBin] = }, {h1.axes[0][iHTBin][0] = },  \t {h1[iHTBin, 'QCD_Incl_Remnant'] = },  {h1[iHTBin, 'QCD_bGen'] = },  {h1[iHTBin, 'QCD_bEnrich'] = } ")

    




    # variable width binned histogram    
    HTBinning = [ 100, 200, 300, 500, 700, 1000, 1500, 2000]
    data1 = np.random.uniform(low=0, high=4000, size=(1000,))
    
    h2 = hist.Hist(hist.axis.Variable(HTBinning, name='xaxis', label='HT [GeV]') )
    print(f"\n\n {data1 = }")
    h2.fill( data1 )


    hC2 = histC.Hist(
        "HT",
        histC.Bin("HT", "HT", 100, 0,2000)
    )
    hC2.fill( HT=data1 )

    hC2_1 = histC.Hist(
        "HT",
        histC.Bin("HT", "HT", HTBinning)
    )
    hC2_1.fill( HT=data1 )

    

    fOut = uproot.recreate("tmp.root")
    fOut['test/htest1'] = h2
    fOut['test/hCtest1'] = hC2.to_hist()
    fOut['test/hCtest1_1'] = hC2_1.to_hist()
    fOut.close()

    

    print(f"\n\n {TestEnum.Test2 = }")
    tmp1 = TestEnum.Test2
    print(f" {tmp1 = },  {tmp1.value = }, {tmp1.name = }, ") 
    print(f" {(tmp1 == TestEnum.Test1) =}")
    print(f" {(tmp1 == TestEnum.Test2) =}")

    print(f"\n\n {TestStrEnum.Test2 = }")
    tmp1 = TestStrEnum.Test2
    print(f" {tmp1 = },  {tmp1.value = }, {tmp1.name = }, ")
    print(f" {(tmp1 == TestStrEnum.Test1) =}")
    print(f" {(tmp1 == TestStrEnum.Test2) =}")


    
