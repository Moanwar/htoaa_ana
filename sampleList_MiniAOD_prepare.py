import os
import subprocess
import json
from collections import OrderedDict as OD
from copy import deepcopy
import argparse
import glob



from htoaa_Settings import *
from DASQueryHelper import getDASDatasetFiles, getDASParentDataset
from SamplesAndCrosssection_NanoAOD_2018 import list_datasetAndXs_2018

printLevel = 0
sXS     = "xs"
sNameSp = "nameSp"


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
    parser.add_argument('-era',                 dest='era', type=str, default=Era_2018, choices=[Era_2016, Era_2017, Era_2018], required=False)
    parser.add_argument('-updateCrossSections', action='store_true', default=False, help='update cross-sections only')
    #parser.add_argument('-addSkimmedNanoAOD',   action='store_true', default=False, help='add skimmed NanoAOD files to the existing sample list')
    args=parser.parse_args()

    era                 = args.era
    updateCrossSections = args.updateCrossSections
    #addSkimmedNanoAOD   = args.addSkimmedNanoAOD
    print(f"era: {era}")
    print(f"{updateCrossSections = }")
    #print(f"{addSkimmedNanoAOD = }")

    list_datasetAndXs = None
    sFileSamplesInfo_toUse = None
    if era == Era_2018:
        list_datasetAndXs = list_datasetAndXs_2018

        
    sFileSamplesInfo_toUse = sFileSamplesInfo[era]
    sFileSamplesInfo_toUse = sFileSamplesInfo_toUse.replace('.json', 'miniAOD_v0.json')

    for datasetName, datasetDetails in list_datasetAndXs.items():
        print(f" { datasetName }")
        datasetParentName = getDASParentDataset(datasetName)

        nEventsTotal, nFiles, files = getDASDatasetFiles(datasetParentName)


        break