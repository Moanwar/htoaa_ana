'''
 Script to prepare Samples_2018UL.json by collecting information with 'DAQ queries'.
 To run:
        1] Set voms proxy:
            voms-proxy-init -voms cms -rfc

        2] python3 samplesList_prepare.py -era 2018
              Add -updateCrossSections to update cross-section values in the existing Samples_2018UL.json

        3] Script produces 'Samples_2018UL_v0.json' output.
           cp Samples_2018UL_v0.json Samples_2018UL.json   if you are satistied with Samples_2018UL_v0.json
'''

# List of cross-section used in this macro also with refereces can be found at 
# https://docs.google.com/spreadsheets/d/1LQDKBWGTdsT1uBumq9tz6RwRTCAVK2G0zmC6qaiqLkw/edit?usp=sharing


import os
import subprocess
import json
from collections import OrderedDict as OD
from copy import deepcopy
import argparse
import glob



from htoaa_Settings import *
from DASQueryHelper import getDASDatasetFiles
from SamplesAndCrosssection_NanoAOD_2018 import list_datasetAndXs_2018

printLevel = 0
sXS     = "xs"
sNameSp = "nameSp"



sNDatasets                = "nDatasets"
sDataset                   = "dataset"
sNanoAOD_nFiles            = "nanoAOD_nFiles"
sNanoAOD                   = "nanoAOD"
sCross_section             = "cross_section"
sNEvents                   = "nEvents"
sSumEvents                 = "sumEvents"
sSkimmedNanoAOD_nFiles      = "skimmedNanoAOD_nFiles"
sSkimmedNanoAOD             = "skimmedNanoAOD"
sNFiles                    = "nFiles"
sampleDetail_dict_template = OD([
    (sCross_section,  -1.),
    (sNEvents,         0),
    (sSumEvents,       0),
    (sNDatasets,       0),
    (sDataset,        []),
    (sNanoAOD_nFiles,  0),
    (sNanoAOD,        []),
    #(sSkimmedNanoAOD_nFiles,  0),
    #(sSkimmedNanoAOD,        []),    
])

    

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
    sFileSamplesInfo_toUse = sFileSamplesInfo_toUse.replace('.json', '_v0.json')

    isSamplesListFromScratch = False
    samples_details = None       
    '''
    if updateCrossSections or addSkimmedNanoAOD:
        # update cross sections
        with open(sFileSamplesInfo[era]) as fSamplesInfo:
            samples_details = json.load(fSamplesInfo)
        print(f"samples_details.keys(): {samples_details.keys()}")
    else:
        samples_details = OD()
    '''
    if os.path.exists( sFileSamplesInfo[era] ):
        with open(sFileSamplesInfo[era]) as fSamplesInfo:
            samples_details = json.load(fSamplesInfo)
        print(f"samples_details.keys(): {samples_details.keys()}")
    else:
        samples_details = OD()
        isSamplesListFromScratch = True


    # Reset sNDatasets, sDataset, sSkimmedNanoAOD etc
    # Should not reset when running with updateCrossSections
    if not updateCrossSections:
        #for sampleName_ in samples_details:
        for iSample, (datasetName, datasetDetails) in enumerate(list_datasetAndXs.items()):
            datasetName_parts            = datasetName.split('/')
            sampleName                   = datasetName_parts[1]
            isMC                         = datasetName_parts[-1] == 'NANOAODSIM'
            if not isMC:
                # for data sample
                sampleName_part2 = (datasetName_parts[2]).split('-')[0] # 'Run2018A-UL2018_MiniAODv2_NanoAODv9-v2'
                sampleName = '%s_%s' % (sampleName, sampleName_part2)  # JetHT_Run2018A

            if sampleName not in samples_details: 
                samples_details[sampleName] = deepcopy(sampleDetail_dict_template)
                '''                 
                if isSamplesListFromScratch:
                    samples_details[sampleName] = deepcopy(sampleDetail_dict_template)
                else:
                    # when adding a new sample to existing samples.json, add in order as in list_datasetAndXs_2018
                    samples_details_tmp_ = samples_details.items()
                    samples_details_tmp_.insert(iSample, (sampleName, deepcopy(sampleDetail_dict_template)))
                    samples_details = samples_details_tmp_'''
            else: # reset sNDatasets, sDataset
                samples_details[sampleName][sNDatasets] = 0
                samples_details[sampleName][sDataset  ] = []
                samples_details[sampleName][sNEvents       ]         = 0
                samples_details[sampleName][sNanoAOD_nFiles]         = 0
                samples_details[sampleName][sNanoAOD       ]         = []


            samples_details[sampleName][sSkimmedNanoAOD] = OD()
            for skimName_ in sPathSkimmedNanoAODs[era]:
                samples_details[sampleName][sSkimmedNanoAOD]['%s_%s' % (skimName_, sNFiles)] = 0
                samples_details[sampleName][sSkimmedNanoAOD][skimName_                     ] = []

            # temperary fix
            #if "skimmedNanoAOD_nFiles" in samples_details[sampleName_]:
            #    samples_details[sampleName_].pop("skimmedNanoAOD_nFiles", None)


    # Now calculate..
    for datasetName, datasetDetails in list_datasetAndXs.items():
        datasetName_parts            = datasetName.split('/')
        sampleName                   = datasetName_parts[1]
        isMC                         = datasetName_parts[-1] == 'NANOAODSIM'
        #print(f"{datasetName = },  {isMC = }")

        if not isMC:
            # for data sample
            sampleName_part2 = (datasetName_parts[2]).split('-')[0] # 'Run2018A-UL2018_MiniAODv2_NanoAODv9-v2'
            sampleName = '%s_%s' % (sampleName, sampleName_part2)  # JetHT_Run2018A

        if updateCrossSections:
            if isMC:
                if sampleName in samples_details:
                    samples_details[sampleName][sCross_section] = datasetDetails[sXS]
                else:
                    print(f"Running updateCrossSections mode, but {sampleName} is not in samples_details \t **** ERROR **** \nTerminating...")
                    exit(0)
            
            continue


        # Add path to skimmed NanoAOD samples
        sDatasetExt = ''
        for sTmp1_ in (datasetName_parts[2]).split('_'): # RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2
            for sTmp2_ in sTmp1_.split('-'):
                if 'ext' in sTmp2_:
                    sDatasetExt = '_%s' % sTmp2_ # '_ext1'
        sSampleNameDir_used = '%s%s' % (datasetName_parts[1], sDatasetExt)

        sSampleTagDir_used = datasetName_parts[2]
        sSampleTagDir_used = sSampleTagDir_used.replace('_NanoAODv9', '')
        #if datasetName_parts[1] == 'SingleMuon':
        sSampleTagDir_used = sSampleTagDir_used.replace('v1', 'v*') # SingleMuon, EGamma  had DatasetTag v2/3 in MiniAOD and v1 in NanoAOD. So use wildcard charester *

        sDataType = 'MC' if isMC else 'Data'

        # Loop over multiple skim version we have
        for skimName_ in sPathSkimmedNanoAODs[era]:
            sPathSkimmedNanoAODs_toUse = sPathSkimmedNanoAODs[era][skimName_][sDataType] 
            # /eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/$SAMPLENAME/r1/PNet_*.root           
            # /eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/PNet_v1_2023_10_06/$SAMPLETAG/$SAMPLENAME/r*/PNet_*.root
            sPathSkimmedNanoAODs_toUse = sPathSkimmedNanoAODs_toUse.replace('$SAMPLENAME', sSampleNameDir_used)
            sPathSkimmedNanoAODs_toUse = sPathSkimmedNanoAODs_toUse.replace('$SAMPLETAG',  sSampleTagDir_used)

            sSkimmedNanoAODs = []
            if "*" in sPathSkimmedNanoAODs_toUse: sSkimmedNanoAODs.extend( glob.glob(sPathSkimmedNanoAODs_toUse) )
            else:                                 sSkimmedNanoAODs.append( sPathSkimmedNanoAODs_toUse )
            print(f"{sPathSkimmedNanoAODs_toUse = }, {sSkimmedNanoAODs = }")
            
            '''if sSkimmedNanoAOD_nFiles not in samples_details[sampleName]:
                samples_details[sampleName][sSkimmedNanoAOD_nFiles]  = len(sSkimmedNanoAODs)
                samples_details[sampleName][sSkimmedNanoAOD]         = sSkimmedNanoAODs
            else:
                samples_details[sampleName][sSkimmedNanoAOD_nFiles] += len(sSkimmedNanoAODs)
                samples_details[sampleName][sSkimmedNanoAOD].extend(   sSkimmedNanoAODs )'''
            samples_details[sampleName][sSkimmedNanoAOD]['%s_%s' % (skimName_, sNFiles)] += len(   sSkimmedNanoAODs )
            samples_details[sampleName][sSkimmedNanoAOD][skimName_                     ].extend(   sSkimmedNanoAODs )

        #if addSkimmedNanoAOD:
        #    continue
        

        #if sampleName not in samples_details:
        #    samples_details[sampleName] = deepcopy(sampleDetail_dict_template)

        samples_details[sampleName][sNDatasets]         += 1
        samples_details[sampleName][sDataset  ].append(    datasetName )
        
        if datasetName_parts[-1] == 'NANOAODSIM':
            # for MC sample
            samples_details[sampleName][sCross_section] = datasetDetails[sXS]
        else:
            # for data sample
            if sCross_section in samples_details[sampleName]:
                del samples_details[sampleName][sCross_section]
            if sSumEvents in samples_details[sampleName]:
                del samples_details[sampleName][sSumEvents]
        
        nEventsTotal, nFiles, files = getDASDatasetFiles(datasetName)
        samples_details[sampleName][sNEvents       ]         += nEventsTotal
        samples_details[sampleName][sNanoAOD_nFiles]         += nFiles
        samples_details[sampleName][sNanoAOD       ].extend(    files        )
        

    # Running the samples_prepare.py on exising Sample.json changes order of samples. 
    # Hence, to forcefully follow order of samples as in list_datasetAndXs, 
    # order samples_details keys to follow order in list_datasetAndXs
    samples_details_inOrder = OD()
    for datasetName, datasetDetails in list_datasetAndXs.items():
        datasetName_parts            = datasetName.split('/')
        sampleName                   = datasetName_parts[1]
        isMC                         = datasetName_parts[-1] == 'NANOAODSIM'
        if not isMC:
            # for data sample
            sampleName_part2 = (datasetName_parts[2]).split('-')[0] # 'Run2018A-UL2018_MiniAODv2_NanoAODv9-v2'
            sampleName = '%s_%s' % (sampleName, sampleName_part2)  # JetHT_Run2018A

        samples_details_inOrder[sampleName] = samples_details[sampleName]
    # Replace samples_details with samples_details_inOrder
    samples_details = samples_details_inOrder 


    if printLevel >= 0:
        print("\n\nsamples:: \n",json.dumps(samples_details, indent=4))

    with open(sFileSamplesInfo_toUse, "w") as fSampleInfo:
        json.dump(samples_details, fSampleInfo, indent=4)

        print(f"\n\nNew sample list wrote to {sFileSamplesInfo_toUse}")
        print(f"Copy {sFileSamplesInfo_toUse} to {sFileSamplesInfo[era]}   if you are satistied with {sFileSamplesInfo_toUse}")
        print(f"i.e.   cp {sFileSamplesInfo_toUse} {sFileSamplesInfo[era]}")
    
    

    


