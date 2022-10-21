import subprocess
import json
from collections import OrderedDict as OD
from copy import deepcopy
import argparse

from htoaa_Settings import *


printLevel = 0

sXS= "xs"
list_datasetAndXs_2018 = OD([

    ## QCD_HTxtoy_PSWeights
    ("/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 187300000.0}),
    ("/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 23590000.0}),
    ("/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1555000.0}),
    ("/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 324500.0}),
    ("/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 30310.0}),

    ("/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM", {sXS: 6444.0}),
    ("/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1127.0}),
    ("/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 109.8}),
    ("/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 21.98}),
    
    ("/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM", {sXS: 1.0}),

    #("", {sXS: }),    
])


sDataset = "dataset"
sNanoAOD_nFiles = "nanoAOD_nFiles"
sNanoAOD = "nanoAOD"
sCross_section = "cross_section"
sNEvents = "nEvents"
sSumEvents = "sumEvents"
sampleDetail_dict_template = OD([
    (sDataset, ""),
    (sNanoAOD_nFiles, -1),
    (sNanoAOD, []),
    (sCross_section, -1.),
    (sNEvents, -1),
    (sSumEvents, -1),
])

def getDatasetFiles(dataset):
    cmd1 = ['bash','-c', 'dasgoclient --query="file dataset=%s" --format=json'%(dataset)]
    if printLevel >= 10:
        print(f"cmd1: {cmd1}")
    output = subprocess.check_output(cmd1) # output in bytes
    output = output.decode("utf-8") # convert output in bytes to string
    output = json.loads(output) # convert output in 'string' to 'dict'
    nFiles = output['nresults']
    files  = []
    nEventsTotal = 0
    
    if nFiles != len(output['data']):
        print(f"nFiles != len(output['data']... something is wrong.. \t\t **** ERROR ****")
        exit(0)
        
    for iFile in range(nFiles):
        if len(output['data'][iFile]['file']) != 1:
            print(f"len(output['data'][iFile]['file']) != 1: Assumption of a single entry list 'output['data'][iFile]['file']' seems wrong... need to follow up.. \t\t **** ERROR **** ")
            exit(0)
            
        file_LFN = output['data'][iFile]['file'][0]['name']
        nEvents  = output['data'][iFile]['file'][0]['nevents']
        if printLevel >= 5:
            print(f"file_LFN: {file_LFN}, nEvents ({type(nEvents)}): {nEvents}, nEventsTotal: {nEventsTotal}  {output['data'][iFile]['file'][0]}")
            
        files.append( file_LFN )
        nEventsTotal += nEvents

    if printLevel >= 3:
        print(f"\ndataset: {dataset}, nEventsTotal: {nEventsTotal}, nFiles: {nFiles}, files: {files}")
    return nEventsTotal, nFiles, files
    

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='htoaa analysis wrapper')
    parser.add_argument('-era', dest='era', type=str, default=Era_2018, choices=[Era_2016, Era_2017, Era_2018], required=False)
    args=parser.parse_args()

    era           = args.era
    print(f"era: {era}")


    list_datasetAndXs = None
    sFileSamplesInfo_toUse = None
    if era == Era_2018:
        list_datasetAndXs = list_datasetAndXs_2018
        sFileSamplesInfo_toUse = sFileSamplesInfo[era]

    sFileSamplesInfo_toUse = sFileSamplesInfo_toUse.replace('.json', '_v0.json')

        
    samples_details = OD()
    for datasetName, datasetDetails in list_datasetAndXs.items():
        sampleDetails_dict = deepcopy(sampleDetail_dict_template)
        print(f"sampleDetails_dict_0: {sampleDetails_dict}")
        sampleName = datasetName.split('/')[1]
        sampleDetails_dict[sDataset] = datasetName
        sampleDetails_dict[sCross_section] = datasetDetails[sXS]

        nEventsTotal, nFiles, files = getDatasetFiles(datasetName)
        sampleDetails_dict[sNEvents] = nEventsTotal
        sampleDetails_dict[sNanoAOD_nFiles] = nFiles
        sampleDetails_dict[sNanoAOD] = files

        samples_details[sampleName] = sampleDetails_dict

    if printLevel >= 0:
        print("\n\nsamples:: \n",json.dumps(samples_details, indent=4))

    with open(sFileSamplesInfo_toUse, "w") as fSampleInfo:
        json.dump(samples_details, fSampleInfo, indent=4)

        print(f"\n\nNew sample list wrote to {sFileSamplesInfo_toUse}")
    
    

    


