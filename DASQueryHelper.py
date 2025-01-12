import subprocess
import json

printLevel = 0

def getDASDatasetFiles(dataset):
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


def getDASParentDataset(dataset):    
    cmd1 = ['bash','-c', 'dasgoclient --query="parent dataset=%s" --format=json'%(dataset)]
    if printLevel >= 10:
        print(f"cmd1: {cmd1}")
    output = subprocess.check_output(cmd1) # output in bytes
    output = output.decode("utf-8") # convert output in bytes to string
    output = json.loads(output) # convert output in 'string' to 'dict'

    print(f"getDASParentDataset({dataset}): {output = }")
    print(f"\n{output.keys() = }")
    print(f"\n{output['nresults'] = }")
    print(f"\n{output['data'] = }")
    print(f"\n{output['data'][0].keys() = }")
    print(f"\n{output['data'][0]['parent'][0] = }")
    print(f"\n{output['data'][0]['parent'][0]['name'] = }")

    nResults = output['nresults']
    files  = []
    nEventsTotal = 0
    
    if nResults != len(output['data']):
        print(f"nFiles != len(output['data']... something is wrong.. \t\t **** ERROR ****")
        exit(0)
        
    for iResult in range(nResults):
        if len(output['data'][iResult]['parent']) != 1:
            print(f"len(output['data'][iResult]['parent']) != 1: Assumption of a single entry list 'output['data'][iFile]['parent']' seems wrong... need to follow up.. \t\t **** ERROR **** ")
            exit(0)
            
        parent_dataset = output['data'][iResult]['parent'][0]['name']

    if printLevel >= 3:
        print(f"\ndataset: {dataset},  parent ({nResults}) {parent_dataset}, ")
    return parent_dataset

