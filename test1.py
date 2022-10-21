import subprocess
import json

cmd1 = 'dasgoclient --query=\"file dataset=\/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8\/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1\/NANOAODSIM\"' 
print(f"cmd1: {cmd1}")
print(f"cmd1.split() ({type(cmd1.split())}): {cmd1.split()}")

output, error = None, None

#process = subprocess.Popen(cmd1.split(), stdout=subprocess.PIPE)
#process = subprocess.Popen(['dasgoclient', '--query="file dataset=/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM"'], stdout=subprocess.PIPE)
#process = subprocess.Popen(['dasgoclient --query="file dataset=/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM"'], stdout=subprocess.PIPE)
#output, error = process.communicate()
print(f"output: {output}, error: {error}")


output = subprocess.check_output(['bash','-c', 'dasgoclient --query="file dataset=/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" --format=json'])
print(f"_2output: {output}, error: {error}")
output = output.decode("utf-8")
print(f"\n\n1 {type(output)}: {output}")

output = json.loads(output)
print(f"\n\n2 {type(output)}: {output}")
print(f"\n\nSiddh3 output.keys(): {output.keys()} ")
print(f"\n\nSiddh4 output['nresults']: {output['nresults']}")
print(f"\n\nSiddh5 output['data'] len(): {len(output['data'])}")
print(f"\n\nSiddh5 output['data'][0] ({type(output['data'][0])}):  {output['data'][0]} ")
print(f"\n\nSiddh5 output['data'][0].keys() ({type(output['data'][0])}):  {output['data'][0].keys()} ")
print(f"\n\nSiddh5 output['data'][0]['file'] ({type(output['data'][0]['file'])}) ({len(output['data'][0]['file'])}):  {output['data'][0]['file']} ")
print(f"\n\nSiddh5 output['data'][0]['file'][0] ({type(output['data'][0]['file'][0])}) :  {output['data'][0]['file'][0]} ")
print(f"\n\nSiddh5 output['data'][0]['file'][0].keys() ({type(output['data'][0]['file'][0])}) :  {output['data'][0]['file'][0].keys()} ")

print(f"\n\nSiddh5 output['data'][0]['file'][0]['name'] ({type(output['data'][0]['file'][0]['name'])}) :  {output['data'][0]['file'][0]['name']} ")
print(f"\n\nSiddh5 output['data'][0]['file'][0]['nevents'] ({type(output['data'][0]['file'][0]['nevents'])}) :  {output['data'][0]['file'][0]['nevents']} ")

