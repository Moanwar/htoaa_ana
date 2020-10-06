import pickle
import numpy as np
import pandas as pd
from analib import PhysObj, Event
from info import trainVars, allVars, cutVars, cutDict, weightDict
from data_manager import processData
import collections
import matplotlib.pyplot as plt
import sys

## accept argument as filename
import argparse
parser = argparse.ArgumentParser(description='fileName for BDT')
parser.add_argument('-f', dest='fileName',
                    help='file name to be processed(without ".root")')
#parser.add_argument('-lf', dest='listFileNames',
#                    help='list of file names to be process (without ".root"')
args = parser.parse_args()

if not (args.fileName): #or args.listFileNames):
    parser.error('please provide name of file to be processed')
    #parser.error('please provide a file name or list of file names to be processed')

## prepare testdata for prediction
#fileName = 'GGH_HPT'
data = pd.DataFrame()

if args.fileName:
    data = processData(args.fileName)
###### !!! TODO:: make the list work properly
# elif args.listFileNames:
#     for fileName in args.listFileNames:
#         tmpData = processData(fileName)
#         data = data.append(tmpData, ignore_index=True, sort = False)


## drop all columns and rows that all nan, then fill nan->0
data = data.dropna(axis = 1, how = 'all') 
data = data.dropna(how = 'all')
data = data.fillna(0)



## load model from file
loaded_model = pickle.load(open('XGB_classifier_8Var.pkl', 'rb'))

## make predictions for test data
## i'm going to cry the dataframe needs to have the same column order
prediction = loaded_model.predict_proba(data[trainVars])
BDTScore = prediction[:,1]
fig, ax = plt.subplots(figsize=(8,8))
ax.hist(BDTScore, bins=20, histtype='step',stacked=True, density=True)
ax.set_title('BDT score')
ax.set_xlabel('BDT Score')
fig.savefig("loadedModel/BDTScore.png")
plt.clf()


## print out how many things predicted as signal and background
prediction_binary = loaded_model.predict(data[trainVars])
count = collections.Counter(prediction_binary)
print('signal-background cutoff: 0.50')
print('signal count: {}'.format(count[1]))
print('background count: {}'.format(count[0]))


## distribution plots

for colName in trainVars:
    #hist_params = {'density': Fase, 'histtype': 'bar', 'fill': True , 'lw':3, 'alpha' : 0.4}
    nbins = 40
    dist = plt.hist(data[colName].values, nbins, weights=data['weights'],
                    density=True, histtype='bar')
    plt.xlabel(colName)
    plt.savefig("loadedModel/distributions_%s.png" % colName)
    plt.clf()

