import pickle
import numpy as np
import pandas as pd
from analib import PhysObj, Event
from info import trainVars, allVars, cutVars, cutDict, weightDict
from data_manager import processData
import collections
import matplotlib.pyplot as plt


## prepare testdata for prediction 
fileName = 'GGH_HPT'
data = processData(fileName) 

## drop all columns and rows that all nan, then fill nan->0
data = data.dropna(axis = 1, how = 'all') 
data = data.dropna(how = 'all')
data = data.fillna(0)



## load model from file
loaded_model = pickle.load(open('XGB_classifier_8Var.pkl', 'rb'))

## make predictions for test data
prediction = loaded_model.predict(data.iloc[:,:-2])
count = collections.Counter(prediction)

## am i also making plots of how many signal and background we found 
## for each region of bin? 
## what is my bin even??????? 
## uh this feels wrong. 
hist = plt.hist(prediction, 2)
plt.title("prediction of number of signal (1) and background(0) events")
plt.savefig("loadedModel/prediction.png")
plt.clf()


## get column names (without the weight, target)
colNames = list(data.columns)
colNames = colNames[:-2]

## distribution plots
for colName in colNames: 
    #hist_params = {'density': Fase, 'histtype': 'bar', 'fill': True , 'lw':3, 'alpha' : 0.4}
    nbins = 8
    dist = plt.hist(data[colName].values, nbins)
    # min_valueS, max_valueS = np.percentile(dataSig[colName], [0.0, 99])        
    # min_valueB, max_valueB = np.percentile(dataBg[colName], [0.0, 99])
    # range_local = (min(min_valueS,min_valueB),  max(max_valueS,max_valueB))
    # valuesS, binsS, _ = plt.hist(
    #     dataSig[colName].values,
    #     range = range_local,
    #     bins = nbins, edgecolor='b', color='b',
    #     label = "Signal", **hist_params
    #     )   
    # to_ymax = max(valuesS)
    # to_ymin = min(valuesS)
    # valuesB, binsB, _ = plt.hist(
    #     dataBg[colName].values,
    #     range = range_local,
    #     bins = nbins, edgecolor='g', color='g',
    #     label = "Background", **hist_params
    #     )
    # to_ymax2 = max(valuesB)
    # to_ymax  = max([to_ymax2, to_ymax])
    # to_ymin2 = min(valuesB)
    # to_ymin  = max([to_ymin2, to_ymin])
    # plt.ylim(ymin=to_ymin*0.1, ymax=to_ymax*1.2)
    # plt.legend(loc='best')
    plt.xlabel(colName)
    plt.savefig("loadedModel/distributions_%s.png" % colName)
    plt.clf()

