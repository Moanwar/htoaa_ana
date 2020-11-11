#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 16:00:41 2020

@author: si_sutantawibul1
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xgboost as xgb
import pickle
import random
import os

from dataManager import processData, ggHPath, BGenPaths, bEnrPaths, allVars, trainVars

from sklearn.metrics import roc_curve, auc, accuracy_score
from sklearn.model_selection import train_test_split


from optparse import OptionParser
parser = OptionParser()
parser.add_option("--ntrees", type="int", dest="ntrees", help="hyp", default = 1000)
parser.add_option("--treeDeph", type="int", dest="treeDeph", help="hyp", default = 2)
parser.add_option("--lr", type="float", dest="lr", help="hyp", default = 0.05)
parser.add_option("--mcw", type="float", dest="mcw", help="hyp", default = 1)
parser.add_option("--doXML", action="store_true", dest="doXML", help="Do save not write the xml file", default=True)
(options, args) = parser.parse_args()

hyppar="ntrees_"+str(options.ntrees)+"_deph_"+str(options.treeDeph)+"_mcw_"+str(options.mcw)+"_lr_"+str(options.lr)
print(hyppar)


data = pd.DataFrame()
data = data.append(processData(ggHPath, 'ggH'), ignore_index=True, sort=False)
#data = data.append(processData(BGenPath, 'BGen'), ignore_index=True, sort = False)
#data = data.append(processData(bEnrPath, 'bEnr'), ignore_index=True, sort = False)

for BGenPath in BGenPaths:
    data = data.append(processData(BGenPath, 'BGen'), ignore_index=True, sort=False)
for bEnrPath in bEnrPaths:
    data = data.append(processData(bEnrPath, 'bEnr'), ignore_index=True, sort=False)

pickle.dump(data, open('data.pkl', 'wb'))
#data = pickle.load(open('data.pkl', 'rb'))


## normalizing the weights?? why do we have to do this? how do we do this?
data.loc[data['target']==0, ['final_weights']] *= 100000/data.loc[data['target']==0]['final_weights'].sum()
data.loc[data['target']==1, ['final_weights']] *= 100000/data.loc[data['target']==1]['final_weights'].sum()

## 
dataSig = data.loc[data.target == 1]
dataBg = data.loc[data.target == 0]

print('Signal event count: ' + str(len(dataSig.index)))
print('Background event count: ' + str(len(dataBg.index)))

## drop events with NaN weights - for safety
data.dropna(subset=['final_weights'],inplace = True)
data.fillna(0)

## split data into training and testing
#randInt = random.randint(0,100)
randInt = 7
print("random int: " + str(randInt))
trainData, testData = train_test_split(data, random_state=randInt)

## training 
cls = xgb.XGBClassifier(
    n_estimators = options.ntrees,
    max_depth = options.treeDeph,
    min_child_weight = options.mcw,
    learning_rate = options.lr,
    )
cls.fit(trainData[trainVars], trainData['target'], sample_weight=(trainData['final_weights']))

print ("XGBoost trained")


## get info for ROC
train_proba = cls.predict_proba(trainData[trainVars])
fpr, tpr, thresholds = roc_curve(trainData['target'], train_proba[:,1])
train_auc = auc(fpr, tpr)
print("XGBoost train set auc - {}".format(train_auc))

test_proba = cls.predict_proba(testData[trainVars])
fprt, tprt, thresholds = roc_curve(testData['target'], test_proba[:,1])
test_auc = auc(fprt, tprt)
print("XGBoost test set auc - {}".format(test_auc))
fig, ax = plt.subplots()

prediction = cls.predict(testData[trainVars])
accuracy = accuracy_score(testData['target'], prediction)
print("XGBoost test accuracy - {}".format(accuracy))


## put the train and test auroc in a file so i can see all the stuff
bdtscorefile = open('bdtScores.txt', 'a')
bdtscorefile.write('\n')
bdtscorefile.write(hyppar)
bdtscorefile.write('\n')
bdtscorefile.write('Train: {} '.format(str(train_auc)))
bdtscorefile.write('\n')
bdtscorefile.write('Test: {} '.format(str(test_auc)))
bdtscorefile.write('\n')
bdtscorefile.write('Diff: {} '.format(str(test_auc-train_auc)))
bdtscorefile.write('\n')


## draw them rocs
fig, ax = plt.subplots(figsize=(8, 8))
train_auc = auc(fpr, tpr)
ax.plot(fpr, tpr, lw=1, color='g',label='XGB train (area = %0.5f)'%(train_auc))
ax.plot(fprt, tprt, lw=1, ls='--',color='g',label='XGB test (area = %0.5f)'%(test_auc) )
ax.set_ylim([0.0,1.0])
ax.set_xlim([0.0,1.0])
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')
ax.legend(loc="lower right")
ax.grid()
ax.set_title(hyppar)
fig.savefig("plots/roc_{}.png".format(hyppar))
plt.clf()


## add the BDT scores to the df so manip is easier. hopefully
trainData.loc[:,'BDTScore'] = train_proba[:,1]
testData.loc[:, 'BDTScore'] = test_proba[:,1]


## making bdt score figs babey
fig, ax = plt.subplots(figsize=(8,8))
ax.hist(trainData.BDTScore.loc[trainData.target == 1], weights=trainData.final_weights.loc[trainData.target == 1],  bins=20, ls = '--', histtype='step',stacked=True, label='train signal', density= True)
ax.hist(trainData.BDTScore.loc[trainData.target == 0], weights=trainData.final_weights.loc[trainData.target == 0], bins=20, ls = '--', histtype='step',stacked=True, label='train background', density=True)
ax.hist(testData.BDTScore.loc[testData.target == 1], weights=testData.final_weights.loc[testData.target == 1], bins=20, histtype='step', stacked=True, label='test signal', fill=False, density=True)
ax.hist(testData.BDTScore.loc[testData.target == 0], weights=testData.final_weights.loc[testData.target == 0], bins=20, histtype='step', stacked=True, label='test background', fill=False, density=True)
ax.legend(loc='lower right')
ax.set_title('BDT score')
ax.set_xlabel('BDT Score')
fig.savefig("plots/BDT_score_{}.png".format(hyppar))
plt.clf()


## distributions of things yeah
for colName in allVars:
    hist_params = {'density': True, 'histtype': 'bar', 'fill': True , 'lw':3, 'alpha' : 0.4}
    nbins = 40
    min_valueS, max_valueS = np.percentile(dataSig[colName], [0, 99.8])
    min_valueB, max_valueB = np.percentile(dataBg[colName], [0, 99.8])
    range_local = (min(min_valueS,min_valueB),  max(max_valueS,max_valueB))
    valuesS, binsS, _ = plt.hist(
        dataSig[colName].values,
        range = range_local,
        bins = nbins, edgecolor='b', color='b',
        label = "Signal", **hist_params
        )
    to_ymax = max(valuesS)
    to_ymin = min(valuesS)
    valuesB, binsB, _ = plt.hist(
        dataBg[colName].values,
        range = range_local,
        bins = nbins, edgecolor='g', color='g',
        label = "Background", **hist_params,
        weights = dataBg['final_weights']
        )
    to_ymax2 = max(valuesB)
    to_ymax  = max([to_ymax2, to_ymax])
    to_ymin2 = min(valuesB)
    to_ymin  = max([to_ymin2, to_ymin])
    plt.ylim(ymin=to_ymin*0.1, ymax=to_ymax*1.2)
    plt.legend(loc='best')

    plt.xlabel(colName)
    plt.savefig("distributions/dist_{}".format(colName))
    plt.clf()


## save model to pickle
pklpath="XGB_classifier_BothBg"
if options.doXML==True :
    pickle.dump(cls, open(pklpath+".pkl", 'wb'))
    file = open(pklpath+"pkl.log","w")
    file.write(str(trainVars)+"\n")
    file.close()












