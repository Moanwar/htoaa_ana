#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:43:51 2020

@author: si_sutantawibul1
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xgboost as xgb
import pickle
import random
import os
from sklearn.model_selection import GridSearchCV
import time

from dataManager import processData, ggHPath, BGenPaths, bEnrPaths, allVars, trainVars

from sklearn.metrics import roc_curve, auc, accuracy_score
from sklearn.model_selection import train_test_split

ntreeslist = [400, 800, 1000]#, 1500, 2000, 2500, 3000]
treedepthlist = [2,3,4,6]#,8,10]
mcwlist = [1,2,3]#,4,5]
lrlist = [0.005, 0.01, 0.05]#, 0.1, 0.2]
earlystoplist = [100]#, 200, 300, 400, 500]

data = pickle.load(open('data.pkl', 'rb'))
## normalizing the weights
data.loc[data['target']==0, ['final_weights']] *= 100000/data.loc[data['target']==0]['final_weights'].sum()
data.loc[data['target']==1, ['final_weights']] *= 100000/data.loc[data['target']==1]['final_weights'].sum()

## result list will be of the form:
## [test AUROC, AUROC diff, ntrees, treedpth, mcw, lr, earlystop, time/round]
results = np.empty([0,8])

for ntrees in ntreeslist:
    for treedepth in treedepthlist:
        for mcw in mcwlist:
            for lr in lrlist:
                for earlystop in earlystoplist:
                    store = np.empty([0,2])
                    start = time.time()
                    for i in range(5):
                        randInt = random.randint(0,100)
                        print("random int: " + str(randInt))
                        trainData, testData = train_test_split(data, random_state=randInt, test_size=0.4)
                        
                        cls = xgb.XGBClassifier(
                            n_estimators = ntrees,
                            max_depth = treedepth,
                            min_child_weight = mcw,
                            learning_rate = lr,
                            )
                        eval_set = [(testData[trainVars], testData['target'])]
                        
                        cls.fit(trainData[trainVars], trainData['target'],
                                sample_weight=(trainData['final_weights']),
                                early_stopping_rounds=earlystop, eval_metric="auc",
                                eval_set = eval_set,
                                sample_weight_eval_set=[testData['final_weights']],
                                verbose=1)
                        train_proba = cls.predict_proba(trainData[trainVars])
                        fpr, tpr, thresholds = roc_curve(trainData['target'], train_proba[:,1])
                        train_auc = auc(fpr, tpr)
                        print("XGBoost train set auc - {}".format(train_auc))
    
                        test_proba = cls.predict_proba(testData[trainVars])
                        fprt, tprt, thresholds = roc_curve(testData['target'], test_proba[:,1])
                        test_auc = auc(fprt, tprt)
                        print("XGBoost test set auc - {}".format(test_auc))
                        store = np.append(store, np.array([[train_auc, test_auc]]), axis=0)
                    avetrain = np.mean(store, axis=0)[0]
                    avetest = np.mean(store, axis=0)[1]
                    stop = time.time()
                    results = np.append(results,
                                        np.array([[avetest, abs(avetest-avetrain), ntrees,
                                                   treedepth, mcw, lr, earlystop, stop-start]]), axis=0)
                    #results.append([avetest, avetest-avetrain, ntrees,
                    #               treedepth, mcw, lr, earlystop, stop-start])



## need to also pickle the result array, so if i fucked up plotting, at least
## the results array is out there
pickle.dump(results, open('gridSearch/gridSearchResults.pkl', 'wb'))
plt.scatter(results[:,0], results[:, 1])
plt.xlabel('test AUC')
plt.ylabel('diff AUC')
plt.savefig('gridSearch/diffVstest.png')
