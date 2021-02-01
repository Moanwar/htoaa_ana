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
from sklearn.model_selection import GridSearchCV

from dataManager import processData, ggHPath, BGenPaths, bEnrPaths, allVars, trainVars, disc, TTJetsPaths, WJetsPaths, ZJetsPaths#, BGenPath, bEnrPath

from sklearn.metrics import roc_curve, auc, accuracy_score
from sklearn.model_selection import train_test_split


from optparse import OptionParser
parser = OptionParser()
parser.add_option("--ntrees", type="int", dest="ntrees", help="hyp", default = 400)
parser.add_option("--treeDeph", type="int", dest="treeDeph", help="hyp", default = 3)
parser.add_option("--lr", type="float", dest="lr", help="hyp", default = 0.01)
parser.add_option("--mcw", type="float", dest="mcw", help="hyp", default = 3)
parser.add_option("--doXML", action="store_true", dest="doXML", help="Do save not write the xml file", default=True)
parser.add_option("--HypOpt", action="store_true", dest="HypOpt", help="If you call this will not do plots with repport", default=False)
parser.add_option("--dist", action='store_true', dest='dist', default=False)
parser.add_option("--randInt", type='int', dest='randInt', default=1)
parser.add_option('--train', dest='train', default=True)
(options, args) = parser.parse_args()

hyppar= "ntrees_"+str(options.ntrees)+"_deph_"+str(options.treeDeph)+"_mcw_"+str(options.mcw)+"_lr_"+str(options.lr)
#hyppar = 'ggH_unweighted ' + ' randInt '+str(options.randInt)
print(hyppar)

if disc == None:
    condition = ''
elif disc == 'h':
    condition = 'high disc'
elif disc == 'm':
    condition = 'high+medium disc'

if options.HypOpt:
    test_size = 0.4
else:
    test_size = 0.4

##############
## commend this out after the first round of running this script
## it will dump the processed datafram to a pickle. Next time, don't have to
## reload the data from ROOT, can just open pickle
## read determines if you want to load data from a root or pikl
root = True
if root:
    data = pd.DataFrame()
    data = data.append(processData(ggHPath, 'ggH'), ignore_index=True, sort=False)
    #data = data.append(processData(BGenPath, 'BGen'), ignore_index=True, sort = False)
    #data = data.append(processData(bEnrPath, 'bEnr'), ignore_index=True, sort = False)

    BGenData = pd.DataFrame()
    bEnrData = pd.DataFrame()
    TTJetsData = pd.DataFrame()
    WJetsData = pd.DataFrame()
    ZJetsData = pd.DataFrame()

    ## getting the BGen, bEnr data into DataFrame format
    for BGenPath in BGenPaths:
        BGenData = BGenData.append(processData(BGenPath, 'BGen'), ignore_index=True, sort=False)
    for bEnrPath in bEnrPaths:
        bEnrData = bEnrData.append(processData(bEnrPath, 'bEnr'), ignore_index=True, sort=False)
    for TTJetsPath in TTJetsPaths:
        TTJetsData = TTJetsData.append(processData(TTJetsPath, 'TTJets'), ignore_index=True, sort=False)
    for WJetsPath in WJetsPaths:
        WJetsData = WJetsData.append(processData(WJetsPath, 'WJets'), ignore_index=True, sort=False)
    for ZJetsPath in ZJetsPaths:
        ZJetsData = ZJetsData.append(processData(ZJetsPath, 'ZJets'), ignore_index=True, sort=False)


    ## reshape them so they are the same size
    bEnrData.final_weights = bEnrData.final_weights*np.sum(bEnrData.final_weights)/np.sum(BGenData.final_weights)

    data = data.append(BGenData, ignore_index=True, sort=False)
    data = data.append(bEnrData, ignore_index=True, sort=False)
    data = data.append(TTJetsData, ignore_index=True, sort=False)
    data = data.append(WJetsData, ignore_index=True, sort=False)
    data = data.append(ZJetsData, ignore_index=True, sort=False)
    pickle.dump(data, open('data.pkl', 'wb'))

##############

if not root:
    data = pickle.load(open('data.pkl', 'rb'))


## this is when you want to be training a new BDT
if options.train:
    ## normalizing the weights
    data.loc[data['target']==0, ['final_weights']] *= 100000/data.loc[data['target']==0]['final_weights'].sum()
    data.loc[data['target']==1, ['final_weights']] *= 100000/data.loc[data['target']==1]['final_weights'].sum()


    dataSig = data.loc[data.target == 1]
    dataBg = data.loc[data.target == 0]

    print('Signal event count: ' + str(len(dataSig.index)))
    print('Background event count: ' + str(len(dataBg.index)))

    ## drop events with NaN weights - for safety
    data.dropna(subset=['final_weights'],inplace = True)
    data.fillna(0)

    print('after drop nan weights: {}'.format(data.shape))

    ## split data into training and testing
    # randInt = random.randint(0,100)
    randInt = options.randInt #7
    print("random int: " + str(randInt))
    trainData, testData = train_test_split(data, random_state=randInt, test_size=test_size)



    #############################

    ## grid search for best parameter. This doesnt work. why
    ## for finding the best hyperparameter yeah ##
    ## this was giving me warning of parameters not being used
    if options.HypOpt:
        param_grid = {
        	'n_estimators': list(range(100, 800, 50)),
        	'min_child_weight': list(range(1,10)),
        	'max_depth': [2,4],
        	'learning_rate': [0.01, 0.05, 0.1]
        }
        scoring = "roc_auc"
        early_stopping_rounds = 200
        cv=3
        cls = xgb.XGBClassifier()

        fit_params = { "eval_set" : [(testData[trainVars],testData['target'])],
                       "eval_metric" : ["auc"],
                       "early_stopping_rounds" : [early_stopping_rounds],
                       'sample_weight': testData['final_weights'].to_numpy()}

        gs = GridSearchCV(cls, fit_params, param_grid=param_grid, scoring=scoring, cv = cv,
                          verbose = 0, n_jobs=None)


        gs.fit(trainData[trainVars], trainData['target'])

        for i, param in enumerate(gs.cv_results_["params"]):
            print("params : {} \n    cv auc = {}  +- {} ".format(param,gs.cv_results_["mean_test_score"][i],gs.cv_results_["std_test_score"][i]))

        print("best parameters",gs.best_params_)
        print("best score",gs.best_score_)
        file = open("plots/XGB_HyperParameterGridSearch_GSCV.log","w")
        file.write(
            (trainVars)+"\n"+
     	    "best parameters"+str(gs.best_params_) + "\n"+
     	    "best score"+str(gs.best_score_)+ "\n"
        )
        for i, param in enumerate(gs.cv_results_["params"]):
            file.write("params : {} \n    cv auc = {}  +- {} {}".format(param,gs.cv_results_["mean_test_score"][i],gs.cv_results_["std_test_score"][i]," \n"))
        file.close()



    #############################

    ## training

    cls = xgb.XGBClassifier(
        n_estimators = options.ntrees,
        max_depth = options.treeDeph,
        min_child_weight = options.mcw,
        learning_rate = options.lr,
        random_state=options.randInt,
        n_jobs=1
        )
    eval_set = [(testData[trainVars], testData['target'])]

    cls.fit(trainData[trainVars], trainData['target'],
            sample_weight=(trainData['final_weights']),
            early_stopping_rounds=100, eval_metric="auc",
            eval_set = eval_set,
            sample_weight_eval_set=[testData['final_weights']],
            verbose=0)

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


    ## save model to pickle
    pklpath="Htoaa_BDT" + condition
    if options.doXML:
        pickle.dump(cls, open(pklpath+".pkl", 'wb'))
        file = open(pklpath+"pkl.log","w")
        file.write(str(trainVars)+"\n")
        file.close()
        print('made XML file')

    if disc == 'h':
        dest = 'plots/high'
    if disc == 'm':
        dest = 'plots/high+medium'
    else:
        dest = 'plots'



    ## draw them rocs
    fig, ax = plt.subplots(figsize=(8,8))
    train_auc = auc(fpr, tpr)
    ax.plot(fpr, tpr, lw=1, color='g',label='XGB train (area = %0.5f)'%(train_auc))
    ax.plot(fprt, tprt, lw=1, ls='--',color='g',label='XGB test (area = %0.5f)'%(test_auc) )
    ax.set_ylim([0.0,1.0])
    ax.set_xlim([0.0,1.0])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend(loc="lower right")
    ax.grid()
    ax.set_title(hyppar + ' ' + condition)
    fig.savefig("{}/roc_{} {}.png".format(dest, hyppar, condition))
    plt.show()
    plt.clf()


    ## add the BDT scores to the df so manip is easier. hopefully
    trainData = trainData.assign(BDTScore = train_proba[:,1])
    testData = testData.assign(BDTScore = test_proba[:,1])

    ## making bdt score figs babey
    density = True
    fig, ax = plt.subplots(figsize=(8,8))
    ax.hist(trainData.BDTScore.loc[trainData.target == 1], weights=trainData.final_weights.loc[trainData.target == 1],  bins=20, ls = '--', histtype='step',stacked=True, label='train signal', density=density)
    ax.hist(trainData.BDTScore.loc[trainData.target == 0], weights=trainData.final_weights.loc[trainData.target == 0], bins=20, ls = '--', histtype='step',stacked=True, label='train background', density=density)
    ax.hist(testData.BDTScore.loc[testData.target == 1], weights=testData.final_weights.loc[testData.target == 1], bins=20, histtype='step', stacked=True, label='test signal', fill=False, density=density)
    ax.hist(testData.BDTScore.loc[testData.target == 0], weights=testData.final_weights.loc[testData.target == 0], bins=20, histtype='step', stacked=True, label='test background', fill=False, density=density)
    ax.legend(loc='lower right')
    ax.set_title('BDT score {} {}'.format(hyppar, condition))
    ax.set_xlabel('BDT Score')
    fig.savefig("{}/BDTScore_{} {}.png".format(dest, hyppar, condition))
    plt.show()
    plt.clf()


    ################### sensitivity plots ####################
    sortedData = testData.loc[testData.target==1].sort_values(by='BDTScore', axis=0, kind='mergesort')
    histHeight = sortedData.final_weights.sum()/10
    edges = [0]
    cumuSum = 0
    edgesloc = []
    for i in range(sortedData.BDTScore.size):
        cumuSum = cumuSum + sortedData.final_weights.iloc[i]
        if cumuSum > histHeight:
            edges.append(sortedData.BDTScore.iloc[i-1])
            cumuSum = sortedData.final_weights.iloc[i]
            edgesloc.append(i-1)
    edges.pop()
    edges.append(sortedData.BDTScore.iloc[-1])

    #edges.append(sortedData.BDTScore.iloc[sortedData.BDTScore.size-1])
    fig, ax = plt.subplots(figsize=(8,8))
    sighist = ax.hist(sortedData.BDTScore, bins=edges,
                      weights=sortedData.final_weights, histtype='step',
                      label='signal', log=True)
    bghist = ax.hist(testData.loc[testData.target==0].BDTScore, bins=edges,
                     weights=testData.loc[testData.target==0].final_weights,
                     histtype='step', label='background', log=True)
    sighistval = sighist[0]/np.sum(sighist[0])
    bghistval = bghist[0]/np.sum(bghist[0])
    sensitivity = np.sqrt(np.sum(np.power(sighistval,2)/bghistval))
    ax.legend(loc='best')

    ax.set_title(f'sensitivity = {sensitivity:.4f}')
    fig.savefig("{}/Sensitivity_{} {}.png".format(dest, hyppar, condition))
    plt.show()
    plt.clf()

    ##########################################################




    ############ feature importance plot #############
    #fig, ax = plt.subplots()
    #f_score_dict = cls.get_booster().get_fscore()
    #print("f_score_dict: {}".format(f_score_dict))
    #f_score_dict = {trainVars[int(k[1:])] : v for k,v in f_score_dict.items()}
    #feat_imp = pd.Series(f_score_dict).sort_values(ascending=True)
    #feat_imp.plot(kind='barh', title='Feature Importances')
    #fig.tight_layout()
    #plt.show()
    #fig.savefig("plots/%s_XGB_importance.png" % hyppar)


    ## distributions of things yeah
    if options.dist:
        for colName in dataSig.columns:
            hist_params = {'density': True, 'histtype': 'bar', 'fill': True , 'lw':3, 'alpha' : 0.4}
            nbins = 40
            minval = 0
            maxval = 99.8
            min_valueS, max_valueS = np.percentile(dataSig[colName], [minval, maxval])
            min_valueB, max_valueB = np.percentile(dataBg[colName], [minval, maxval])
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
            plt.savefig("distributions/dist_{}.png".format(colName))
            plt.clf()

        ## lhe ht plots by bg type
        # trainData.BDTScore.loc[trainData.target == 1]
        hist_params = {'density': True, 'histtype': 'step', 'bins': 40}
        plt.hist(TTJetsData.LHE_HT, weights=TTJetsData.final_weights, label='TTBar', color='b', **hist_params)
        plt.hist(ZJetsData.LHE_HT, weights=ZJetsData.final_weights, label='ZJets', color='r', **hist_params)
        plt.hist(WJetsData.LHE_HT, weights=WJetsData.final_weights, label='WJets', color='g', **hist_params)
        plt.hist(BGenData.LHE_HT, weights=BGenData.LHE_weights, label='BGen', color='k', **hist_params)
        plt.hist(bEnrData.LHE_HT, weights=bEnrData.LHE_weights, label='bEnr', color='m', **hist_params)
        plt.legend(loc='best')
        plt.savefig("{}/LHE HT_{} {}.png".format(dest, hyppar, condition))
        plt.xlabel('LHE_HT')
        #plt.show()
        plt.clf()

else:
    print('we fc')
    ## this is if you want to use existing BDT to test
    ## or should i be doing that in datavs mc. no i should do it here too






# ## plotting signal with bdtscore > 0.5 to see what's up
# highBDTScore = testData[testData.BDTScore > 0.5]
# nbins = 20
# for var in allVars:
#     if 'pt' in var:
#         plt.hist2d(highBDTScore.BDTScore, highBDTScore[var], bins=nbins, range=[[0.5,1],[200,600]])
#     else:
#         plt.hist2d(highBDTScore.BDTScore, highBDTScore[var], bins=nbins)
#     plt.xlabel('BDTScore unweighted')
#     plt.ylabel(var)
#     plt.title(var)
#     plt.show()

#for colName in allVars:













