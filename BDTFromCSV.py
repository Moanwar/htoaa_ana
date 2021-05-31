import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, accuracy_score
import random

CSVpath = '' 
test_size = 20
ntrees = 400
treeDepth = 3
mcw = 3
lr = 0.01
trainVars = [    
    'FatJet_pt',
    'FatJet_eta',
    'FatJet_mass',
    'FatJet_btagCSVV2',
    'FatJet_btagDeepB',
    'FatJet_msoftdrop',
    'FatJet_btagDDBvL',
    'FatJet_deepTagMD_H4qvsQCD',
    'FatJet_n2b1',
    'SubJet_mass(1)',
    'SubJet_mass(2)',
    'SubJet_tau1(1)',
    'FatJet_nSV'
    ]


data = pd.read_csv(CSVpath, header=0)
dataSig = data.loc[data.isSig == 1]
dataBg = data.loc[data.isSig == 0]

print('Signal event count: ' + str(len(dataSig.index)))
print('Background event count: ' + str(len(dataBg.index)))


## split data into training and testing
# randInt = random.randint(0,100)
randInt = 1 #random.randint(0,1000)
print("random int: " + str(randInt))
trainData, testData = train_test_split(data, random_state=randInt, test_size=test_size)

cls = xgb.XGBClassifier(
    n_estimators = ntrees,
    max_depth = treeDepth,
    min_child_weight = mcw,
    learning_rate = lr,
    random_state=randInt,
    n_jobs=1
    )
eval_set = [(testData[trainVars], testData['isSig'])]

cls.fit(trainData[trainVars], trainData['isSig'],
           sample_weight=(trainData['final_weights']),
           early_stopping_rounds=100, eval_metric="auc",
           eval_set = eval_set,
           sample_weight_eval_set=[testData['final_weights']],
           verbose=0)

print ("XGBoost trained")


## get info for ROC
train_proba = cls.predict_proba(trainData[trainVars])
fpr, tpr, thresholds = roc_curve(trainData['isSig'], train_proba[:,1])
train_auc = auc(fpr, tpr)
print("XGBoost train set auc - {}".format(train_auc))

test_proba = cls.predict_proba(testData[trainVars])
fprt, tprt, thresholds = roc_curve(testData['isSig'], test_proba[:,1])
test_auc = auc(fprt, tprt)
print("XGBoost test set auc - {}".format(test_auc))
fig, ax = plt.subplots()

prediction = cls.predict(testData[trainVars])
accuracy = accuracy_score(testData['isSig'], prediction)
print("XGBoost test accuracy - {}".format(accuracy))



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
ax.set_title('ROC')
#fig.savefig("{}/roc_{} {}.png".format(dest, hyppar, condition))
plt.show()
plt.clf()

## add the BDT scores to the df so manip is easier. hopefully
trainData = trainData.assign(BDTScore = train_proba[:,1])
testData = testData.assign(BDTScore = test_proba[:,1])

## making bdt score figs babey
density = True
fig, ax = plt.subplots(figsize=(8,8))
ax.hist(trainData.BDTScore.loc[trainData.isSig == 1], weights=trainData.final_weights.loc[trainData.isSig == 1],  bins=20, ls = '--', histtype='step',stacked=True, label='train signal', density=density)
ax.hist(trainData.BDTScore.loc[trainData.isSig == 0], weights=trainData.final_weights.loc[trainData.isSig == 0], bins=20, ls = '--', histtype='step',stacked=True, label='train background', density=density)
ax.hist(testData.BDTScore.loc[testData.isSig == 1], weights=testData.final_weights.loc[testData.isSig == 1], bins=20, histtype='step', stacked=True, label='test signal', fill=False, density=density)
ax.hist(testData.BDTScore.loc[testData.isSig == 0], weights=testData.final_weights.loc[testData.isSig == 0], bins=20, histtype='step', stacked=True, label='test background', fill=False, density=density)
ax.legend(loc='lower right')
ax.set_title('BDT score')
ax.set_xlabel('BDT Score')
#fig.savefig("{}/BDTScore_{} {}.png".format(dest, hyppar, condition))
plt.show()
plt.clf()