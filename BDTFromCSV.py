import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, accuracy_score
import random

CSVpathBg = 'CSV/TTBarMuEleFat_Bkg.csv'
CSVpathSig = 'CSV/TTBarMuEleFat_Sig.csv'
 
test_size = 0.2
ntrees = 400
treeDepth = 3
mcw = 3
lr = 0.01
# trainVars = [    
#     'FatJet_pt',
#     'FatJet_eta',
#     'FatJet_mass',
#     'FatJet_btagCSVV2',
#     'FatJet_btagDeepB',
#     'FatJet_msoftdrop',
#     'FatJet_btagDDBvL',
#     'FatJet_deepTagMD_H4qvsQCD',
#     'FatJet_n2b1',
#     'SubJet_mass(1)',
#     'SubJet_mass(2)',
#     'SubJet_tau1(1)',
#     'FatJet_nSV'
#     ]
trainVars=['msoft_fat', 'mass_fat', 'deepB_max_jet',
       'flavB_max_jet', 'pt_llJ', 'pt_llvJ', 'MT_llvJ', 'mass_llJ',
       'mass_llvJ', 'mass_max_lep_fat', 'mass_min_lep_fat', 'dR_max_lep_fat',
       'dR_min_lep_fat', 'dEta_max_lep_fat', 'dEta_min_lep_fat', 'dR_ll_fat',
       'dEta_ll_fat']


# data = pd.read_csv(CSVpathbg, header=0)
# dataSig = data.loc[data.isSignal == 1]
# dataBg = data.loc[data.isSignal == 0]
colunstripped = [' isSignal', ' msoft_fat', ' mass_fat', ' deepB_max_jet', ' flavB_max_jet', 
        ' pt_llJ',' pt_llvJ', ' MT_llvJ', ' mass_llJ', ' mass_llvJ', ' mass_max_lep_fat', 
        ' mass_min_lep_fat', ' dR_max_lep_fat', ' dR_min_lep_fat',  ' dEta_max_lep_fat',
        ' dEta_min_lep_fat', ' dR_ll_fat', ' dEta_ll_fat']
colstripped = [x.strip() for x in colunstripped]
cols = dict(zip(colunstripped, colstripped))
dataBg = pd.read_csv(CSVpathBg, header=0)
dataBg.rename(columns=cols, inplace=True)
dataBg = dataBg.sample(frac=0.6)


dataSig = pd.read_csv(CSVpathSig, header=0)
dataSig.rename(columns=cols, inplace=True)
data = dataSig.append(dataBg, ignore_index=True)
#%%
print('Signal event count: ' + str(len(dataSig.index)))
print('Background event count: ' + str(len(dataBg.index)))


## split data into training and testing
# randInt = random.randint(0,100)
randInt = random.randint(0,1000)
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
eval_set = [(testData[trainVars], testData['isSignal'])]

cls.fit(trainData[trainVars], trainData['isSignal'],
           #sample_weight=(trainData['final_weights']),
           #early_stopping_rounds=100, 
           eval_metric="auc",
           eval_set = eval_set,
           #sample_weight_eval_set=[testData['final_weights']],
           verbose=0)

print ("XGBoost trained")


## get info for ROC
train_proba = cls.predict_proba(trainData[trainVars])
fpr, tpr, thresholds = roc_curve(trainData['isSignal'], train_proba[:,1])
train_auc = auc(fpr, tpr)
print("XGBoost train set auc - {}".format(train_auc))

test_proba = cls.predict_proba(testData[trainVars])
fprt, tprt, thresholds = roc_curve(testData['isSignal'], test_proba[:,1])
test_auc = auc(fprt, tprt)
print("XGBoost test set auc - {}".format(test_auc))
fig, ax = plt.subplots()

print(f"diff - {test_auc-train_auc}")


prediction = cls.predict(testData[trainVars])
accuracy = accuracy_score(testData['isSignal'], prediction)
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
ax.hist(trainData.BDTScore.loc[trainData.isSignal == 1],  bins=20, ls = '--', histtype='step',stacked=True, label='train signal', density=density)
ax.hist(trainData.BDTScore.loc[trainData.isSignal == 0], bins=20, ls = '--', histtype='step',stacked=True, label='train background', density=density)
ax.hist(testData.BDTScore.loc[testData.isSignal == 1], bins=20, histtype='step', stacked=True, label='test signal', fill=False, density=density)
ax.hist(testData.BDTScore.loc[testData.isSignal == 0],  bins=20, histtype='step', stacked=True, label='test background', fill=False, density=density)
ax.legend(loc='best')
ax.set_title('BDT score')
ax.set_xlabel('BDT Score')
#fig.savefig("{}/BDTScore_{} {}.png".format(dest, hyppar, condition))
plt.show()
plt.clf()