import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, accuracy_score
import random
import pickle


CSVpathBg = 'CSV/TTBarMuEleFat_Bkg.csv'
CSVpathSig = 'CSV/TTBarMuEleFat_Sig.csv'
 
test_size = 0.2

n_estimators = [400, 800, 1200, 1600]
max_depth = [2, 4, 6, 8]
min_child_weight = [2,4, 6,8, 10]
learning_rate = [0.0001, 0.001, 0.01]



params =  { 
    'n_estimators' : 400,
    'max_depth' : 2,
    'min_child_weight' : 2,
    'learning_rate' : 0.02,
    'objective' : 'binary:logistic', 
    'eval_metric' : 'auc',
    'reg_alpha' : 1,
    }
num_boost_round = 999
trainVars=['msoft_fat', 
           'mass_fat', 
           'deepB_max_jet',
           'flavB_max_jet', 
           'pt_llJ', 
           'pt_llvJ', 
           'MT_llvJ', 
           'mass_llJ',
           'mass_llvJ', 
           'mass_max_lep_fat', 
           'mass_min_lep_fat', 
           'dR_max_lep_fat',
           'dR_min_lep_fat', 
           'dEta_max_lep_fat', 
           'dEta_min_lep_fat', 
           'dR_ll_fat',
           'dEta_ll_fat']

dataBg = pd.read_csv(CSVpathBg, header=0)
dataBg.columns = dataBg.columns.str.replace(' ', '')
dataBg = dataBg.sample(frac=0.6)


dataSig = pd.read_csv(CSVpathSig, header=0)
dataSig.columns = dataSig.columns.str.replace(' ', '')
data = dataSig.append(dataBg, ignore_index=True)

print('Signal event count: ' + str(len(dataSig.index)))
print('Background event count: ' + str(len(dataBg.index)))


## split data into training and testing
randInt = 576#random.randint(0,1000) #576
print("random int: " + str(randInt))
trainData, testData = train_test_split(data, random_state=randInt, test_size=test_size)


diff = []
testauc = []
paramslist = []

# n_estimators = [400, 800, 1200, 1600]
# max_depth = [2, 4, 6, 8]
# min_child_weight = [2,4, 6,8, 10]
# learning_rate = [0.0001, 0.001, 0.01]



# for i in n_estimators:
#     for j in max_depth:
#         for k in min_child_weight:
#             for l in learning_rate:
cls = xgb.XGBClassifier(
    random_state=randInt,
    **params,
    )
eval_set = [(testData[trainVars], testData['isSignal'])]

cls.fit(trainData[trainVars], trainData['isSignal'],
           early_stopping_rounds=100, 
           eval_metric="auc",
           eval_set = eval_set,
           verbose=1)

'''dtrain = xgb.DMatrix(trainData[trainVars], label=trainData['isSignal'])
cv_results = xgb.cv(
    params,
    dtrain=dtrain,
    num_boost_round=num_boost_round,
    seed=42,
    nfold=5,
    metrics={'auc'},
    early_stopping_rounds=100,
)
print(cv_results)'''

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



plotdir = 'csvpickles/plots'
paramcap = f'ntrees{params["n_estimators"]}_depth{params["max_depth"]}_mcw{params["min_child_weight"]}_lr{params["learning_rate"]}_Int{randInt}'
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
ax.set_title(f'ROC {paramcap}')
#%%
#fig.savefig("{}/roc_{} {}.png".format(dest, hyppar, condition))
plt.savefig(f'{plotdir}/roc_{paramcap}.png')
plt.show()
plt.clf()

## add the BDT scores to the df so manip is easier. hopefully
trainData = trainData.assign(BDTScore = train_proba[:,1])
testData = testData.assign(BDTScore = test_proba[:,1])

## making bdt score figs babey
density = True
fig, ax = plt.subplots(figsize=(8,8))
bins = np.linspace(0,1,21)
ax.hist(trainData.BDTScore.loc[trainData.isSignal == 1],  bins=bins, ls = '--', histtype='step',stacked=True, label='train signal', density=density)
ax.hist(trainData.BDTScore.loc[trainData.isSignal == 0], bins=bins, ls = '--', histtype='step',stacked=True, label='train background', density=density)
ax.hist(testData.BDTScore.loc[testData.isSignal == 1], bins=bins, histtype='step', stacked=True, label='test signal', fill=False, density=density)
ax.hist(testData.BDTScore.loc[testData.isSignal == 0],  bins=bins, histtype='step', stacked=True, label='test background', fill=False, density=density)
ax.legend(loc='best')
ax.set_title('BDT score')
ax.set_xlabel('BDT Score')
plt.savefig(f'{plotdir}/bdtscore_{paramcap}.png')
#fig.savefig("{}/BDTScore_{} {}.png".format(dest, hyppar, condition))
plt.show()
plt.clf()

## feature importance 
xgb.plot_importance(cls)
plt.savefig(f'{plotdir}/featureImportance_{paramcap}.png')
plt.show()
       
## save model
pickle.dump(cls, open(f'csvpickles/model_ntrees{params["n_estimators"]}_depth{params["max_depth"]}_mcw{params["min_child_weight"]}_lr{params["learning_rate"]}','wb'))         
                
#                 testauc.append(test_auc)
#                 diff.append(np.abs(train_auc-test_auc))
#                 paramslist.append((i,j,k,l))
                

# results = pd.DataFrame(testauc, columns = ['testauc'])
# results = results.assign(diff=diff)
# results = results.assign(params=paramslist)

# #%%
# #results = pickle.load(open('csvpickls/results.pkl','rb'))
# plt.scatter(results['diff'], results['testauc'])
# plt.xlabel('|testauc - train auc|')
# plt.ylabel('test auc')
