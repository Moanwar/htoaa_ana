# htoaa BDT

## Running Trained BDT

`htoaa_predict.py` will create plots of BDTScore + distribution from the root file it analyses. Variables that get turned into distribution plots can be found in `info.py` under the variable trainVars. Ability to process multiple root files at once will be added in the future.

**required packages**: xgboost, pandas, numpy, munch, uproot, sklearn

1. clone the repo
2. inside the main htoaa folder, run `mkdir loadedModel`
3. to use the trained BDT to predict another root file, (\<path to rootfile\> must not include '.root')
  
   ```python htoaa_predict.py -f <path to rootfile>``` 



## Description of files
**analib.py**: library used to make cuts on data and propagate the cuts for all events

**data_manager.py**: contains function that processes the data from root file into form that can be fed into the xgboost BDT

**info.py**: contains 1) names of signal/background files. 2) corresponding weights. 3) variables used to train, cut data, and associated values for cutting

**htoaa_BDT2.py**: trains the BDT. Makes roc plots, distribution plots, and BDT score plots. Need to make folders called 'plots', and 'distributions' before running. 

**htoaa_predict.py**: run this to use the trained model to predict data. Requires the path to root file as an argument

**XGB_classifier_8Var.pkl**: the pickle file of the Trained BDT. This one corresponds to ntrees = 1000, depth = 2, lr = 0.05, min child weight = 1

