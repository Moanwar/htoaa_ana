# htoaa

analib.py: library used to make cuts on data and propagate the cuts for all events

data_manager.py: contains function that processes the data from root file into form that can be fed into the xgboost BDT

info.py: contains 1) names of signal/background files. 2) corresponding weights. 3) variables used to train, cut data, and associated values for cutting

htoaa_BDT2.py: the BDT. Makes roc plots, distribution plots, and importance plots for each variable. Need to make folders called 'plots', and 'distributions' before running. 

htoaa_predict.py: run this to use the trained model to predict data. Change file name to your data file before running. 

XGB_classifier_8Var.pkl: the pickle file of the Trained BDT. This one corresponds to ntrees = 1000, depth = 2, lr = 0.05, min child weight = 1

To run the BDT pickle: 
1. must have these packages installed: xgboost, pandas, numpy, munch, uproot, sklearn
2. go into htoaa_predict.py and change `'GGH_HPT'` to the name of the datafile you want to test (without the .root extension)
3. run `python htoaa_predict.py`
