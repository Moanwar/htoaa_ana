universe = vanilla 
arguments = $(dataset)
executable = condor_exec_vbf.sh
getenv = TRUE 
log = logs/analyis_pnet_$(dataset)_condor_skimmed.log 
output = logs/analyis_pnet_$(dataset)_condor_skimmed.out 
error = logs/analyis_pnet_$(dataset)_condor_skimmed.error 
transfer_input_files = htoaa_Analysis_dataOnly_VBF_v13.py, list_samples_2018.json, list_samples_2017.json, list_samples_Post2016.json, list_samples_Pre2016.json 
transfer_output_files = "" 
request_memory = 8GB
notification = never 
should_transfer_files = YES 
when_to_transfer_output = ON_EXIT 
+JobFlavour = "tomorrow" 

Queue dataset from DatasetList_skimmed.txt
