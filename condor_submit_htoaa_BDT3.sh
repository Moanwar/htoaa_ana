universe = vanilla 
executable = condor_exec_htoaa_BDT3.sh 
getenv = TRUE 
log = htoaa_BDT3.log 
output = htoaa_BDT3.out 
error = htoaa_BDT3.error 
notification = never 
should_transfer_files = YES 
when_to_transfer_output = ON_EXIT 
+JobFlavour = "espresso" 
queue 
