# htoaa Analysis Framework using Coffea

## Seting up environment using conda
1. Install 'conda' on your system/lxplus. Suggestion: Follow https://docs.anaconda.com/anaconda/install/linux/

2. Once conda is installed, set up conda environment named ' myCondaEnv' to install the required software libraries/packages using
   ```
   conda env create --name myCondaEnv  --file environment_lxplus.yml 
   ```
3. Activate conda envirnment every time you login into your system/lxplus
   ```
   conda activate myCondaEnv 
   ```


## Set proxy when you open a new terminal on lxplus or after every few days
```
voms-proxy-init --rfc --voms cms -valid 192:00 --out ~/x509up_u108989
```

## Analysis steps
Analysis is an iterative process, and different corrections are needed to take into account. The different corrections considered so far are listed below with details on obtaining those corrections. If you want to run the analysis for the first iteration or intermediate corrections are already calculated then you can jump to 'Running analysis' sub-section.

### Prepare samples' list with all details (Samples_2018UL.csv)
Add/update dataset name and cross-sections for all samples in 'samplesList_prepare.py'.\
To run:
```
python3 samplesList_prepare.py -era <era>
cp Samples_2018UL_v0.json Samples_2018UL.json
```
Add command line option '-updateCrossSections' to update the cross-section value in the existing sample list.


### Calculate sumEvents
Calculate sumEvents as the number of events with positive generator weight minus the number of events with negative weight.\ 
Run:
```
python3 htoaa_Wrapper.py -analyze countSumEventsInSample.py -era <era> -v <version name> 
```
It runs jobs in HT Condor, and produces final output hadd root file.
Set path of the output hadd root file in 'sFIn_analysis_stage1_dict' variable in 'samplesList_update.py' and run
```
python3 samplesList_update.py -era <era>
cp Samples_2018UL_updated.json Samples_2018UL.json
```
This updates samples' full detail list with sumEvents. Now samples' full detail list (Samples_2018UL.csv) is ready for the further analysis. Now on the previous analysis steps need not to run again.



### Calculate HT reweight SFs for QCD_bGenFilter sample
In scripts/CopyHistos.py, update 'sIpFile' to point to the latest analyze_htoaa_stage1.root file. \
Run to copy the required histograms into a seperate file stored in data/:
```
cd scripts/
python3 CopyHistos.py
```
Calculate HT reweights for QCD_bGen sample:
```
cd corrections/
python3 cal_HTRewght_QCDbGen.py -era <era>
```
Copy "HTRewgt" SFs printed on the terminal into "Corrections" variable in htoaa_Settings.py file.



### Calculate lumi-scale for MC samples stitch with reweighting overlapping phase space
Calculate lumi-scale for MC samples stitch with reweighting overlapping phase space, instead of removing overlapping phase space.\
In htoaa_Samples.py, set 'QCDInclMode = 2' to use QCD_Incl_madgraph-PSWeights sample for QCD_Incl.\
Run:
```
python3 calculateLumiscale.py -era <era>
```
 

### Calculate luminosity for triggers used in the analysis
HLT paths for Run2 can be found at https://twiki.cern.ch/twiki/bin/viewauth/CMS/HLTPathsRunIIList. \
Recommendations for Run2 luminosityis at https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2. However, luminosity for the triggers used in the analysis can be computed using brilcalc tool (https://twiki.cern.ch/twiki/bin/view/CMS/BrilcalcQuickStart) by running the following commands on lxplus. However, always refer https://twiki.cern.ch/twiki/bin/view/CMS/BrilcalcQuickStart for the latest recommendations. 
```
source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env

brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i <Golden JSON file> --hltpath "<your trigger paths>_v*" -o output.csv
```

Provide Golden JSON file recommended for the data era and HLT path of your interest. \
For e.g. luminosity for "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 path in 2018 data can be calculated using
```
brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_v*" -o output_brilcalc_314472-325175_UL18_HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_v.csv
```


### Running analysis
'htoaa_Analysis_Example.py' is a central analysis macro to run the GGF htoaa analysis. It takes input and output file name and other settings in input through config.json file. Command to run on individual config.json file:
```
python3 htoaa_Analysis_Example.py <config>.json
```
Histograms stored in the output root file have the following naming convention:
```
evt/<sampleCategory from <config>.json>/<histogram name>_<systematics>
```

'htoaa_Wrapper.py' prepares the config.json files for differnt data/MC samples and run with htoaa_Analysis_Example.py in parallel HT Condor jobs.
Command to run htoaa analysis macro on data and MC samples:
```
python3 htoaa_Wrapper.py -analyze htoaa_Analysis_Example.py -era <era> -run_mode condor -v <version name>  -xrdcpIpAftNResub 0
```
Append '-server tifr' to the previous command to run on TIFR server. \
Append '-ntuples SkimmedNanoAOD -nFilesPerJob 1' to run on the new skimmed NanoAOD files.


### Data and MC stack plots
Run 'scripts/PlotHistos1D_DataVsMC.ipynb' (in jupyter-notebook) to make Data-vs-MC comparison stack plots. A list of histograms for the stack plots are set in 'scripts/HistogramListForPlottingDataVsMC_Analysis_Example.py' file, which get imported inside 'scripts/PlotHistos1D_DataVsMC.ipynb' file.



## Description of files
**htoaa_Settings.py**: File listing main analysis options.

**htoaa_CommonTools.py**: Common functions used to run the analysis

**htoaa_Samples.py**: List of data and MC samples to run for the analysis. This list is stored information in python-dictionary format. Keys of the dictionary represent the 'sample category' and values of the dictionary list sample dataset physics name. 

**Samples_2018UL.csv**: Dictionary of all data and MC samples with NanoAOD files, cross-section etc listed for each sample.

**samplesList_prepare.py**: Prepare 'Samples_2018UL.csv' from needed samples using CMS-DAS (Data Aggregation System) tool.


## Persistent screen session on lxplus
Instructions from reference https://hsf-training.github.io/analysis-essentials/shell-extras/persistent-screen.html worked in early 2023, but somehow not working now.

Currently working instructions are from https://frankenthal.dev/post/ssh_kerberos_keytabs_macos/

### Setting up password-less kerberos token
Store encrypted passwork for lxplus in '~/keytab' file by running the following commands on lxplus terminal:
```
$ ktutil 
addent -password -p user@CERN.CH -k 3 -e arcfour-hmac-md5
(type your password)
wkt ~/keytab
quit
``` 
Replace 'user' with your lxplus username.

### Making use of the keytab

Use the keytab file to authenticate to kinit.

```
kinit -kt ~/keytab user@CERN.CH
```
If the last command runs without any error then passwordless kerberos token is generated successfully.

### Using k5reauth to automatically refresh your kerberos token
Start screen session using
```
k5reauth -f -i 3600 -p <user> -k  ~/keytab  -- screen -S <session name>
```
