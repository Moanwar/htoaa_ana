# htoaa Analysis Framework using Coffea

## Seting environment using conda
1. Install 'conda' on your system/lxplus. Suggestion: Follow https://docs.anaconda.com/anaconda/install/linux/
2. Once conda is installed, set up conda environment named ' myCondaEnv' to install the required software libraries/packages using
   ```
   conda env create --name myCondaEnv -f environment_lxplus.yml
   ```
3. Activate conda envirnment every time you login into your system/lxplus
   ```
   conda env create --name myCondaEnv -f environment_lxplus.yml
   ```


## Running analysis
'htoaa_Analysis_wCoffea.py' is a central analysis macro to run the GGF htoaa analysis. It takes input and output file name and other settings in input through config.json file. Command to run on individual config.json file:
```
python3 htoaa_Analysis_wCoffea.py config.json
```

'htoaa_Wrapper.py' prepares the config.json files for differnt data/MC samples and run with htoaa_Analysis_wCoffea.py as parallel HT Condor jobs.
Command to run htoaa analysis macro on data and MC samples:
```
python3 htoaa_Wrapper.py -era <era> -run_mode condor -v <version name>
```


## Description of files
**htoaa_Settings.py**: File listing main analysis options.

**htoaa_CommonTools.py**: Common functions used to run the analysis

**htoaa_Samples.py**: List of data and MC samples to run for the analysis.

**Samples_2018UL.csv**: Dictionary of all data and MC samples with NanoAOD files, cross-section etc listed for each sample.

**prepare_SamplesList.py**: Prepare 'Samples_2018UL.csv' from needed samples using CMS-DAS (Data Aggregation System) tool.
