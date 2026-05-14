#! /usr/bin/env python
## Create ROOT TTrees with VBF events passing X4b selections

import os
import sys
import json
import argparse
import urllib.request
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)
R.gErrorIgnoreLevel = R.kWarning

# Default year if not found in dataset name
DEFAULT_YEAR = '2018'

# JSON file mapping for input files
JSON_FILES = {
    '2016preVFP': 'list_samples_Pre2016.json',
    '2016postVFP': 'list_samples_Post2016.json',
    '2017': 'list_samples_2017.json',
    '2018': 'list_samples_2018.json',
}

# Golden JSON URLs
GOLDEN_JSON_URLS = {
    '2016preVFP': 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt',
    '2016postVFP': 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt',
    '2017': 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',
    '2018': 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
}

# Define VBF categories
VBF_CATEGORIES = ['VBFjjLoPtLo', 'VBFjjLoPtHi', 'VBFjjHiPtLo', 'VBFjjHiPtHi']
X4B_SELECTIONS = ['X4bSR', 'X4bSB', 'X4bSB_only']


def extract_year_from_dataset(dataset_name):
    """Extract year from dataset name. Returns DEFAULT_YEAR if not found."""
    if '2016' in dataset_name:
        if 'HIPM' in dataset_name or 'Pre' in dataset_name or 'pre' in dataset_name:
            return '2016preVFP'
        else:
            return '2016postVFP'
    elif '2017' in dataset_name:
        return '2017'
    elif '2018' in dataset_name:
        return '2018'
    else:
        print(f"Warning: Could not determine year from '{dataset_name}'. Using default {DEFAULT_YEAR}")
        return DEFAULT_YEAR


def load_golden_json(era):
    """Load golden JSON for a given era from CERN. Exits if download fails."""
    url = GOLDEN_JSON_URLS.get(era)
    if not url:
        print(f"Error: No golden JSON URL defined for era {era}")
        sys.exit(1)
    
    print(f"Loading golden JSON from: {url}")
    try:
        with urllib.request.urlopen(url) as response:
            content = response.read().decode('utf-8')
            # Parse JSON format (run: [[lumi_start, lumi_end], ...])
            golden_json = json.loads(content)
            print(f"✓ Golden JSON loaded successfully for {era}")
            return golden_json
    except Exception as e:
        print(f"Error: Failed to load golden JSON for {era}")
        print(f"Details: {e}")
        sys.exit(1)


def is_good_lumi(golden_json, run, lumi):
    """Check if a given run/lumi is in the golden JSON."""
    if not golden_json:
        return True  # Should not happen as we exit on failure
    
    run_str = str(run)
    if run_str not in golden_json:
        return False
    
    for lumi_range in golden_json[run_str]:
        if lumi_range[0] <= lumi <= lumi_range[1]:
            return True
    
    return False


def get_file_list(dataset_name, year):
    """Get list of ROOT files for a dataset from the appropriate JSON file."""
    json_file = JSON_FILES.get(year)
    if not json_file:
        print(f"Error: No JSON file defined for year {year}")
        return []
    
    if not os.path.exists(json_file):
        print(f"Error: JSON file {json_file} not found in current directory")
        return []
    
    with open(json_file, 'r') as f:
        filesets = json.load(f)
    
    if dataset_name not in filesets:
        print(f"Error: Dataset '{dataset_name}' not found in {json_file}")
        print(f"Available datasets: {list(filesets.keys())}")
        return []
    
    return filesets[dataset_name]


def extract_primary_dataset_from_path(input_files):
    """Extract primary dataset type from file path. Returns 'JetHT', 'BTagCSV', or None."""
    for file_path in input_files:
        if 'JetHT' in file_path:
            return 'JetHT'
        elif 'BTagCSV' in file_path:
            return 'BTagCSV'
        elif 'SingleMuon' in file_path:
            return 'SingleMuon'
        elif 'SingleElectron' in file_path:
            return 'SingleElectron'
        elif 'MET' in file_path:
            return 'MET'
    return None


def check_triggers_with_priority(chain, year, current_dataset_type, verbose=False):
    """
    Check triggers with priority handling to avoid double-counting across datasets.
    Checks all triggers from original code, gracefully handles missing branches.
    
    Priority order: JetHT > BTagCSV
    """
    
    # Define trigger bits for each primary dataset (full HLT names as they appear in file)
    # Keep ALL triggers from original code for each era
    if '2018' in year:
        triggers = {
            'JetHT': [
                'HLT_AK8PFJet500', 'HLT_PFHT1050', 'HLT_PFJet500',
                'HLT_AK8PFHT800_TrimMass50', 'HLT_AK8PFJet400_TrimMass30'
            ],
            'BTagCSV': [
                'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4',
                'HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71',
                'HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2',
                'HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1',
                'HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5'
            ]
        }
    elif '2017' in year:
        triggers = {
            'JetHT': [
                'HLT_AK8PFJet500', 'HLT_PFHT1050', 'HLT_PFJet500',
                'HLT_AK8PFHT800_TrimMass50', 'HLT_AK8PFJet400_TrimMass30',
                'HLT_AK8PFJet360_TrimMass30'
            ],
            'BTagCSV': [
                'HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2',
                'HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2',
                'HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5',
                'HLT_AK8PFHT750_TrimMass50', 'HLT_AK8PFJet380_TrimMass30',
                'HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33',
                'HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0'
            ]
        }
    elif '2016' in year:
        triggers = {
            'JetHT': [
                'HLT_PFJet450', 'HLT_DiCentralPFJet430',
                'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5', 'HLT_PFHT750_4JetPt50',
                'HLT_PFHT800', 'HLT_PFHT900', 'HLT_AK8PFJet360_TrimMass30',
                'HLT_AK8PFJet450', 'HLT_AK8PFHT650_TrimR0p1PT0p03Mass50',
                'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50'
            ],
            'BTagCSV': [
                'HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20',
                'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20',
                'HLT_PFHT400_SixJet30_DoubleBTagCSV_p056',
                'HLT_PFHT450_SixJet40_BTagCSV_p056',
                'HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20',
                'HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6',
                'HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160',
                'HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6',
                'HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172',
                'HLT_DoubleJet90_Double30_TripleBTagCSV_p08',
                'HLT_QuadJet45_TripleBTagCSV_p087',
                'HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460',
                'HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500',
                'HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200',
                'HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240'
            ]
        }
    else:
        return False, {'keep': False, 'reason': f'Unknown year: {year}'}
    
    # Check which triggers fired (gracefully handle missing branches)
    fired = {'JetHT': False, 'BTagCSV': False}
    
    for pd_type in ['JetHT', 'BTagCSV']:
        for hlt_name in triggers[pd_type]:
            # Check if branch exists before accessing
            if hasattr(chain, hlt_name):
                if getattr(chain, hlt_name):
                    fired[pd_type] = True
                    if verbose:
                        print(f"  Fired {pd_type} trigger: {hlt_name}")
                    break  # One trigger is enough for this PD type
            # If branch doesn't exist, skip it (it's not available in this dataset)
    
    # Priority logic: JetHT has priority over BTagCSV
    if fired['JetHT']:
        if current_dataset_type == 'JetHT':
            return True, {'keep': True, 'priority': 'JetHT', 'fired': fired}
        else:
            return False, {'keep': False, 'priority': 'JetHT', 'fired': fired, 'reason': 'Event belongs to higher priority JetHT'}
    
    elif fired['BTagCSV']:
        if current_dataset_type == 'BTagCSV':
            return True, {'keep': True, 'priority': 'BTagCSV', 'fired': fired}
        else:
            return False, {'keep': False, 'priority': 'BTagCSV', 'fired': fired, 'reason': 'Event only fires BTagCSV, not JetHT'}
    
    else:
        # No relevant triggers fired
        return False, {'keep': False, 'priority': None, 'fired': fired, 'reason': 'No triggers fired'}
    
def setup_output_tree(output_file, tree_name):
    """Create a TTree with all branches."""
    tree = R.TTree(tree_name, tree_name)
    
    # Event identification branches
    run = R.vector('int')()
    luminosityBlock = R.vector('unsigned int')()
    event = R.vector('unsigned long long')()
    
    # Fat jet branches
    fatjet_pt = R.vector('float')()
    fatjet_eta = R.vector('float')()
    fatjet_phi = R.vector('float')()
    fatjet_mass = R.vector('float')()
    
    # PNet mass branches
    pnet_massH_v2b = R.vector('float')()
    pnet_34massAa = R.vector('float')()
    
    # VBF jet 1 branches
    vbfjet1_pt = R.vector('float')()
    vbfjet1_eta = R.vector('float')()
    vbfjet1_phi = R.vector('float')()
    vbfjet1_mass = R.vector('float')()
    
    # VBF jet 2 branches
    vbfjet2_pt = R.vector('float')()
    vbfjet2_eta = R.vector('float')()
    vbfjet2_phi = R.vector('float')()
    vbfjet2_mass = R.vector('float')()
    
    # VBF dijet variables
    vbfjj_deta = R.vector('float')()
    vbfjj_mass = R.vector('float')()
    
    # X4b score
    x4b_score = R.vector('float')()
    
    # Trigger info
    trigger_fired = R.vector('int')()
    trigger_priority = R.vector('int')()
    
    # Branch setup
    tree.Branch('run', run)
    tree.Branch('luminosityBlock', luminosityBlock)
    tree.Branch('event', event)
    tree.Branch('fatjet_pt', fatjet_pt)
    tree.Branch('fatjet_eta', fatjet_eta)
    tree.Branch('fatjet_phi', fatjet_phi)
    tree.Branch('fatjet_mass', fatjet_mass)
    tree.Branch('vbfjet1_pt', vbfjet1_pt)
    tree.Branch('vbfjet1_eta', vbfjet1_eta)
    tree.Branch('vbfjet1_phi', vbfjet1_phi)
    tree.Branch('vbfjet1_mass', vbfjet1_mass)
    tree.Branch('vbfjet2_pt', vbfjet2_pt)
    tree.Branch('vbfjet2_eta', vbfjet2_eta)
    tree.Branch('vbfjet2_phi', vbfjet2_phi)
    tree.Branch('vbfjet2_mass', vbfjet2_mass)
    tree.Branch('vbfjj_deta', vbfjj_deta)
    tree.Branch('vbfjj_mass', vbfjj_mass)
    tree.Branch('x4b_score', x4b_score)
    tree.Branch('pnet_massH_v2b', pnet_massH_v2b)
    tree.Branch('pnet_34massAa', pnet_34massAa)
    tree.Branch('trigger_fired', trigger_fired)
    tree.Branch('trigger_priority', trigger_priority)
    
    vectors = {
        'run': run, 'luminosityBlock': luminosityBlock, 'event': event,
        'fatjet_pt': fatjet_pt, 'fatjet_eta': fatjet_eta, 'fatjet_phi': fatjet_phi, 'fatjet_mass': fatjet_mass,
        'vbfjet1_pt': vbfjet1_pt, 'vbfjet1_eta': vbfjet1_eta, 'vbfjet1_phi': vbfjet1_phi, 'vbfjet1_mass': vbfjet1_mass,
        'vbfjet2_pt': vbfjet2_pt, 'vbfjet2_eta': vbfjet2_eta, 'vbfjet2_phi': vbfjet2_phi, 'vbfjet2_mass': vbfjet2_mass,
        'vbfjj_deta': vbfjj_deta, 'vbfjj_mass': vbfjj_mass, 'x4b_score': x4b_score,
        'pnet_massH_v2b': pnet_massH_v2b, 'pnet_34massAa': pnet_34massAa,
        'trigger_fired': trigger_fired, 'trigger_priority': trigger_priority
    }
    
    return tree, vectors


def fill_tree(tree, vectors, values):
    """Fill a tree with values."""
    for key, val in values.items():
        vec = vectors[key]
        vec.clear()
        vec.push_back(val)
    tree.Fill()


def process_dataset(dataset_name, output_name=None):
    """Process a single dataset and create output ROOT file."""
    
    # Extract year from dataset name
    year = extract_year_from_dataset(dataset_name)
    print(f"\n{'='*70}")
    print(f"Processing: {dataset_name}")
    print(f"Year: {year}")
    
    # Load golden JSON for this era (for data)
    print("Loading golden JSON certification...")
    golden_json = load_golden_json(year)
    
    # Get list of input files
    input_files = get_file_list(dataset_name, year)
    if not input_files:
        print(f"Error: No input files found for {dataset_name}")
        return False
    
    print(f"Found {len(input_files)} input files")
    
    # Determine primary dataset type from input files
    primary_dataset = extract_primary_dataset_from_path(input_files)
    if primary_dataset not in ['JetHT', 'BTagCSV']:
        print(f"Warning: Unknown primary dataset type '{primary_dataset}' for {dataset_name}")
        print("Trigger priority handling will be skipped - using simple OR logic")
        use_priority = False
    else:
        use_priority = True
        print(f"Primary dataset: {primary_dataset}")
        print(f"Using priority handling: JetHT > BTagCSV")
    
    # Create output ROOT file
    if output_name is None:
        output_name = f"{dataset_name}_output.root"
    
    output_file = R.TFile(output_name, "RECREATE")
    print(f"Output file: {output_name}")
    
    # Create trees for each VBF category and X4b selection
    trees = {}
    vectors_dict = {}
    
    for vbf_cat in VBF_CATEGORIES:
        for x4b_sel in X4B_SELECTIONS:
            tree_name = f"{vbf_cat}_{x4b_sel}"
            tree, vectors = setup_output_tree(output_file, tree_name)
            trees[tree_name] = tree
            vectors_dict[tree_name] = vectors
    
    # Build TChain
    chain = R.TChain("Events")
    for f in input_files:
        if os.path.exists(f):
            chain.Add(f)
            print(f"  Added: {os.path.basename(f)}")
        else:
            print(f"  Warning: File not found: {f}")
    
    nEntries = chain.GetEntries()
    print(f"\nTotal events to process: {nEntries}")
    
    # Counters
    counters = {
        'total': 0,
        'pass_golden': 0,
        'pass_trigger': 0,
        'pass_fatH': 0,
        'pass_VBF': 0,
    }
    
    for vbf_cat in VBF_CATEGORIES:
        for x4b_sel in X4B_SELECTIONS:
            counters[f'{vbf_cat}_{x4b_sel}'] = 0
    
    # Main event loop
    for iEvt in range(nEntries):
        if iEvt % 10000 == 0:
            print(f"Processing event {iEvt} / {nEntries}")
        
        chain.GetEntry(iEvt)
        counters['total'] += 1
        
        # 1. Golden JSON filter (data only)
        run_val = getattr(chain, 'run', 0)
        ls_val = getattr(chain, 'luminosityBlock', 0)
        if not is_good_lumi(golden_json, run_val, ls_val):
            continue
        counters['pass_golden'] += 1
        
        # 2. nPV
        if chain.PV_npvsGood <= 0:
            continue
        
        # 3. Noise filters
        noise = (chain.Flag_goodVertices and chain.Flag_globalSuperTightHalo2016Filter and
                 chain.Flag_HBHENoiseFilter and chain.Flag_HBHENoiseIsoFilter and
                 chain.Flag_eeBadScFilter and chain.Flag_BadPFMuonFilter and
                 chain.Flag_BadPFMuonDzFilter and chain.Flag_EcalDeadCellTriggerPrimitiveFilter)
        if not noise:
            continue
        
        # 4. Trigger requirement with priority handling
        if use_priority:
            trig_pass, trig_info = check_triggers_with_priority(chain, year, primary_dataset, verbose=False)
            if not trig_pass:
                continue
            counters['pass_trigger'] += 1
            trigger_fired_code = 1
            trigger_priority_code = 1 if trig_info.get('priority') == 'JetHT' else 2
        else:
            trig_pass = False
            if hasattr(chain, 'HLT_AK8PFJet500') and chain.HLT_AK8PFJet500:
                trig_pass = True
            elif hasattr(chain, 'HLT_PFHT1050') and chain.HLT_PFHT1050:
                trig_pass = True
            elif hasattr(chain, 'HLT_PFJet500') and chain.HLT_PFJet500:
                trig_pass = True
            
            if not trig_pass:
                continue
            counters['pass_trigger'] += 1
            trigger_fired_code = 1
            trigger_priority_code = 0
        
        # 5. Find leading Higgs candidate (AK8 jet)
        xFatH = -1
        xFatH_X4b = -999.0
        xFatTops = []
        xFatWZs = []
        
        nFatJet = getattr(chain, 'nFatJet', 0)
        for iFat in range(nFatJet):
            fatjet_pt = getattr(chain, 'FatJet_pt', [0])[iFat] if hasattr(chain, 'FatJet_pt') else 0
            if fatjet_pt <= 250: continue
            
            fatjet_eta = getattr(chain, 'FatJet_eta', [0])[iFat] if hasattr(chain, 'FatJet_eta') else 0
            if abs(fatjet_eta) >= 2.4: continue
            
            fatjet_jetId = getattr(chain, 'FatJet_jetId', [0])[iFat] if hasattr(chain, 'FatJet_jetId') else 0
            if fatjet_jetId != 6: continue
            
            fatjet_msoftdrop = getattr(chain, 'FatJet_msoftdrop', [0])[iFat] if hasattr(chain, 'FatJet_msoftdrop') else 0
            if fatjet_msoftdrop <= 20: continue
            
            # Top tagging (for veto)
            if hasattr(chain, 'FatJet_particleNet_TvsQCD'):
                top_score = chain.FatJet_particleNet_TvsQCD[iFat]
                if ('2018' in year and top_score > 0.970) or \
                   ('2017' in year and top_score > 0.970) or \
                   ('2016' in year and top_score > 0.957):
                    #if fatjet_pt > 300:
                    xFatTops.append(iFat)
            
            # W/Z tagging (for veto)
            if hasattr(chain, 'FatJet_particleNet_WZvsQCD'):
                wz_score = chain.FatJet_particleNet_WZvsQCD[iFat]
                if ('2018' in year and wz_score > 0.9873) or \
                   ('2017' in year and wz_score > 0.9858) or \
                   ('2016' in year and wz_score > 0.9843):
                    xFatWZs.append(iFat)
            
            # Get X4b score for Higgs
            fatjet_Xbb = getattr(chain, 'FatJet_particleNetMD_XbbvsQCD', [0])[iFat] if hasattr(chain, 'FatJet_particleNetMD_XbbvsQCD') else 0
            if fatjet_Xbb <= 0.75: continue
            
            # Try different possible branch names for X4b score
            x4b_score = -999
            if hasattr(chain, 'FatJet_PNet_X4b_v2a_Haa4b_score'):
                score_a = getattr(chain, 'FatJet_PNet_X4b_v2a_Haa4b_score', [0])[iFat]
                score_b = getattr(chain, 'FatJet_PNet_X4b_v2b_Haa4b_score', [0])[iFat]
                x4b_score = 0.5 * (score_a + score_b)
            elif hasattr(chain, 'FatJet_X4b_score'):
                x4b_score = getattr(chain, 'FatJet_X4b_score', [0])[iFat]
            
            if x4b_score > xFatH_X4b:
                xFatH = iFat
                xFatH_X4b = x4b_score
        
        if xFatH < 0:
            continue
        
        # Get PNet mass values
        pnet_massH_v2b_val = -999.0
        pnet_34massAa_val = -999.0
        if hasattr(chain, 'FatJet_PNet_massH_v2b'):
            pnet_massH_v2b_val = getattr(chain, 'FatJet_PNet_massH_v2b', [0])[xFatH]
        if hasattr(chain, 'FatJet_PNet_34massAa'):
            pnet_34massAa_val = getattr(chain, 'FatJet_PNet_34massAa', [0])[xFatH]
        
        counters['pass_fatH'] += 1
        
        # 6. Top and WZ veto
        xFatTop = -1
        xFatWZ = -1
        if len(xFatTops) > 0:
            xFatTop = xFatTops[0] if xFatTops[0] != xFatH else (xFatTops[1] if len(xFatTops) > 1 else -1)
        if len(xFatWZs) > 0:
            xFatWZ = xFatWZs[0] if xFatWZs[0] != xFatH else (xFatWZs[1] if len(xFatWZs) > 1 else -1)
        
        if xFatTop >= 0 or xFatWZ >= 0:
            continue
        
        # 7. Muon veto
        hasMuon = False
        nMuon = getattr(chain, 'nMuon', 0)
        for iMu in range(nMuon):
            mu_pt = getattr(chain, 'Muon_pt', [0])[iMu] if hasattr(chain, 'Muon_pt') else 0
            if mu_pt <= 26.0: continue
            mu_eta = getattr(chain, 'Muon_eta', [0])[iMu] if hasattr(chain, 'Muon_eta') else 0
            if abs(mu_eta) >= 2.4: continue
            hasMuon = True
            break
        if hasMuon:
            continue
        
        # 8. Electron veto
        hasElectron = False
        nElectron = getattr(chain, 'nElectron', 0)
        for iEle in range(nElectron):
            ele_pt = getattr(chain, 'Electron_pt', [0])[iEle] if hasattr(chain, 'Electron_pt') else 0
            if ele_pt <= 30.0: continue
            ele_eta = getattr(chain, 'Electron_eta', [0])[iEle] if hasattr(chain, 'Electron_eta') else 0
            if abs(ele_eta) >= 2.5: continue
            hasElectron = True
            break
        if hasElectron:
            continue
        
        # 9. MET veto
        met_pt = getattr(chain, 'MET_pt', 0)
        if met_pt > 200:
            continue
        
        # 10. Find VBF jets (light-flavor AK4 jets)
        vFatH_pt = getattr(chain, 'FatJet_pt', [0])[xFatH] if hasattr(chain, 'FatJet_pt') else 0
        vFatH_eta = getattr(chain, 'FatJet_eta', [0])[xFatH] if hasattr(chain, 'FatJet_eta') else 0
        vFatH_phi = getattr(chain, 'FatJet_phi', [0])[xFatH] if hasattr(chain, 'FatJet_phi') else 0
        vFatH_mass = getattr(chain, 'FatJet_mass', [0])[xFatH] if hasattr(chain, 'FatJet_mass') else 0
        
        vFatH_vec = R.TLorentzVector()
        vFatH_vec.SetPtEtaPhiM(vFatH_pt, vFatH_eta, vFatH_phi, vFatH_mass)
        
        btagWPM = 0.2783 if '2018' in year else (0.3040 if '2017' in year else 0.2489)
        vJetsLF = []
        hasBJet = False
        
        nJet = getattr(chain, 'nJet', 0)
        for iJet in range(nJet):
            jet_pt = getattr(chain, 'Jet_pt', [0])[iJet] if hasattr(chain, 'Jet_pt') else 0
            if jet_pt <= 30: continue
            
            jet_eta = getattr(chain, 'Jet_eta', [0])[iJet] if hasattr(chain, 'Jet_eta') else 0
            jet_phi = getattr(chain, 'Jet_phi', [0])[iJet] if hasattr(chain, 'Jet_phi') else 0
            jet_mass = getattr(chain, 'Jet_mass', [0])[iJet] if hasattr(chain, 'Jet_mass') else 0
            
            vJet = R.TLorentzVector()
            vJet.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_mass)
            
            if vJet.DeltaR(vFatH_vec) <= 0.8: continue
            
            # Check if b-tagged
            jet_btag = getattr(chain, 'Jet_btagDeepFlavB', [0])[iJet] if hasattr(chain, 'Jet_btagDeepFlavB') else 0
            isBJet = (abs(jet_eta) < 2.4 and jet_btag > btagWPM)
            if isBJet:
                hasBJet = True
                break
            else:
                vJetsLF.append(vJet)
        
        if hasBJet:
            continue
        
        if len(vJetsLF) < 2:
            continue
        
        # Sort by pT and take two leading
        vJetsLF.sort(key=lambda j: j.Pt(), reverse=True)
        q1, q2 = vJetsLF[0], vJetsLF[1]
        
        # VBF cuts
        deta = abs(q1.Eta() - q2.Eta())
        mjj = (q1 + q2).M()
        
        if deta <= 2.2 or mjj <= 450:
            continue
        
        counters['pass_VBF'] += 1
        
        # Determine VBF sub-category
        fatjet_pt_val = vFatH_pt
        is_VBFjjHi = (deta > 3.0 and mjj > 900)
        
        if is_VBFjjHi:
            if fatjet_pt_val < 400:
                vbf_cat = 'VBFjjHiPtLo'
            else:
                vbf_cat = 'VBFjjHiPtHi'
        else:
            if fatjet_pt_val < 400:
                vbf_cat = 'VBFjjLoPtLo'
            else:
                vbf_cat = 'VBFjjLoPtHi'
        
        # Prepare values for trees
        evt_val = getattr(chain, 'event', 0)
        
        vbf1_pt = q1.Pt()
        vbf1_eta = q1.Eta()
        vbf1_phi = q1.Phi()
        vbf1_mass = q1.M()
        
        vbf2_pt = q2.Pt()
        vbf2_eta = q2.Eta()
        vbf2_phi = q2.Phi()
        vbf2_mass = q2.M()
        
        score = xFatH_X4b
        
        common_values = {
            'run': run_val, 'luminosityBlock': ls_val, 'event': evt_val,
            'fatjet_pt': fatjet_pt_val, 'fatjet_eta': vFatH_eta, 'fatjet_phi': vFatH_phi, 'fatjet_mass': vFatH_mass,
            'vbfjet1_pt': vbf1_pt, 'vbfjet1_eta': vbf1_eta, 'vbfjet1_phi': vbf1_phi, 'vbfjet1_mass': vbf1_mass,
            'vbfjet2_pt': vbf2_pt, 'vbfjet2_eta': vbf2_eta, 'vbfjet2_phi': vbf2_phi, 'vbfjet2_mass': vbf2_mass,
            'vbfjj_deta': deta, 'vbfjj_mass': mjj, 'x4b_score': score,
            'pnet_massH_v2b': pnet_massH_v2b_val,
            'pnet_34massAa': pnet_34massAa_val,
            'trigger_fired': trigger_fired_code,
            'trigger_priority': trigger_priority_code
        }
        
        # Fill trees based on X4b score
        if score > 0.96:
            tree_name = f"{vbf_cat}_X4bSR"
            fill_tree(trees[tree_name], vectors_dict[tree_name], common_values)
            counters[tree_name] += 1
        
        if score > 0.84:
            tree_name = f"{vbf_cat}_X4bSB"
            fill_tree(trees[tree_name], vectors_dict[tree_name], common_values)
            counters[tree_name] += 1
        
        if 0.84 < score <= 0.96:
            tree_name = f"{vbf_cat}_X4bSB_only"
            fill_tree(trees[tree_name], vectors_dict[tree_name], common_values)
            counters[tree_name] += 1
    
    output_file.Write()
    output_file.Close()
    
    print(f"\n{'='*70}")
    print(f"Summary for {dataset_name}:")
    print(f"{'='*70}")
    print(f"Total events processed: {counters['total']}")
    print(f"Pass trigger: {counters['pass_trigger']}")
    print(f"Pass fatH presel: {counters['pass_fatH']}")
    print(f"Pass VBF selection: {counters['pass_VBF']}")
    
    print(f"\nTree counts:")
    for vbf_cat in VBF_CATEGORIES:
        print(f"\n  {vbf_cat}:")
        for x4b_sel in X4B_SELECTIONS:
            tree_name = f"{vbf_cat}_{x4b_sel}"
            print(f"    {x4b_sel}: {counters[tree_name]}")
    
    print(f"\nOutput saved to: {output_name}")
    
    return True


def main():
    parser = argparse.ArgumentParser(description='Process VBF events and create TTrees with X4b selections')
    parser.add_argument('datasets', type=str, nargs='+', help='List of dataset names to process')
    parser.add_argument('--output', '-o', type=str, default=None, help='Output ROOT file name (default: {dataset_name}_output.root)')
    
    args = parser.parse_args()
    
    print(f"\n{'='*70}")
    print(f"Processing {len(args.datasets)} dataset(s)")
    print(f"{'='*70}")
    
    for dataset in args.datasets:
        process_dataset(dataset, args.output)
    
    print(f"\n{'='*70}")
    print("All done!")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
