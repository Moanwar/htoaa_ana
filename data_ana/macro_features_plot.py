import ROOT
import os
import json
from collections import defaultdict

# Disable batch mode if you want to see plots interactively
ROOT.gROOT.SetBatch(False)

# Define your files and channels
files = {
    '2016': 'output_2016.root',
    '2017': 'output_2017.root', 
    '2018': 'output_2018.root',
    'run2': 'output_run2.root'
}

channels = ['VBFjjHiPtLo', 'VBFjjHiPtHi', 'VBFjjLoPtLo', 'VBFjjLoPtHi']
tree_types = ['X4bSR', 'X4bSB_only']  # Note: SB_only not SB
variables = ['x4b_score', 'vbfjj_mass', 'vbfjj_deta', 'fatjet_mass', 'fatjet_pt']

# Mass regions: {name: (mH_min, mH_max, mA_min, mA_max)}
mass_regions = {
    'main': (110, 140, 37, 45),
    'mA_upper': (110, 140, 45, 49),
    'mA_lower': (110, 140, 33, 37),
    'mH_lower': (100, 110, 37, 45),
    'mH_upper': (140, 160, 37, 45)
}

# Variable ranges for plots (adjust as needed)
var_ranges = {
    'x4b_score': (0, 1, 50),
    'vbfjj_mass': (0, 2000, 50),  # Adjust based on your physics
    'vbfjj_deta': (0, 10, 50),
    'fatjet_mass': (0, 500, 50),
    'fatjet_pt': (0, 1000, 50)
}

def apply_mass_cuts(entry, mH_min, mH_max, mA_min, mA_max):
    """
    Apply mass region cuts.
    ASSUMING: pnet_massH_v2b is mH, pnet_34massAa is mA
    AND: these are vectors - taking first element
    """
    try:
        mH = entry.pnet_massH_v2b[0] if len(entry.pnet_massH_v2b) > 0 else -999
        mA = entry.pnet_34massAa[0] if len(entry.pnet_34massAa) > 0 else -999
        
        if (mH_min <= mH < mH_max) and (mA_min <= mA < mA_max):
            return True
    except:
        pass
    return False

def process_file(filename, year, output_dir='plots'):
    """Process a single ROOT file"""
    f = ROOT.TFile(filename, 'READ')
    event_counts = defaultdict(int)
    
    for channel in channels:
        for tree_type in tree_types:
            tree_name = f"{channel}_{tree_type}"
            tree = f.Get(tree_name)
            
            if not tree:
                print(f"Warning: {tree_name} not found in {filename}")
                continue
            
            # Create output directory structure
            year_dir = os.path.join(output_dir, year, channel, tree_type)
            os.makedirs(year_dir, exist_ok=True)
            
            # Loop over mass regions
            for region_name, (mH_min, mH_max, mA_min, mA_max) in mass_regions.items():
                region_dir = os.path.join(year_dir, region_name)
                os.makedirs(region_dir, exist_ok=True)
                
                # Create histograms for this region
                histograms = {}
                for var in variables:
                    hist_name = f"{year}_{channel}_{tree_type}_{region_name}_{var}"
                    hist = ROOT.TH1F(hist_name, hist_name, 
                                    var_ranges[var][2], var_ranges[var][0], var_ranges[var][1])
                    histograms[var] = hist
                
                # Fill histograms
                entry_count = 0
                for entry in tree:
                    if apply_mass_cuts(entry, mH_min, mH_max, mA_min, mA_max):
                        entry_count += 1
                        for var in variables:
                            try:
                                # Assuming all variables are vectors, take first element
                                value = getattr(entry, var)[0] if len(getattr(entry, var)) > 0 else -999
                                if value != -999:
                                    histograms[var].Fill(value)
                            except:
                                pass
                
                # Store event count
                key = f"{year}/{channel}/{tree_type}/{region_name}"
                event_counts[key] = entry_count
                
                # Save and plot histograms
                for var, hist in histograms.items():
                    if hist.GetEntries() > 0:
                        canvas = ROOT.TCanvas(f"c_{hist.GetName()}", "", 800, 600)
                        hist.Draw()
                        hist.SetLineColor(1)
                        hist.SetTitle(f"{year} {channel} {tree_type} {region_name};{var};Events")
                        
                        # Save as root and png
                        output_root = os.path.join(region_dir, f"{var}.root")
                        output_png = os.path.join(region_dir, f"{var}.png")
                        
                        f_out = ROOT.TFile(output_root, 'RECREATE')
                        hist.Write()
                        f_out.Close()
                        
                        canvas.SaveAs(output_png)
                        canvas.Close()
                
                print(f"  {year}/{channel}/{tree_type}/{region_name}: {entry_count} events")
    
    f.Close()
    return event_counts

def main():
    # Process each file
    all_counts = {}
    for year, filename in files.items():
        if not os.path.exists(filename):
            print(f"File {filename} not found, skipping...")
            continue
        print(f"\nProcessing {year} ({filename})...")
        counts = process_file(filename, year)
        all_counts.update(counts)
    
    # Print summary of all event counts
    print("\n" + "="*80)
    print("SUMMARY OF EVENT COUNTS:")
    print("="*80)
    for key, count in sorted(all_counts.items()):
        print(f"{key}: {count} events")
    
    # Save counts to JSON
    with open('event_counts.json', 'w') as f:
        json.dump(all_counts, f, indent=2)
    
    print("\nPlots saved in 'plots/' directory")
    print("Event counts saved in 'event_counts.json'")

if __name__ == "__main__":
    main()
