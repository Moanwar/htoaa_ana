import ROOT
import os
from collections import defaultdict

# Disable ROOT's signal handling and enable batch mode
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input files
files = {
    '2016': 'output_2016.root',
    '2017': 'output_2017.root',
    '2018': 'output_2018.root',
    'run2': 'output_run2.root'
}

# Channels
channels = ["VBFjjHiPtLo", "VBFjjHiPtHi", "VBFjjLoPtLo", "VBFjjLoPtHi"]

tree_types = ['X4bSR', 'X4bSB_only']

# Variables to plot (1D)
variables = ['x4b_score', 'vbfjj_mass', 'vbfjj_deta', 'fatjet_mass', 'fatjet_pt']
var_labels = {
    'x4b_score': 'X4b Score',
    'vbfjj_mass': 'VBFjj Mass (GeV)',
    'vbfjj_deta': 'VBFjj #Delta#eta',
    'fatjet_mass': 'FatJet Mass (GeV)',
    'fatjet_pt': 'FatJet p_{T} (GeV)'
}

var_ranges = {
    'x4b_score': (0, 1, 50),
    'vbfjj_mass': (0, 2000, 50),
    'vbfjj_deta': (0, 10, 50),
    'fatjet_mass': (0, 500, 50),
    'fatjet_pt': (0, 1000, 50)
}

# Mass regions
mH_min, mH_max = 110, 140
mA_min, mA_max = 37, 45

# Sideband regions
mH_bin_width = 5
mA_bin_width = 2

mH_lower_min, mH_lower_max = 115, 135
mH_upper_min, mH_upper_max = 120, 130
mA_left_min = mA_min - 2 * mA_bin_width
mA_left_max = mA_min
mA_right_min = mA_max
mA_right_max = mA_max + 2 * mA_bin_width

# Define all regions
mass_regions = {
    'Signal': (mH_min, mH_max, mA_min, mA_max),
    'Lower_mH_Sideband': (mH_lower_min, mH_lower_max, mA_min, mA_max),
    'Upper_mH_Sideband': (mH_upper_min, mH_upper_max, mA_min, mA_max),
    'Left_mA_Sideband': (mH_min, mH_max, mA_left_min, mA_left_max),
    'Right_mA_Sideband': (mH_min, mH_max, mA_right_min, mA_right_max)
}

# Colors (same as your macro)
colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kOrange+7]

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def apply_mass_cuts(entry, mH_min, mH_max, mA_min, mA_max):
    """Apply mass region cuts"""
    try:
        # Handle both scalar and array cases
        mH = entry.pnet_massH_v2b
        mA = entry.pnet_massAa
        
        # If it's an array with length > 0, take first element
        if hasattr(mH, '__len__') and len(mH) > 0:
            mH = mH[0]
        if hasattr(mA, '__len__') and len(mA) > 0:
            mA = mA[0]
        
        # If it's None or invalid, skip
        if mH is None or mA is None:
            return False
            
        if (mH_min <= mH < mH_max) and (mA_min <= mA < mA_max):
            return True
    except:
        pass
    return False

def get_variable_value(entry, var):
    """Extract variable value"""
    try:
        val = getattr(entry, var)
        # Handle both scalar and array cases
        if hasattr(val, '__len__') and len(val) > 0:
            return val[0]
        return val if val is not None else -999
    except:
        return -999

def get_region_entries(tree, mH_min, mH_max, mA_min, mA_max):
    """Count entries in a specific mass region"""
    if not tree:
        return 0
    
    count = 0
    for entry in tree:
        if apply_mass_cuts(entry, mH_min, mH_max, mA_min, mA_max):
            count += 1
    return count

def create_1D_plot(histograms_dict, year, tree_type, region_name, variable, region_label):
    """Create a 1D plot with multiple channels overlaid"""
    if not histograms_dict:
        return False
    
    c = ROOT.TCanvas(f"c_{year}_{tree_type}_{region_name}_{variable}", 
                     f"{year} {tree_type} {region_name} {variable}", 800, 600)
    
    legend = ROOT.TLegend(0.65, 0.65, 0.9, 0.9)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    
    max_y = 0
    first = True
    
    for i, channel in enumerate(channels):
        if channel not in histograms_dict:
            continue
            
        hist = histograms_dict[channel]
        entries = hist.GetEntries()
        
        hist.SetLineColor(colors[i % len(colors)])
        hist.SetLineWidth(2)
        
        if first:
            hist.Draw("HIST")
            hist.GetXaxis().SetTitle(var_labels[variable])
            hist.GetYaxis().SetTitle("Events")
            hist.SetTitle(f"{year} {tree_type} - {region_label}")
            first = False
        else:
            hist.Draw("HIST SAME")
        
        if hist.GetMaximum() > max_y:
            max_y = hist.GetMaximum()
        
        display_name = channel.replace('VBFjj', '')
        legend.AddEntry(hist, f"{display_name} (N = {int(entries)})", "l")
    
    if histograms_dict:
        for hist in histograms_dict.values():
            hist.GetYaxis().SetRangeUser(0, max_y * 1.3)
    
    legend.Draw()
    
    c.Update()
    c.SaveAs(f"plots/{year}_{tree_type}_{region_name}_{variable}.png")
    c.SaveAs(f"plots/{year}_{tree_type}_{region_name}_{variable}.pdf")
    c.Close()
    
    return True

# ============================================================================
# MAIN PROCESSING
# ============================================================================

def main():
    # Create output directory
    os.system("mkdir -p plots")
    
    # Set style
    ROOT.gStyle.SetOptStat(0)
    
    print("\n" + "="*100)
    print("LOADING DATA")
    print("="*100)
    
    # Dictionary to store histograms directly (avoid keeping tree references)
    # Instead of storing trees, we'll fill histograms immediately
    event_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    
    # Process each file one by one and close it immediately
    for year, filename in files.items():
        if not os.path.exists(filename):
            print(f"Warning: {filename} not found, skipping...")
            continue
        
        print(f"\nProcessing {year} from {filename}...")
        f = ROOT.TFile.Open(filename, "READ")
        
        if not f or f.IsZombie():
            print(f"  Could not open {filename}")
            continue
        
        # Process each tree type and channel
        for tree_type in tree_types:
            for channel in channels:
                tree_name = f"{channel}_{tree_type}"
                tree = f.Get(tree_name)
                
                if not tree or tree.IsZombie():
                    print(f"  Warning: {tree_name} not found")
                    continue
                
                print(f"  Processing {year}/{tree_type}/{channel}")
                
                # Count events in each region and store for later printing
                for region_name, (mH_lo, mH_hi, mA_lo, mA_hi) in mass_regions.items():
                    count = 0
                    for entry in tree:
                        if apply_mass_cuts(entry, mH_lo, mH_hi, mA_lo, mA_hi):
                            count += 1
                    event_counts[year][tree_type][region_name][channel] = count
        
        # Close file immediately after processing
        f.Close()
        print(f"  Closed {filename}")
    
    # ========================================================================
    # PRINT EVENT COUNTS (matching your macro style)
    # ========================================================================
    
    print("\n" + "="*100)
    print("EVENT COUNTS SUMMARY")
    print("="*100)
    
    year_cats = ["run2", "2016", "2017", "2018"]
    
    # Print by year
    for tree_type in tree_types:
        print(f"\n{'='*100}")
        print(f"TREE TYPE: {tree_type}")
        print(f"{'='*100}")
        
        for year_cat in year_cats:
            if year_cat not in event_counts:
                print(f"\nWarning: {year_cat} not found in data")
                continue
            
            print(f"\n{'='*100}")
            print(f"YEAR: {year_cat.upper()}")
            print(f"{'='*100}")
            
            for region_name, (mH_lo, mH_hi, mA_lo, mA_hi) in mass_regions.items():
                print(f"\n--- {region_name} (mH: {mH_lo}-{mH_hi}, mA: {mA_lo}-{mA_hi}) ---")
                for channel in channels:
                    count = event_counts[year_cat][tree_type][region_name].get(channel, 0)
                    print(f"  {channel:20s}: {count:10.0f} events")
        
        # Print by channel
        print(f"\n{'='*100}")
        print(f"BY CHANNEL - {tree_type}")
        print(f"{'='*100}")
        
        for channel in channels:
            print(f"\n{'='*100}")
            print(f"CHANNEL: {channel}")
            print(f"{'='*100}")
            
            for region_name, (mH_lo, mH_hi, mA_lo, mA_hi) in mass_regions.items():
                print(f"\n--- {region_name} (mH: {mH_lo}-{mH_hi}, mA: {mA_lo}-{mA_hi}) ---")
                for year_cat in year_cats:
                    count = event_counts[year_cat][tree_type][region_name].get(channel, 0)
                    print(f"  {year_cat:10s}: {count:10.0f} events")
    
    # ========================================================================
    # CREATE 1D PLOTS (need to re-open files to fill histograms)
    # ========================================================================
    
    print("\n" + "="*100)
    print("CREATING 1D PLOTS")
    print("="*100)
    
    region_labels = {
        'Signal': f'Signal (mH: {mH_min}-{mH_max}, mA: {mA_min}-{mA_max})',
        'Lower_mH_Sideband': f'Lower mH Sideband (mH: {mH_lower_min}-{mH_lower_max}, mA: {mA_min}-{mA_max})',
        'Upper_mH_Sideband': f'Upper mH Sideband (mH: {mH_upper_min}-{mH_upper_max}, mA: {mA_min}-{mA_max})',
        'Left_mA_Sideband': f'Left mA Sideband (mH: {mH_min}-{mH_max}, mA: {mA_left_min}-{mA_left_max})',
        'Right_mA_Sideband': f'Right mA Sideband (mH: {mH_min}-{mH_max}, mA: {mA_right_min}-{mA_right_max})'
    }
    
    # Re-open files to create histograms
    for year, filename in files.items():
        if not os.path.exists(filename) or year not in year_cats:
            continue
        
        print(f"\nCreating plots for {year}...")
        f = ROOT.TFile.Open(filename, "READ")
        
        if not f or f.IsZombie():
            continue
        
        for tree_type in tree_types:
            for region_name, (mH_lo, mH_hi, mA_lo, mA_hi) in mass_regions.items():
                for variable in variables:
                    histograms = {}
                    
                    for channel in channels:
                        tree_name = f"{channel}_{tree_type}"
                        tree = f.Get(tree_name)
                        
                        if not tree or tree.IsZombie():
                            continue
                        
                        # Create histogram
                        min_val, max_val, n_bins = var_ranges[variable]
                        hist_name = f"{year}_{channel}_{tree_type}_{region_name}_{variable}"
                        hist = ROOT.TH1F(hist_name, hist_name, n_bins, min_val, max_val)
                        
                        # Fill histogram
                        for entry in tree:
                            if apply_mass_cuts(entry, mH_lo, mH_hi, mA_lo, mA_hi):
                                value = get_variable_value(entry, variable)
                                if value != -999 and min_val <= value <= max_val:
                                    hist.Fill(value)
                        
                        if hist.GetEntries() > 0:
                            histograms[channel] = hist
                    
                    # Create plot
                    if histograms:
                        create_1D_plot(histograms, year, tree_type, region_name, variable, region_labels[region_name])
                        print(f"  Created: {year}_{tree_type}_{region_name}_{variable}")
        
        f.Close()
    
    print("\n" + "="*100)
    print(f"Done! All plots saved in 'plots/' directory")
    print("="*100)

if __name__ == "__main__":
    main()
