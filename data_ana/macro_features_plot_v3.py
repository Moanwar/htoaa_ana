import ROOT
import os
from collections import defaultdict

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

# Channels (similar to sub_cats in your macro)
channels = ["VBFjjHiPtLo", "VBFjjHiPtHi", "VBFjjLoPtLo", "VBFjjLoPtHi"]
channel_labels = {
    'VBFjjHiPtLo': 'HiPtLo',
    'VBFjjHiPtHi': 'HiPtHi',
    'VBFjjLoPtLo': 'LoPtLo',
    'VBFjjLoPtHi': 'LoPtHi'
}

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

# Mass regions (using your exact definitions)
mH_min, mH_max = 110, 140
mA_min, mA_max = 37, 45

# Sideband regions as in your macro
mH_bin_width = 5
mA_bin_width = 2

mH_lower_min, mH_lower_max = 115, 135  # Lower mH sideband
mH_upper_min, mH_upper_max = 120, 130  # Upper mH sideband
mA_left_min = mA_min - 2 * mA_bin_width  # 33 GeV
mA_left_max = mA_min  # 37 GeV
mA_right_min = mA_max  # 45 GeV
mA_right_max = mA_max + 2 * mA_bin_width  # 49 GeV

# Define all regions
mass_regions = {
    'Signal': (mH_min, mH_max, mA_min, mA_max),
    'Lower_mH_Sideband': (mH_lower_min, mH_lower_max, mA_min, mA_max),
    'Upper_mH_Sideband': (mH_upper_min, mH_upper_max, mA_min, mA_max),
    'Left_mA_Sideband': (mH_min, mH_max, mA_left_min, mA_left_max),
    'Right_mA_Sideband': (mH_min, mH_max, mA_right_min, mA_right_max)
}

# Colors for plotting
colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kOrange+7]
markers = [20, 21, 22, 23]

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def apply_mass_cuts(entry, mH_min, mH_max, mA_min, mA_max):
    """Apply mass region cuts using pnet_massH_v2b and pnet_34massAa"""
    try:
        mH = entry.pnet_massH_v2b[0] if len(entry.pnet_massH_v2b) > 0 else -999
        mA = entry.pnet_34massAa[0] if len(entry.pnet_34massAa) > 0 else -999
        
        if (mH_min <= mH < mH_max) and (mA_min <= mA < mA_max):
            return True
    except:
        pass
    return False

def get_variable_value(entry, var):
    """Extract variable value (assuming vector, take first element)"""
    try:
        val = getattr(entry, var)
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

def create_1D_plot(histograms_dict, title, filename, x_title, year_cat, region_name):
    """Create a 1D plot with multiple channels overlaid"""
    if not histograms_dict:
        return False
    
    canvas = ROOT.TCanvas(f"c_{filename}", title, 800, 600)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)
    
    legend = ROOT.TLegend(0.65, 0.65, 0.9, 0.9)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    
    max_y = 0
    first = True
    
    for i, (channel, hist) in enumerate(histograms_dict.items()):
        entries = hist.GetEntries()
        hist.SetLineColor(colors[i % len(colors)])
        hist.SetLineWidth(2)
        hist.SetFillStyle(3001)
        hist.SetFillColorAlpha(colors[i % len(colors)], 0.2)
        
        if first:
            hist.Draw("HIST")
            hist.GetXaxis().SetTitle(x_title)
            hist.GetYaxis().SetTitle("Events / bin")
            hist.SetTitle(title)
            first = False
        else:
            hist.Draw("HIST SAME")
        
        if hist.GetMaximum() > max_y:
            max_y = hist.GetMaximum()
        
        legend.AddEntry(hist, f"{channel_labels[channel]} (N = {int(entries)})", "f")
    
    if histograms_dict:
        for hist in histograms_dict.values():
            hist.GetYaxis().SetRangeUser(0, max_y * 1.3)
    
    legend.Draw()
    
    # Add CMS label
    cms_text = ROOT.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextFont(42)
    cms_text.SetTextSize(0.04)
    cms_text.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{Simulation}}")
    
    # Add year and region
    info_text = ROOT.TLatex()
    info_text.SetNDC()
    info_text.SetTextFont(42)
    info_text.SetTextSize(0.035)
    info_text.DrawLatex(0.7, 0.96, f"{year_cat}, {region_name}")
    
    canvas.Update()
    canvas.SaveAs(f"plots/1D_{filename}.png")
    canvas.SaveAs(f"plots/1D_{filename}.pdf")
    canvas.Close()
    
    return True

def create_2D_plot(tree, channel, year, region_name, x_min, x_max, y_min, y_max):
    """Create a 2D histogram for mH vs mA in a specific region"""
    if not tree:
        return False
    
    # Create 2D histogram
    hist_name = f"{year}_{channel}_{region_name}_2D"
    n_bins_x = int((x_max - x_min) / mH_bin_width)
    n_bins_y = int((y_max - y_min) / mA_bin_width)
    
    hist2D = ROOT.TH2D(hist_name, f"{year} {channel_labels[channel]} {region_name};m_{H} (GeV);m_{A} (GeV)",
                       n_bins_x, x_min, x_max, n_bins_y, y_min, y_max)
    
    # Fill histogram
    for entry in tree:
        try:
            mH = entry.pnet_massH_v2b[0] if len(entry.pnet_massH_v2b) > 0 else -999
            mA = entry.pnet_34massAa[0] if len(entry.pnet_34massAa) > 0 else -999
            
            if (x_min <= mH < x_max) and (y_min <= mA < y_max):
                hist2D.Fill(mH, mA)
        except:
            pass
    
    if hist2D.GetEntries() == 0:
        return False
    
    # Create canvas
    canvas = ROOT.TCanvas(f"c_2D_{hist_name}", hist_name, 800, 600)
    canvas.SetRightMargin(0.15)
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.12)
    
    # Draw with COLZ and TEXT to show bin contents
    hist2D.Draw("COLZ TEXT")
    
    # Set palette
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    
    # Add CMS label
    cms_text = ROOT.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextFont(42)
    cms_text.SetTextSize(0.04)
    cms_text.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{Simulation}}")
    
    canvas.Update()
    canvas.SaveAs(f"plots/2D_{year}_{channel}_{region_name}.png")
    canvas.SaveAs(f"plots/2D_{year}_{channel}_{region_name}.pdf")
    canvas.Close()
    
    return True

# ============================================================================
# MAIN PROCESSING
# ============================================================================

def main():
    # Create output directory
    os.system("mkdir -p plots")
    
    # Set style
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    
    # Store all data
    data = {}
    event_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    
    print("\n" + "="*100)
    print("LOADING DATA")
    print("="*100)
    
    # Load data from all files
    for year, filename in files.items():
        if not os.path.exists(filename):
            print(f"Warning: {filename} not found, skipping...")
            continue
        
        print(f"\nLoading {year} from {filename}...")
        f = ROOT.TFile.Open(filename, "READ")
        
        if not f or f.IsZombie():
            print(f"  Could not open {filename}")
            continue
        
        data[year] = {}
        
        for channel in channels:
            for tree_type in tree_types:
                tree_name = f"{channel}_{tree_type}"
                tree = f.Get(tree_name)
                
                if tree and not tree.IsZombie():
                    if tree_type not in data[year]:
                        data[year][tree_type] = {}
                    data[year][tree_type][channel] = tree
                    print(f"  Loaded {year}/{tree_type}/{channel}")
        
        f.Close()
    
    # ========================================================================
    # PRINT EVENT COUNTS (like your macro)
    # ========================================================================
    
    print("\n" + "="*100)
    print("EVENT COUNTS SUMMARY")
    print("="*100)
    
    years_to_show = ["run2", "2016", "2017", "2018"]
    
    for tree_type in tree_types:
        print(f"\n{'='*100}")
        print(f"TREE TYPE: {tree_type}")
        print(f"{'='*100}")
        
        for year in years_to_show:
            if year not in data or tree_type not in data[year]:
                print(f"\nWarning: {year}/{tree_type} not found in data")
                continue
            
            print(f"\n{'='*100}")
            print(f"YEAR: {year.upper()}")
            print(f"{'='*100}")
            
            for region_name, (mH_lo, mH_hi, mA_lo, mA_hi) in mass_regions.items():
                print(f"\n--- {region_name} (mH: {mH_lo}-{mH_hi}, mA: {mA_lo}-{mA_hi}) ---")
                for channel in channels:
                    if channel in data[year][tree_type]:
                        tree = data[year][tree_type][channel]
                        count = get_region_entries(tree, mH_lo, mH_hi, mA_lo, mA_hi)
                        print(f"  {channel:20s}: {count:10.0f} events")
                        event_counts[year][tree_type][region_name][channel] = count
        
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
                for year in years_to_show:
                    if year in data and tree_type in data[year] and channel in data[year][tree_type]:
                        tree = data[year][tree_type][channel]
                        count = get_region_entries(tree, mH_lo, mH_hi, mA_lo, mA_hi)
                        print(f"  {year:10s}: {count:10.0f} events")
    
    # ========================================================================
    # CREATE 1D PLOTS FOR EACH VARIABLE
    # ========================================================================
    
    print("\n" + "="*100)
    print("CREATING 1D PLOTS")
    print("="*100)
    
    for tree_type in tree_types:
        for year in years_to_show:
            if year not in data or tree_type not in data[year]:
                continue
            
            for region_name, (mH_lo, mH_hi, mA_lo, mA_hi) in mass_regions.items():
                for variable in variables:
                    # Create histograms for each channel
                    histograms = {}
                    
                    for channel in channels:
                        if channel not in data[year][tree_type]:
                            continue
                        
                        tree = data[year][tree_type][channel]
                        
                        # Create histogram
                        hist_name = f"{year}_{channel}_{tree_type}_{region_name}_{variable}"
                        n_bins = 50
                        if variable == 'x4b_score':
                            hist = ROOT.TH1F(hist_name, hist_name, n_bins, 0, 1)
                        elif variable == 'vbfjj_mass':
                            hist = ROOT.TH1F(hist_name, hist_name, n_bins, 0, 2000)
                        elif variable == 'vbfjj_deta':
                            hist = ROOT.TH1F(hist_name, hist_name, n_bins, 0, 10)
                        elif variable == 'fatjet_mass':
                            hist = ROOT.TH1F(hist_name, hist_name, n_bins, 0, 500)
                        elif variable == 'fatjet_pt':
                            hist = ROOT.TH1F(hist_name, hist_name, n_bins, 0, 1000)
                        
                        # Fill histogram
                        for entry in tree:
                            if apply_mass_cuts(entry, mH_lo, mH_hi, mA_lo, mA_hi):
                                value = get_variable_value(entry, variable)
                                if value != -999:
                                    hist.Fill(value)
                        
                        if hist.GetEntries() > 0:
                            histograms[channel] = hist
                    
                    # Create plot
                    if histograms:
                        title = f"{year} {tree_type} {region_name}"
                        filename = f"{year}_{tree_type}_{region_name}_{variable}"
                        create_1D_plot(histograms, title, filename, var_labels[variable], year, region_name)
                        print(f"  Created: {filename}")
    
    # ========================================================================
    # CREATE 2D PLOTS (mH vs mA)
    # ========================================================================
    
    print("\n" + "="*100)
    print("CREATING 2D PLOTS (mH vs mA)")
    print("="*100)
    
    for tree_type in tree_types:
        for year in years_to_show:
            if year not in data or tree_type not in data[year]:
                continue
            
            for channel in channels:
                if channel not in data[year][tree_type]:
                    continue
                
                tree = data[year][tree_type][channel]
                
                for region_name, (mH_lo, mH_hi, mA_lo, mA_hi) in mass_regions.items():
                    if create_2D_plot(tree, channel, year, region_name, mH_lo, mH_hi, mA_lo, mA_hi):
                        print(f"  Created: 2D_{year}_{channel}_{region_name}")
    
    # ========================================================================
    # SUMMARY
    # ========================================================================
    
    print("\n" + "="*100)
    print("SUMMARY OF EVENT COUNTS BY REGION")
    print("="*100)
    
    for tree_type in tree_types:
        print(f"\n{tree_type}:")
        for region_name in mass_regions.keys():
            print(f"\n  {region_name}:")
            for year in years_to_show:
                if year in event_counts and tree_type in event_counts[year] and region_name in event_counts[year][tree_type]:
                    total = sum(event_counts[year][tree_type][region_name].values())
                    print(f"    {year:10s}: {total:10.0f} total events")
                    for channel, count in event_counts[year][tree_type][region_name].items():
                        print(f"      {channel_labels[channel]:15s}: {count:10.0f}")
    
    print("\n" + "="*100)
    print(f"Done! All plots saved in 'plots/' directory")
    print("="*100)

if __name__ == "__main__":
    main()
