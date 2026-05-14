import ROOT
import os
import json
from collections import defaultdict

# Disable batch mode if you want to see plots interactively (set to True for batch processing)
ROOT.gROOT.SetBatch(False)

# Define your files and channels
files = {
    '2016': 'output_2016.root',
    '2017': 'output_2017.root', 
    '2018': 'output_2018.root',
    'run2': 'output_run2.root'
}

channels = ['VBFjjHiPtLo', 'VBFjjHiPtHi', 'VBFjjLoPtLo', 'VBFjjLoPtHi']
# Define colors for each channel (ROOT colors)
channel_colors = {
    'VBFjjHiPtLo': ROOT.kRed,
    'VBFjjHiPtHi': ROOT.kBlue,
    'VBFjjLoPtLo': ROOT.kGreen+2,
    'VBFjjLoPtHi': ROOT.kMagenta
}
# Define markers/styles
channel_styles = {
    'VBFjjHiPtLo': 20,  # circle
    'VBFjjHiPtHi': 21,  # square
    'VBFjjLoPtLo': 22,  # triangle
    'VBFjjLoPtHi': 23   # diamond
}

tree_types = ['X4bSR', 'X4bSB_only']
variables = ['x4b_score', 'vbfjj_mass', 'vbfjj_deta', 'fatjet_mass', 'fatjet_pt']

# Variable labels for plotting
var_labels = {
    'x4b_score': 'X4b Score',
    'vbfjj_mass': 'VBFjj Mass (GeV)',
    'vbfjj_deta': 'VBFjj #Delta#eta',
    'fatjet_mass': 'FatJet Mass (GeV)',
    'fatjet_pt': 'FatJet p_{T} (GeV)'
}

# Mass regions: {name: (mH_min, mH_max, mA_min, mA_max)}
mass_regions = {
    'main': (110, 140, 37, 45),
    'mA_upper': (110, 140, 45, 49),
    'mA_lower': (110, 140, 33, 37),
    'mH_lower': (100, 110, 37, 45),
    'mH_upper': (140, 160, 37, 45)
}

# Region labels for plotting
region_labels = {
    'main': 'mH: 110-140, mA: 37-45',
    'mA_upper': 'mH: 110-140, mA: 45-49',
    'mA_lower': 'mH: 110-140, mA: 33-37',
    'mH_lower': 'mH: 100-110, mA: 37-45',
    'mH_upper': 'mH: 140-160, mA: 37-45'
}

# Variable ranges for plots (min, max, nbins) - adjust as needed
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
    Using pnet_massH_v2b for mH and pnet_34massAa for mA
    Taking first element of vectors
    """
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

def make_plot_with_channels(histograms_dict, year, tree_type, region_name, variable, output_dir):
    """Create a plot with all channels overlaid"""
    
    canvas_name = f"{year}_{tree_type}_{region_name}_{variable}"
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 900, 700)
    
    # Set pad for better visibility
    pad = ROOT.TPad("pad", "pad", 0, 0, 1, 1)
    pad.SetLeftMargin(0.15)
    pad.SetRightMargin(0.15)
    pad.SetTopMargin(0.1)
    pad.SetBottomMargin(0.12)
    pad.Draw()
    pad.cd()
    
    # Create legend
    legend = ROOT.TLegend(0.70, 0.70, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    
    # Draw first histogram to set the frame
    first_channel = list(histograms_dict.keys())[0]
    first_hist = histograms_dict[first_channel]
    first_hist.SetTitle(f"{year} {tree_type} {region_labels[region_name]};{var_labels[variable]};Events / bin")
    first_hist.SetLineColor(channel_colors[first_channel])
    first_hist.SetLineWidth(2)
    first_hist.SetFillStyle(3001)  # Hatched
    first_hist.SetFillColorAlpha(channel_colors[first_channel], 0.3)
    first_hist.Draw("HIST")
    
    # Add event count to legend for first channel
    entries = first_hist.GetEntries()
    legend.AddEntry(first_hist, f"{first_channel} (N={int(entries)})", "f")
    
    # Draw remaining histograms
    for channel, hist in histograms_dict.items():
        if channel == first_channel:
            continue
        hist.SetLineColor(channel_colors[channel])
        hist.SetLineWidth(2)
        hist.SetFillStyle(3001)
        hist.SetFillColorAlpha(channel_colors[channel], 0.3)
        hist.Draw("HIST SAME")
        
        entries = hist.GetEntries()
        legend.AddEntry(hist, f"{channel} (N={int(entries)})", "f")
    
    # Draw legend
    legend.Draw()
    
    # Add CMS-like header (optional - adjust as needed)
    cms_text = ROOT.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextFont(42)
    cms_text.SetTextSize(0.04)
    cms_text.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{Simulation}}")
    
    # Add year info
    year_text = ROOT.TLatex()
    year_text.SetNDC()
    year_text.SetTextFont(42)
    year_text.SetTextSize(0.035)
    year_text.DrawLatex(0.85, 0.96, f"{year}")
    
    # Save plot
    output_file = os.path.join(output_dir, f"{variable}.png")
    canvas.SaveAs(output_file)
    output_root = os.path.join(output_dir, f"{variable}.root")
    canvas.SaveAs(output_root)
    
    canvas.Close()
    return canvas

def process_file(filename, year, output_base_dir='plots'):
    """Process a single ROOT file"""
    f = ROOT.TFile(filename, 'READ')
    event_counts = defaultdict(int)
    
    # For each tree_type and region and variable, collect histograms per channel
    for tree_type in tree_types:
        for region_name in mass_regions.keys():
            for variable in variables:
                # Create output directory
                output_dir = os.path.join(output_base_dir, year, tree_type, region_name)
                os.makedirs(output_dir, exist_ok=True)
                
                # Dictionary to store histograms for this combination
                histograms = {}
                
                # Process each channel
                for channel in channels:
                    tree_name = f"{channel}_{tree_type}"
                    tree = f.Get(tree_name)
                    
                    if not tree:
                        print(f"Warning: {tree_name} not found in {filename}")
                        continue
                    
                    # Create histogram for this channel
                    hist_name = f"{year}_{channel}_{tree_type}_{region_name}_{variable}"
                    hist = ROOT.TH1F(hist_name, hist_name,
                                   var_ranges[variable][2], 
                                   var_ranges[variable][0], 
                                   var_ranges[variable][1])
                    
                    # Fill histogram with cuts
                    entry_count = 0
                    mH_min, mH_max, mA_min, mA_max = mass_regions[region_name]
                    
                    for entry in tree:
                        if apply_mass_cuts(entry, mH_min, mH_max, mA_min, mA_max):
                            entry_count += 1
                            value = get_variable_value(entry, variable)
                            if value != -999:
                                hist.Fill(value)
                    
                    # Store histogram if it has entries
                    if hist.GetEntries() > 0:
                        histograms[channel] = hist
                    
                    # Store event count
                    key = f"{year}/{channel}/{tree_type}/{region_name}"
                    event_counts[key] = entry_count
                
                # Make plot with all channels together
                if histograms:
                    print(f"  Plotting {year}/{tree_type}/{region_name}/{variable}...")
                    make_plot_with_channels(histograms, year, tree_type, region_name, variable, output_dir)
                else:
                    print(f"  No entries for {year}/{tree_type}/{region_name}/{variable}")
    
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
    
    # Organize by year, then tree_type, then region
    summary = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for key, count in sorted(all_counts.items()):
        year, channel, tree_type, region = key.split('/')
        summary[year][tree_type][region][channel] = count
    
    for year in sorted(summary.keys()):
        print(f"\n{year}:")
        for tree_type in summary[year].keys():
            print(f"  {tree_type}:")
            for region in summary[year][tree_type].keys():
                print(f"    {region_labels[region]}:")
                for channel in channels:
                    count = summary[year][tree_type][region].get(channel, 0)
                    print(f"      {channel}: {count} events")
    
    # Save counts to JSON
    with open('event_counts.json', 'w') as f:
        json.dump(all_counts, f, indent=2)
    
    print("\n" + "="*80)
    print(f"Plots saved in 'plots/' directory")
    print("Event counts saved in 'event_counts.json'")
    print("="*80)

if __name__ == "__main__":
    main()
