import ROOT
import os

# Define paths and categories
base_path = "/afs/cern.ch/work/m/moanwar/public/hto2ato4b/2DAlphabetfiles_VBF_sys_v2/20251021_DataMC"
years = ["2016preVFP", "2016postVFP", "2017", "2018"]
sub_cats = ["VBFHiPTHi", "VBFHiPTLo", "VBFLoPTHi", "VBFLoPTLo"]

# Region of interest
mH_min, mH_max = 110, 140
mA_min, mA_max = 37, 45

# Bin widths
mH_bin_width = 5  # GeV
mA_bin_width = 2  # GeV

# Calculate sideband regions (exactly 2 bins adjacent, not including the ROI)
# For mA: left sideband (2 bins to the left)
mA_left_min = mA_min - 2 * mA_bin_width  # 37 - 4 = 33 GeV
mA_left_max = mA_min  # 37 GeV

# For mA: right sideband (2 bins to the right)
mA_right_min = mA_max  # 45 GeV
mA_right_max = mA_max + 2 * mA_bin_width  # 45 + 4 = 49 GeV

# For mH: left sideband (100-110 GeV, which is 2 bins of 5 GeV)
mH_left_min = 100
mH_left_max = 110

# For mH: right sideband (140-150? or 140-160? Using 140-150 for 2 bins)
mH_right_min = mH_max  # 140 GeV
mH_right_max = mH_max + 4 * mH_bin_width  # 140 + 20 = 160 GeV

print(f"\nSideband definitions:")
print(f"mA left sideband: {mA_left_min}-{mA_left_max} GeV")
print(f"mA right sideband: {mA_right_min}-{mA_right_max} GeV")
print(f"mH left sideband: {mH_left_min}-{mH_left_max} GeV")
print(f"mH right sideband: {mH_right_min}-{mH_right_max} GeV")

# Dictionary to store all histograms
# Structure: data[year_cat][sub_cat] = 2D histogram
data = {}

print("\nLoading histograms...")

# First, load all histograms without cloning
for year in years:
    for sub_cat in sub_cats:
        file_path = f"{base_path}/{year}/VBFjj/2DAlphabet_inputFiles/{sub_cat}/{sub_cat}_Data_{year}.root"
        
        print(f"Trying: {file_path}")
        
        # Open file
        f = ROOT.TFile.Open(file_path)
        if not f or f.IsZombie():
            print(f"  Warning: Could not open {file_path}")
            continue
            
        # Get histogram
        hist_name = f"{sub_cat}_Data_{year}_pnet_4a_WP40_Pass_Nom"
        hist = f.Get(hist_name)
        
        if hist and not hist.IsZombie():
            # Store the histogram and keep the file open? No - better to clone properly
            # But let's make a proper clone that ROOT will manage
            year_key = year
            if year_key not in data:
                data[year_key] = {}
            
            # Proper cloning with new name
            cloned_hist = hist.Clone(f"{sub_cat}_{year}_cloned")
            cloned_hist.SetDirectory(0)  # Detach from file
            data[year_key][sub_cat] = cloned_hist
            print(f"  Loaded: {year_key} - {sub_cat}")
        else:
            print(f"  Warning: Could not find {hist_name} in file")
        f.Close()

# Combine 2016preVFP and 2016postVFP into "2016"
if "2016preVFP" in data or "2016postVFP" in data:
    data["2016"] = {}
    for sub_cat in sub_cats:
        # Create a new histogram for the sum
        summed_hist = None
        for era in ["2016preVFP", "2016postVFP"]:
            if era in data and sub_cat in data[era]:
                if summed_hist is None:
                    # Create a copy of the first histogram
                    summed_hist = data[era][sub_cat].Clone(f"2016_{sub_cat}_sum")
                    summed_hist.SetDirectory(0)
                else:
                    # Add the next one
                    summed_hist.Add(data[era][sub_cat])
        
        if summed_hist:
            data["2016"][sub_cat] = summed_hist
            print(f"Combined 2016 - {sub_cat}")

# Create "Run2" by summing all years
data["Run2"] = {}
for sub_cat in sub_cats:
    summed_hist = None
    # Sum 2016 (combined), 2017, 2018
    for year in ["2016", "2017", "2018"]:
        if year in data and sub_cat in data[year]:
            if summed_hist is None:
                summed_hist = data[year][sub_cat].Clone(f"Run2_{sub_cat}_sum")
                summed_hist.SetDirectory(0)
            else:
                summed_hist.Add(data[year][sub_cat])
    
    if summed_hist:
        data["Run2"][sub_cat] = summed_hist
        print(f"Combined Run2 - {sub_cat}")

# Check what we have loaded
print("\nLoaded data summary:")
for year in data.keys():
    print(f"  {year}: {list(data[year].keys())}")

def get_region_integral(hist2D, x_min, x_max, y_min, y_max):
    """Get integral of 2D histogram in specified region"""
    if not hist2D or hist2D.IsZombie():
        return 0
    
    x_bin_min = hist2D.GetXaxis().FindBin(x_min)
    x_bin_max = hist2D.GetXaxis().FindBin(x_max - 0.001)
    y_bin_min = hist2D.GetYaxis().FindBin(y_min)
    y_bin_max = hist2D.GetYaxis().FindBin(y_max - 0.001)
    
    # Check bounds
    if (x_bin_min < 1 or x_bin_max > hist2D.GetNbinsX() or 
        y_bin_min < 1 or y_bin_max > hist2D.GetNbinsY()):
        return 0
    
    integral = hist2D.Integral(x_bin_min, x_bin_max, y_bin_min, y_bin_max)
    return integral

def get_1D_projection(hist2D, proj_axis='x'):
    """
    Create 1D projection from 2D histogram for the region mH:110-140, mA:37-45
    proj_axis='x' for mH projection (sum over mA range)
    proj_axis='y' for mA projection (sum over mH range)
    """
    if not hist2D or hist2D.IsZombie():
        return None
    
    # Find bins corresponding to ROI
    x_bins = hist2D.GetXaxis().FindBin(mH_min), hist2D.GetXaxis().FindBin(mH_max-0.001)
    y_bins = hist2D.GetYaxis().FindBin(mA_min), hist2D.GetYaxis().FindBin(mA_max-0.001)
    
    # Ensure bin indices are valid
    if x_bins[0] < 1 or x_bins[1] > hist2D.GetNbinsX() or y_bins[0] < 1 or y_bins[1] > hist2D.GetNbinsY():
        print(f"  Warning: ROI bins out of range: x=({x_bins[0]},{x_bins[1]}), y=({y_bins[0]},{y_bins[1]})")
        return None
    
    try:
        if proj_axis == 'x':
            # Project onto x-axis (mH) for the specified y range (mA 37-45)
            proj = hist2D.ProjectionX(f"{hist2D.GetName()}_projX", y_bins[0], y_bins[1])
            proj.GetXaxis().SetTitle("m_{H} (GeV)")
            proj.GetYaxis().SetTitle("Events")
            proj.SetTitle(f"mH distribution (mA: {mA_min}-{mA_max} GeV)")
            
            # Also restrict x-axis to mH range 110-140
            proj.GetXaxis().SetRange(x_bins[0], x_bins[1])
            
        else:
            # Project onto y-axis (mA) for the specified x range (mH 110-140)
            proj = hist2D.ProjectionY(f"{hist2D.GetName()}_projY", x_bins[0], x_bins[1])
            proj.GetXaxis().SetTitle("m_{A} (GeV)")
            proj.GetYaxis().SetTitle("Events")
            proj.SetTitle(f"mA distribution (mH: {mH_min}-{mH_max} GeV)")
            
            # Also restrict x-axis to mA range 37-45
            proj.GetXaxis().SetRange(y_bins[0], y_bins[1])
        
        proj.SetDirectory(0)  # Detach from any file
        return proj
    except Exception as e:
        print(f"  Warning: Failed to project histogram {hist2D.GetName()}: {e}")
        return None

def create_2D_plot(hist2D, title, filename):
    """Create a 2D heatmap plot with bin contents shown"""
    if not hist2D or hist2D.IsZombie():
        return False
    
    # Clone the histogram to avoid modifying original
    hist_clone = hist2D.Clone(f"{hist2D.GetName()}_clone")
    hist_clone.SetDirectory(0)
    
    # Create canvas
    c = ROOT.TCanvas(f"c_{filename}", title, 800, 600)
    c.SetRightMargin(0.15)
    
    # Draw with COLZ option (color palette with z-axis)
    hist_clone.Draw("COLZ")
    
    # Set titles
    hist_clone.GetXaxis().SetTitle("m_{H} (GeV)")
    hist_clone.GetYaxis().SetTitle("m_{A} (GeV)")
    hist_clone.SetTitle(title)
    
    # Set palette
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    
    # Update and save
    c.Update()
    c.Modified()
    c.SaveAs(f"plots/2D_{filename}.png")
    c.SaveAs(f"plots/2D_{filename}.pdf")
    
    return True

# Create output directory
os.system("mkdir -p plots")

# Set style
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kViridis)

# ============================================================================
# PRINT EVENT COUNTS TO SCREEN
# ============================================================================
print("\n" + "="*100)
print("EVENT COUNTS SUMMARY")
print("="*100)

# Group 1: By year-category (Run2, 2016, 2017, 2018)
year_cats = ["Run2", "2016", "2017", "2018"]
colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kOrange+7]

for year_cat in year_cats:
    if year_cat not in data:
        print(f"Warning: {year_cat} not found in data")
        continue
    
    print(f"\n{'='*100}")
    print(f"YEAR: {year_cat}")
    print(f"{'='*100}")
    
    # Signal region
    print(f"\n--- SIGNAL REGION (mH: {mH_min}-{mH_max}, mA: {mA_min}-{mA_max}) ---")
    for sub_cat in sub_cats:
        if sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_min, mA_max)
            print(f"  {sub_cat:20s}: {total:10.0f} events")
    
    # mA sidebands
    print(f"\n--- mA SIDEBANDS (mH: {mH_min}-{mH_max}) ---")
    print(f"  Left sideband (mA: {mA_left_min}-{mA_left_max}):")
    for sub_cat in sub_cats:
        if sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_left_min, mA_left_max)
            print(f"    {sub_cat:18s}: {total:10.0f} events")
    
    print(f"  Right sideband (mA: {mA_right_min}-{mA_right_max}):")
    for sub_cat in sub_cats:
        if sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_right_min, mA_right_max)
            print(f"    {sub_cat:18s}: {total:10.0f} events")
    
    # mH sidebands
    print(f"\n--- mH SIDEBANDS (mA: {mA_min}-{mA_max}) ---")
    print(f"  Left sideband (mH: {mH_left_min}-{mH_left_max}):")
    for sub_cat in sub_cats:
        if sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_left_min, mH_left_max, mA_min, mA_max)
            print(f"    {sub_cat:18s}: {total:10.0f} events")
    
    print(f"  Right sideband (mH: {mH_right_min}-{mH_right_max}):")
    for sub_cat in sub_cats:
        if sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_right_min, mH_right_max, mA_min, mA_max)
            print(f"    {sub_cat:18s}: {total:10.0f} events")

# Group 2: By sub-category
print(f"\n{'='*100}")
print("BY SUB-CATEGORY")
print("="*100)

for sub_cat in sub_cats:
    print(f"\n{'='*100}")
    print(f"SUB-CATEGORY: {sub_cat}")
    print(f"{'='*100}")
    
    # Signal region
    print(f"\n--- SIGNAL REGION (mH: {mH_min}-{mH_max}, mA: {mA_min}-{mA_max}) ---")
    for year_cat in year_cats:
        if year_cat in data and sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_min, mA_max)
            print(f"  {year_cat:10s}: {total:10.0f} events")
    
    # mA sidebands
    print(f"\n--- mA SIDEBANDS (mH: {mH_min}-{mH_max}) ---")
    print(f"  Left sideband (mA: {mA_left_min}-{mA_left_max}):")
    for year_cat in year_cats:
        if year_cat in data and sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_left_min, mA_left_max)
            print(f"    {year_cat:8s}: {total:10.0f} events")
    
    print(f"  Right sideband (mA: {mA_right_min}-{mA_right_max}):")
    for year_cat in year_cats:
        if year_cat in data and sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_right_min, mA_right_max)
            print(f"    {year_cat:8s}: {total:10.0f} events")
    
    # mH sidebands
    print(f"\n--- mH SIDEBANDS (mA: {mA_min}-{mA_max}) ---")
    print(f"  Left sideband (mH: {mH_left_min}-{mH_left_max}):")
    for year_cat in year_cats:
        if year_cat in data and sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_left_min, mH_left_max, mA_min, mA_max)
            print(f"    {year_cat:8s}: {total:10.0f} events")
    
    print(f"  Right sideband (mH: {mH_right_min}-{mH_right_max}):")
    for year_cat in year_cats:
        if year_cat in data and sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_right_min, mH_right_max, mA_min, mA_max)
            print(f"    {year_cat:8s}: {total:10.0f} events")

# ============================================================================
# CREATE 1D PLOTS (mH and mA projections)
# ============================================================================
print("\n" + "="*100)
print("CREATING 1D PLOTS")
print("="*100)

# Group 1: By year-category (Run2, 2016, 2017, 2018)
for year_cat in year_cats:
    if year_cat not in data:
        print(f"Warning: {year_cat} not found in data")
        continue
    
    # Create mH plot for this year-category
    c_mH = ROOT.TCanvas(f"c_mH_{year_cat}", f"{year_cat} - mH Distribution", 800, 600)
    legend_mH = ROOT.TLegend(0.65, 0.65, 0.9, 0.9)
    
    max_y = 0
    projections = []
    
    for i, sub_cat in enumerate(sub_cats):
        if sub_cat in data[year_cat]:
            proj = get_1D_projection(data[year_cat][sub_cat], 'x')
            if proj:
                total_events = proj.Integral()
                proj.SetLineColor(colors[i])
                proj.SetLineWidth(2)
                proj.SetMarkerColor(colors[i])
                proj.SetMarkerStyle(20)
                proj.SetMarkerSize(0.8)
                projections.append((proj, total_events))
                
                if proj.GetMaximum() > max_y:
                    max_y = proj.GetMaximum()
    
    for idx, (proj, total_events) in enumerate(projections):
        if idx == 0:
            proj.Draw("HIST")
            proj.GetYaxis().SetRangeUser(0, max_y*1.3)
        else:
            proj.Draw("HIST SAME")
        legend_mH.AddEntry(proj, f"{sub_cats[idx].replace('VBF', '')} (N = {total_events:.0f})", "l")
    
    if projections:
        legend_mH.Draw()
        c_mH.Update()
        c_mH.Modified()
        c_mH.SaveAs(f"plots/mH_{year_cat}.png")
        c_mH.SaveAs(f"plots/mH_{year_cat}.pdf")
        print(f"Created mH plot for {year_cat}")
    
    # Create mA plot for this year-category
    c_mA = ROOT.TCanvas(f"c_mA_{year_cat}", f"{year_cat} - mA Distribution", 800, 600)
    legend_mA = ROOT.TLegend(0.65, 0.65, 0.9, 0.9)
    
    max_y = 0
    projections = []
    
    for i, sub_cat in enumerate(sub_cats):
        if sub_cat in data[year_cat]:
            proj = get_1D_projection(data[year_cat][sub_cat], 'y')
            if proj:
                total_events = proj.Integral()
                proj.SetLineColor(colors[i])
                proj.SetLineWidth(2)
                proj.SetMarkerColor(colors[i])
                proj.SetMarkerStyle(20)
                proj.SetMarkerSize(0.8)
                projections.append((proj, total_events))
                
                if proj.GetMaximum() > max_y:
                    max_y = proj.GetMaximum()
    
    for idx, (proj, total_events) in enumerate(projections):
        if idx == 0:
            proj.Draw("HIST")
            proj.GetYaxis().SetRangeUser(0, max_y*1.3)
        else:
            proj.Draw("HIST SAME")
        legend_mA.AddEntry(proj, f"{sub_cats[idx].replace('VBF', '')} (N = {total_events:.0f})", "l")
    
    if projections:
        legend_mA.Draw()
        c_mA.Update()
        c_mA.Modified()
        c_mA.SaveAs(f"plots/mA_{year_cat}.png")
        c_mA.SaveAs(f"plots/mA_{year_cat}.pdf")
        print(f"Created mA plot for {year_cat}")

# Group 2: By sub-category
for sub_cat in sub_cats:
    # Create mH plot for this sub-category
    c_mH = ROOT.TCanvas(f"c_mH_{sub_cat}", f"{sub_cat} - mH Distribution", 800, 600)
    legend_mH = ROOT.TLegend(0.65, 0.65, 0.9, 0.9)
    
    max_y = 0
    projections = []
    
    for i, year_cat in enumerate(year_cats):
        if year_cat in data and sub_cat in data[year_cat]:
            proj = get_1D_projection(data[year_cat][sub_cat], 'x')
            if proj:
                total_events = proj.Integral()
                proj.SetLineColor(colors[i])
                proj.SetLineWidth(2)
                proj.SetMarkerColor(colors[i])
                proj.SetMarkerStyle(20)
                proj.SetMarkerSize(0.8)
                projections.append((proj, total_events, year_cat))
                
                if proj.GetMaximum() > max_y:
                    max_y = proj.GetMaximum()
    
    for idx, (proj, total_events, year_cat) in enumerate(projections):
        if idx == 0:
            proj.Draw("HIST")
            proj.GetYaxis().SetRangeUser(0, max_y*1.3)
        else:
            proj.Draw("HIST SAME")
        legend_mH.AddEntry(proj, f"{year_cat} (N = {total_events:.0f})", "l")
    
    if projections:
        legend_mH.Draw()
        c_mH.Update()
        c_mH.Modified()
        c_mH.SaveAs(f"plots/mH_{sub_cat}.png")
        c_mH.SaveAs(f"plots/mH_{sub_cat}.pdf")
        print(f"Created mH plot for {sub_cat}")
    
    # Create mA plot for this sub-category
    c_mA = ROOT.TCanvas(f"c_mA_{sub_cat}", f"{sub_cat} - mA Distribution", 800, 600)
    legend_mA = ROOT.TLegend(0.65, 0.65, 0.9, 0.9)
    
    max_y = 0
    projections = []
    
    for i, year_cat in enumerate(year_cats):
        if year_cat in data and sub_cat in data[year_cat]:
            proj = get_1D_projection(data[year_cat][sub_cat], 'y')
            if proj:
                total_events = proj.Integral()
                proj.SetLineColor(colors[i])
                proj.SetLineWidth(2)
                proj.SetMarkerColor(colors[i])
                proj.SetMarkerStyle(20)
                proj.SetMarkerSize(0.8)
                projections.append((proj, total_events, year_cat))
                
                if proj.GetMaximum() > max_y:
                    max_y = proj.GetMaximum()
    
    for idx, (proj, total_events, year_cat) in enumerate(projections):
        if idx == 0:
            proj.Draw("HIST")
            proj.GetYaxis().SetRangeUser(0, max_y*1.3)
        else:
            proj.Draw("HIST SAME")
        legend_mA.AddEntry(proj, f"{year_cat} (N = {total_events:.0f})", "l")
    
    if projections:
        legend_mA.Draw()
        c_mA.Update()
        c_mA.Modified()
        c_mA.SaveAs(f"plots/mA_{sub_cat}.png")
        c_mA.SaveAs(f"plots/mA_{sub_cat}.pdf")
        print(f"Created mA plot for {sub_cat}")

# ============================================================================
# CREATE 2D PLOTS
# ============================================================================
print("\n" + "="*100)
print("CREATING 2D PLOTS")
print("="*100)

# 2D plots: By year-category showing all sub-categories
for year_cat in year_cats:
    if year_cat not in data:
        continue
    
    print(f"\nCreating 2D plots for {year_cat}...")
    
    for sub_cat in sub_cats:
        if sub_cat in data[year_cat]:
            title = f"{year_cat} - {sub_cat}"
            filename = f"{year_cat}_{sub_cat}"
            if create_2D_plot(data[year_cat][sub_cat], title, filename):
                print(f"  Created 2D plot for {filename}")

# 2D plots: By sub-category showing all years
for sub_cat in sub_cats:
    print(f"\nCreating 2D plots for {sub_cat}...")
    
    for year_cat in year_cats:
        if year_cat in data and sub_cat in data[year_cat]:
            title = f"{sub_cat} - {year_cat}"
            filename = f"{sub_cat}_{year_cat}"
            if create_2D_plot(data[year_cat][sub_cat], title, filename):
                print(f"  Created 2D plot for {filename}")

print("\n" + "="*100)
print("Done! Check the 'plots/' directory for output.")
print("="*100)
