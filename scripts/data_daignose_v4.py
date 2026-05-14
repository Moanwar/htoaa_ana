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

# Calculate sideband regions
# For mH sidebands
mH_lower_min, mH_lower_max = 100, 110  # Lower mH sideband
mH_upper_min, mH_upper_max = 140, 160  # Upper mH sideband

# For mA sidebands (exactly 2 bins adjacent, not including the ROI)
mA_left_min = mA_min - 2 * mA_bin_width  # 33 GeV
mA_left_max = mA_min  # 37 GeV
mA_right_min = mA_max  # 45 GeV
mA_right_max = mA_max + 2 * mA_bin_width  # 49 GeV

print(f"\nRegion definitions:")
print(f"Signal region: mH: {mH_min}-{mH_max}, mA: {mA_min}-{mA_max}")
print(f"Lower mH sideband: mH: {mH_lower_min}-{mH_lower_max}, mA: {mA_min}-{mA_max}")
print(f"Upper mH sideband: mH: {mH_upper_min}-{mH_upper_max}, mA: {mA_min}-{mA_max}")
print(f"Left mA sideband: mH: {mH_min}-{mH_max}, mA: {mA_left_min}-{mA_left_max}")
print(f"Right mA sideband: mH: {mH_min}-{mH_max}, mA: {mA_right_min}-{mA_right_max}")

# Dictionary to store all histograms
data = {}

print("\nLoading histograms...")

# Load all histograms
for year in years:
    for sub_cat in sub_cats:
        file_path = f"{base_path}/{year}/VBFjj/2DAlphabet_inputFiles/{sub_cat}/{sub_cat}_Data_{year}.root"
        
        print(f"Trying: {file_path}")
        
        f = ROOT.TFile.Open(file_path)
        if not f or f.IsZombie():
            print(f"  Warning: Could not open {file_path}")
            continue
            
        hist_name = f"{sub_cat}_Data_{year}_pnet_4a_WP40_Pass_Nom"
        hist = f.Get(hist_name)
        
        if hist and not hist.IsZombie():
            year_key = year
            if year_key not in data:
                data[year_key] = {}
            
            cloned_hist = hist.Clone(f"{sub_cat}_{year}_cloned")
            cloned_hist.SetDirectory(0)
            data[year_key][sub_cat] = cloned_hist
            print(f"  Loaded: {year_key} - {sub_cat}")
        else:
            print(f"  Warning: Could not find {hist_name} in file")
        f.Close()

# Combine 2016preVFP and 2016postVFP into "2016"
if "2016preVFP" in data or "2016postVFP" in data:
    data["2016"] = {}
    for sub_cat in sub_cats:
        summed_hist = None
        for era in ["2016preVFP", "2016postVFP"]:
            if era in data and sub_cat in data[era]:
                if summed_hist is None:
                    summed_hist = data[era][sub_cat].Clone(f"2016_{sub_cat}_sum")
                    summed_hist.SetDirectory(0)
                else:
                    summed_hist.Add(data[era][sub_cat])
        
        if summed_hist:
            data["2016"][sub_cat] = summed_hist
            print(f"Combined 2016 - {sub_cat}")

# Create "Run2" by summing all years
data["Run2"] = {}
for sub_cat in sub_cats:
    summed_hist = None
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

# Check loaded data
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
    
    if (x_bin_min < 1 or x_bin_max > hist2D.GetNbinsX() or 
        y_bin_min < 1 or y_bin_max > hist2D.GetNbinsY()):
        return 0
    
    return hist2D.Integral(x_bin_min, x_bin_max, y_bin_min, y_bin_max)

def create_2D_plot_region(hist2D, title, filename, x_min, x_max, y_min, y_max):
    """Create a 2D plot showing only the specified region with bin contents"""
    if not hist2D or hist2D.IsZombie():
        return False
    
    # Find bin ranges
    x_bin_min = hist2D.GetXaxis().FindBin(x_min)
    x_bin_max = hist2D.GetXaxis().FindBin(x_max - 0.001)
    y_bin_min = hist2D.GetYaxis().FindBin(y_min)
    y_bin_max = hist2D.GetYaxis().FindBin(y_max - 0.001)
    
    if (x_bin_min < 1 or x_bin_max > hist2D.GetNbinsX() or 
        y_bin_min < 1 or y_bin_max > hist2D.GetNbinsY()):
        print(f"  Warning: Region {x_min}-{x_max}, {y_min}-{y_max} out of range for {filename}")
        return False
    
    # Get the bin edges as arrays
    x_edges = []
    for i in range(x_bin_min, x_bin_max + 2):
        x_edges.append(hist2D.GetXaxis().GetBinLowEdge(i))
    
    y_edges = []
    for i in range(y_bin_min, y_bin_max + 2):
        y_edges.append(hist2D.GetYaxis().GetBinLowEdge(i))
    
    # Convert to arrays for ROOT
    x_array = ROOT.vector('double')(x_edges)
    y_array = ROOT.vector('double')(y_edges)
    
    # Create new histogram using array constructors
    n_x_bins = len(x_edges) - 1
    n_y_bins = len(y_edges) - 1
    
    region_hist = ROOT.TH2D(f"{filename}_region", title, 
                            n_x_bins, x_array.data(),
                            n_y_bins, y_array.data())
    
    # Fill with content from original
    for i in range(x_bin_min, x_bin_max + 1):
        for j in range(y_bin_min, y_bin_max + 1):
            content = hist2D.GetBinContent(i, j)
            region_hist.SetBinContent(i - x_bin_min + 1, j - y_bin_min + 1, content)
    
    # Create canvas
    c = ROOT.TCanvas(f"c_{filename}", title, 800, 600)
    c.SetRightMargin(0.15)
    
    # Draw with COLZ TEXT to show bin contents
    region_hist.Draw("COLZ TEXT")
    
    # Set titles
    region_hist.GetXaxis().SetTitle("m_{H} (GeV)")
    region_hist.GetYaxis().SetTitle("m_{A} (GeV)")
    
    # Set palette
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    
    # Update and save
    c.Update()
    c.Modified()
    c.SaveAs(f"plots/2D_{filename}.png")
    c.SaveAs(f"plots/2D_{filename}.pdf")
    
    # Clean up
    del region_hist
    
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

year_cats = ["Run2", "2016", "2017", "2018"]

# Print by year
for year_cat in year_cats:
    if year_cat not in data:
        print(f"Warning: {year_cat} not found in data")
        continue
    
    print(f"\n{'='*100}")
    print(f"YEAR: {year_cat}")
    print(f"{'='*100}")
    
    print(f"\n--- SIGNAL REGION (mH: {mH_min}-{mH_max}, mA: {mA_min}-{mA_max}) ---")
    for sub_cat in sub_cats:
        if sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_min, mA_max)
            print(f"  {sub_cat:20s}: {total:10.0f} events")
    
    print(f"\n--- LOWER mH SIDEBAND (mH: {mH_lower_min}-{mH_lower_max}, mA: {mA_min}-{mA_max}) ---")
    for sub_cat in sub_cats:
        if sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_lower_min, mH_lower_max, mA_min, mA_max)
            print(f"  {sub_cat:20s}: {total:10.0f} events")
    
    print(f"\n--- UPPER mH SIDEBAND (mH: {mH_upper_min}-{mH_upper_max}, mA: {mA_min}-{mA_max}) ---")
    for sub_cat in sub_cats:
        if sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_upper_min, mH_upper_max, mA_min, mA_max)
            print(f"  {sub_cat:20s}: {total:10.0f} events")
    
    print(f"\n--- LEFT mA SIDEBAND (mH: {mH_min}-{mH_max}, mA: {mA_left_min}-{mA_left_max}) ---")
    for sub_cat in sub_cats:
        if sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_left_min, mA_left_max)
            print(f"  {sub_cat:20s}: {total:10.0f} events")
    
    print(f"\n--- RIGHT mA SIDEBAND (mH: {mH_min}-{mH_max}, mA: {mA_right_min}-{mA_right_max}) ---")
    for sub_cat in sub_cats:
        if sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_right_min, mA_right_max)
            print(f"  {sub_cat:20s}: {total:10.0f} events")

# Print by sub-category
print(f"\n{'='*100}")
print("BY SUB-CATEGORY")
print("="*100)

for sub_cat in sub_cats:
    print(f"\n{'='*100}")
    print(f"SUB-CATEGORY: {sub_cat}")
    print(f"{'='*100}")
    
    print(f"\n--- SIGNAL REGION (mH: {mH_min}-{mH_max}, mA: {mA_min}-{mA_max}) ---")
    for year_cat in year_cats:
        if year_cat in data and sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_min, mA_max)
            print(f"  {year_cat:10s}: {total:10.0f} events")
    
    print(f"\n--- LOWER mH SIDEBAND (mH: {mH_lower_min}-{mH_lower_max}, mA: {mA_min}-{mA_max}) ---")
    for year_cat in year_cats:
        if year_cat in data and sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_lower_min, mH_lower_max, mA_min, mA_max)
            print(f"  {year_cat:10s}: {total:10.0f} events")
    
    print(f"\n--- UPPER mH SIDEBAND (mH: {mH_upper_min}-{mH_upper_max}, mA: {mA_min}-{mA_max}) ---")
    for year_cat in year_cats:
        if year_cat in data and sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_upper_min, mH_upper_max, mA_min, mA_max)
            print(f"  {year_cat:10s}: {total:10.0f} events")
    
    print(f"\n--- LEFT mA SIDEBAND (mH: {mH_min}-{mH_max}, mA: {mA_left_min}-{mA_left_max}) ---")
    for year_cat in year_cats:
        if year_cat in data and sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_left_min, mA_left_max)
            print(f"  {year_cat:10s}: {total:10.0f} events")
    
    print(f"\n--- RIGHT mA SIDEBAND (mH: {mH_min}-{mH_max}, mA: {mA_right_min}-{mA_right_max}) ---")
    for year_cat in year_cats:
        if year_cat in data and sub_cat in data[year_cat]:
            total = get_region_integral(data[year_cat][sub_cat], mH_min, mH_max, mA_right_min, mA_right_max)
            print(f"  {year_cat:10s}: {total:10.0f} events")

# ============================================================================
# CREATE 1D PLOTS (mH and mA projections) - Keeping these as before
# ============================================================================
print("\n" + "="*100)
print("CREATING 1D PLOTS")
print("="*100)

colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kOrange+7]

def get_1D_projection(hist2D, proj_axis='x'):
    if not hist2D or hist2D.IsZombie():
        return None
    
    x_bins = hist2D.GetXaxis().FindBin(mH_min), hist2D.GetXaxis().FindBin(mH_max-0.001)
    y_bins = hist2D.GetYaxis().FindBin(mA_min), hist2D.GetYaxis().FindBin(mA_max-0.001)
    
    if x_bins[0] < 1 or x_bins[1] > hist2D.GetNbinsX() or y_bins[0] < 1 or y_bins[1] > hist2D.GetNbinsY():
        return None
    
    try:
        if proj_axis == 'x':
            proj = hist2D.ProjectionX(f"{hist2D.GetName()}_projX", y_bins[0], y_bins[1])
            proj.GetXaxis().SetTitle("m_{H} (GeV)")
            proj.GetYaxis().SetTitle("Events")
            proj.SetTitle(f"mH distribution (mA: {mA_min}-{mA_max} GeV)")
            proj.GetXaxis().SetRange(x_bins[0], x_bins[1])
        else:
            proj = hist2D.ProjectionY(f"{hist2D.GetName()}_projY", x_bins[0], x_bins[1])
            proj.GetXaxis().SetTitle("m_{A} (GeV)")
            proj.GetYaxis().SetTitle("Events")
            proj.SetTitle(f"mA distribution (mH: {mH_min}-{mH_max} GeV)")
            proj.GetXaxis().SetRange(y_bins[0], y_bins[1])
        
        proj.SetDirectory(0)
        return proj
    except:
        return None

# 1D plots by year
for year_cat in year_cats:
    if year_cat not in data:
        continue
    
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
        c_mH.SaveAs(f"plots/mH_{year_cat}.png")
        print(f"Created mH plot for {year_cat}")
    
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
        c_mA.SaveAs(f"plots/mA_{year_cat}.png")
        print(f"Created mA plot for {year_cat}")

# 1D plots by sub-category
for sub_cat in sub_cats:
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
        c_mH.SaveAs(f"plots/mH_{sub_cat}.png")
        print(f"Created mH plot for {sub_cat}")
    
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
        c_mA.SaveAs(f"plots/mA_{sub_cat}.png")
        print(f"Created mA plot for {sub_cat}")

# ============================================================================
# CREATE 2D PLOTS FOR EACH REGION
# ============================================================================
print("\n" + "="*100)
print("CREATING 2D PLOTS")
print("="*100)

# Define all regions to plot
regions = [
    ("Signal", mH_min, mH_max, mA_min, mA_max),
    ("Lower_mH_Sideband", mH_lower_min, mH_lower_max, mA_min, mA_max),
    ("Upper_mH_Sideband", mH_upper_min, mH_upper_max, mA_min, mA_max),
    ("Left_mA_Sideband", mH_min, mH_max, mA_left_min, mA_left_max),
    ("Right_mA_Sideband", mH_min, mH_max, mA_right_min, mA_right_max)
]

# Create 2D plots for each year and sub-category
for year_cat in year_cats:
    if year_cat not in data:
        continue
    
    print(f"\nCreating 2D plots for {year_cat}...")
    
    for sub_cat in sub_cats:
        if sub_cat not in data[year_cat]:
            continue
        
        for region_name, x_min, x_max, y_min, y_max in regions:
            title = f"{year_cat} - {sub_cat} - {region_name}"
            filename = f"{year_cat}_{sub_cat}_{region_name}"
            
            if create_2D_plot_region(data[year_cat][sub_cat], title, filename, x_min, x_max, y_min, y_max):
                print(f"  Created: {filename}")

print("\n" + "="*100)
print("Done! Check the 'plots/' directory for output.")
print("="*100)
