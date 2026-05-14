import ROOT
import os

# Define paths and categories
base_path = "/afs/cern.ch/work/m/moanwar/public/hto2ato4b/2DAlphabetfiles_VBF_sys_v2/20251021_DataMC"
years = ["2016preVFP", "2016postVFP", "2017", "2018"]
sub_cats = ["VBFHiPTHi", "VBFHiPTLo", "VBFLoPTHi", "VBFLoPTLo"]

# Region of interest
mH_min, mH_max = 110, 140
mA_min, mA_max = 37, 45

# Dictionary to store all histograms
# Structure: data[year_cat][sub_cat] = 2D histogram
data = {}

print("Loading histograms...")

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

# Function to create 1D projection from 2D histogram in ROI
def get_1D_projectionv(hist2D, proj_axis='x'):
    """
    Create 1D projection from 2D histogram for the region mH:110-140, mA:37-45
    proj_axis='x' for mH projection (sum over mA)
    proj_axis='y' for mA projection (sum over mH)
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
            proj = hist2D.ProjectionX(f"{hist2D.GetName()}_projX", y_bins[0], y_bins[1])
            proj.GetXaxis().SetTitle("m_{H} (GeV)")
            proj.GetYaxis().SetTitle("Events")
            proj.SetTitle(f"mH distribution (mA: {mA_min}-{mA_max} GeV)")
        else:
            proj = hist2D.ProjectionY(f"{hist2D.GetName()}_projY", x_bins[0], x_bins[1])
            proj.GetXaxis().SetTitle("m_{A} (GeV)")
            proj.GetYaxis().SetTitle("Events")
            proj.SetTitle(f"mA distribution (mH: {mH_min}-{mH_max} GeV)")
        
        proj.SetDirectory(0)  # Detach from any file
        return proj
    except:
        print(f"  Warning: Failed to project histogram {hist2D.GetName()}")
        return None

# Create output directory
os.system("mkdir -p plots")

# Set style
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kViridis)

# Group 1: By year-category (Run2, 2016, 2017, 2018)
year_cats = ["Run2", "2016", "2017", "2018"]
colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kOrange+7]

for year_cat in year_cats:
    if year_cat not in data:
        print(f"Warning: {year_cat} not found in data")
        continue
    
    # Create mH plot for this year-category
    c_mH = ROOT.TCanvas(f"c_mH_{year_cat}", f"{year_cat} - mH Distribution", 800, 600)
    legend_mH = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    
    first = True
    max_y = 0
    projections = []
    
    for i, sub_cat in enumerate(sub_cats):
        if sub_cat in data[year_cat]:
            proj = get_1D_projection(data[year_cat][sub_cat], 'x')
            if proj:
                proj.SetLineColor(colors[i])
                proj.SetLineWidth(2)
                proj.SetMarkerColor(colors[i])
                proj.SetMarkerStyle(20)
                proj.SetMarkerSize(0.8)
                projections.append(proj)
                
                if proj.GetMaximum() > max_y:
                    max_y = proj.GetMaximum()
    
    # Draw all projections
    for i, proj in enumerate(projections):
        if i == 0:
            #proj.Draw("HIST E")
            proj.Draw("HIST")

            proj.GetYaxis().SetRangeUser(0, max_y*1.3)
        else:
            #proj.Draw("HIST E SAME")
            proj.Draw("HIST SAME")
        legend_mH.AddEntry(proj, sub_cats[i].replace("VBF", ""), "l")
    
    if projections:  # Only draw legend if we have projections
        legend_mH.Draw()
        c_mH.Update()
        c_mH.SaveAs(f"plots/mH_{year_cat}.png")
        print(f"Created mH plot for {year_cat}")
    
    # Create mA plot for this year-category
    c_mA = ROOT.TCanvas(f"c_mA_{year_cat}", f"{year_cat} - mA Distribution", 800, 600)
    legend_mA = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    
    first = True
    max_y = 0
    projections = []
    
    for i, sub_cat in enumerate(sub_cats):
        if sub_cat in data[year_cat]:
            proj = get_1D_projection(data[year_cat][sub_cat], 'y')
            if proj:
                proj.SetLineColor(colors[i])
                proj.SetLineWidth(2)
                proj.SetMarkerColor(colors[i])
                proj.SetMarkerStyle(20)
                proj.SetMarkerSize(0.8)
                projections.append(proj)
                
                if proj.GetMaximum() > max_y:
                    max_y = proj.GetMaximum()
    
    for i, proj in enumerate(projections):
        if i == 0:
            #proj.Draw("HIST E")
            proj.Draw("HIST")
            proj.GetYaxis().SetRangeUser(0, max_y*1.3)
        else:
            #proj.Draw("HIST E SAME")
            proj.Draw("HIST SAME")
        legend_mA.AddEntry(proj, sub_cats[i].replace("VBF", ""), "l")
    
    if projections:
        legend_mA.Draw()
        c_mA.Update()
        c_mA.SaveAs(f"plots/mA_{year_cat}.png")
        print(f"Created mA plot for {year_cat}")

# Group 2: By sub-category
for sub_cat in sub_cats:
    # Create mH plot for this sub-category
    c_mH = ROOT.TCanvas(f"c_mH_{sub_cat}", f"{sub_cat} - mH Distribution", 800, 600)
    legend_mH = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    
    first = True
    max_y = 0
    projections = []
    
    for i, year_cat in enumerate(year_cats):
        if year_cat in data and sub_cat in data[year_cat]:
            proj = get_1D_projection(data[year_cat][sub_cat], 'x')
            if proj:
                proj.SetLineColor(colors[i])
                proj.SetLineWidth(2)
                proj.SetMarkerColor(colors[i])
                proj.SetMarkerStyle(20)
                proj.SetMarkerSize(0.8)
                projections.append(proj)
                
                if proj.GetMaximum() > max_y:
                    max_y = proj.GetMaximum()
    
    for i, proj in enumerate(projections):
        if i == 0:
            #proj.Draw("HIST E")
            proj.Draw("HIST")
            proj.GetYaxis().SetRangeUser(0, max_y*1.3)
        else:
            #proj.Draw("HIST E SAME")
            proj.Draw("HIST SAME")
        legend_mH.AddEntry(proj, year_cats[i], "l")
    
    if projections:
        legend_mH.Draw()
        c_mH.Update()
        c_mH.SaveAs(f"plots/mH_{sub_cat}.png")
        print(f"Created mH plot for {sub_cat}")
    
    # Create mA plot for this sub-category
    c_mA = ROOT.TCanvas(f"c_mA_{sub_cat}", f"{sub_cat} - mA Distribution", 800, 600)
    legend_mA = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    
    first = True
    max_y = 0
    projections = []
    
    for i, year_cat in enumerate(year_cats):
        if year_cat in data and sub_cat in data[year_cat]:
            proj = get_1D_projection(data[year_cat][sub_cat], 'y')
            if proj:
                proj.SetLineColor(colors[i])
                proj.SetLineWidth(2)
                proj.SetMarkerColor(colors[i])
                proj.SetMarkerStyle(20)
                proj.SetMarkerSize(0.8)
                projections.append(proj)
                
                if proj.GetMaximum() > max_y:
                    max_y = proj.GetMaximum()
    
    for i, proj in enumerate(projections):
        if i == 0:
            #proj.Draw("HIST E")
            proj.Draw("HIST")
            proj.GetYaxis().SetRangeUser(0, max_y*1.3)
        else:
            #proj.Draw("HIST E SAME")
            proj.Draw("HIST SAME")
        legend_mA.AddEntry(proj, year_cats[i], "l")
    
    if projections:
        legend_mA.Draw()
        c_mA.Update()
        c_mA.SaveAs(f"plots/mA_{sub_cat}.png")
        print(f"Created mA plot for {sub_cat}")

print("\nDone! Check the 'plots/' directory for output.")
