import uproot
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_systematics_2D(rootfile, outdir="plots", histname="hLeadingFatJetPNet_massH_v2b_vs_massA34a_VBFIncl_Xto4bv2_SRWP40_Nom"):
    f = uproot.open(rootfile)

    #if "evt" in f.keys():
    # e.g. f["evt/VBFHtoaato4b_mA_40"]
    #sample_dir = f["evt/ggHtoaato4b_mA_30"]
    sample_dir = f["evt/VBFHtoaato4b_mA_30"]

    #else:
    # No directory: histograms live at the top level
    #sample_dir = f
        
    # Now you can safely do:
    # --- Nominal histogram (2D) ---
    h2_nom = sample_dir[histname]
    nom_values_2d, x_edges, y_edges = h2_nom.to_numpy()
    
    # Bin centers
    x_centers = 0.5 * (x_edges[1:] + x_edges[:-1])
    y_centers = 0.5 * (y_edges[1:] + y_edges[:-1])

    # Find all systematics with same prefix
    prefix = histname.replace("Nom", "")
    syst_hists = [k for k in sample_dir.keys() if k.startswith(prefix) and not k.endswith("Nom;1")]

    os.makedirs(outdir, exist_ok=True)

    for syst in syst_hists:
        if not ("Up" in syst or "Down" in syst):
            continue

        syst_tag = syst.replace(";1", "")
        syst_base = syst_tag[:-2] if syst_tag.endswith(("Up","Down")) else None
        if syst_base is None:
            continue

        h2_up = sample_dir.get(syst_base + "Up")
        h2_dn = sample_dir.get(syst_base + "Down")
        if h2_up is None or h2_dn is None:
            continue

        up_values_2d = h2_up.to_numpy()[0]
        dn_values_2d = h2_dn.to_numpy()[0]

        # --- X-axis projection (massH) ---
        nom_x = np.sum(nom_values_2d, axis=1)
        up_x  = np.sum(up_values_2d, axis=1)
        dn_x  = np.sum(dn_values_2d, axis=1)

        fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios":[3,1]}, sharex=True)
        ax.step(x_centers, nom_x, where="mid", label="Nominal", color="black")
        ax.step(x_centers, up_x, where="mid", label="Up", color="red")
        ax.step(x_centers, dn_x, where="mid", label="Down", color="blue")
        ax.set_ylabel("Events")
        ax.legend()
        ratio_up = np.divide(up_x, nom_x, out=np.ones_like(up_x), where=nom_x!=0)
        ratio_dn = np.divide(dn_x, nom_x, out=np.ones_like(dn_x), where=nom_x!=0)
        rax.axhline(1.0, color="black", linestyle="--")
        rax.step(x_centers, ratio_up, where="mid", color="red")
        rax.step(x_centers, ratio_dn, where="mid", color="blue")
        rax.set_ylabel("Ratio")
        rax.set_xlabel("Leading FatJet PNet massH [GeV]")
        syst_label = syst_base.replace(prefix, "")
        plt.tight_layout()
        outname_x = os.path.join(outdir, f"{syst_label}_massH.png")
        plt.savefig(outname_x)
        plt.close()
        print(f"Saved {outname_x}")

        # --- Y-axis projection (massA) ---
        nom_y = np.sum(nom_values_2d, axis=0)
        up_y  = np.sum(up_values_2d, axis=0)
        dn_y  = np.sum(dn_values_2d, axis=0)

        fig, (ay, ray) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios":[3,1]}, sharex=True)
        ay.step(y_centers, nom_y, where="mid", label="Nominal", color="black")
        ay.step(y_centers, up_y, where="mid", label="Up", color="red")
        ay.step(y_centers, dn_y, where="mid", label="Down", color="blue")
        ay.set_ylabel("Events")
        ay.legend()
        ratio_up_y = np.divide(up_y, nom_y, out=np.ones_like(up_y), where=nom_y!=0)
        ratio_dn_y = np.divide(dn_y, nom_y, out=np.ones_like(dn_y), where=nom_y!=0)
        ray.axhline(1.0, color="black", linestyle="--")
        ray.step(y_centers, ratio_up_y, where="mid", color="red")
        ray.step(y_centers, ratio_dn_y, where="mid", color="blue")
        ray.set_ylabel("Ratio")
        ray.set_xlabel("Leading FatJet PNet massA [GeV]")
        syst_label = syst_base.replace(prefix, "")
        plt.tight_layout()
        outname_y = os.path.join(outdir, f"{syst_label}_massA.png")
        plt.savefig(outname_y)
        plt.close()
        print(f"Saved {outname_y}")
plot_systematics_2D(
    rootfile="analyze_htoaa_SUSY_VBFH_HToAATo4B_Pt150_M-30_TuneCP5_13TeV_madgraph_pythia8_0_0.root",
    outdir="plots_massH_1D_massA_1D",
    histname="hLeadingFatJetPNet_massH_v2b_vs_massA34a_VBFIncl_Xto4bv2_SRWP40_Nom"
)
#hLeadingFatJetPNet_massH_v2b_vs_massA34a_VBFIncl_Xto4bv2_SRWP40_Nom
#hLeadingFatJetPNet_massH_v2b_vs_massA34a_gg0lHi_Xto4bv2_SRWP40_Nom

#"/afs/cern.ch/work/m/moanwar/public/hto2ato4b/2DAlphabetfiles_VBF_sys_v2/20251021_DataMC/2018/VBFjj/2DAlphabet_inputFiles/VBFHiPTHi/VBFHiPTHi_VBFHtoaato4b_mA_50_2018.root"
#    histname="VBFHiPTHi_VBFHtoaato4b_mA_50_2018_pnet_34a_WP40_Pass_Nom"
#"/afs/cern.ch/work/m/moanwar/private/siddhesh86/analysis/20251021_DataMC/2018/VBFjj/analyze_htoaa_stage1.root"
