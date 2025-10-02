import uproot
import matplotlib.pyplot as plt
import numpy as np
import os
def plot_systematics(rootfile, outdir="plots", histname="hLeadingFatJetPt_VBFIncl_Nom"):
    f = uproot.open(rootfile)

    # Navigate into dirs
    sample_dir = f["evt/VBFHtoaato4b_mA_40"]

    # Nominal histogram
    h_nom = sample_dir[histname]
    nom_values, nom_edges = h_nom.to_numpy()
    bin_centers = 0.5 * (nom_edges[1:] + nom_edges[:-1])

    # Find all histograms with same prefix
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

        # fetch from sample_dir, not f
        h_up = sample_dir.get(syst_base + "Up")
        h_dn = sample_dir.get(syst_base + "Down")
        if h_up is None or h_dn is None:
            continue

        up_values, _ = h_up.to_numpy()
        dn_values, _ = h_dn.to_numpy()

        # --- Plot ---
        fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios":[3,1]}, sharex=True)

        # Main pad
        ax.step(bin_centers, nom_values, where="mid", label="Nominal", color="black")
        ax.step(bin_centers, up_values, where="mid", label="Up", color="red")
        ax.step(bin_centers, dn_values, where="mid", label="Down", color="blue")
        ax.set_ylabel("Events")
        ax.legend()

        # Ratio pad
        ratio_up = np.divide(up_values, nom_values, out=np.ones_like(up_values), where=nom_values!=0)
        ratio_dn = np.divide(dn_values, nom_values, out=np.ones_like(dn_values), where=nom_values!=0)
        rax.axhline(1.0, color="black", linestyle="--")
        rax.step(bin_centers, ratio_up, where="mid", color="red")
        rax.step(bin_centers, ratio_dn, where="mid", color="blue")
        rax.set_ylabel("Ratio")
        rax.set_xlabel("Leading FatJet $p_T$ [GeV]")

        # Save
        syst_label = syst_base.replace(prefix, "")
        plt.tight_layout()
        outname = os.path.join(outdir, f"{syst_label}.png")
        plt.savefig(outname)
        plt.close()
        print(f"Saved {outname}")


# Example usage
plot_systematics(
    rootfile="analyze_htoaa_SUSY_VBFH_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8_0_0.root",
    outdir="plots_VBFHtoaato4b_mA_40",
    histname="hLeadingFatJetPt_VBFIncl_Nom"
)
