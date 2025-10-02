import uproot
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_systematics_2D(rootfile, outdir="plots", histname="hLeadingFatJetPNet_massH_v2b_vs_massA34a_VBFIncl_Xto4bv2_SRWP40_Nom"):
    f = uproot.open(rootfile)
    sample_dir = f["evt/VBFHtoaato4b_mA_40"]

    # --- Nominal histogram (2D) ---
    h2_nom = sample_dir[histname]
    nom_values_2d, x_edges, y_edges = h2_nom.to_numpy()
    # Project onto x-axis (massH)
    nom_values = np.sum(nom_values_2d, axis=1)
    bin_centers = 0.5 * (x_edges[1:] + x_edges[:-1])

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

        up_values = np.sum(h2_up.to_numpy()[0], axis=1)
        dn_values = np.sum(h2_dn.to_numpy()[0], axis=1)

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
        rax.set_xlabel("Leading FatJet PNet massH [GeV]")

        # Save
        syst_label = syst_base.replace(prefix, "")
        plt.tight_layout()
        outname = os.path.join(outdir, f"{syst_label}.png")
        plt.savefig(outname)
        plt.close()
        print(f"Saved {outname}")


# Example usage
plot_systematics_2D(
    rootfile="analyze_htoaa_SUSY_VBFH_HToAATo4B_Pt150_M-40_TuneCP5_13TeV_madgraph_pythia8_0_0.root",
    outdir="plots_massH_1D",
    histname="hLeadingFatJetPNet_massH_v2b_vs_massA34a_VBFIncl_Xto4bv2_SRWP40_Nom"
)
