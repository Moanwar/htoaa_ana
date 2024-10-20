import uproot
import awkward as ak
import vector
vector.register_awkward()
from matplotlib import pyplot as plt
import mplhep as hep
hep.style.use("CMS")
import math
import itertools
import os
from typing import Dict, List, Tuple
from numpy.typing import ArrayLike
import numpy as np
import correctionlib
import awkward as ak
import fastjet
#from vector import LorentzVector
from coffea.nanoevents.methods import vector
from coffea import nanoevents
from coffea import processor
#from coffea.nanoevents.methods import candidate
from coffea.analysis_tools import Weights, PackedSelection
from hist import Hist
ak.behavior.update(vector.behavior)
import sys
sys.path.append('/afs/cern.ch/work/m/moanwar/private/myAnalysis_htoaa/htoaa_ana/LundReweighting/utils')
#from LundReweighter import *
from LundReweighting.utils.Utils import *
from Lp_Htoaa4b_tools import (selectMETFilters)
from Lp_Htoaa4b_vars import *

from coffea.nanoevents import NanoEventsFactory
#root_file_path = "/eos/user/m/moanwar/htoaa/analysis/VBF_PFoutput/HToAATo4B_15/PNet_v1_skim_Hto4b_0p8_VBFH_M-15.root"
#root_file_path = "/eos/user/m/moanwar/htoaa/analysis/VBF_PFoutput/HToAATo4B_30/PNet_v4_skim_Hto4b_0p8_VBFH_M-30.root"
root_file_path = "/eos/user/m/moanwar/htoaa/analysis/VBF_PFoutput/HToAATo4B_55/PNet_v1_skim_Hto4b_0p8_VBFH_M-55.root"
events = nanoevents.NanoEventsFactory.from_root(
    root_file_path,
    schemaclass=nanoevents.NanoAODSchema,
).events()

# Initialize Event Data
Era, isMC, sTrgSelection = '2018', True, 'Trg_Combo_AK4AK8Jet_HT'
frac=1
events = events[:int(frac * len(events))]  # Use full event fraction

# Filter desired fields
desired_fields = ["Electron", "Muon", "FatJet", "genWeight", "GenPart", "FatJetPFCands", "PFCands", "PV", "Flag", "HLT", "L1"]
filtered_events = events[desired_fields]

# Pad array to target length with a given value
def pad_val(arr, target, value, axis=0, to_numpy=True, clip=True):
    padded_arr = ak.fill_none(ak.pad_none(arr, target, axis=axis, clip=clip), value, axis=axis)
    return padded_arr.to_numpy() if to_numpy else padded_arr

# Add selections and update cutflow
def add_selection(name, sel, selection, cutflow=None, isData=False, signGenWeights=None):
    selection.add(name, sel)
    if cutflow is not None:
        cutflow[name] = np.sum(selection.all(*selection.names)) if isData else np.sum(signGenWeights[selection.all(*selection.names)])

# Initialize MC-related variables
isData = False
signGenWeights = None if isData else np.sign(filtered_events["genWeight"])
n_events = len(filtered_events) if isData else int(np.sum(signGenWeights))
selection = PackedSelection()  # Initialize selection object
cutflow = {"all": n_events}     # Store the total number of events

# Preselection Cuts
preselection_cut_0 = pad_val(filtered_events.PV.npvsGood >= 1, len(events), False)
preselection_cut_1 = pad_val(selectMETFilters(filtered_events.Flag, Era, isMC), len(events), False)
add_selection("Good PV and MET Filters", preselection_cut_0 & preselection_cut_1, selection, cutflow, isData, signGenWeights)

# FatJet Cuts
fatjets = filtered_events.FatJet
leadingFatJet = ak.firsts(fatjets)
preselection_cut_vals = {"pt_Thsh": 250, "eta_Thsh": 2.4, "min_msd": 50, "max_msd": 9999, "FatJetZHbb_Thsh": 0.7}

# Combine FatJet pT, eta, jetId, and mass soft drop cuts
preselection_fatjet = pad_val(
    (leadingFatJet.pt > preselection_cut_vals["pt_Thsh"]) &
    (abs(leadingFatJet.eta) < preselection_cut_vals["eta_Thsh"]) &
    (leadingFatJet.jetId == int(JetIDs.tightIDPassingLeptonVeto)) &
    (leadingFatJet.msoftdrop > preselection_cut_vals["min_msd"]) &
    (leadingFatJet.msoftdrop < preselection_cut_vals["max_msd"]), 
    len(events), False
)
add_selection("FatJet Cuts", preselection_fatjet.astype(bool), selection, cutflow, isData, signGenWeights)

# Trigger Selection
mask_Trgs = np.full(len(events), False)
for HLTName, L1TList in Triggers_perEra[Era][sTrgSelection].items():
    mask_HLT = events.HLT[HLTName.replace('HLT_', '')]
    mask_L1Ts = np.any([events.L1[L1.replace('L1_', '')] for L1 in L1TList], axis=0)
    mask_Trg_i = (mask_HLT & mask_L1Ts)
    mask_Trgs  = (mask_Trgs | mask_Trg_i) 

preselection_cut_trig = pad_val(mask_Trgs, len(events), False)
add_selection("Trigger Selection", preselection_cut_trig.astype(bool), selection, cutflow, isData, signGenWeights)

# Handle invalid deepTagMD values for ZHbb and ZHcc
ZHbb_tag = ak.where(leadingFatJet.deepTagMD_ZHbbvsQCD >= 0, leadingFatJet.deepTagMD_ZHbbvsQCD, 0)
ZHcc_tag = ak.where(leadingFatJet.deepTagMD_ZHccvsQCD >= 0, leadingFatJet.deepTagMD_ZHccvsQCD, 0)

# Calculate corrected ZHbb tagging score
corrected_ZHbb_tag = ZHbb_tag * (1 - ZHcc_tag) / (1 - (ZHbb_tag * ZHcc_tag))

# Apply selection based on ZHbb tagging threshold
preselection_ZHbb = pad_val(corrected_ZHbb_tag > preselection_cut_vals["FatJetZHbb_Thsh"], len(events), False)

# Add selection to cutflow
add_selection("leadingFatJetZHbb", preselection_ZHbb.astype(bool), selection, cutflow, isData, signGenWeights)

# Define and filter tight muons and electrons  
muons = ak.with_field(filtered_events.Muon, 0, "flavor")
electrons = ak.with_field(filtered_events.Electron, 1, "flavor")
tight_muons = (muons.pt > 10) & (abs(muons.eta) < 2.4) & ((muons.mediumPromptId) | (muons.highPtId > 0)) & (muons.miniPFRelIso_all <= 0.1) & (abs(muons.dxy) < 0.2) & (abs(muons.dz) < 0.5)
tight_electrons = (electrons.pt > 10) & ((abs(electrons.eta) < 1.44) | (abs(electrons.eta) > 1.57)) & (abs(electrons.eta) < 2.5) & (electrons.mvaFall17V2noIso_WP90)
# Count leptons and apply cut
n_leptons = ak.sum(tight_muons, axis=1) + ak.sum(tight_electrons, axis=1)
add_selection("no isolated tight leptons", pad_val(n_leptons == 0, len(filtered_events), False).astype(bool), selection, cutflow, isData, signGenWeights)

print(cutflow)
### Gen Inormation ###

# Constants for Particle IDs and Flags
HIGGS_PDGID = 25
A_BOSON_PDGID = 36
b_PDGID = 5
GEN_FLAGS = ["fromHardProcess", "isLastCopy"]
FILL_NONE_VALUE = -99999
PAD_VAL = -99999

# Skim variables for Gen Higgs and A boson 4-vectors
skim_vars = {"eta": "Eta", "phi": "Phi", "mass": "Mass", "pt": "Pt", "pdgId": "pdgId"}

# Select Higgs bosons
higgs = filtered_events.GenPart[
    (abs(filtered_events.GenPart.pdgId) == HIGGS_PDGID) & filtered_events.GenPart.hasFlags(GEN_FLAGS)
]
# Store Higgs 4-vector information
GenHiggsVars = {f"GenHiggs{key}": higgs[var].to_numpy() for var, key in skim_vars.items()}
GenHiggsVars["GenHiggsChildren"] = abs(higgs.children.pdgId[:, :, 0]).to_numpy()

# Check for H -> AA decay where A bosons exist
is_A_boson = abs(higgs.children.pdgId) == A_BOSON_PDGID
has_HtoAA = ak.sum(ak.flatten(is_A_boson, axis=2), axis=1) == 2
add_selection("HtoAA", has_HtoAA, selection, cutflow, isData, signGenWeights)

# Extract A boson children (H -> AA -> ...)
AA_bosons = ak.flatten(higgs.children[is_A_boson], axis=2)
GenAAVars = {f"GenAA{key}": AA_bosons[var].to_numpy() for var, key in skim_vars.items()}

# Check if both A bosons decay to 2 b-quarks (H -> AA -> 4b)
AA_children = AA_bosons.children
is_b_quark = abs(AA_children.pdgId) == b_PDGID
decays_to_4b = ak.all(ak.all(is_b_quark, axis=2), axis=1)
add_selection("AA decays to 4b", decays_to_4b, selection, cutflow, isData, signGenWeights)

# Apply selections and display filtered events and cutflow
events_after_cut = filtered_events[selection.all(*selection.names)]

#find the Higgs jet idx, dR(gen_Higgs, jet)<0.8 & dR(gen_A_1, jet)<0.8 & dR(gen_A_2, jet)<0.8

# Extract Gen Higgs particles after previous cuts
higgs = events_after_cut.GenPart[
    (abs(events_after_cut.GenPart.pdgId) == HIGGS_PDGID) & events_after_cut.GenPart.hasFlags(GEN_FLAGS)
]
GenHiggsVars = {f"GenHiggs{key}": higgs[var].to_numpy() for var, key in skim_vars.items()}

# Check for A boson (H -> AA decay)
higgs_children = higgs.children
is_A_boson = abs(higgs_children.pdgId) == A_BOSON_PDGID
A_bosons = ak.flatten(higgs_children[is_A_boson], axis=2)

# Store A boson information for future use
GenAAVars = {f"GenAA{key}": pad_val(A_bosons[var], 2, FILL_NONE_VALUE, axis=1) for var, key in skim_vars.items()}
GenA1Vars = {f"GenAA1{key}": pad_val(A_bosons[var], 2, FILL_NONE_VALUE, axis=1)[:, 0] for var, key in skim_vars.items()}
GenA2Vars = {f"GenAA2{key}": pad_val(A_bosons[var], 2, FILL_NONE_VALUE, axis=1)[:, 1] for var, key in skim_vars.items()}

# Extract and match FatJet with Gen Higgs (H -> AA -> 4b) within deltaR < 0.8
fatjets = events_after_cut.FatJet
HAA = ak.pad_none(higgs, 1, axis=1, clip=True)[:, 0]  # Pad to keep the first Higgs particle
AAdr = fatjets[:, :1].delta_r(HAA)  # Calculate deltaR between Higgs and our leading FatJet
match_dR = 0.6
HAA_match = AAdr <= match_dR  # Check if the FatJet matches within dR threshold

# Store matched FatJet information
HAAJets = ak.pad_none(fatjets[HAA_match], 1, axis=1)[:, 0]  # Get the first matched FatJet
filled_pt = pad_val(HAAJets.pt, len(HAAJets.pt), FILL_NONE_VALUE)

# Apply final selection on matched HAA Jets
selection_has_HAAjet = filled_pt > 0
events_final = events_after_cut[selection_has_HAAjet]

# ==================== Gen-Higgs and AA Information After Final Cuts ===================

# Collect Gen Higgs Information
higgs = events_final.GenPart[
    (abs(events_final.GenPart.pdgId) == HIGGS_PDGID) * events_final.GenPart.hasFlags(GEN_FLAGS)
]
GenHiggsVars = {
    f"GenHiggs{key}": higgs[var].to_numpy() for (var, key) in skim_vars.items()
}
higgs_children = higgs.children

# Collect Gen-Higgs children AA information, even if it might not be needed for matching
is_AA = abs(higgs_children.pdgId) == A_BOSON_PDGID
AA = ak.flatten(higgs_children[is_AA], axis=2)

# Save GenAA variables for Higgs children (A1, A2)
GenAAVars = {
    f"GenAA{key}": AA[var].to_numpy() for (var, key) in skim_vars.items()
}
GenA1Vars = {
    f"GenAA{key}": AA[var][:, 0].to_numpy() for (var, key) in skim_vars.items()
}
GenA2Vars = {
    f"GenAA{key}": AA[var][:, 1].to_numpy() for (var, key) in skim_vars.items()
}

# ============================ FatJet Information ============================

# Collect FatJet information
fatjets = events_final.FatJet
HAA = ak.pad_none(higgs, 1, axis=1, clip=True)[:, 0]  # Pad and clip Gen Higgs array
AAdr = fatjets[:, :1].delta_r(HAA)  # Calculate delta R between leadingfatjet and Higgs
HAA_match = AAdr <= match_dR  # Match FatJets to Higgs within delta R
HAAJets = ak.pad_none(fatjets[HAA_match], 1, axis=1)[:, 0]  # Handle None values in HAAJets array

# ====================== 4q (AA Children) Information =========================

# Collect AA children (4 quarks) information
AA_children = AA.children
quarks = abs(AA_children.pdgId) == b_PDGID  # Check if AA children are b-quarks

# Save Gen4q variables (for the 4 quarks) and pad where necessary
Gen4qVars = {
    f"Gen4q{key}": ak.to_numpy(
        ak.fill_none(
            ak.pad_none(ak.pad_none(AA_children[var], 2, axis=1, clip=True), 2, axis=2, clip=True),
            FILL_NONE_VALUE,
        )
    )
    for (var, key) in skim_vars.items()
}

# ======================= Matching Jets and Quarks =========================

# Stack the FatJet's pt, eta, phi, mass into a 4-vector array for HAA jets
higgs_jet_4vec = np.array(np.stack(
    (np.array(HAAJets.pt), np.array(HAAJets.eta), np.array(HAAJets.phi), np.array(HAAJets.mass)), axis=1
))

# Reshape Gen4q variables for matching
eta = Gen4qVars["Gen4qEta"].reshape(Gen4qVars["Gen4qEta"].shape[0], -1)
phi = Gen4qVars["Gen4qPhi"].reshape(Gen4qVars["Gen4qPhi"].shape[0], -1)
pdI = Gen4qVars["Gen4qpdgId"].reshape(Gen4qVars["Gen4qpdgId"].shape[0], -1)

# Create arrays with eta-phi and pdgId information for quark matching
gen_parts_eta_phi_HAA = np.array(np.dstack((eta, phi)))
gen_parts_pdg_ids_HAA = np.array(np.dstack((pdI))).T

# Function to count the number of quarks inside AK8 jets (within a given delta R)
def count_quarks_in_jets(jet_4vec, gen_parts_eta_phi, delta_r_cut=0.8):
    num_jets = len(jet_4vec)
    num_quarks_in_jets = np.zeros(num_jets, dtype=int)

    for i in range(num_jets):
        jet_eta, jet_phi = jet_4vec[i][1], jet_4vec[i][2]

        quark_eta_phi = gen_parts_eta_phi[i]
        quark_eta, quark_phi = quark_eta_phi[:, 0], quark_eta_phi[:, 1]

        delta_eta = jet_eta - quark_eta
        delta_phi = jet_phi - quark_phi

        delta_r_squared = delta_eta**2 + delta_phi**2
        quarks_in_jet = np.sqrt(delta_r_squared) < delta_r_cut

        num_quarks_in_jets[i] = np.sum(quarks_in_jet)

    return num_quarks_in_jets

# Apply the quark-counting function to HAA jets and store the result
result = count_quarks_in_jets(higgs_jet_4vec, gen_parts_eta_phi_HAA)

# ==================== Final Cut Selection ====================

H3q4q_cut = result >= 3

# Get HAA jet index and FatJetPFCands 4-vectors
HAA_match = AAdr <= match_dR
HAA_jet_idx = np.argmax(pad_val(HAA_match, 3, False, 1, True), axis=1)

HAA_FatJetPFCands = events_final.FatJetPFCands.jetIdx == HAA_jet_idx
HAA_FatJetPFCands_pFCandsIdx = events_final.FatJetPFCands.pFCandsIdx[HAA_FatJetPFCands]

# Collect PFCands 4-vector
pt_array, eta_array, phi_array, mass_array = (ak.Array(events_final.PFCands[var]) for var in ["pt", "eta", "phi", "mass"])
selected_pt = pt_array[HAA_FatJetPFCands_pFCandsIdx]

selected_arrays = [pad_val(selected_pt, 150, 0, 1, True) for selected_pt in [selected_pt, eta_array[HAA_FatJetPFCands_pFCandsIdx],
    phi_array[HAA_FatJetPFCands_pFCandsIdx], mass_array[HAA_FatJetPFCands_pFCandsIdx]]]

selected_pt_padded, selected_eta_padded, selected_phi_padded, selected_mass_padded = selected_arrays

# Calculate PFCands 4-vector
pf_cands_px = selected_pt_padded * np.cos(selected_phi_padded)
pf_cands_py = selected_pt_padded * np.sin(selected_phi_padded)
pf_cands_pz = selected_pt_padded * np.sinh(selected_eta_padded)
pf_cands_E = np.sqrt(pf_cands_px**2 + pf_cands_py**2 + pf_cands_pz**2 + selected_mass_padded**2)
pf_cands_pxpypzE = np.dstack((pf_cands_px, pf_cands_py, pf_cands_pz, pf_cands_E))

# Final cut selection using H4q3q
#gen_parts_pdg_ids = np.array(Gen4qVars["Gen4qpdgId"])
gen_parts_eta_phi_HAA = gen_parts_eta_phi_HAA[H3q4q_cut]
gen_parts_pdg_id_HAA = gen_parts_pdg_ids_HAA[H3q4q_cut]

# Apply cuts to pf_cands and higgs_jet_4vec
pf_cands_pxpypzE = pf_cands_pxpypzE[H3q4q_cut]
higgs_jet_4vec = higgs_jet_4vec[H3q4q_cut]

# Print shapes for debugging
print(cutflow)


print(higgs_jet_4vec[0])
nom_weights = np.ones(len(pf_cands_pxpypzE))
f_ratio_name = '../LundReweighting/data//ratio_2018.root'
f_ratio = ROOT.TFile.Open(f_ratio_name)
LP_rw = LundReweighter(f_ratio = f_ratio)
LP_weights = LP_rw.get_all_weights(pf_cands_pxpypzE, gen_parts_eta_phi_HAA, higgs_jet_4vec, gen_parts_pdg_ids = gen_parts_pdg_id_HAA)
LP_rw.print_unmatched()
for key in LP_weights.keys():
    if('nom' in key or 'up' in key or 'down' in key):
        if(isinstance(LP_weights[key], np.ndarray)) : LP_weights[key] *= nom_weights

print("Bad match frac %.2f" % np.mean(LP_weights['bad_match']))
print("Reclustered bad match frac %.2f" % np.mean(LP_weights['reclust_still_bad_match']))

LPweight_nom  = LP_weights['nom']
LPweight_up   = LP_weights['sys_up']
LPweight_down = LP_weights['sys_down']

print(" LPweight_nom ",LPweight_nom)
print(" LPweight_up ",LPweight_up)
print(" LPweight_down ",LPweight_down)

import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep
import boost_histogram as bh
from cycler import cycler
import numpy as np

# Font style setup
use_helvet = False
if use_helvet:
    CMShelvet = hep.style.CMS
    CMShelvet['font.sans-serif'] = ['Helvetica', 'Arial']
    plt.style.use(CMShelvet)
else:
    plt.style.use(hep.style.CMS)

# Figure setup
plt.figure(figsize=(10, 10))
ax = plt.gca()
plt.grid()
hep.cms.label(data=False, year="2018", ax=ax, fontname='sans-serif')

# Define histogram parameters
nbins, x_min, x_max = 20, 0, 230.0
bins = np.linspace(x_min, x_max, nbins + 1)

# Create histograms
hist_before = bh.Histogram(bh.axis.Regular(nbins, x_min, x_max), storage=bh.storage.Weight())
hist_before.fill(higgs_jet_4vec[:, 3])  # No weights
hist_before_value = hist_before.view().value
hist_before_err = np.sqrt(hist_before.view().variance)

hist_after = bh.Histogram(bh.axis.Regular(nbins, x_min, x_max), storage=bh.storage.Weight())
hist_after.fill(higgs_jet_4vec[:, 3], weight=LPweight_nom)  # Nominal weights
hist_after_value = hist_after.view().value
hist_after_err = np.sqrt(hist_after.view().variance)

# Plot histograms
hep.histplot(hist_before_value, bins=bins, yerr=hist_before_err, label='Before LP reweighting', lw=2, edges=False, histtype="step")
hep.histplot(hist_after_value, bins=bins, yerr=hist_after_err, label='After LP reweighting', lw=2, edges=False, histtype="step")

# Add labels and legend
plt.legend(loc='upper left', frameon=False, fontsize=20)
y_min, y_max = plt.gca().get_ylim()
plt.text(0.08, 0.83 * y_max, "VBF M-A 15 GeV", fontsize=20)

# X and Y axis labels
plt.xlabel(r'Leading FatJet Mass [GeV]', fontsize=20, ha='right', x=1)
plt.ylabel('Events [un-weighted]', fontsize=20, ha='right', y=1)

# Adjust tick size
plt.xticks(size=14)
plt.yticks(size=14)

# Save plot
plt.savefig(f"Jet_Mass_weighted.pdf", bbox_inches='tight')
