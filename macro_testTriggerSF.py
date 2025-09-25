import json
import matplotlib.pyplot as plt
import numpy as np

# Load the JSON
with open("data/correction/mc/TrgEffSF/Hadronic/GGF_triggerSF_soup_1D_Inc_pT_2018.json") as f:
    data = json.load(f)

# Navigate into correction content
correction = data["corrections"][0]  # only one in your file
contents = {entry["key"]: entry["value"] for entry in correction["data"]["content"]}

# Extract pT bin edges
pt_edges = contents["nominal"]["edges"][0]  # one list inside
pt_centers = 0.5 * (np.array(pt_edges[:-1]) + np.array(pt_edges[1:]))

# Extract values
nominal = np.array(contents["nominal"]["content"])
stat_up = np.array(contents["stat_up"]["content"])
stat_dn = np.array(contents["stat_dn"]["content"])

# Plot
plt.figure(figsize=(8,6))
plt.step(pt_centers, nominal, where="mid", label="Nominal", color="black")
plt.step(pt_centers, stat_up, where="mid", label="Stat Up", color="red", linestyle="--")
plt.step(pt_centers, stat_dn, where="mid", label="Stat Down", color="blue", linestyle="--")

plt.xlabel("FatJet pT [GeV]")
plt.ylabel("Trigger SF")
plt.title("Trigger scale factor (2018)")
plt.legend()
plt.grid(True, linestyle=":")
plt.tight_layout()
plt.show()
