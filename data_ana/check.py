#! /usr/bin/env python
import ROOT as R

# Replace with one of your actual file paths
test_file = "/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/PNet_v2_2024_11_22/JetHT/r1_Run2018A/PNet_v1_Skim_0_0.root"  # or the full path to your file

f = R.TFile.Open(test_file)
tree = f.Get("Events")

# Print all trigger-related branches
print("Trigger-related branches in your file:")
branches = [branch.GetName() for branch in tree.GetListOfBranches()]
trigger_branches = [b for b in branches if 'HLT' in b or 'hlt' in b or 'L1' in b]

for i, branch in enumerate(sorted(trigger_branches)[:50]):
    print(f"  {branch}")

# Also check if any of our expected triggers exist
expected_triggers = [
    'HLT_AK8PFJet500', 'HLT_PFHT1050', 'HLT_PFJet500',
    'HLT_AK8PFHT800_TrimMass50', 'HLT_AK8PFJet400_TrimMass30',
    'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4'
]

print("\nChecking for expected triggers:")
for trig in expected_triggers:
    branch_name = trig.replace('HLT_', '')
    if branch_name in trigger_branches:
        print(f"  ✓ Found: {trig} (branch: {branch_name})")
    else:
        print(f"  ✗ Not found: {trig}")

f.Close()
