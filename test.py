import uproot
import awkward as ak

arr = [[10], [12, 3], [33, 22, 2], [2]]

print("arr: {}".format(arr))
print("arr: {}".format(ak.flatten(arr)))
