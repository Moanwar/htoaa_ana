'''
Convert json file to csv format
Command: python3 scripts/convert_json_to_csv.py <json file>.json

'''

import sys
#import json
#import csv
import pandas as pd



#print("sys.argv: {}".format(sys.argv))
if len(sys.argv) != 2:
    print("ERROR: Give json file in argument. \n")
    exit(0)

sInFile = sys.argv[1]
sOutFile = sInFile.replace('.json', '.csv')
print("Convert: %s --> %s" % (sInFile, sOutFile))


data = pd.read_json(sInFile, orient='index')
data.to_csv(sOutFile)

print("data: {}".format(data))

