import os, sys

sys.path.append( os.path.abspath('../') )
print(f"{os.path.abspath('../') = }")

from htoaa_Settings import *

# Luminosities_perTrigger, Triggers_perEra

#sTriggerSet = 'Trg_Combo_AK4AK8Jet_HT_VBF';  Datasets = ['JetHT', 'BTagCSV']
sTriggerSet = 'Trg_Combo_MET';  Datasets = ['MET']

triggers_set = set()

for era in Triggers_perEra:
    #print(f"{era = }, {Triggers_perEra[era][sTriggerSet] = }")
    #Triggers_perEra[era][sTriggerSet]
    for dataset in Datasets:
        if dataset not in Triggers_perEra[era][sTriggerSet]: continue
        triggers_set.update( list(Triggers_perEra[era][sTriggerSet][dataset].keys()) )

print(f" triggers_set {len(triggers_set)}: {triggers_set} \n\n")

print("Latex table::")
for trg in triggers_set:
    sTrg = '\\verb|%s|' % (trg)
    sLine = '%-75s ' % (sTrg)

    for era in [Era_2016preVFP, Era_2016postVFP, Era_2017, Era_2018]:
        if trg in Luminosities_perTrigger[era]:
            sLine += ' & %4.1f' % (Luminosities_perTrigger[era][trg][0])
        else:
            sLine += ' & %4s' % ('-')
    sLine += ' \\\\'
    print(sLine)