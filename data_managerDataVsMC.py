import uproot
import pandas as pd
from analib import PhysObj, Event
#from info import allVars, cutVars, cutDict, weightDict
import sys
import os
import pickle 

dataPath = '/home/chosila/Projects/htoaa/data/2018D_Parked_promptD-v1_200218_214714_Skim_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root'
ggHPath = '/home/chosila/Projects/htoaa/MC/nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId_999k.root'
#ggHPath = 'GGH_HPT.root'
BGenPath = '/home/chosila/Projects/htoaa/MC/QCD_BGen_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root'
bEnrPath = '/home/chosila/Projects/htoaa/MC/QCD_bEnriched_nFat1_doubB_0p8_deepB_Med_massH_90_200_msoft_90_200_pT_240_Mu_pT_6_IP_2_softId.root'


BGenWeight = [1, 0.259, 0.0515, 0.01666, 0.00905, 0.003594, 0.001401]
bEnrWeight =[ 1, 0.33, 0.034, 0.034, 0.0024, 0.00024, 0.00044]

trainVars = [
    'FatJet_pt', 
    'FatJet_eta', 
    'FatJet_mass', 
    'FatJet_btagCSVV2', 
    'FatJet_btagDeepB', 
    'FatJet_msoftdrop', 
    'FatJet_btagDDBvL',
    'FatJet_deepTagMD_H4qvsQCD'
    ]
cutVars = [
    'FatJet_btagDDBvL', 
    'FatJet_btagDeepB', 
    'FatJet_mass', 
    'FatJet_msoftdrop', 
    'FatJet_pt',
    'Muon_pt'
    ]
plotVars = [
    'PV_npvs',
    'PV_npvsGood',
    'Muon_pt',
    'Muon_eta',
    'Muon_ip3d'
    ]
extraVars = ['LHE_HT',
             'Muon_dxy',
             'Muon_dxyErr']

allVars = list(set(cutVars + trainVars + extraVars + plotVars))

cutValues = [0.8, 0.4184, 90, 90, 240, 6]

cutDict = dict(zip(cutVars, cutValues))

muonR = pickle.load(open('muontensor/MuonRtensor.p', 'rb'))
muonL = pickle.load(open('muontensor/MuonLtensor.p', 'rb'))

ptkeys = list(muonL.keys())
ipkeys = list(muonL[ptkeys[1]].keys())


def processData (filePath, tag):
    ## open file, get events
    fileName, fileExtension = os.path.splitext(filePath)
    
    print(fileName)

    if fileExtension != '.root':
        print('this program only supports .root  files')
        sys.exit()

    f = uproot.open(fileName + '.root')
    events = f.get('Events')
    ## make PhysObj of the event
    data = PhysObj('data_' + fileName)
    

    if tag == 'data':
        allVars.remove('LHE_HT')
    for var in allVars: 
        data[var] = pd.DataFrame(events.array(var))
        ## makes eta positive only
        if 'eta' in var: 
            data[var] = data[var].abs()
    
    ## make event object
    ev = Event(data)
    
    ## apply cuts
    for cutVar in cutVars:
        data.cut(data[cutVar] > cutDict[cutVar])
        if (cutVar == 'FatJet_mass' or cutVar == 'FatJet_msoftdrop'):
            data.cut(data[cutVar] < 200)
    data.cut(data['Muon_eta'] < 2.4)
    data.cut(data['Muon_ip3d'] < 0.5)
    data.cut((data['Muon_dxy']/data['Muon_dxyErr']).abs() > 2)

    ## sync Events
    ev.sync()

    if not (data.FatJet_pt.empty):
        ## get the max pt index for each event
        ## then just loop through the PhysicsObjs and extract the 
        ## ones that matches that index. Not 

        colidx = data['FatJet_pt'].idxmax(axis = 1).to_numpy()
        rowidx = list(range(len(colidx)))
        maxPtData = pd.DataFrame()

        for var in allVars:
            npArr = data[var].to_numpy()
            maxPtData[var] = npArr[rowidx, colidx]

        
        ## LHE weights
        #if tag == 'data':
        #    maxPtData['LHE_weights'] = 1
        if tag == 'ggH':
            maxPtData['LHE_weights'] = 43920/999000
        elif tag == 'BGen':
            maxPtData.loc[(maxPtData['LHE_HT']>200) & (maxPtData['LHE_HT']<300),
                          'LHE_weights'] = BGenWeight[0]
            maxPtData.loc[(maxPtData['LHE_HT']>300) & (maxPtData['LHE_HT']<500),
                          'LHE_weights'] = BGenWeight[1]
            maxPtData.loc[(maxPtData['LHE_HT']>500) & (maxPtData['LHE_HT']<700),
                          'LHE_weights'] = BGenWeight[2]
            maxPtData.loc[(maxPtData['LHE_HT']>700) & (maxPtData['LHE_HT']<1000),
                          'LHE_weights'] = BGenWeight[3]
            maxPtData.loc[(maxPtData['LHE_HT']>1000) & (maxPtData['LHE_HT']<1500),
                          'LHE_weights'] = BGenWeight[4]
            maxPtData.loc[(maxPtData['LHE_HT']>1500) & (maxPtData['LHE_HT']<2000),
                          'LHE_weights'] = BGenWeight[5]
            maxPtData.loc[maxPtData['LHE_HT']>2000,
                          'LHE_weights'] = BGenWeight[6]
            print(BGenWeight[6])

        elif tag == 'bEnr':
            maxPtData.loc[(maxPtData['LHE_HT']>200) & (maxPtData['LHE_HT']<300),
                          'LHE_weights'] = bEnrWeight[0]
            maxPtData.loc[(maxPtData['LHE_HT']>300) & (maxPtData['LHE_HT']<500),
                          'LHE_weights'] = bEnrWeight[1]
            maxPtData.loc[(maxPtData['LHE_HT']>500) & (maxPtData['LHE_HT']<700),
                          'LHE_weights'] = bEnrWeight[2]
            maxPtData.loc[(maxPtData['LHE_HT']>700) & (maxPtData['LHE_HT']<1000),
                          'LHE_weights'] = bEnrWeight[3]
            maxPtData.loc[(maxPtData['LHE_HT']>1000) & (maxPtData['LHE_HT']<1500),
                          'LHE_weights'] = bEnrWeight[4]
            maxPtData.loc[(maxPtData['LHE_HT']>1500) & (maxPtData['LHE_HT']<2000),
                          'LHE_weights'] = bEnrWeight[5]
            maxPtData.loc[maxPtData['LHE_HT']>2000,
                          'LHE_weights'] = bEnrWeight[6]

        ## Muon weights 
        if tag != 'data': 
            for ptIdx, pt in enumerate(ptkeys):
                if pt == 'meta':
                    continue
                for ipIdx, ip in enumerate(ipkeys): 
                    maxPtData.loc[(maxPtData.Muon_pt > pt) & (maxPtData < ptkeys[ptIdx+1]) &
                                  (maxPtData.Muon_ip3d > ip) & (maxPtData.Muon_ip3d < ipkeys[ipIdx +1]) &
                                  (maxPtData.Muon_eta < 1.5), 'muon_weights'] = muonR[pt][ip]['L']
                    maxPtData.loc[(maxPtData.Muon_pt > pt) & (maxPtData < ptkeys[ptIdx+1]) &
                                  (maxPtData.Muon_ip3d > ip) & (maxPtData.Muon_ip3d < ipkeys[ipIdx +1]) &
                                  (maxPtData.Muon_eta > 1.5), 'muon_weights'] = muonR[pt][ip]['H']


        #### !!!!! TODO !!!! #####
        ## lumi weights


    else:
        maxPtData = pd.DataFrame()

    return maxPtData
#     except:
#         print('in except')
#         return maxPtData
#         
        
    
        
        ## what if instead, cut it using max???
        ## that didn't work lol
        # data.cut(data.FatJet_pt.max())
        # print(data.FatJet_pt)
        
        
        
        ## rename columns 
        # jetNums = list(range(1, 9)) # for naming the columns
        # wideData = pd.DataFrame()
        # colNames = []
        
        # for var in allVars: 
        #     colValues = [var + "_" + str(i) for i in jetNums]
        #     colNames = colNames + colValues
        #     colDict = dict(list(enumerate(colValues)))
        #     data[var] = data[var].rename(columns = colDict)
            
            ## slicing so only 4 jets.
            ## heck some of them only have 3
            # if (len(data[var].columns) > 3):
            #     dfToAppend = data[var].iloc[:, [0,1,2,3]]
            # else:
            #     dfToAppend = data[var].iloc[:, [0,1,2]]
            
            ## slice down to only 3 jets
            # dfToAppend = data[var].iloc[:, [0, 1, 2]]
            
            # if var == allVars[0]: 
            #     wideData = wideData.append(dfToAppend)
            # else: 
            #     wideData = wideData.join(dfToAppend, sort = False)
            # if var == allVars[0]:
            #     wideData = wideData.append(data[var])
            # else:
            #     wideData = wideData.join(data[var], sort=False)
           
        ## add info about whether it is a signal or bg, and add weight
        ## idk if i need to know what process each one is but here we go
        ##wideData['process'] = fileName
        
        ## right now I am just cutting the data out. might try to use weights later
        # wideData['weights'] = weightDict[fileName]
        #wideData = wideData.sample(frac = weightDict[fileName])
        # if fileName == 'GGH_HPT':
        #     wideData['target'] = 1
        # else: 
        #     wideData['target'] = 0
    
        # return wideData

        
        
# def main():
#     processData('GGH_HPT')

# if __name__ == "__main__":
#     main()
    

#def insertWeight(maxPtData, minHT, maxHT):


    
    
    
    
    
