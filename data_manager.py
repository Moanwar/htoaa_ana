import uproot
import pandas as pd
from analib import PhysObj, Event
from info import allVars, cutVars, cutDict, weightDict

def processData (fileName): 
    ## open file, get events
    f = uproot.open(fileName + '.root')
    events = f.get('Events')
    ## make PhysObj of the event
    data = PhysObj('data_' + fileName)
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
            
    ## sync Events
    ev.sync()


    if not (data.FatJet_pt.empty):
        ## get the max pt index for each event
        ## then just loop through the PhysicsObjs and extract the 
        ## ones that matches that index. Not 
        colidx = data['FatJet_pt'].idxmax(axis = 1).to_numpy()
        #print(colidx)
        rowidx = list(range(len(colidx)))
        maxPtData = pd.DataFrame()
        for var in allVars: 
            npArr = data[var].to_numpy()
            maxPtData[var] = npArr[rowidx, colidx]        
            
        maxPtData['weights'] = weightDict[fileName]
        if fileName == 'GGH_HPT':
            maxPtData['target'] = 1
        else: 
            maxPtData['target'] = 0
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
    

    
    
    
    
    
    
    
