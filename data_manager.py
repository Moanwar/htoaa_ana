import uprootimport numpy as npimport pandas as pdfrom analib import PhysObj, Eventfrom info import trainVars, allVars, cutVars, cutDict, weightDictdef processData (fileName):     ## open file, get events    f = uproot.open(fileName + '.root')    events = f.get('Events')    ## make PhysObj of the event    data = PhysObj('data_' + fileName)    for var in allVars:         data[var] = pd.DataFrame(events.array(var))        if 'eta' in var:             data[var] = data[var].abs()            ## make event object    ev = Event(data)        ## apply cuts    for cutVar in cutVars:         data.cut(data[cutVar] > cutDict[cutVar])        ## sync Events    ev.sync()        ## rename columns     jetNums = list(range(1, data.FatJet_pt.shape[1]+1)) # for naming the columns    wideData = pd.DataFrame()    colNames = []        for var in allVars:         colValues = [var + "_" + str(i) for i in jetNums]        colNames = colNames + colValues        colDict = dict(list(enumerate(colValues)))        data[var] = data[var].rename(columns = colDict)                if var == allVars[0]:            wideData = wideData.append(data[var])        else:            wideData = wideData.join(data[var], sort=False)           ## add info about whether it is a signal or bg, and add weight    ## idk if i need to know what process each one is but here we go    wideData['process'] = fileName        ## right now I am just cutting the data out. might try to use weights later    wideData['weights'] = weightDict[fileName]    #wideData = wideData.sample(frac = weightDict[fileName])    if fileName == 'GGH_HPT':        wideData['target'] = 1    else:         wideData['target'] = 0    wideData = wideData.dropna(how = 'all')    return wideData, colNames                                