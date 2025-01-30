import argparse

import ROOT as R
import ctypes


def printHistogramBinCotent(h, x):
    if not h: return (-1, 0)

    xBin        = h.FindBin(x)
    BinContent  = h.GetBinContent(xBin)
    eBinContent = h.GetBinError(xBin)
    #print(f"  {xBin}: {x}, {BinContent} +- {eBinContent}")
    return (BinContent, eBinContent)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='printEventYields')
    parser.add_argument('-ipDir',    dest='sIpFileDir',           type=str,   default=None,  required=True)
    parser.add_argument('-category', dest='category',           type=str,   default=None,  required=True)
    args=parser.parse_args()
    print("args: {}".format(args))

    sIpFileDir         = args.sIpFileDir
    category         = args.category
    

    #Categories = ['gg0lIncl']
    Categories = [category]
    Processes = ['Data', 'ggHtoaato4b_mA_15', 'ggHtoaato4b_mA_30', 'ggHtoaato4b_mA_55', ]
    Years = ['2018']
    mass = 'mass'
    WPs = ['WP40', 'WP60', 'WP80', ]
    SignalRegions = ['Pass', 'Fail']
    Systematics = 'Nom'

    for cat in Categories:
        for proc in Processes:
            for year in Years:
                sIpFile = '%s/%s_%s_%s.root' % (sIpFileDir, cat,proc,year)
                fIn = R.TFile(sIpFile)
                print("File %s:" % (sIpFile))

                print("WP  : \t Pass, \t\t\t Fail")
                for wp in WPs:
                    sPrint = '%s:\t' % (wp)
                    for reg in SignalRegions:
                        sHisto = '%s_%s_%s_%s_%s_%s_%s' % (cat,proc,year, mass,wp,reg,Systematics)
                        h   = fIn.Get(sHisto)
                        h.StatOverflows()

                        nEventsUnwgt = h.GetEntries()
                        errNEvents = ctypes.c_double(0)
                        nEvents = h.IntegralAndError(1, h.GetNbinsX(), 1, h.GetNbinsY(), errNEvents)

                        sPrint += '%.1f +- %.1f,\t' % (nEvents, errNEvents.value)
                        #sPrint += '%.1f +- %.1f (%d)' % (nEvents, errNEvents.value, nEventsUnwgt)
                        #print(f"{h.GetEntries() = }, {h.GetEffectiveEntries() = }, {h.Integral() = }, {h.GetSumOfWeights() = }, ")
                    
                    print('%s' % sPrint)
                print("")

                fIn.Close()


