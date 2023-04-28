import argparse

import ROOT as R


def printHistogramBinCotent(h, x):
    if not h: return (-1, 0)

    xBin        = h.FindBin(x)
    BinContent  = h.GetBinContent(xBin)
    eBinContent = h.GetBinError(xBin)
    #print(f"  {xBin}: {x}, {BinContent} +- {eBinContent}")
    return (BinContent, eBinContent)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='printHisto1DBinContent')
    parser.add_argument('-ipFile',    dest='sIpFile',           type=str,   default=None,  required=True)
    parser.add_argument('-histogram', dest='sHistoName',        type=str,   default=None,  required=True)
    parser.add_argument('-x',         dest='x_forBinContent',   type=float, default=None,  required=True)
    args=parser.parse_args()
    print("args: {}".format(args))

    sIpFile         = args.sIpFile
    sHistoName      = args.sHistoName
    x_forBinContent = args.x_forBinContent

    
    
    sHistogramsList = [
        "evt/Data/%s"                       % (sHistoName.replace('central', 'noweight')),
        "evt/QCD_bEnrich/%s"                % (sHistoName),
        "evt/QCD_bGen/%s"                   % (sHistoName),
        "evt/QCDIncl/%s"                    % (sHistoName),
        "evt/QCDIncl_PSWeight/%s"           % (sHistoName),
        "evt/TTJets/%s"                     % (sHistoName),
        "evt/ZJets/%s"                      % (sHistoName),
        "evt/WJets/%s"                      % (sHistoName),
        "evt/SUSY_GluGluH_01J_HToAATo4B/%s" % (sHistoName),
    ]

    fIn = R.TFile(sIpFile)
    print(f"Input file: {sIpFile}")
    
    for sHisto in sHistogramsList:
        #print(f"{sHisto = }")
        h   = fIn.Get(sHisto)
        BinContent, eBinContent = printHistogramBinCotent(h, x_forBinContent)
        print("\t %-50s: \t x: %6g, \t n: %10d  +-  %10d " % (sHisto, x_forBinContent, BinContent, eBinContent))
