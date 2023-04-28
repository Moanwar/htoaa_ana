

from collections import OrderedDict as OD

import ROOT as R



if __name__ == "__main__":

    sIpFile = "/home/siddhesh/Work/CMS/htoaa/analysis/20230419_forHTSamplesStitch/2018/analyze_htoaa_stage1.root"
    sOpFile = "/home/siddhesh/Work/CMS/htoaa/htoaa/data/correction/mc/HTSamplesStitch/LHE_HT_2018.root"

    sHistos_list = OD([
        # (<histogram name in input file>,   <histogram name in output file. If empty, use input file histogram name>)
        ("evt/QCD_bEnrich/hCutFlow_central",           ""),
        ("evt/QCD_bEnrich/hCutFlowWeighted_central",   ""),
        ("evt/QCD_bEnrich/hGenLHE_HT_all_central",     ""),
        
        ("evt/QCD_bGen/hCutFlow_central",              ""),
        ("evt/QCD_bGen/hCutFlowWeighted_central",      ""),
        ("evt/QCD_bGen/hGenLHE_HT_all_central",        ""),
        
        ("evt/QCDIncl/hCutFlow_central",               ""),
        ("evt/QCDIncl/hCutFlowWeighted_central",       ""),
        ("evt/QCDIncl/hGenLHE_HT_all_central",         ""),
        
        ("evt/QCDIncl_PSWeight/hCutFlow_central",               ""),
        ("evt/QCDIncl_PSWeight/hCutFlowWeighted_central",       ""),
        ("evt/QCDIncl_PSWeight/hGenLHE_HT_all_central",         ""),
        
        ("evt/ZJets/hCutFlow_central",               ""),
        ("evt/ZJets/hCutFlowWeighted_central",       ""),
        ("evt/ZJets/hGenLHE_HT_all_central",         ""),
        
        ("evt/WJets/hCutFlow_central",               ""),
        ("evt/WJets/hCutFlowWeighted_central",       ""),
        ("evt/WJets/hGenLHE_HT_all_central",         ""),
        
    ])


    fIpFile = R.TFile(sIpFile)
    print(f"Input file: {sIpFile}")
    if not fIpFile.IsOpen():
        print(f"Could not open input file {sIpFile}")
        exit(0)

    fOpFile = R.TFile(sOpFile, "RECREATE")
    print(f"Output file: {sOpFile}")
    
    for sIpHisto, sOpHisto in sHistos_list.items():
        if not sOpHisto:
             sOpHisto = sIpHisto
        print(f"sIpHisto: {sIpHisto}, \t sOpHisto: {sOpHisto}.")

        h = fIpFile.Get(sIpHisto)

        fOpFile.cd()
        if '/' in sOpHisto:
            sOpHisto_parts = sOpHisto.split('/')
            dirName = '/'.join(sOpHisto_parts[:-1])
            #print(f"sOpHisto_parts: {sOpHisto_parts} ")
            #print(f"{dirName = }")

            #print(f"{fOpFile.GetList() = }")
            #print(f"{fOpFile.GetListOfKeys() = }")
            #print(f"{fOpFile.FindObject(dirName) = }")
            
            dir1 = fOpFile.mkdir( '/'.join(sOpHisto_parts[:-1]) )
            
            #dir1.cd()
            fOpFile.cd(dirName)
            h.Write()

    
    fOpFile.Close()
    
    print(f"\nCopied histograms from \n{sIpFile} \nto \n{sOpFile}")
