

from collections import OrderedDict as OD

import ROOT as R



if __name__ == "__main__":

    '''
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

    '''


    sIpFile = "/home/siddhesh/Work/CMS/htoaa/analysis/20230519_StitchOverlapPhSpRemoval/2018/analyze_htoaa_stage1.root"
    sOpFile = "/home/siddhesh/Work/CMS/htoaa/htoaa/data/correction/mc/HTSamplesStitch/LHE_HT_2018.root"

    sHistos_list = OD([
        # (<histogram name in input file>,   <histogram name in output file. If empty, use input file histogram name>)
        ("evt/QCD_bEnrich/hCutFlow_central",                     ""),
        ("evt/QCD_bEnrich/hCutFlowWeighted_central",             ""),
        ("evt/QCD_bEnrich/hGenLHE_HT_all_central",               ""),
        ("evt/QCD_bEnrich/hGenLHE_HT_SelQCDbEnrich_central",     ""),
        
        ("evt/QCD_bGen/hCutFlow_central",                    ""),
        ("evt/QCD_bGen/hCutFlowWeighted_central",            ""),
        ("evt/QCD_bGen/hGenLHE_HT_all_central",              ""),
        ("evt/QCD_bGen/hGenLHE_HT_SelQCDbGen_central",       ""),
        
        ("evt/QCD_Incl/hCutFlow_central",                     ""),
        ("evt/QCD_Incl/hCutFlowWeighted_central",             ""),
        ("evt/QCD_Incl/hGenLHE_HT_all_central",               ""),
        ("evt/QCD_Incl/hGenLHE_HT_SelQCDbEnrich_central",     ""),
        ("evt/QCD_Incl/hGenLHE_HT_SelQCDbGen_central",        ""),
        
        ("evt/QCD_Incl_PSWeight/hCutFlow_central",                     ""),
        ("evt/QCD_Incl_PSWeight/hCutFlowWeighted_central",             ""),
        ("evt/QCD_Incl_PSWeight/hGenLHE_HT_all_central",               ""),
        ("evt/QCD_Incl_PSWeight/hGenLHE_HT_SelQCDbEnrich_central",     ""),
        ("evt/QCD_Incl_PSWeight/hGenLHE_HT_SelQCDbGen_central",        ""),
        
        ("evt/TTJets_NLO/hCutFlow_central",               ""),
        ("evt/TTJets_NLO/hCutFlowWeighted_central",       ""),
        ("evt/TTJets_NLO/hGenLHE_HT_all_central",         ""),
        
        ("evt/TTJets_Incl/hCutFlow_central",               ""),
        ("evt/TTJets_Incl/hCutFlowWeighted_central",       ""),
        ("evt/TTJets_Incl/hGenLHE_HT_all_central",         ""),
        
        ("evt/TTJets_HT/hCutFlow_central",               ""),
        ("evt/TTJets_HT/hCutFlowWeighted_central",       ""),
        ("evt/TTJets_HT/hGenLHE_HT_all_central",         ""),
        
        ("evt/ZJetsToQQ_HT/hCutFlow_central",               ""),
        ("evt/ZJetsToQQ_HT/hCutFlowWeighted_central",       ""),
        ("evt/ZJetsToQQ_HT/hGenLHE_HT_all_central",         ""),
        
        ("evt/WJetsToQQ_HT/hCutFlow_central",               ""),
        ("evt/WJetsToQQ_HT/hCutFlowWeighted_central",       ""),
        ("evt/WJetsToQQ_HT/hGenLHE_HT_all_central",         ""),
        
        ("evt/WJetsToLNu_Incl/hCutFlow_central",               ""),
        ("evt/WJetsToLNu_Incl/hCutFlowWeighted_central",       ""),
        ("evt/WJetsToLNu_Incl/hGenLHE_HT_all_central",         ""),
        
        ("evt/WJetsToLNu_HT/hCutFlow_central",               ""),
        ("evt/WJetsToLNu_HT/hCutFlowWeighted_central",       ""),
        ("evt/WJetsToLNu_HT/hGenLHE_HT_all_central",         ""),
        
        ("evt/W1JetsToLNu/hCutFlow_central",               ""),
        ("evt/W1JetsToLNu/hCutFlowWeighted_central",       ""),
        ("evt/W1JetsToLNu/hGenLHE_HT_all_central",         ""),
        
        ("evt/W2JetsToLNu/hCutFlow_central",               ""),
        ("evt/W2JetsToLNu/hCutFlowWeighted_central",       ""),
        ("evt/W2JetsToLNu/hGenLHE_HT_all_central",         ""),
        
        ("evt/W3JetsToLNu/hCutFlow_central",               ""),
        ("evt/W3JetsToLNu/hCutFlowWeighted_central",       ""),
        ("evt/W3JetsToLNu/hGenLHE_HT_all_central",         ""),
        
        ("evt/W4JetsToLNu/hCutFlow_central",               ""),
        ("evt/W4JetsToLNu/hCutFlowWeighted_central",       ""),
        ("evt/W4JetsToLNu/hGenLHE_HT_all_central",         ""),
        
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
