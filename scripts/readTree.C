
#include <glob.h>
#include <string>
#include <iostream>

std::vector<std::string> glob(const char *pattern) {
    glob_t g;
    glob(pattern, GLOB_TILDE, nullptr, &g); // one should ensure glob returns 0!
    std::vector<std::string> filelist;
    filelist.reserve(g.gl_pathc);
    for (size_t i = 0; i < g.gl_pathc; ++i) {
        filelist.emplace_back(g.gl_pathv[i]);
    }
    globfree(&g);
    return filelist;
}

TTree * getTTreeFromTFile(const char* sFile, const char* sTree) {
    TFile *f = new TFile(sFile);
    TTree *tr = (TTree*)f->Get(sTree);
    return tr;
}

TTree * getTTreeFromTChain(const char* sFile, const char* sTree) {
    TChain *tr = new TChain(sTree);
    for (const auto &filename : glob(sFile)) {
       tr->Add(filename.c_str());
    }
    return tr;
}


void readTrees() {
    //TChain ch("Events");
    //ch.Add("/eos/cms/store/group/phys_susy/HToaaTo4b/GluGluH_MINNLO_NANOGEN/glugluH_minnlo_2/Run3Summer22_nanogen/250608_105711/0000/wmLHEGS_glugluH_minnlo_*.root");
    TTree *Events = getTTreeFromTChain("/eos/cms/store/group/phys_susy/HToaaTo4b/GluGluH_MINNLO_NANOGEN/glugluH_minnlo_2/Run3Summer22_nanogen/250608_105711/0000/wmLHEGS_glugluH_minnlo_*.root", "Events");


    TFile *fOut               = new TFile("/eos/cms/store/user/ssawant/htoaa/analysis/HiggsPtRewgts/2018/GenHiggsPt_GGH_NNLO.root", "recreate");

    TH1F *hGenHiggsPt_Nom     = new TH1F("hGenHiggsPt_Nom", "", 2000,            0,               2000);
    TH1F *hGenHiggsLog2Pt_Nom = new TH1F("hGenHiggsLog2Pt_Nom", "", 200,  TMath::Log2(1),    TMath::Log2(2000));

    Events->Draw("GenPart_pt>>hGenHiggsPt_Nom",                  "(abs(GenPart_pdgId)==25) && (GenPart_status==62) * (1 - 2*(genWeight < 0))");
    Events->Draw("TMath::Log2(GenPart_pt)>>hGenHiggsLog2Pt_Nom", "(abs(GenPart_pdgId)==25) && (GenPart_status==62) * (1 - 2*(genWeight < 0))");
    
    fOut->cd();
    hGenHiggsPt_Nom->Write();
    hGenHiggsLog2Pt_Nom->Write();
    fOut->Close();

    std::cout << "readTrees done" << std::endl;
}
