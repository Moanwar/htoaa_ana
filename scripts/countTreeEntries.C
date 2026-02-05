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
 
//void countTreesEntries(const char* sFile) {
//void countTreesEntries(TString sFile) {
//void countTreesEntries(char* sFile) {
void countTreesEntries() {
    const char* sFile = "/eos/cms/store/group/phys_susy/HToaaTo4b/GluGluH_MINNLO_NANOGEN/glugluH_minnlo_2/Run3Summer22_nanogen/250608_105711/0000/wmLHEGS_glugluH_minnlo_*.root";
    std::cout << sFile << ": " << std::endl;
    /*
    //TChain ch("Events");
    //ch.Add("/eos/cms/store/group/phys_susy/HToaaTo4b/GluGluH_MINNLO_NANOGEN/glugluH_minnlo_2/Run3Summer22_nanogen/250608_105711/0000/wmLHEGS_glugluH_minnlo_*.root");
    //TTree *Events = getTTreeFromTChain("/eos/cms/store/group/phys_susy/HToaaTo4b/GluGluH_MINNLO_NANOGEN/glugluH_minnlo_2/Run3Summer22_nanogen/250608_105711/0000/wmLHEGS_glugluH_minnlo_*.root", "Events");
    TTree *Events = getTTreeFromTChain(sFile, "Events");

    
    std::cout << sFile << ": " << Events->GetEntries() << std::endl;
    */

}