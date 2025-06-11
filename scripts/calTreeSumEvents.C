#include <string>

void calTreeSumEvents(string sNtuple, string sTreeName="Events") {
    std::cout << "sNtuple: " << sNtuple.c_str() << ", sTreeName: " << sTreeName.c_str() << std::endl;
    TChain ch(sTreeName.c_str());
    ch.Add(sNtuple.c_str());

    long long nTotalEvents = ch.GetEntries();
    long long nNegEvents = ch.Draw("nGenPart", "genWeight<0", "goff");
    std::cout << "nTotalEvents: " << nTotalEvents << ", nNegEvents: " << nNegEvents << std::endl;
    std::cout << "* (" << nTotalEvents << " - 2*" << nNegEvents << ")/" << nTotalEvents << "  = " << double(nTotalEvents - 2*nNegEvents)/nTotalEvents << std::endl;
}