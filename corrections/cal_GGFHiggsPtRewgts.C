/*
Macro to derive and test GEN-Higgs pT reweighting for GGHToAATo4B madgraph sample w.r.t. GGHTo2B MINLO sample

To run:
    root -l
    .L cal_GGFHiggsPtRewgts.C
    cal_GGFHiggsPtRewgts_run();
    cal_GGFHiggsPtRewgts_fit(); 
*/






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


void cal_GGFHiggsPtRewgts() {
    TFile *fGGHToBBIncl = new TFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8/hadd/GluGluHToBB_M-125_NNLOPS_NanoAODv9.root");
    TFile *fGGHToBBPt   = new TFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/GluGluHToBB_Pt-200ToInf_M-125_TuneCP5_MINLO_13TeV-powheg-pythia8/hadd/GluGluHToBB_Pt-200ToInf_M-125_NanoAODv9.root");
    TFile *fGGHToAATo4BPt = new TFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-All_TuneCP5_13TeV_madgraph_pythia8/hadd/GluGluH_01J_HToAATo4B_Pt150_M-All_NanoAODv9.root");

    //double xsGGHToBBIncl = 48.61 * 0.582;
    //double xsGGHToBBPt   = 0.2740; 
    //double xsGGHToAATo4BPt = 48.61 * 0.057;
    // scale GGHToBB and GGHToAATo4B to 48.61 pb cross section
    double xsGGHToBBIncl = 48.61 ;
    double xsGGHToBBPt   = 0.2740 / 0.572 / 0.582; // 0.2740; 
    double xsGGHToAATo4BPt =  48.61 * 0.057 / 1.38; 

    TTree *trGGHToBBIncl = (TTree*)fGGHToBBIncl->Get("Events");
    long nEventsGenWgtPos_GGHToBBIncl = trGGHToBBIncl->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHToBBIncl = trGGHToBBIncl->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHToBBIncl       = nEventsGenWgtPos_GGHToBBIncl - nEventsGenWgtNeg_GGHToBBIncl;
    std::cout << "GGHToBBIncl nEvents: " << trGGHToBBIncl->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHToBBIncl: " << nEventsGenWgtPos_GGHToBBIncl;
    std::cout << ", nEventsGenWgtNeg_GGHToBBIncl: " << nEventsGenWgtNeg_GGHToBBIncl;
    std::cout << ", nEventsTot_GGHToBBIncl: " << nEventsTot_GGHToBBIncl << "\n";

    TTree *trGGHToBBPt = (TTree*)fGGHToBBPt->Get("Events");
    long nEventsGenWgtPos_GGHToBBPt = trGGHToBBPt->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHToBBPt = trGGHToBBPt->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHToBBPt       = nEventsGenWgtPos_GGHToBBPt - nEventsGenWgtNeg_GGHToBBPt;
    std::cout << "GGHToBBPt nEvents: " << trGGHToBBPt->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHToBBPt: " << nEventsGenWgtPos_GGHToBBPt;
    std::cout << ", nEventsGenWgtNeg_GGHToBBPt: " << nEventsGenWgtNeg_GGHToBBPt;
    std::cout << ", nEventsTot_GGHToBBPt: " << nEventsTot_GGHToBBPt << "\n";

    TTree *trGGHToAATo4BPt = (TTree*)fGGHToAATo4BPt->Get("Events");
    long nEventsGenWgtPos_GGHToAATo4BPt = trGGHToAATo4BPt->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHToAATo4BPt = trGGHToAATo4BPt->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHToAATo4BPt       = nEventsGenWgtPos_GGHToAATo4BPt - nEventsGenWgtNeg_GGHToAATo4BPt;
    std::cout << "GGHToAATo4BPt nEvents: " << trGGHToAATo4BPt->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHToAATo4BPt: " << nEventsGenWgtPos_GGHToAATo4BPt;
    std::cout << ", nEventsGenWgtNeg_GGHToAATo4BPt: " << nEventsGenWgtNeg_GGHToAATo4BPt;
    std::cout << ", nEventsTot_GGHToAATo4BPt: " << nEventsTot_GGHToAATo4BPt << "\n";

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas *c1 = new TCanvas("c1", "c1", 600,500);
    c1->SetLogy();
    c1->SetGrid();
    c1->cd();
      
    TString sCut_GGHToBBIncl = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld)", xsGGHToBBIncl, nEventsTot_GGHToBBIncl);
    std::cout << "sCut_GGHToBBIncl: " <<  sCut_GGHToBBIncl.Data() << "\n";
    TH1D *hGenHiggsPt_GGHToBBIncl = new TH1D("hGenHiggsPt_GGHToBBIncl", "hGenHiggsPt_GGHToBBIncl", 150, 0, 1500); 
    hGenHiggsPt_GGHToBBIncl->SetMarkerStyle(20);
    hGenHiggsPt_GGHToBBIncl->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToBBIncl->SetLineColor(kBlack);
    hGenHiggsPt_GGHToBBIncl->SetMarkerColor(kBlack);
    hGenHiggsPt_GGHToBBIncl->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHToBBIncl->GetYaxis()->SetTitle("No. of events");
    trGGHToBBIncl->Draw("GenPart_pt >> hGenHiggsPt_GGHToBBIncl", sCut_GGHToBBIncl.Data());

    
    TString sCut_GGHToBBPt = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld)", xsGGHToBBPt, nEventsTot_GGHToBBPt);
    std::cout << "sCut_GGHToBBPt: " <<  sCut_GGHToBBPt.Data() << "\n";
    TH1D *hGenHiggsPt_GGHToBBPt   = new TH1D("hGenHiggsPt_GGHToBBPt", "hGenHiggsPt_GGHToBBPt", 150, 0, 1500);
    hGenHiggsPt_GGHToBBPt->SetMarkerStyle(20);
    hGenHiggsPt_GGHToBBPt->SetMarkerSize(0.5);   
    hGenHiggsPt_GGHToBBPt->SetLineColor(kBlue);
    hGenHiggsPt_GGHToBBPt->SetMarkerColor(kBlue);
    hGenHiggsPt_GGHToBBPt->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    trGGHToBBPt->Draw("GenPart_pt >> hGenHiggsPt_GGHToBBPt", sCut_GGHToBBPt.Data(), "same");

 
    TString sCut_GGHToAATo4BPt_0 = Form("( (GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld) )", xsGGHToAATo4BPt, nEventsTot_GGHToAATo4BPt);
    std::cout << "sCut_GGHToAATo4BPt_0: " <<  sCut_GGHToAATo4BPt_0.Data() << "\n";
    TH1D *hGenHiggsPt_GGHToAATo4BPt_0   = new TH1D("hGenHiggsPt_GGHToAATo4BPt_0", "hGenHiggsPt_GGHToAATo4BPt_0", 150, 0, 1500);
    hGenHiggsPt_GGHToAATo4BPt_0->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BPt_0->SetMarkerSize(0.5);     
    hGenHiggsPt_GGHToAATo4BPt_0->SetLineColor(kCyan);
    hGenHiggsPt_GGHToAATo4BPt_0->SetMarkerColor(kCyan);
    trGGHToAATo4BPt->Draw("GenPart_pt >> hGenHiggsPt_GGHToAATo4BPt_0", sCut_GGHToAATo4BPt_0.Data(), "same");


    TString sCut_GGHToAATo4BPt = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld) * (min(max(3.9 - (0.4 * log2(GenPart_pt)), 0.1), 10))", xsGGHToAATo4BPt, nEventsTot_GGHToAATo4BPt);
    std::cout << "sCut_GGHToAATo4BPt: " <<  sCut_GGHToAATo4BPt.Data() << "\n";
    TH1D *hGenHiggsPt_GGHToAATo4BPt   = new TH1D("hGenHiggsPt_GGHToAATo4BPt", "hGenHiggsPt_GGHToAATo4BPt", 150, 0, 1500);
    hGenHiggsPt_GGHToAATo4BPt->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BPt->SetMarkerSize(0.5);     
    hGenHiggsPt_GGHToAATo4BPt->SetLineColor(kRed);
    hGenHiggsPt_GGHToAATo4BPt->SetMarkerColor(kRed);
    trGGHToAATo4BPt->Draw("GenPart_pt >> hGenHiggsPt_GGHToAATo4BPt", sCut_GGHToAATo4BPt.Data(), "same");

    TLegend *leg1 = new TLegend(0.5,0.75,0.99,0.99);
    leg1->AddEntry(hGenHiggsPt_GGHToBBIncl, "GGToHTo2B inclusive", "lep");
    leg1->AddEntry(hGenHiggsPt_GGHToBBPt, "GGToHTo2B pT(H)>200", "lep");
    leg1->AddEntry(hGenHiggsPt_GGHToAATo4BPt_0, "GGToHToAATo4B pT(H)>150", "lep");
    leg1->AddEntry(hGenHiggsPt_GGHToAATo4BPt, "GGToHToAATo4B pT(H)>150 reweighted", "lep");
    leg1->Draw();


    TCanvas *c2 = new TCanvas("c2", "c2", 600,500);
    c2->SetLogy();
    c2->SetGrid();
    c2->cd();
    TH1D *hGenHiggsPt_GGHToBBIncl_cl   = (TH1D*)hGenHiggsPt_GGHToBBIncl->Clone(Form("%s_clone", hGenHiggsPt_GGHToBBIncl->GetName()));
    TH1D *hGenHiggsPt_GGHToBBPt_cl     = (TH1D*)hGenHiggsPt_GGHToBBPt->Clone(Form("%s_clone", hGenHiggsPt_GGHToBBPt->GetName()));
    TH1D *hGenHiggsPt_GGHToAATo4BPt_cl = (TH1D*)hGenHiggsPt_GGHToAATo4BPt->Clone(Form("%s_clone", hGenHiggsPt_GGHToAATo4BPt->GetName()));
    double PtPoint = 400;
    int binPtPoint = hGenHiggsPt_GGHToBBPt->FindBin(PtPoint);
    double sf_GGHToBBIncl = hGenHiggsPt_GGHToBBPt->GetBinContent(binPtPoint) / hGenHiggsPt_GGHToBBIncl->GetBinContent(binPtPoint);
    double sf_GGHToAATo4BPt = hGenHiggsPt_GGHToBBPt->GetBinContent(binPtPoint) / hGenHiggsPt_GGHToAATo4BPt->GetBinContent(binPtPoint);
    std::cout << "PtPoint for normalization: " << PtPoint << ", binPtPoint: " << binPtPoint << ", sf_GGHToBBIncl: " << sf_GGHToBBIncl << ", sf_GGHToAATo4BIncl: " << sf_GGHToAATo4BPt << "\n";
    hGenHiggsPt_GGHToBBIncl_cl->Scale(sf_GGHToBBIncl);
    hGenHiggsPt_GGHToAATo4BPt_cl->Scale(sf_GGHToAATo4BPt);
    hGenHiggsPt_GGHToBBIncl_cl->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHToBBIncl_cl->GetYaxis()->SetTitle("No. of events");
    hGenHiggsPt_GGHToBBIncl_cl->Draw();
    hGenHiggsPt_GGHToBBPt_cl->Draw("same");
    hGenHiggsPt_GGHToAATo4BPt_cl->Draw("same");

    TLegend *leg2 = new TLegend(0.5,0.75,0.99,0.99);
    leg2->AddEntry(hGenHiggsPt_GGHToBBIncl_cl, "GGToHTo2B inclusive", "lep");
    leg2->AddEntry(hGenHiggsPt_GGHToBBPt_cl, "GGToHTo2B pT(H)>200", "lep");
    leg2->AddEntry(hGenHiggsPt_GGHToAATo4BPt_cl, "GGToHToAATo4B pT(H)>150 reweighted", "lep");
    leg2->Draw();



    TCanvas *c3 = new TCanvas("c3", "c3", 600,500);
    c3->SetGrid();
    c3->cd();

    TH1D *hRatio_GGHToAATo4B = (TH1D*)hGenHiggsPt_GGHToBBPt_cl->Clone("hRatio_GGHToAATo4B");
    hRatio_GGHToAATo4B->SetMarkerStyle(20);
    hRatio_GGHToAATo4B->SetMarkerSize(0.5);   
    hRatio_GGHToAATo4B->SetLineColor(kRed);
    hRatio_GGHToAATo4B->SetMarkerColor(kRed);
    hRatio_GGHToAATo4B->GetYaxis()->SetRangeUser(0, 2);
    hRatio_GGHToAATo4B->Divide(hGenHiggsPt_GGHToAATo4BPt_cl, hGenHiggsPt_GGHToBBPt_cl);
    hRatio_GGHToAATo4B->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hRatio_GGHToAATo4B->GetYaxis()->SetTitle("No. of events");
    hRatio_GGHToAATo4B->Draw();


    //TH1D *hRatio_GGHToBB = (TH1D*)hGenHiggsPt_GGHToBBPt_cl->Clone("hRatio_GGHToBB");
    //hRatio_GGHToBB->SetLineColor(kBlack);
    //hRatio_GGHToBB->Divide(hGenHiggsPt_GGHToBBIncl_cl, hGenHiggsPt_GGHToBBPt_cl);
    //hRatio_GGHToBB->Draw("same");

    TLegend *leg3 = new TLegend(0.5,0.75,0.99,0.99);
    //leg3->AddEntry(hRatio_GGHToBB, "GGToHTo2B inclusive", "lep");
    leg3->AddEntry(hRatio_GGHToAATo4B, "GGToHToAATo4B pT(H)>150 reweighted", "lep");
    leg3->Draw();



    TCanvas *c4 = new TCanvas("c4", "c4", 600,500);
    c4->SetGrid();
    c4->cd();

    TH1D *hRatio_GGHToAATo4B_1 = (TH1D*)hGenHiggsPt_GGHToAATo4BPt->Clone("hRatio_GGHToAATo4B_1");
    hRatio_GGHToAATo4B_1->SetMarkerStyle(20);
    hRatio_GGHToAATo4B_1->SetMarkerSize(0.5);
    hRatio_GGHToAATo4B_1->SetLineColor(kRed);
    hRatio_GGHToAATo4B_1->SetMarkerColor(kRed);
    hRatio_GGHToAATo4B_1->GetYaxis()->SetRangeUser(0, 2);
    hRatio_GGHToAATo4B_1->Divide(hGenHiggsPt_GGHToAATo4BPt, hGenHiggsPt_GGHToBBPt);
    hRatio_GGHToAATo4B_1->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hRatio_GGHToAATo4B_1->GetYaxis()->SetTitle("No. of events");
    hRatio_GGHToAATo4B_1->Draw();


    TH1D *hRatio_GGHToBB_1 = (TH1D*)hGenHiggsPt_GGHToBBIncl->Clone("hRatio_GGHToBB_1");
    hRatio_GGHToBB_1->SetMarkerStyle(20);
    hRatio_GGHToBB_1->SetMarkerSize(0.5);
    hRatio_GGHToBB_1->SetLineColor(kBlack);
    hRatio_GGHToBB_1->SetMarkerColor(kBlack);
    hRatio_GGHToBB_1->Divide(hGenHiggsPt_GGHToBBIncl, hGenHiggsPt_GGHToBBPt);
    hRatio_GGHToBB_1->Draw("same");

    TLegend *leg4 = new TLegend(0.5,0.75,0.99,0.99);    
    leg4->AddEntry(hRatio_GGHToAATo4B_1, "GGToHToAATo4B pT(H)>150 reweighted", "lep");
    leg4->AddEntry(hRatio_GGHToBB_1, "GGToHTo2B inclusive", "lep");
    leg4->Draw();


    c1->SaveAs("hGenHPt_GGHTo2B_GGHToAATo4B_overlay.root");
    c2->SaveAs("hGenHPt_GGHTo2B_GGHToAATo4B_overlay_scaled.root");
    c3->SaveAs("hGenHPt_GGHTo2B_GGHToAATo4B_ratio_wScaled.root");
    c4->SaveAs("hGenHPt_GGHTo2B_GGHToAATo4B_ratio_woScaled.root");

}


void cal_GGFHiggsPtRewgts_check1(int mA = 20) {
    //int mA = 20;
    //TFile *fGGHToAATo4BIncl = new TFile("/eos/cms/store/user/ssawant/NanoAOD/SUSY_GluGluH_01J_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/hadded/hadded_NanoAOD.root");
    //TFile *fGGHToAATo4BPt   = new TFile(Form("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-%d_TuneCP5_13TeV_madgraph_pythia8/r1/PNet_v1.root",mA));

    double xsGGHToAATo4BIncl = 48.61;
    double xsGGHToAATo4BPt   = 48.61 * 0.057; 
    double xsGGHToAATo4BPt_v1 = 48.61 * 0.057 / 1.38; 

    TTree *trGGHToAATo4BIncl = getTTreeFromTChain(Form("/eos/cms/store/user/ssawant/NanoAOD/SUSY_GluGluH_01J_HToAATo4B_M-%d_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/*.root", mA), "Events");
    TTree *trGGHToAATo4BPt   = getTTreeFromTFile(Form("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-%d_TuneCP5_13TeV_madgraph_pythia8/r1/PNet_v1.root",mA), "Events");

    //TTree *trGGHToAATo4BIncl = (TTree*)fGGHToAATo4BIncl->Get("Events");
    long nEventsGenWgtPos_GGHToAATo4BIncl = trGGHToAATo4BIncl->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHToAATo4BIncl = trGGHToAATo4BIncl->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHToAATo4BIncl       = nEventsGenWgtPos_GGHToAATo4BIncl - nEventsGenWgtNeg_GGHToAATo4BIncl;
    std::cout << "GGHToAATo4BIncl nEvents: " << trGGHToAATo4BIncl->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHToAATo4BIncl: " << nEventsGenWgtPos_GGHToAATo4BIncl;
    std::cout << ", nEventsGenWgtNeg_GGHToAATo4BIncl: " << nEventsGenWgtNeg_GGHToAATo4BIncl;
    std::cout << ", nEventsTot_GGHToAATo4BIncl: " << nEventsTot_GGHToAATo4BIncl << "\n";

    //TTree *trGGHToAATo4BPt = (TTree*)fGGHToAATo4BPt->Get("Events");
    long nEventsGenWgtPos_GGHToAATo4BPt = trGGHToAATo4BPt->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHToAATo4BPt = trGGHToAATo4BPt->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHToAATo4BPt       = nEventsGenWgtPos_GGHToAATo4BPt - nEventsGenWgtNeg_GGHToAATo4BPt;
    std::cout << "GGHToAATo4BPt nEvents: " << trGGHToAATo4BPt->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHToAATo4BPt: " << nEventsGenWgtPos_GGHToAATo4BPt;
    std::cout << ", nEventsGenWgtNeg_GGHToAATo4BPt: " << nEventsGenWgtNeg_GGHToAATo4BPt;
    std::cout << ", nEventsTot_GGHToAATo4BPt: " << nEventsTot_GGHToAATo4BPt << "\n";


    TCanvas *c1 = new TCanvas("c1", "c1", 600,500);
    c1->SetLogy();
    c1->SetGrid();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);        
    c1->cd();
      
    TString sCut_GGHToAATo4BIncl = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld)", xsGGHToAATo4BIncl, nEventsTot_GGHToAATo4BIncl);
    std::cout << "sCut_GGHToAATo4BIncl: " <<  sCut_GGHToAATo4BIncl.Data() << "\n";
    trGGHToAATo4BIncl->SetLineColor(kRed);
    TH1D *hGenHiggsPt_GGHToAATo4BIncl = new TH1D("hGenHiggsPt_GGHToAATo4BIncl", "hGenHiggsPt_GGHToAATo4BIncl", 2000, 0, 2000); 
    hGenHiggsPt_GGHToAATo4BIncl->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BIncl->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BIncl->SetLineColor(kRed);
    hGenHiggsPt_GGHToAATo4BIncl->SetMarkerColor(kRed);    
    hGenHiggsPt_GGHToAATo4BIncl->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHToAATo4BIncl->GetYaxis()->SetTitle("No. of events");
    trGGHToAATo4BIncl->Draw("GenPart_pt >> hGenHiggsPt_GGHToAATo4BIncl", sCut_GGHToAATo4BIncl.Data());

    
    TString sCut_GGHToAATo4BPt = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld)", xsGGHToAATo4BPt, nEventsTot_GGHToAATo4BPt);
    std::cout << "sCut_GGHToAATo4BPt: " <<  sCut_GGHToAATo4BPt.Data() << "\n";
    trGGHToAATo4BPt->SetLineColor(kBlue);
    TH1D *hGenHiggsPt_GGHToAATo4BPt   = new TH1D("hGenHiggsPt_GGHToAATo4BPt", "hGenHiggsPt_GGHToAATo4BPt", 2000, 0, 2000);
    hGenHiggsPt_GGHToAATo4BPt->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BPt->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BPt->SetLineColor(kBlue);
    hGenHiggsPt_GGHToAATo4BPt->SetMarkerColor(kBlue);
    hGenHiggsPt_GGHToAATo4BPt->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHToAATo4BPt->GetYaxis()->SetTitle("No. of events");    
    trGGHToAATo4BPt->Draw("GenPart_pt >> hGenHiggsPt_GGHToAATo4BPt", sCut_GGHToAATo4BPt.Data(), "same");

    
    TString sCut_GGHToAATo4BPt_v1 = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld)", xsGGHToAATo4BPt_v1, nEventsTot_GGHToAATo4BPt);
    std::cout << "sCut_GGHToAATo4BPt_v1: " <<  sCut_GGHToAATo4BPt_v1.Data() << "\n";
    TH1D *hGenHiggsPt_GGHToAATo4BPt_v1   = new TH1D("hGenHiggsPt_GGHToAATo4BPt_v1", "hGenHiggsPt_GGHToAATo4BPt_v1", 2000, 0, 2000);
    hGenHiggsPt_GGHToAATo4BPt_v1->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BPt_v1->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BPt_v1->SetLineColor(kMagenta);
    hGenHiggsPt_GGHToAATo4BPt_v1->SetMarkerColor(kMagenta);
    hGenHiggsPt_GGHToAATo4BPt_v1->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHToAATo4BPt_v1->GetYaxis()->SetTitle("No. of events");    
    trGGHToAATo4BPt->Draw("GenPart_pt >> hGenHiggsPt_GGHToAATo4BPt_v1", sCut_GGHToAATo4BPt_v1.Data(), "same");


    TLegend *leg1 = new TLegend(0.5,0.75,0.99,0.99);
    leg1->AddEntry(hGenHiggsPt_GGHToAATo4BIncl, Form("GGToHToAATo4B inclusive, #sigma: %g",xsGGHToAATo4BIncl), "lep");
    leg1->AddEntry(hGenHiggsPt_GGHToAATo4BPt, Form("GGToHToAATo4B pT(H)>150 GeV, #sigma: %g",xsGGHToAATo4BPt), "lep");
    leg1->AddEntry(hGenHiggsPt_GGHToAATo4BPt_v1, Form("GGToHToAATo4B pT(H)>150 GeV, #sigma: %g",xsGGHToAATo4BPt_v1), "lep");
    leg1->Draw();

    double ptMin_forIntegral = 300;
    int binPtMin_forIntegral = hGenHiggsPt_GGHToAATo4BIncl->FindBin(ptMin_forIntegral);
    int nBins = hGenHiggsPt_GGHToAATo4BIncl->GetNbinsX();
    std::cout << "nEvents with pT > " << ptMin_forIntegral;
    std::cout <<  ", hGenHiggsPt_GGHToAATo4BIncl: " << hGenHiggsPt_GGHToAATo4BIncl->Integral(binPtMin_forIntegral, nBins);
    std::cout <<  ", hGenHiggsPt_GGHToAATo4BPt: " << hGenHiggsPt_GGHToAATo4BPt->Integral(binPtMin_forIntegral, nBins) << "\n";


    TCanvas *c2 = new TCanvas("c2", "c2", 600,500);
    //c2->SetLogy();
    c2->SetGrid();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);     
    c2->cd();


    TH1D *hRatio_GGHToAATo4B = (TH1D*)hGenHiggsPt_GGHToAATo4BPt->Clone("hRatio_GGHToAATo4B");
    hRatio_GGHToAATo4B->SetLineColor(kBlue);
    hRatio_GGHToAATo4B->SetMarkerColor(kBlue);
    hRatio_GGHToAATo4B->SetMarkerStyle(20);
    hRatio_GGHToAATo4B->SetMarkerSize(0.5);
    hRatio_GGHToAATo4B->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hRatio_GGHToAATo4B->GetYaxis()->SetTitle("Exclusive / Inclusive");    
    hRatio_GGHToAATo4B->Divide(hGenHiggsPt_GGHToAATo4BPt, hGenHiggsPt_GGHToAATo4BIncl);
    hRatio_GGHToAATo4B->Draw();

    TH1D *hRatio_GGHToAATo4B_v1 = (TH1D*)hGenHiggsPt_GGHToAATo4BPt_v1->Clone("hRatio_GGHToAATo4B_v1");
    hRatio_GGHToAATo4B_v1->SetLineColor(kMagenta);
    hRatio_GGHToAATo4B_v1->SetMarkerColor(kMagenta);
    hRatio_GGHToAATo4B_v1->SetMarkerStyle(20);
    hRatio_GGHToAATo4B_v1->SetMarkerSize(0.5);
    hRatio_GGHToAATo4B_v1->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hRatio_GGHToAATo4B_v1->GetYaxis()->SetTitle("Exclusive / Inclusive");    
    hRatio_GGHToAATo4B_v1->Divide(hGenHiggsPt_GGHToAATo4BPt_v1, hGenHiggsPt_GGHToAATo4BIncl);
    hRatio_GGHToAATo4B_v1->Draw("same");    

    TLegend *leg2 = new TLegend(0.5,0.75,0.99,0.99);
    leg2->AddEntry(hRatio_GGHToAATo4B, Form("GGToHToAATo4B pT(H)>150 GeV, #sigma: %g",xsGGHToAATo4BPt), "lep");
    leg2->AddEntry(hRatio_GGHToAATo4B_v1, Form("GGToHToAATo4B pT(H)>150 GeV, #sigma: %g",xsGGHToAATo4BPt_v1), "lep");
    leg2->Draw();    

    c1->SaveAs(Form("hGenHPt_GGHToAATo4B_overlay_mA%d.root",mA));
    c2->SaveAs(Form("hGenHPt_GGHToAATo4B_ratio_mA%d.root",mA));
}




void cal_GGFHiggsPtRewgts_check2() {
    TFile *fGGHTo2BIncl = new TFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8/hadd/GluGluHToBB_M-125_NNLOPS_NanoAODv9.root");
    TFile *fGGHTo2BPt   = new TFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/GluGluHToBB_Pt-200ToInf_M-125_TuneCP5_MINLO_13TeV-powheg-pythia8/hadd/GluGluHToBB_Pt-200ToInf_M-125_NanoAODv9.root");

    double xsGGHTo2BIncl = 48.61 * 0.582;
    double xsGGHTo2BPt   = 0.2740; 
    double xsGGHTo2BPt_v1 = 0.2740 / 0.572; // 0.2740 / 0.333; 

    TTree *trGGHTo2BIncl = (TTree*)fGGHTo2BIncl->Get("Events");
    long nEventsGenWgtPos_GGHTo2BIncl = trGGHTo2BIncl->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHTo2BIncl = trGGHTo2BIncl->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHTo2BIncl       = nEventsGenWgtPos_GGHTo2BIncl - nEventsGenWgtNeg_GGHTo2BIncl;
    std::cout << "GGHTo2BIncl nEvents: " << trGGHTo2BIncl->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHTo2BIncl: " << nEventsGenWgtPos_GGHTo2BIncl;
    std::cout << ", nEventsGenWgtNeg_GGHTo2BIncl: " << nEventsGenWgtNeg_GGHTo2BIncl;
    std::cout << ", nEventsTot_GGHTo2BIncl: " << nEventsTot_GGHTo2BIncl << "\n";

    TTree *trGGHTo2BPt = (TTree*)fGGHTo2BPt->Get("Events");
    long nEventsGenWgtPos_GGHTo2BPt = trGGHTo2BPt->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHTo2BPt = trGGHTo2BPt->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHTo2BPt       = nEventsGenWgtPos_GGHTo2BPt - nEventsGenWgtNeg_GGHTo2BPt;
    std::cout << "GGHTo2BPt nEvents: " << trGGHTo2BPt->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHTo2BPt: " << nEventsGenWgtPos_GGHTo2BPt;
    std::cout << ", nEventsGenWgtNeg_GGHTo2BPt: " << nEventsGenWgtNeg_GGHTo2BPt;
    std::cout << ", nEventsTot_GGHTo2BPt: " << nEventsTot_GGHTo2BPt << "\n";


    TCanvas *c1 = new TCanvas("c1", "c1", 600,500);
    c1->SetLogy();
    c1->SetGrid();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);        
    c1->cd();
      
    TString sCut_GGHTo2BIncl = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld)", xsGGHTo2BIncl, nEventsTot_GGHTo2BIncl);
    std::cout << "sCut_GGHTo2BIncl: " <<  sCut_GGHTo2BIncl.Data() << "\n";
    trGGHTo2BIncl->SetLineColor(kRed);
    TH1D *hGenHiggsPt_GGHTo2BIncl = new TH1D("hGenHiggsPt_GGHTo2BIncl", "hGenHiggsPt_GGHTo2BIncl", 150, 0, 1500); 
    hGenHiggsPt_GGHTo2BIncl->SetMarkerStyle(20);
    hGenHiggsPt_GGHTo2BIncl->SetMarkerSize(0.5);
    hGenHiggsPt_GGHTo2BIncl->SetLineColor(kRed);
    hGenHiggsPt_GGHTo2BIncl->SetMarkerColor(kRed);    
    hGenHiggsPt_GGHTo2BIncl->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHTo2BIncl->GetYaxis()->SetTitle("No. of events");
    trGGHTo2BIncl->Draw("GenPart_pt >> hGenHiggsPt_GGHTo2BIncl", sCut_GGHTo2BIncl.Data());

    
    TString sCut_GGHTo2BPt = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld)", xsGGHTo2BPt, nEventsTot_GGHTo2BPt);
    std::cout << "sCut_GGHTo2BPt: " <<  sCut_GGHTo2BPt.Data() << "\n";
    trGGHTo2BPt->SetLineColor(kBlue);
    TH1D *hGenHiggsPt_GGHTo2BPt   = new TH1D("hGenHiggsPt_GGHTo2BPt", "hGenHiggsPt_GGHTo2BPt", 150, 0, 1500);
    hGenHiggsPt_GGHTo2BPt->SetMarkerStyle(20);
    hGenHiggsPt_GGHTo2BPt->SetMarkerSize(0.5);
    hGenHiggsPt_GGHTo2BPt->SetLineColor(kBlue);
    hGenHiggsPt_GGHTo2BPt->SetMarkerColor(kBlue);
    hGenHiggsPt_GGHTo2BPt->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHTo2BPt->GetYaxis()->SetTitle("No. of events");    
    trGGHTo2BPt->Draw("GenPart_pt >> hGenHiggsPt_GGHTo2BPt", sCut_GGHTo2BPt.Data(), "same");

    
    TString sCut_GGHTo2BPt_v1 = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld)", xsGGHTo2BPt_v1, nEventsTot_GGHTo2BPt);
    std::cout << "sCut_GGHTo2BPt_v1: " <<  sCut_GGHTo2BPt_v1.Data() << "\n";
    TH1D *hGenHiggsPt_GGHTo2BPt_v1   = new TH1D("hGenHiggsPt_GGHTo2BPt_v1", "hGenHiggsPt_GGHTo2BPt_v1", 150, 0, 1500);
    hGenHiggsPt_GGHTo2BPt_v1->SetMarkerStyle(20);
    hGenHiggsPt_GGHTo2BPt_v1->SetMarkerSize(0.5);
    hGenHiggsPt_GGHTo2BPt_v1->SetLineColor(kMagenta);
    hGenHiggsPt_GGHTo2BPt_v1->SetMarkerColor(kMagenta);
    hGenHiggsPt_GGHTo2BPt_v1->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHTo2BPt_v1->GetYaxis()->SetTitle("No. of events");    
    trGGHTo2BPt->Draw("GenPart_pt >> hGenHiggsPt_GGHTo2BPt_v1", sCut_GGHTo2BPt_v1.Data(), "same");
    

    TLegend *leg1 = new TLegend(0.5,0.75,0.99,0.99);
    leg1->AddEntry(hGenHiggsPt_GGHTo2BIncl, Form("GGToHToAATo4B inclusive, #sigma: %g",xsGGHTo2BIncl), "lep");
    leg1->AddEntry(hGenHiggsPt_GGHTo2BPt, Form("GGToHToAATo4B pT(H)>150 GeV, #sigma: %g",xsGGHTo2BPt), "lep");
    leg1->AddEntry(hGenHiggsPt_GGHTo2BPt_v1, Form("GGToHToAATo4B pT(H)>150 GeV, #sigma: %g",xsGGHTo2BPt_v1), "lep");
    leg1->Draw();

    double ptMin_forIntegral = 350;
    int binPtMin_forIntegral = hGenHiggsPt_GGHTo2BIncl->FindBin(ptMin_forIntegral);
    int nBins = hGenHiggsPt_GGHTo2BIncl->GetNbinsX();
    std::cout << "nEvents with pT > " << ptMin_forIntegral;
    std::cout <<  ", hGenHiggsPt_GGHTo2BIncl: " << hGenHiggsPt_GGHTo2BIncl->Integral(binPtMin_forIntegral, nBins);
    std::cout <<  ", hGenHiggsPt_GGHTo2BPt: " << hGenHiggsPt_GGHTo2BPt->Integral(binPtMin_forIntegral, nBins) << "\n";


    TCanvas *c2 = new TCanvas("c2", "c2", 600,500);
    //c2->SetLogy();
    c2->SetGrid();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);     
    c2->cd();


    TH1D *hRatio_GGHTo2B = (TH1D*)hGenHiggsPt_GGHTo2BPt->Clone("hRatio_GGHTo2B");
    hRatio_GGHTo2B->SetLineColor(kBlue);
    hRatio_GGHTo2B->SetMarkerColor(kBlue);
    hRatio_GGHTo2B->SetMarkerStyle(20);
    hRatio_GGHTo2B->SetMarkerSize(0.5);
    hRatio_GGHTo2B->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hRatio_GGHTo2B->GetYaxis()->SetTitle("Exclusive / Inclusive");  
    hRatio_GGHTo2B->GetYaxis()->SetRangeUser(0, 1.7);  
    hRatio_GGHTo2B->Divide(hGenHiggsPt_GGHTo2BPt, hGenHiggsPt_GGHTo2BIncl);
    hRatio_GGHTo2B->Draw();

    
    TH1D *hRatio_GGHTo2B_v1 = (TH1D*)hGenHiggsPt_GGHTo2BPt_v1->Clone("hRatio_GGHTo2B_v1");
    hRatio_GGHTo2B_v1->SetLineColor(kMagenta);
    hRatio_GGHTo2B_v1->SetMarkerColor(kMagenta);
    hRatio_GGHTo2B_v1->SetMarkerStyle(20);
    hRatio_GGHTo2B_v1->SetMarkerSize(0.5);
    hRatio_GGHTo2B_v1->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hRatio_GGHTo2B_v1->GetYaxis()->SetTitle("Exclusive / Inclusive");    
    hRatio_GGHTo2B_v1->Divide(hGenHiggsPt_GGHTo2BPt_v1, hGenHiggsPt_GGHTo2BIncl);
    hRatio_GGHTo2B_v1->Draw("same");    
    

    TLegend *leg2 = new TLegend(0.5,0.75,0.99,0.99);
    leg2->AddEntry(hRatio_GGHTo2B, Form("GGToHToAATo4B pT(H)>150 GeV, #sigma: %g",xsGGHTo2BPt), "lep");
    leg2->AddEntry(hRatio_GGHTo2B_v1, Form("GGToHToAATo4B pT(H)>150 GeV, #sigma: %g",xsGGHTo2BPt_v1), "lep");
    leg2->Draw();    

    c1->SaveAs("hGenHPt_GGHTo2B_overlay.root");
    c2->SaveAs("hGenHPt_GGHTo2B_ratio.root");
}


void cal_GGFHiggsPtRewgts_run() {
    //TFile *fGGHToBBIncl = new TFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8/hadd/GluGluHToBB_M-125_NNLOPS_NanoAODv9.root");
    //TFile *fGGHToBBPt   = new TFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/GluGluHToBB_Pt-200ToInf_M-125_TuneCP5_MINLO_13TeV-powheg-pythia8/hadd/GluGluHToBB_Pt-200ToInf_M-125_NanoAODv9.root");
    //TFile *fGGHToAATo4BPt = new TFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-All_TuneCP5_13TeV_madgraph_pythia8/hadd/GluGluH_01J_HToAATo4B_Pt150_M-All_NanoAODv9.root");

    //double xsGGHToBBIncl = 48.61 * 0.582;
    //double xsGGHToBBPt   = 0.2740; 
    //double xsGGHToAATo4BPt = 48.61 * 0.057;
    // scale GGHToBB and GGHToAATo4B to 48.61 pb cross section
    double xsGGHToBBIncl     = 1 ;
    double xsGGHToBBPt       = 1; // 0.2740; 
    double xsGGHToAATo4BIncl = 1;
    double xsGGHToAATo4BPt   =  1; 

    //TTree *trGGHToBBIncl = (TTree*)fGGHToBBIncl->Get("Events");
    TTree *trGGHToBBIncl     = getTTreeFromTFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/GluGluHToBB_M-125_TuneCP5_MINLO_NNLOPS_13TeV-powheg-pythia8/hadd/GluGluHToBB_M-125_NNLOPS_NanoAODv9.root", "Events");
    TTree *trGGHToBBPt       = getTTreeFromTFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/GluGluHToBB_Pt-200ToInf_M-125_TuneCP5_MINLO_13TeV-powheg-pythia8/hadd/GluGluHToBB_Pt-200ToInf_M-125_NanoAODv9.root", "Events");  
    TTree *trGGHToAATo4BIncl = getTTreeFromTChain("/eos/cms/store/user/ssawant/NanoAOD/SUSY_GluGluH_01J_HToAATo4B_M-*/*/NANOAODSIM/*.root", "Events");
    TTree *trGGHToAATo4BPt   = getTTreeFromTFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-All_TuneCP5_13TeV_madgraph_pythia8/hadd/GluGluH_01J_HToAATo4B_Pt150_M-All_NanoAODv9.root", "Events");  
    

    long nEventsGenWgtPos_GGHToBBIncl = trGGHToBBIncl->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHToBBIncl = trGGHToBBIncl->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHToBBIncl       = nEventsGenWgtPos_GGHToBBIncl - nEventsGenWgtNeg_GGHToBBIncl;
    std::cout << "GGHToBBIncl nEvents: " << trGGHToBBIncl->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHToBBIncl: " << nEventsGenWgtPos_GGHToBBIncl;
    std::cout << ", nEventsGenWgtNeg_GGHToBBIncl: " << nEventsGenWgtNeg_GGHToBBIncl;
    std::cout << ", nEventsTot_GGHToBBIncl: " << nEventsTot_GGHToBBIncl << "\n";

    //TTree *trGGHToBBPt = (TTree*)fGGHToBBPt->Get("Events");
    long nEventsGenWgtPos_GGHToBBPt = trGGHToBBPt->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHToBBPt = trGGHToBBPt->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHToBBPt       = nEventsGenWgtPos_GGHToBBPt - nEventsGenWgtNeg_GGHToBBPt;
    std::cout << "GGHToBBPt nEvents: " << trGGHToBBPt->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHToBBPt: " << nEventsGenWgtPos_GGHToBBPt;
    std::cout << ", nEventsGenWgtNeg_GGHToBBPt: " << nEventsGenWgtNeg_GGHToBBPt;
    std::cout << ", nEventsTot_GGHToBBPt: " << nEventsTot_GGHToBBPt << "\n";

    long nEventsGenWgtPos_GGHToAATo4BIncl = trGGHToAATo4BIncl->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHToAATo4BIncl = trGGHToAATo4BIncl->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHToAATo4BIncl       = nEventsGenWgtPos_GGHToAATo4BIncl - nEventsGenWgtNeg_GGHToAATo4BIncl;
    std::cout << "GGHToAATo4BIncl nEvents: " << trGGHToAATo4BIncl->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHToAATo4BIncl: " << nEventsGenWgtPos_GGHToAATo4BIncl;
    std::cout << ", nEventsGenWgtNeg_GGHToAATo4BIncl: " << nEventsGenWgtNeg_GGHToAATo4BIncl;
    std::cout << ", nEventsTot_GGHToAATo4BIncl: " << nEventsTot_GGHToAATo4BIncl << "\n";

    //TTree *trGGHToAATo4BPt = (TTree*)fGGHToAATo4BPt->Get("Events");
    long nEventsGenWgtPos_GGHToAATo4BPt = trGGHToAATo4BPt->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHToAATo4BPt = trGGHToAATo4BPt->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHToAATo4BPt       = nEventsGenWgtPos_GGHToAATo4BPt - nEventsGenWgtNeg_GGHToAATo4BPt;
    std::cout << "GGHToAATo4BPt nEvents: " << trGGHToAATo4BPt->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHToAATo4BPt: " << nEventsGenWgtPos_GGHToAATo4BPt;
    std::cout << ", nEventsGenWgtNeg_GGHToAATo4BPt: " << nEventsGenWgtNeg_GGHToAATo4BPt;
    std::cout << ", nEventsTot_GGHToAATo4BPt: " << nEventsTot_GGHToAATo4BPt << "\n";

    /*
    No. of events in the samples: prints
    GGHToBBIncl nEvents: 14861428
    nEventsGenWgtPos_GGHToBBIncl: 14336537, nEventsGenWgtNeg_GGHToBBIncl: 524891, nEventsTot_GGHToBBIncl: 13811646
    GGHToBBPt nEvents: 497000
    nEventsGenWgtPos_GGHToBBPt: 496887, nEventsGenWgtNeg_GGHToBBPt: 113, nEventsTot_GGHToBBPt: 496774
    GGHToAATo4BIncl nEvents: 4394177
    nEventsGenWgtPos_GGHToAATo4BIncl: 4394177, nEventsGenWgtNeg_GGHToAATo4BIncl: 0, nEventsTot_GGHToAATo4BIncl: 4394177
    GGHToAATo4BPt nEvents: 5221010
    nEventsGenWgtPos_GGHToAATo4BPt: 5221010, nEventsGenWgtNeg_GGHToAATo4BPt: 0, nEventsTot_GGHToAATo4BPt: 5221010    
    */

    int nBinsPt = 500;    double PtMin = 0;    double PtMax = 2000;
    int nBinsLog2Pt = 500;    double Log2PtMin = log2(1);    double Log2PtMax = log2(2000);
    int nBinsLog10Pt = 500;    double Log10PtMin = log10(1);    double Log10PtMax = log10(2000);

    // GGHToBBIncl ----------------------------------------------------------------------------------------------------------------------------------------------
    TString sCut_GGHToBBIncl = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld)", xsGGHToBBIncl, nEventsTot_GGHToBBIncl);
    std::cout << "sCut_GGHToBBIncl: " <<  sCut_GGHToBBIncl.Data() << "\n";

    TH1D *hGenHiggsPt_GGHToBBIncl   = new TH1D("hGenHiggsPt_GGHToBBIncl", "hGenHiggsPt_GGHToBBIncl", nBinsPt, PtMin, PtMax);
    hGenHiggsPt_GGHToBBIncl->SetMarkerStyle(20);
    hGenHiggsPt_GGHToBBIncl->SetMarkerSize(0.5);   
    hGenHiggsPt_GGHToBBIncl->SetLineColor(kBlue);
    hGenHiggsPt_GGHToBBIncl->SetMarkerColor(kBlue);
    hGenHiggsPt_GGHToBBIncl->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    trGGHToBBIncl->Draw("GenPart_pt >> hGenHiggsPt_GGHToBBIncl", sCut_GGHToBBIncl.Data(), "goff");

    std::cout << "sCut_GGHToBBIncl: " <<  sCut_GGHToBBIncl.Data() << ", log2Pt \n";
    TH1D *hGenHiggsLog2Pt_GGHToBBIncl   = new TH1D("hGenHiggsLog2Pt_GGHToBBIncl", "hGenHiggsLog2Pt_GGHToBBIncl", nBinsLog2Pt, Log2PtMin, Log2PtMax);
    hGenHiggsLog2Pt_GGHToBBIncl->SetMarkerStyle(20);
    hGenHiggsLog2Pt_GGHToBBIncl->SetMarkerSize(0.5);   
    hGenHiggsLog2Pt_GGHToBBIncl->SetLineColor(kBlue);
    hGenHiggsLog2Pt_GGHToBBIncl->SetMarkerColor(kBlue);
    hGenHiggsLog2Pt_GGHToBBIncl->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    trGGHToBBIncl->Draw("log2(GenPart_pt) >> hGenHiggsLog2Pt_GGHToBBIncl", sCut_GGHToBBIncl.Data(), "goff");

    std::cout << "sCut_GGHToBBIncl: " <<  sCut_GGHToBBIncl.Data() << ", log10Pt \n";

    TH1D *hGenHiggsLog10Pt_GGHToBBIncl   = new TH1D("hGenHiggsLog10Pt_GGHToBBIncl", "hGenHiggsLog10Pt_GGHToBBIncl", nBinsLog10Pt, Log10PtMin, Log10PtMax);
    hGenHiggsLog10Pt_GGHToBBIncl->SetMarkerStyle(20);
    hGenHiggsLog10Pt_GGHToBBIncl->SetMarkerSize(0.5);   
    hGenHiggsLog10Pt_GGHToBBIncl->SetLineColor(kBlue);
    hGenHiggsLog10Pt_GGHToBBIncl->SetMarkerColor(kBlue);
    hGenHiggsLog10Pt_GGHToBBIncl->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    trGGHToBBIncl->Draw("log10(GenPart_pt) >> hGenHiggsLog10Pt_GGHToBBIncl", sCut_GGHToBBIncl.Data(), "goff");


    // GGHToBBPt ------------------------------------------------------------------------------------------------------------------------------------------
    TString sCut_GGHToBBPt = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld)", xsGGHToBBPt, nEventsTot_GGHToBBPt);
    std::cout << "sCut_GGHToBBPt: " <<  sCut_GGHToBBPt.Data() << "\n";

    TH1D *hGenHiggsPt_GGHToBBPt   = new TH1D("hGenHiggsPt_GGHToBBPt", "hGenHiggsPt_GGHToBBPt", nBinsPt, PtMin, PtMax);
    hGenHiggsPt_GGHToBBPt->SetMarkerStyle(20);
    hGenHiggsPt_GGHToBBPt->SetMarkerSize(0.5);   
    hGenHiggsPt_GGHToBBPt->SetLineColor(kBlue);
    hGenHiggsPt_GGHToBBPt->SetMarkerColor(kBlue);
    hGenHiggsPt_GGHToBBPt->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    trGGHToBBPt->Draw("GenPart_pt >> hGenHiggsPt_GGHToBBPt", sCut_GGHToBBPt.Data(), "goff");

    std::cout << "sCut_GGHToBBPt: " <<  sCut_GGHToBBPt.Data() << ", log2Pt \n";
    TH1D *hGenHiggsLog2Pt_GGHToBBPt   = new TH1D("hGenHiggsLog2Pt_GGHToBBPt", "hGenHiggsLog2Pt_GGHToBBPt", nBinsLog2Pt, Log2PtMin, Log2PtMax);
    hGenHiggsLog2Pt_GGHToBBPt->SetMarkerStyle(20);
    hGenHiggsLog2Pt_GGHToBBPt->SetMarkerSize(0.5);   
    hGenHiggsLog2Pt_GGHToBBPt->SetLineColor(kBlue);
    hGenHiggsLog2Pt_GGHToBBPt->SetMarkerColor(kBlue);
    hGenHiggsLog2Pt_GGHToBBPt->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    trGGHToBBPt->Draw("log2(GenPart_pt) >> hGenHiggsLog2Pt_GGHToBBPt", sCut_GGHToBBPt.Data(), "goff");

    std::cout << "sCut_GGHToBBPt: " <<  sCut_GGHToBBPt.Data() << ", log10Pt \n";

    TH1D *hGenHiggsLog10Pt_GGHToBBPt   = new TH1D("hGenHiggsLog10Pt_GGHToBBPt", "hGenHiggsLog10Pt_GGHToBBPt", nBinsLog10Pt, Log10PtMin, Log10PtMax);
    hGenHiggsLog10Pt_GGHToBBPt->SetMarkerStyle(20);
    hGenHiggsLog10Pt_GGHToBBPt->SetMarkerSize(0.5);   
    hGenHiggsLog10Pt_GGHToBBPt->SetLineColor(kBlue);
    hGenHiggsLog10Pt_GGHToBBPt->SetMarkerColor(kBlue);
    hGenHiggsLog10Pt_GGHToBBPt->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    trGGHToBBPt->Draw("log10(GenPart_pt) >> hGenHiggsLog10Pt_GGHToBBPt", sCut_GGHToBBPt.Data(), "goff");


    // GGHToAATo4BIncl ------------------------------------------------------------------------------------------------------------------------------------------------------------
    TString sCut_GGHToAATo4BIncl_0 = Form("( (GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld) )", xsGGHToAATo4BIncl, nEventsTot_GGHToAATo4BIncl);
    std::cout << "sCut_GGHToAATo4BIncl_0: " <<  sCut_GGHToAATo4BIncl_0.Data() << "\n";

    TH1D *hGenHiggsPt_GGHToAATo4BIncl_0   = new TH1D("hGenHiggsPt_GGHToAATo4BIncl_0", "hGenHiggsPt_GGHToAATo4BIncl_0", nBinsPt, PtMin, PtMax);
    hGenHiggsPt_GGHToAATo4BIncl_0->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BIncl_0->SetMarkerSize(0.5);     
    hGenHiggsPt_GGHToAATo4BIncl_0->SetLineColor(kRed);
    hGenHiggsPt_GGHToAATo4BIncl_0->SetMarkerColor(kRed);
    trGGHToAATo4BIncl->Draw("GenPart_pt >> hGenHiggsPt_GGHToAATo4BIncl_0", sCut_GGHToAATo4BIncl_0.Data(), "goff");

    std::cout << "sCut_GGHToAATo4BIncl_0: " <<  sCut_GGHToAATo4BIncl_0.Data() << ", log2Pt\n";

    TH1D *hGenHiggsLog2Pt_GGHToAATo4BIncl_0   = new TH1D("hGenHiggsLog2Pt_GGHToAATo4BIncl_0", "hGenHiggsLog2Pt_GGHToAATo4BIncl_0", nBinsLog2Pt, Log2PtMin, Log2PtMax);
    hGenHiggsLog2Pt_GGHToAATo4BIncl_0->SetMarkerStyle(20);
    hGenHiggsLog2Pt_GGHToAATo4BIncl_0->SetMarkerSize(0.5);     
    hGenHiggsLog2Pt_GGHToAATo4BIncl_0->SetLineColor(kRed);
    hGenHiggsLog2Pt_GGHToAATo4BIncl_0->SetMarkerColor(kRed);
    trGGHToAATo4BIncl->Draw("log2(GenPart_pt) >> hGenHiggsLog2Pt_GGHToAATo4BIncl_0", sCut_GGHToAATo4BIncl_0.Data(), "goff");

    std::cout << "sCut_GGHToAATo4BIncl_0: " <<  sCut_GGHToAATo4BIncl_0.Data() << ", log10Pt\n";

    TH1D *hGenHiggsLog10Pt_GGHToAATo4BIncl_0   = new TH1D("hGenHiggsLog10Pt_GGHToAATo4BIncl_0", "hGenHiggsLog10Pt_GGHToAATo4BIncl_0", nBinsLog10Pt, Log10PtMin, Log10PtMax);
    hGenHiggsLog10Pt_GGHToAATo4BIncl_0->SetMarkerStyle(20);
    hGenHiggsLog10Pt_GGHToAATo4BIncl_0->SetMarkerSize(0.5);     
    hGenHiggsLog10Pt_GGHToAATo4BIncl_0->SetLineColor(kRed);
    hGenHiggsLog10Pt_GGHToAATo4BIncl_0->SetMarkerColor(kRed);
    trGGHToAATo4BIncl->Draw("log10(GenPart_pt) >> hGenHiggsLog10Pt_GGHToAATo4BIncl_0", sCut_GGHToAATo4BIncl_0.Data(), "goff");


    // GGHToAATo4BIncl Reweighted --------------------------------------------------------------------------------------------------------------------------------------------------------
    TString sCut_GGHToAATo4BIncl_Rewgted = Form("( (GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld) * min(max(1.45849 + -0.00400668*GenPart_pt + 4.02577e-06*pow(GenPart_pt, 2) + -1.38804e-09*pow(GenPart_pt, 3), 0.09), 1.02) )", xsGGHToAATo4BIncl, nEventsTot_GGHToAATo4BIncl);
    std::cout << "sCut_GGHToAATo4BIncl_Rewgted: " <<  sCut_GGHToAATo4BIncl_Rewgted.Data() << "\n";

    TH1D *hGenHiggsPt_GGHToAATo4BIncl_Rewgted   = new TH1D("hGenHiggsPt_GGHToAATo4BIncl_Rewgted", "hGenHiggsPt_GGHToAATo4BIncl_Rewgted", nBinsPt, PtMin, PtMax);
    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->SetMarkerSize(0.5);     
    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->SetLineColor(kRed);
    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->SetMarkerColor(kRed);
    trGGHToAATo4BIncl->Draw("GenPart_pt >> hGenHiggsPt_GGHToAATo4BIncl_Rewgted", sCut_GGHToAATo4BIncl_Rewgted.Data(), "goff");

    std::cout << "sCut_GGHToAATo4BIncl_Rewgted: " <<  sCut_GGHToAATo4BIncl_Rewgted.Data() << ", log2Pt\n";

    TH1D *hGenHiggsLog2Pt_GGHToAATo4BIncl_Rewgted   = new TH1D("hGenHiggsLog2Pt_GGHToAATo4BIncl_Rewgted", "hGenHiggsLog2Pt_GGHToAATo4BIncl_Rewgted", nBinsLog2Pt, Log2PtMin, Log2PtMax);
    hGenHiggsLog2Pt_GGHToAATo4BIncl_Rewgted->SetMarkerStyle(20);
    hGenHiggsLog2Pt_GGHToAATo4BIncl_Rewgted->SetMarkerSize(0.5);     
    hGenHiggsLog2Pt_GGHToAATo4BIncl_Rewgted->SetLineColor(kRed);
    hGenHiggsLog2Pt_GGHToAATo4BIncl_Rewgted->SetMarkerColor(kRed);
    trGGHToAATo4BIncl->Draw("log2(GenPart_pt) >> hGenHiggsLog2Pt_GGHToAATo4BIncl_Rewgted", sCut_GGHToAATo4BIncl_Rewgted.Data(), "goff");

    std::cout << "sCut_GGHToAATo4BIncl_Rewgted: " <<  sCut_GGHToAATo4BIncl_Rewgted.Data() << ", log10Pt\n";

    TH1D *hGenHiggsLog10Pt_GGHToAATo4BIncl_Rewgted   = new TH1D("hGenHiggsLog10Pt_GGHToAATo4BIncl_Rewgted", "hGenHiggsLog10Pt_GGHToAATo4BIncl_Rewgted", nBinsLog10Pt, Log10PtMin, Log10PtMax);
    hGenHiggsLog10Pt_GGHToAATo4BIncl_Rewgted->SetMarkerStyle(20);
    hGenHiggsLog10Pt_GGHToAATo4BIncl_Rewgted->SetMarkerSize(0.5);     
    hGenHiggsLog10Pt_GGHToAATo4BIncl_Rewgted->SetLineColor(kRed);
    hGenHiggsLog10Pt_GGHToAATo4BIncl_Rewgted->SetMarkerColor(kRed);
    trGGHToAATo4BIncl->Draw("log10(GenPart_pt) >> hGenHiggsLog10Pt_GGHToAATo4BIncl_Rewgted", sCut_GGHToAATo4BIncl_Rewgted.Data(), "goff");



    // GGHToAATo4BPt -----------------------------------------------------------------------------------------------------------------------------------------------------------
    TString sCut_GGHToAATo4BPt_0 = Form("( (GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld) )", xsGGHToAATo4BPt, nEventsTot_GGHToAATo4BPt);
    std::cout << "sCut_GGHToAATo4BPt_0: " <<  sCut_GGHToAATo4BPt_0.Data() << "\n";

    TH1D *hGenHiggsPt_GGHToAATo4BPt_0   = new TH1D("hGenHiggsPt_GGHToAATo4BPt_0", "hGenHiggsPt_GGHToAATo4BPt_0", nBinsPt, PtMin, PtMax);
    hGenHiggsPt_GGHToAATo4BPt_0->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BPt_0->SetMarkerSize(0.5);     
    hGenHiggsPt_GGHToAATo4BPt_0->SetLineColor(kRed);
    hGenHiggsPt_GGHToAATo4BPt_0->SetMarkerColor(kRed);
    trGGHToAATo4BPt->Draw("GenPart_pt >> hGenHiggsPt_GGHToAATo4BPt_0", sCut_GGHToAATo4BPt_0.Data(), "goff");

    std::cout << "sCut_GGHToAATo4BPt_0: " <<  sCut_GGHToAATo4BPt_0.Data() << ", log2Pt\n";

    TH1D *hGenHiggsLog2Pt_GGHToAATo4BPt_0   = new TH1D("hGenHiggsLog2Pt_GGHToAATo4BPt_0", "hGenHiggsLog2Pt_GGHToAATo4BPt_0", nBinsLog2Pt, Log2PtMin, Log2PtMax);
    hGenHiggsLog2Pt_GGHToAATo4BPt_0->SetMarkerStyle(20);
    hGenHiggsLog2Pt_GGHToAATo4BPt_0->SetMarkerSize(0.5);     
    hGenHiggsLog2Pt_GGHToAATo4BPt_0->SetLineColor(kRed);
    hGenHiggsLog2Pt_GGHToAATo4BPt_0->SetMarkerColor(kRed);
    trGGHToAATo4BPt->Draw("log2(GenPart_pt) >> hGenHiggsLog2Pt_GGHToAATo4BPt_0", sCut_GGHToAATo4BPt_0.Data(), "goff");

    std::cout << "sCut_GGHToAATo4BPt_0: " <<  sCut_GGHToAATo4BPt_0.Data() << ", log10Pt\n";

    TH1D *hGenHiggsLog10Pt_GGHToAATo4BPt_0   = new TH1D("hGenHiggsLog10Pt_GGHToAATo4BPt_0", "hGenHiggsLog10Pt_GGHToAATo4BPt_0", nBinsLog10Pt, Log10PtMin, Log10PtMax);
    hGenHiggsLog10Pt_GGHToAATo4BPt_0->SetMarkerStyle(20);
    hGenHiggsLog10Pt_GGHToAATo4BPt_0->SetMarkerSize(0.5);     
    hGenHiggsLog10Pt_GGHToAATo4BPt_0->SetLineColor(kRed);
    hGenHiggsLog10Pt_GGHToAATo4BPt_0->SetMarkerColor(kRed);
    trGGHToAATo4BPt->Draw("log10(GenPart_pt) >> hGenHiggsLog10Pt_GGHToAATo4BPt_0", sCut_GGHToAATo4BPt_0.Data(), "goff");


    // GGHToAATo4BPt reweighted -------------------------------------------------------------------------------------------------------------------------------------------------
    TString sCut_GGHToAATo4BPt_Rewgted = Form("( (GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld) * min(max(1.45849 + -0.00400668*GenPart_pt + 4.02577e-06*pow(GenPart_pt, 2) + -1.38804e-09*pow(GenPart_pt, 3), 0.09), 1.02) )", xsGGHToAATo4BPt, nEventsTot_GGHToAATo4BPt);
    std::cout << "sCut_GGHToAATo4BPt_Rewgted: " <<  sCut_GGHToAATo4BPt_Rewgted.Data() << "\n";

    TH1D *hGenHiggsPt_GGHToAATo4BPt_Rewgted   = new TH1D("hGenHiggsPt_GGHToAATo4BPt_Rewgted", "hGenHiggsPt_GGHToAATo4BPt_Rewgted", nBinsPt, PtMin, PtMax);
    hGenHiggsPt_GGHToAATo4BPt_Rewgted->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BPt_Rewgted->SetMarkerSize(0.5);     
    hGenHiggsPt_GGHToAATo4BPt_Rewgted->SetLineColor(kBlue);
    hGenHiggsPt_GGHToAATo4BPt_Rewgted->SetMarkerColor(kBlue);
    trGGHToAATo4BPt->Draw("GenPart_pt >> hGenHiggsPt_GGHToAATo4BPt_Rewgted", sCut_GGHToAATo4BPt_Rewgted.Data(), "goff");

    std::cout << "sCut_GGHToAATo4BPt_Rewgted: " <<  sCut_GGHToAATo4BPt_Rewgted.Data() << ", log2Pt\n";

    TH1D *hGenHiggsLog2Pt_GGHToAATo4BPt_Rewgted   = new TH1D("hGenHiggsLog2Pt_GGHToAATo4BPt_Rewgted", "hGenHiggsLog2Pt_GGHToAATo4BPt_Rewgted", nBinsLog2Pt, Log2PtMin, Log2PtMax);
    hGenHiggsLog2Pt_GGHToAATo4BPt_Rewgted->SetMarkerStyle(20);
    hGenHiggsLog2Pt_GGHToAATo4BPt_Rewgted->SetMarkerSize(0.5);     
    hGenHiggsLog2Pt_GGHToAATo4BPt_Rewgted->SetLineColor(kBlue);
    hGenHiggsLog2Pt_GGHToAATo4BPt_Rewgted->SetMarkerColor(kBlue);
    trGGHToAATo4BPt->Draw("log2(GenPart_pt) >> hGenHiggsLog2Pt_GGHToAATo4BPt_Rewgted", sCut_GGHToAATo4BPt_Rewgted.Data(), "goff");

    std::cout << "sCut_GGHToAATo4BPt_Rewgted: " <<  sCut_GGHToAATo4BPt_Rewgted.Data() << ", log10Pt\n";

    TH1D *hGenHiggsLog10Pt_GGHToAATo4BPt_Rewgted   = new TH1D("hGenHiggsLog10Pt_GGHToAATo4BPt_Rewgted", "hGenHiggsLog10Pt_GGHToAATo4BPt_Rewgted", nBinsLog10Pt, Log10PtMin, Log10PtMax);
    hGenHiggsLog10Pt_GGHToAATo4BPt_Rewgted->SetMarkerStyle(20);
    hGenHiggsLog10Pt_GGHToAATo4BPt_Rewgted->SetMarkerSize(0.5);     
    hGenHiggsLog10Pt_GGHToAATo4BPt_Rewgted->SetLineColor(kBlue);
    hGenHiggsLog10Pt_GGHToAATo4BPt_Rewgted->SetMarkerColor(kBlue);
    trGGHToAATo4BPt->Draw("log10(GenPart_pt) >> hGenHiggsLog10Pt_GGHToAATo4BPt_Rewgted", sCut_GGHToAATo4BPt_Rewgted.Data(), "goff");




    TFile *fOut=new TFile("hGenHiggsPt_GGHToBB_GGHToAATo4B_xsUnity_4GeVBins.root", "recreate");
    fOut->cd();

    hGenHiggsPt_GGHToBBIncl->Write();
    hGenHiggsLog2Pt_GGHToBBIncl->Write();
    hGenHiggsLog10Pt_GGHToBBIncl->Write();

    hGenHiggsPt_GGHToBBPt->Write();
    hGenHiggsLog2Pt_GGHToBBPt->Write();
    hGenHiggsLog10Pt_GGHToBBPt->Write();

    hGenHiggsPt_GGHToAATo4BIncl_0->Write();
    hGenHiggsLog2Pt_GGHToAATo4BIncl_0->Write();
    hGenHiggsLog10Pt_GGHToAATo4BIncl_0->Write();

    hGenHiggsPt_GGHToAATo4BPt_0->Write();
    hGenHiggsLog2Pt_GGHToAATo4BPt_0->Write();
    hGenHiggsLog10Pt_GGHToAATo4BPt_0->Write();

    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->Write();
    hGenHiggsLog2Pt_GGHToAATo4BIncl_Rewgted->Write();
    hGenHiggsLog10Pt_GGHToAATo4BIncl_Rewgted->Write();

    hGenHiggsPt_GGHToAATo4BPt_Rewgted->Write();
    hGenHiggsLog2Pt_GGHToAATo4BPt_Rewgted->Write();
    hGenHiggsLog10Pt_GGHToAATo4BPt_Rewgted->Write();


    fOut->Close();
    std::cout << "Wrote output file" << "\n";

}


void cal_GGFHiggsPtRewgts_fit() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0); 

    TFile *fIn=new TFile("hGenHiggsPt_GGHToBB_GGHToAATo4B_xsUnity_4GeVBins.root");

    TH1D *hGenHiggsPt_GGHTo2BIncl             = (TH1D*)fIn->Get("hGenHiggsPt_GGHToBBIncl");
    TH1D *hGenHiggsPt_GGHTo2BPt               = (TH1D*)fIn->Get("hGenHiggsPt_GGHToBBPt");
    TH1D *hGenHiggsPt_GGHToAATo4BIncl_0       = (TH1D*)fIn->Get("hGenHiggsPt_GGHToAATo4BIncl_0");
    TH1D *hGenHiggsPt_GGHToAATo4BIncl_Rewgted       = (TH1D*)fIn->Get("hGenHiggsPt_GGHToAATo4BIncl_Rewgted");
    TH1D *hGenHiggsPt_GGHToAATo4BPt_0         = (TH1D*)fIn->Get("hGenHiggsPt_GGHToAATo4BPt_0");
    TH1D *hGenHiggsPt_GGHToAATo4BPt_Rewgted   = (TH1D*)fIn->Get("hGenHiggsPt_GGHToAATo4BPt_Rewgted");

    double ptThreshold     = 300;
    int    bin_ptThreshold = hGenHiggsPt_GGHTo2BIncl->FindBin(ptThreshold);
    int    nTotalBins      = hGenHiggsPt_GGHTo2BIncl->GetNbinsX();

    double scale_GGHTo2BHighPt = hGenHiggsPt_GGHTo2BIncl->Integral(bin_ptThreshold, nTotalBins) / hGenHiggsPt_GGHTo2BPt->Integral(bin_ptThreshold, nTotalBins);
    std::cout << "hGenHiggsPt_GGHTo2BPt->Integral(bin_ptThreshold, nTotalBins): " << hGenHiggsPt_GGHTo2BPt->Integral(bin_ptThreshold, nTotalBins);
    std::cout << ", hGenHiggsPt_GGHTo2BIncl->Integral(bin_ptThreshold, nTotalBins): " << hGenHiggsPt_GGHTo2BIncl->Integral(bin_ptThreshold, nTotalBins) ;
    std::cout << "scale_GGHTo2BHighPt: " << scale_GGHTo2BHighPt << "\n";
    TH1D *hGenHiggsPt_GGHTo2BCombined = (TH1D*)hGenHiggsPt_GGHTo2BIncl->Clone("hGenHiggsPt_GGHTo2BCombined");
    for (int i=bin_ptThreshold; i<=nTotalBins; i++){
        double y_ = hGenHiggsPt_GGHTo2BPt->GetBinContent(i);
        double ey_ = hGenHiggsPt_GGHTo2BPt->GetBinError(i);
        hGenHiggsPt_GGHTo2BCombined->SetBinContent(i,  y_ * scale_GGHTo2BHighPt);
        hGenHiggsPt_GGHTo2BCombined->SetBinError(  i, ey_ * scale_GGHTo2BHighPt);
    }

    double scale_GGHToAATo4BHighPt = hGenHiggsPt_GGHToAATo4BIncl_0->Integral(bin_ptThreshold, nTotalBins) / hGenHiggsPt_GGHToAATo4BPt_0->Integral(bin_ptThreshold, nTotalBins);
    std::cout << "hGenHiggsPt_GGHToAATo4BPt_0->Integral(bin_ptThreshold, nTotalBins): " << hGenHiggsPt_GGHToAATo4BPt_0->Integral(bin_ptThreshold, nTotalBins);
    std::cout << ", hGenHiggsPt_GGHToAATo4BIncl_0->Integral(bin_ptThreshold, nTotalBins): " << hGenHiggsPt_GGHToAATo4BIncl_0->Integral(bin_ptThreshold, nTotalBins) ;
    std::cout << "scale_GGHToAATo4BHighPt: " << scale_GGHToAATo4BHighPt << "\n";
    TH1D *hGenHiggsPt_GGHToAATo4BCombined_0 = (TH1D*)hGenHiggsPt_GGHToAATo4BIncl_0->Clone("hGenHiggsPt_GGHToAATo4BCombined_0");
    for (int i=bin_ptThreshold; i<=nTotalBins; i++){
        double y_ = hGenHiggsPt_GGHToAATo4BPt_0->GetBinContent(i);
        double ey_ = hGenHiggsPt_GGHToAATo4BPt_0->GetBinError(i);
        hGenHiggsPt_GGHToAATo4BCombined_0->SetBinContent(i,  y_ * scale_GGHToAATo4BHighPt);
        hGenHiggsPt_GGHToAATo4BCombined_0->SetBinError(  i, ey_ * scale_GGHToAATo4BHighPt);
    }


    scale_GGHToAATo4BHighPt = hGenHiggsPt_GGHToAATo4BIncl_Rewgted->Integral(bin_ptThreshold, nTotalBins) / hGenHiggsPt_GGHToAATo4BPt_Rewgted->Integral(bin_ptThreshold, nTotalBins);
    std::cout << "hGenHiggsPt_GGHToAATo4BPt_Rewgted->Integral(bin_ptThreshold, nTotalBins): " << hGenHiggsPt_GGHToAATo4BPt_Rewgted->Integral(bin_ptThreshold, nTotalBins);
    std::cout << ", hGenHiggsPt_GGHToAATo4BIncl_Rewgted->Integral(bin_ptThreshold, nTotalBins): " << hGenHiggsPt_GGHToAATo4BIncl_Rewgted->Integral(bin_ptThreshold, nTotalBins) ;
    std::cout << "scale_GGHToAATo4BHighPt: " << scale_GGHToAATo4BHighPt << "\n";
    TH1D *hGenHiggsPt_GGHToAATo4BCombined_Rewgted = (TH1D*)hGenHiggsPt_GGHToAATo4BIncl_Rewgted->Clone("hGenHiggsPt_GGHToAATo4BCombined_Rewgted");
    for (int i=bin_ptThreshold; i<=nTotalBins; i++){
        double y_ = hGenHiggsPt_GGHToAATo4BPt_Rewgted->GetBinContent(i);
        double ey_ = hGenHiggsPt_GGHToAATo4BPt_Rewgted->GetBinError(i);
        hGenHiggsPt_GGHToAATo4BCombined_Rewgted->SetBinContent(i,  y_ * scale_GGHToAATo4BHighPt);
        hGenHiggsPt_GGHToAATo4BCombined_Rewgted->SetBinError(  i, ey_ * scale_GGHToAATo4BHighPt);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 600,500);
    c1->SetLogy();
    c1->SetGrid();    
    c1->cd();

    hGenHiggsPt_GGHTo2BIncl->SetMarkerStyle(20);
    hGenHiggsPt_GGHTo2BIncl->SetMarkerSize(0.5);
    hGenHiggsPt_GGHTo2BIncl->SetLineColor(kBlack);
    hGenHiggsPt_GGHTo2BIncl->SetMarkerColor(kBlack);    
    hGenHiggsPt_GGHTo2BIncl->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHTo2BIncl->GetYaxis()->SetTitle("No. of events");
    hGenHiggsPt_GGHTo2BIncl->Draw();

    hGenHiggsPt_GGHTo2BPt->SetMarkerStyle(20);
    hGenHiggsPt_GGHTo2BPt->SetMarkerSize(0.5);
    hGenHiggsPt_GGHTo2BPt->SetLineColor(kBlue);
    hGenHiggsPt_GGHTo2BPt->SetMarkerColor(kBlue);  
    hGenHiggsPt_GGHTo2BPt->Draw("same");   

    hGenHiggsPt_GGHTo2BCombined->SetMarkerStyle(20);
    hGenHiggsPt_GGHTo2BCombined->SetMarkerSize(0.5);
    hGenHiggsPt_GGHTo2BCombined->SetLineColor(kRed);
    hGenHiggsPt_GGHTo2BCombined->SetMarkerColor(kRed);  
    hGenHiggsPt_GGHTo2BCombined->Draw("same");  

    TLegend *leg1 = new TLegend(0.5,0.75,0.99,0.99);
    leg1->AddEntry(hGenHiggsPt_GGHTo2BIncl, "GGToHTo2B inclusive, #sigma=1 ", "lep");
    leg1->AddEntry(hGenHiggsPt_GGHTo2BPt, "GGToHTo2B pT(H)>200, #sigma=1", "lep");
    leg1->AddEntry(hGenHiggsPt_GGHTo2BCombined, "GGToHTo2B combined, #sigma=1", "lep");
    leg1->Draw();    

    TCanvas *c2 = new TCanvas("c2", "c2", 600,500);
    //c2->SetLogy();
    c2->SetGrid();    
    c2->cd();

    TH1D *hRatio_GenHiggsPt_GGHTo2B_InclOverComb = (TH1D*)hGenHiggsPt_GGHTo2BCombined->Clone("hRatio_GenHiggsPt_GGHTo2B_InclOverComb");
    hRatio_GenHiggsPt_GGHTo2B_InclOverComb->Divide(hGenHiggsPt_GGHTo2BIncl, hGenHiggsPt_GGHTo2BCombined);
    int nRebin = 4;
    hRatio_GenHiggsPt_GGHTo2B_InclOverComb->Rebin(nRebin); hRatio_GenHiggsPt_GGHTo2B_InclOverComb->Scale(double(1./nRebin));

    hRatio_GenHiggsPt_GGHTo2B_InclOverComb->SetMarkerStyle(20);
    hRatio_GenHiggsPt_GGHTo2B_InclOverComb->SetMarkerSize(0.5);
    hRatio_GenHiggsPt_GGHTo2B_InclOverComb->SetLineColor(kBlack);
    hRatio_GenHiggsPt_GGHTo2B_InclOverComb->SetMarkerColor(kBlack);    
    hRatio_GenHiggsPt_GGHTo2B_InclOverComb->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hRatio_GenHiggsPt_GGHTo2B_InclOverComb->GetYaxis()->SetTitle("GGHTo2B Incl/Combined"); 
    hRatio_GenHiggsPt_GGHTo2B_InclOverComb->GetYaxis()->SetRangeUser(0,2);
    hRatio_GenHiggsPt_GGHTo2B_InclOverComb->Draw(); 

    TLegend *leg2 = new TLegend(0.5,0.85,0.99,0.99);
    leg2->AddEntry(hRatio_GenHiggsPt_GGHTo2B_InclOverComb, "GGToHTo2B inclusive/combined", "lep");
    leg2->Draw();  



    TCanvas *c3 = new TCanvas("c3", "c3", 600,500);
    c3->SetLogy();
    c3->SetGrid();    
    c3->cd();

    hGenHiggsPt_GGHToAATo4BIncl_0->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BIncl_0->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BIncl_0->SetLineColor(kBlack);
    hGenHiggsPt_GGHToAATo4BIncl_0->SetMarkerColor(kBlack);    
    hGenHiggsPt_GGHToAATo4BIncl_0->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHToAATo4BIncl_0->GetYaxis()->SetTitle("No. of events");
    hGenHiggsPt_GGHToAATo4BIncl_0->Draw();

    hGenHiggsPt_GGHToAATo4BPt_0->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BPt_0->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BPt_0->SetLineColor(kBlue);
    hGenHiggsPt_GGHToAATo4BPt_0->SetMarkerColor(kBlue);  
    hGenHiggsPt_GGHToAATo4BPt_0->Draw("same");   

    hGenHiggsPt_GGHToAATo4BCombined_0->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BCombined_0->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BCombined_0->SetLineColor(kRed);
    hGenHiggsPt_GGHToAATo4BCombined_0->SetMarkerColor(kRed);  
    hGenHiggsPt_GGHToAATo4BCombined_0->Draw("same");  

    TLegend *leg3 = new TLegend(0.5,0.75,0.99,0.99);
    leg3->AddEntry(hGenHiggsPt_GGHToAATo4BIncl_0, "GGToHToAATo4B inclusive, #sigma=1 ", "lep");
    leg3->AddEntry(hGenHiggsPt_GGHToAATo4BPt_0, "GGToHToAATo4B pT(H)>150, #sigma=1", "lep");
    leg3->AddEntry(hGenHiggsPt_GGHToAATo4BCombined_0, "GGToHToAATo4B combined, #sigma=1", "lep");
    leg3->Draw();    

    TCanvas *c4 = new TCanvas("c4", "c4", 600,500);
    //c2->SetLogy();
    c4->SetGrid();    
    c4->cd();

    TH1D *hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb = (TH1D*)hGenHiggsPt_GGHToAATo4BCombined_0->Clone("hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb");
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->Divide(hGenHiggsPt_GGHToAATo4BIncl_0, hGenHiggsPt_GGHToAATo4BCombined_0);
    //int nRebin = 4;
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->Rebin(nRebin); hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->Scale(double(1./nRebin));

    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->SetMarkerStyle(20);
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->SetMarkerSize(0.5);
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->SetLineColor(kBlack);
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->SetMarkerColor(kBlack);    
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->GetYaxis()->SetTitle("GGHToAATo4B Incl/Combined"); 
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->GetYaxis()->SetRangeUser(0,2);
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->Draw(); 

    TLegend *leg4 = new TLegend(0.5,0.85,0.99,0.99);
    leg4->AddEntry(hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb, "GGToHTo2B inclusive/combined", "lep");
    leg4->Draw();  


    TCanvas *c5 = new TCanvas("c5", "c5", 600,500);
    c5->SetLogy();
    c5->SetGrid();    
    c5->cd();

    hGenHiggsPt_GGHTo2BCombined->SetMarkerStyle(20);
    hGenHiggsPt_GGHTo2BCombined->SetMarkerSize(0.5);
    hGenHiggsPt_GGHTo2BCombined->SetLineColor(kBlack);
    hGenHiggsPt_GGHTo2BCombined->SetMarkerColor(kBlack);    
    hGenHiggsPt_GGHTo2BCombined->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHTo2BCombined->GetYaxis()->SetTitle("No. of events");
    hGenHiggsPt_GGHTo2BCombined->Draw();

    hGenHiggsPt_GGHToAATo4BCombined_0->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BCombined_0->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BCombined_0->SetLineColor(kBlue);
    hGenHiggsPt_GGHToAATo4BCombined_0->SetMarkerColor(kBlue);  
    hGenHiggsPt_GGHToAATo4BCombined_0->Draw("same");    

    hGenHiggsPt_GGHToAATo4BCombined_Rewgted->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BCombined_Rewgted->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BCombined_Rewgted->SetLineColor(kRed);
    hGenHiggsPt_GGHToAATo4BCombined_Rewgted->SetMarkerColor(kRed);  
    hGenHiggsPt_GGHToAATo4BCombined_Rewgted->Draw("same");      

    TLegend *leg5 = new TLegend(0.5,0.75,0.99,0.99);
    leg5->AddEntry(hGenHiggsPt_GGHTo2BCombined, "GGToHTo2B combined, #sigma=1 ", "lep");
    leg5->AddEntry(hGenHiggsPt_GGHToAATo4BCombined_0, "GGToHToAATo4B combined, #sigma=1", "lep");
    leg5->AddEntry(hGenHiggsPt_GGHToAATo4BCombined_Rewgted, "GGToHToAATo4B combined - reweighted, #sigma=1", "lep");
    leg5->Draw();      

    TH1D *hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0 = (TH1D*)hGenHiggsPt_GGHTo2BCombined->Clone("hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0");
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0->Divide(hGenHiggsPt_GGHTo2BCombined, hGenHiggsPt_GGHToAATo4BCombined_0);
    nRebin = 4;
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0->Rebin(nRebin); hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0->Scale(double(1./nRebin));

    TH1D *hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_Rewgted = (TH1D*)hGenHiggsPt_GGHTo2BCombined->Clone("hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_Rewgted");
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_Rewgted->Divide(hGenHiggsPt_GGHTo2BCombined, hGenHiggsPt_GGHToAATo4BCombined_Rewgted);
    nRebin = 4;
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_Rewgted->Rebin(nRebin); hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_Rewgted->Scale(double(1./nRebin));

    TH1D *hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit = (TH1D*)hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0->Clone("hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit");

    TCanvas *c6 = new TCanvas("c6", "c6", 600,500);
    //c5->SetLogy();
    c6->SetGrid();    
    c6->cd();

    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit->SetMarkerStyle(20);
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit->SetMarkerSize(0.5);
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit->SetLineColor(kBlack);
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit->SetMarkerColor(kBlack);    
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit->GetYaxis()->SetTitle("GGHTo2B / GGHToAATo4B"); 
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit->GetYaxis()->SetRangeUser(0,1.4);
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit->Draw();     

    TF1 *fGenHiggsPtWgt_fitRegion = new TF1("fGenHiggsPtWgt_fitRegion", "[0] + [1]*x + [2]*pow(x, 2) + [3]*pow(x, 3)", 150,1100);
    fGenHiggsPtWgt_fitRegion->SetLineColor(kRed);
    fGenHiggsPtWgt_fitRegion->SetLineWidth(2);
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit->Fit("fGenHiggsPtWgt_fitRegion", "R");
    /*
    FCN=116.486 FROM MIGRAD    STATUS=CONVERGED     121 CALLS         122 TOTAL
                     EDM=2.29576e-17    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           1.45849e+00   9.84930e-03   5.63601e-06   6.30360e-09
   2  p1          -4.00668e-03   6.85713e-05   1.26226e-08  -3.37749e-05
   3  p2           4.02577e-06   1.41755e-07   2.03692e-11   2.96507e-02
   4  p3          -1.38804e-09   8.74525e-11   7.37155e-13  -2.63216e+00

   Chi2/NDF: 2.08011
    */

    TString sFitFunction_fitRegion = Form("%g + %g*x + %g*pow(x, 2) + %g*pow(x, 3)", \
        fGenHiggsPtWgt_fitRegion->GetParameter(0), 
        fGenHiggsPtWgt_fitRegion->GetParameter(1), 
        fGenHiggsPtWgt_fitRegion->GetParameter(2), 
        fGenHiggsPtWgt_fitRegion->GetParameter(3)
        );
    TString sFitFunction_full = Form("min(max(%s, 0.09), 1.02)", sFitFunction_fitRegion.Data()); // min(max(1.45849 + -0.00400668*x + 4.02577e-06*pow(x, 2) + -1.38804e-09*pow(x, 3), 0.09), 1.02)
    std::cout << "sFitFunction_full: " << sFitFunction_full << "\n";

    TF1 *fGenHiggsPtWgt_full = new TF1("fGenHiggsPtWgt_full", sFitFunction_full, 0, 2000);
    fGenHiggsPtWgt_full->SetLineColor(kBlue);
    fGenHiggsPtWgt_full->SetLineWidth(1);   
    fGenHiggsPtWgt_full->Draw("same");


    TLegend *leg6 = new TLegend(0.,0.9,1.0,1.25);
    //leg6->SetTextSize(0.4);
    leg6->AddEntry(hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0, "GGHTo2B / GGHToAATo4B",        "lep");
    leg6->AddEntry(fGenHiggsPtWgt_fitRegion,                     sFitFunction_fitRegion.Data(),  "l"  );
    leg6->AddEntry(fGenHiggsPtWgt_full,                          sFitFunction_full.Data(),       "l"  );
    leg6->Draw(); 



    TCanvas *c7 = new TCanvas("c7", "c7", 600,500);
    c7->SetLogy();
    c7->SetGrid();    
    c7->cd();

    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->SetLineColor(kBlack);
    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->SetMarkerColor(kBlack);    
    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->GetYaxis()->SetTitle("No. of events");
    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->Draw();

    hGenHiggsPt_GGHToAATo4BPt_Rewgted->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BPt_Rewgted->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BPt_Rewgted->SetLineColor(kBlue);
    hGenHiggsPt_GGHToAATo4BPt_Rewgted->SetMarkerColor(kBlue);  
    hGenHiggsPt_GGHToAATo4BPt_Rewgted->Draw("same");   

    hGenHiggsPt_GGHToAATo4BCombined_Rewgted->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BCombined_Rewgted->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BCombined_Rewgted->SetLineColor(kRed);
    hGenHiggsPt_GGHToAATo4BCombined_Rewgted->SetMarkerColor(kRed);  
    hGenHiggsPt_GGHToAATo4BCombined_Rewgted->Draw("same");  

    TLegend *leg7 = new TLegend(0.5,0.75,0.99,0.99);
    leg7->AddEntry(hGenHiggsPt_GGHToAATo4BIncl_Rewgted, "GGToHToAATo4B inclusive, #sigma=1 ", "lep");
    leg7->AddEntry(hGenHiggsPt_GGHToAATo4BPt_Rewgted, "GGToHToAATo4B pT(H)>150, #sigma=1", "lep");
    leg7->AddEntry(hGenHiggsPt_GGHToAATo4BCombined_Rewgted, "GGToHToAATo4B combined, #sigma=1", "lep");
    leg7->Draw();    

    TCanvas *c8 = new TCanvas("c8", "c8", 600,500);
    //c2->SetLogy();
    c8->SetGrid();    
    c8->cd();

    TH1D *hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted = (TH1D*)hGenHiggsPt_GGHToAATo4BCombined_Rewgted->Clone("hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted");
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted->Divide(hGenHiggsPt_GGHToAATo4BIncl_Rewgted, hGenHiggsPt_GGHToAATo4BCombined_Rewgted);
    //int nRebin = 4;
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted->Rebin(nRebin); hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted->Scale(double(1./nRebin));

    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted->SetMarkerStyle(20);
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted->SetMarkerSize(0.5);
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted->SetLineColor(kBlack);
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted->SetMarkerColor(kBlack);    
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted->GetYaxis()->SetTitle("GGHToAATo4B Incl/Combined"); 
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted->GetYaxis()->SetRangeUser(0,2);
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted->Draw(); 

    TLegend *leg8 = new TLegend(0.5,0.85,0.99,0.99);
    leg8->AddEntry(hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb_Rewgted, "GGToHTo2B inclusive/combined", "lep");
    leg8->Draw();  


    TCanvas *c9 = new TCanvas("c9", "c9", 600,500);
    //c2->SetLogy();
    c9->SetGrid();    
    c9->cd();

    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_Rewgted->Draw();

    TLegend *leg9 = new TLegend(0.5,0.85,0.99,0.99);
    leg9->AddEntry(hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_Rewgted, "GGToHTo2B / GGHToAATo4B-reweighted", "lep");
    leg9->Draw(); 



    TFile *fOut=new TFile("hGenHiggsPt_GGHToAATo4B_ptReweighting.root", "recreate");
    fOut->cd();

    hGenHiggsPt_GGHTo2BIncl->Write(); 
    hGenHiggsPt_GGHTo2BPt->Write();
    hGenHiggsPt_GGHTo2BCombined->Write();
    hRatio_GenHiggsPt_GGHTo2B_InclOverComb->Write();

    hGenHiggsPt_GGHToAATo4BIncl_0->Write();
    hGenHiggsPt_GGHToAATo4BPt_0->Write();
    hGenHiggsPt_GGHToAATo4BCombined_0->Write();
    hRatio_GenHiggsPt_GGHToAATo4B_InclOverComb->Write();

    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0->Write();
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_0_fit->Write();

    hGenHiggsPt_GGHToAATo4BIncl_Rewgted->Write();
    hGenHiggsPt_GGHToAATo4BPt_Rewgted->Write();
    hGenHiggsPt_GGHToAATo4BCombined_Rewgted->Write();
    hRatio_GenHiggsPt_GGHTo2B_over_GGHToAATo4B_Rewgted->Write();


    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
    c6->Write();
    c7->Write();
    c8->Write();
    c9->Write();

    fOut->Close();
    std::cout << "Wrote output file" << "\n";    
}