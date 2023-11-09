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


void cal_GGFHiggsPtRewgts_check1() {
    TFile *fGGHToAATo4BIncl = new TFile("/eos/cms/store/user/ssawant/NanoAOD/SUSY_GluGluH_01J_HToAATo4B_M-12_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/hadded/hadded_NanoAOD.root");
    TFile *fGGHToAATo4BPt   = new TFile("/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-12_TuneCP5_13TeV_madgraph_pythia8/r1/PNet_v1.root");

    double xsGGHToAATo4BIncl = 48.61;
    double xsGGHToAATo4BPt   = 48.61 * 0.057; 
    double xsGGHToAATo4BPt_v1 = 48.61 * 0.057 / 1.38; 

    TTree *trGGHToAATo4BIncl = (TTree*)fGGHToAATo4BIncl->Get("Events");
    long nEventsGenWgtPos_GGHToAATo4BIncl = trGGHToAATo4BIncl->Draw("run", "(genWeight >= 0)", "goff");
    long nEventsGenWgtNeg_GGHToAATo4BIncl = trGGHToAATo4BIncl->Draw("run", "(genWeight < 0)", "goff");
    long nEventsTot_GGHToAATo4BIncl       = nEventsGenWgtPos_GGHToAATo4BIncl - nEventsGenWgtNeg_GGHToAATo4BIncl;
    std::cout << "GGHToAATo4BIncl nEvents: " << trGGHToAATo4BIncl->GetEntries() << "\n";
    std::cout << "nEventsGenWgtPos_GGHToAATo4BIncl: " << nEventsGenWgtPos_GGHToAATo4BIncl;
    std::cout << ", nEventsGenWgtNeg_GGHToAATo4BIncl: " << nEventsGenWgtNeg_GGHToAATo4BIncl;
    std::cout << ", nEventsTot_GGHToAATo4BIncl: " << nEventsTot_GGHToAATo4BIncl << "\n";

    TTree *trGGHToAATo4BPt = (TTree*)fGGHToAATo4BPt->Get("Events");
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
    TH1D *hGenHiggsPt_GGHToAATo4BIncl = new TH1D("hGenHiggsPt_GGHToAATo4BIncl", "hGenHiggsPt_GGHToAATo4BIncl", 150, 0, 1500); 
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
    TH1D *hGenHiggsPt_GGHToAATo4BPt   = new TH1D("hGenHiggsPt_GGHToAATo4BPt", "hGenHiggsPt_GGHToAATo4BPt", 150, 0, 1500);
    hGenHiggsPt_GGHToAATo4BPt->SetMarkerStyle(20);
    hGenHiggsPt_GGHToAATo4BPt->SetMarkerSize(0.5);
    hGenHiggsPt_GGHToAATo4BPt->SetLineColor(kBlue);
    hGenHiggsPt_GGHToAATo4BPt->SetMarkerColor(kBlue);
    hGenHiggsPt_GGHToAATo4BPt->GetXaxis()->SetTitle("GEN Higgs pT [GeV]");
    hGenHiggsPt_GGHToAATo4BPt->GetYaxis()->SetTitle("No. of events");    
    trGGHToAATo4BPt->Draw("GenPart_pt >> hGenHiggsPt_GGHToAATo4BPt", sCut_GGHToAATo4BPt.Data(), "same");

    
    TString sCut_GGHToAATo4BPt_v1 = Form("(GenPart_pdgId == 25 && GenPart_status == 62) * (1 - 2*(genWeight < 0)) * (%f/%ld)", xsGGHToAATo4BPt_v1, nEventsTot_GGHToAATo4BPt);
    std::cout << "sCut_GGHToAATo4BPt_v1: " <<  sCut_GGHToAATo4BPt_v1.Data() << "\n";
    TH1D *hGenHiggsPt_GGHToAATo4BPt_v1   = new TH1D("hGenHiggsPt_GGHToAATo4BPt_v1", "hGenHiggsPt_GGHToAATo4BPt_v1", 150, 0, 1500);
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

    c1->SaveAs("hGenHPt_GGHToAATo4B_overlay.root");
    c2->SaveAs("hGenHPt_GGHToAATo4B_ratio.root");
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