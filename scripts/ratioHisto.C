

void ratioHisto() {
    TFile *f=new TFile("/eos/cms/store/user/ssawant/htoaa/analysis/20250717_DataMC/2018/gg0l/analyze_htoaa_stage1.root");

    TH1D *hQCD_bEnr_Fail = (TH1D*)f->Get("evt/QCD_bEnr/hnAk4JetsCentral_nonbTag_nonoverlaping_leadingFatJet_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom");
    TH1D *hQCD_BGen_Fail = (TH1D*)f->Get("evt/QCD_BGen/hnAk4JetsCentral_nonbTag_nonoverlaping_leadingFatJet_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom");
    TH1D *hQQCD_Incl_Fail = (TH1D*)f->Get("evt/QCD_Incl/hnAk4JetsCentral_nonbTag_nonoverlaping_leadingFatJet_gg0lIncl_Xto4bv2_SBplusSRWP40_Nom");
    TH1D *hQCD_Fail = (TH1D*)hQCD_bEnr_Fail->Clone("hQCD_Fail");
    hQCD_Fail->Add(hQCD_BGen_Fail);
    hQCD_Fail->Add(hQQCD_Incl_Fail);
    hQCD_Fail->SetLineColor(1);
    

    TH1D *hQCD_bEnr_Pass = (TH1D*)f->Get("evt/QCD_bEnr/hnAk4JetsCentral_nonbTag_nonoverlaping_leadingFatJet_gg0lIncl_Xto4bv2_SRWP40_Nom");
    TH1D *hQCD_BGen_Pass = (TH1D*)f->Get("evt/QCD_BGen/hnAk4JetsCentral_nonbTag_nonoverlaping_leadingFatJet_gg0lIncl_Xto4bv2_SRWP40_Nom");
    TH1D *hQQCD_Incl_Pass = (TH1D*)f->Get("evt/QCD_Incl/hnAk4JetsCentral_nonbTag_nonoverlaping_leadingFatJet_gg0lIncl_Xto4bv2_SRWP40_Nom");
    TH1D *hQCD_Pass = (TH1D*)hQCD_bEnr_Pass->Clone("hQCD_Pass");
    hQCD_Pass->Add(hQCD_BGen_Pass);
    hQCD_Pass->Add(hQQCD_Incl_Pass);
    hQCD_Pass->SetLineColor(2);

    TH1D *hQCD_Fail1 = (TH1D*)hQCD_Fail->Clone("hQCD_Fail1");
    hQCD_Fail1->Add(hQCD_Fail, hQCD_Pass, 1, -1);
    hQCD_Fail1->SetLineColor(4);


    TCanvas *c1 = new TCanvas("c1", "c1", 500,400);
    c1->cd();   
    c1->SetLogy();

    hQCD_Fail->Draw();
    hQCD_Pass->Draw("same");
    hQCD_Fail1->Draw("same");

    TH1D *hRatio = (TH1D*)hQCD_Pass->Clone("hRatio");
    hRatio->Divide(hQCD_Pass, hQCD_Fail1);

    TCanvas *c2 = new TCanvas("c2", "c2", 500,400);
    c2->cd(); 

    hRatio->Draw();

    


    //TH1D *hQCD_bEnrich_Fail = (TH1D*)f->Get("evt/QCD_bEnr/hnAk4JetsCentral_nonbTag_nonoverlaping_leadingFatJet_gg0lHi_Xto4bv2_SBplusSRWP40_Nom")
    //std::cout << "h " << hQCD_bEnrich_Fail->GetNbinsX() << std::endl;


}