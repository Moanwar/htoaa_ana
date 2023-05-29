

void checkSFs() {
  
  TFile *f0 = new TFile("LHE_HT_2018.root");
  TH1D *hQCDIncl=(TH1D*)f0->Get("evt/QCDIncl_PSWeight/hGenLHE_HT_all_central");
  hQCDIncl->SetMarkerColor(1);
  hQCDIncl->SetLineColor(1);

  //hQCDIncl->GetYaxis()->SetRangeUser(1e2, 1e13);

  TFile *f1 = new TFile("HTSamplesStitchSF_2018.root");
  TH1D *hQCDbGen=(TH1D*)f1->Get("QCD_bGen/hGenLHE_HT_all_central");
  hQCDbGen->SetMarkerColor(2);
  hQCDbGen->SetLineColor(2);

  TH1D *hQCDbGenCorr=(TH1D*)f1->Get("QCD_bGen/hGenLHE_HT_all_central_wHTSamplesStitch_HTCorrSFFromHighToLowHT");
  hQCDbGenCorr->SetMarkerColor(3);
  hQCDbGenCorr->SetLineColor(3);

  hQCDIncl->Draw();
  hQCDbGen->Draw("same");
  hQCDbGenCorr->Draw("same");

  TH1D *hQCDbGenCorr1=(TH1D*)f1->Get("QCD_bGen/hGenLHE_HT_all_central_wHTSamplesStitch1_HTCorrSFFromHighToLowHT");
  hQCDbGenCorr1->SetMarkerColor(4);
  hQCDbGenCorr1->SetLineColor(4);
  hQCDbGenCorr1->Draw("same");

  TH1D *hQCDbGenCorr2=(TH1D*)f1->Get("QCD_bGen/hGenLHE_HT_all_central_wHTSamplesStitch1_HTCorrSFFromLowToHighHT");
  hQCDbGenCorr2->SetMarkerColor(5);
  hQCDbGenCorr2->SetLineColor(5);
  hQCDbGenCorr2->Draw("same");
  
  //return;
    
  double kHT = 2100;
  
  TH1D *hQCDbGen_1=(TH1D*)hQCDbGen->Clone("hQCDbGen_1");
  hQCDbGen_1->SetMarkerColor(6);
  hQCDbGen_1->SetLineColor(6);
  hQCDbGen_1->Scale( hQCDIncl->GetBinContent(kHT) / hQCDbGen->GetBinContent(kHT) );
  
  TH1D *hQCDbGenCorr_1=(TH1D*)hQCDbGenCorr1->Clone("hQCDbGenCorr_1");
  hQCDbGenCorr_1->SetMarkerColor(7);
  hQCDbGenCorr_1->SetLineColor(7);
  hQCDbGenCorr_1->Scale( hQCDIncl->GetBinContent(kHT) / hQCDbGenCorr1->GetBinContent(kHT) );


  hQCDbGen_1->Draw("same");
  hQCDbGenCorr_1->Draw("same");

  TLegend *leg = new TLegend(0.5, 0.5, 0.9, 0.9);
  leg->AddEntry(hQCDIncl, "QCD Incl", "l");
  leg->AddEntry(hQCDbGen, "QCD bGen", "l");
  leg->AddEntry(hQCDbGenCorr, "QCD bGen HTStitch - step0", "l");
  leg->AddEntry(hQCDbGenCorr1, "QCD bGen HTStitch - step1", "l");
  
  leg->AddEntry(hQCDbGen_1, "QCD bGen - scaled", "l");
  leg->AddEntry(hQCDbGenCorr_1, "QCD bGen HTStitch - step1 - scaled", "l");
  leg->Draw();
  

  
    
}
