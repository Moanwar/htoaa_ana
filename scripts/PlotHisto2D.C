

void PlotHisto2D() {
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.12);
	
  // use large Times-Roman fonts
  gStyle->SetTextFont(132);
  gStyle->SetTextSize(0.05);
	
  gStyle->SetLabelFont(132,"x");
  gStyle->SetLabelFont(132,"y");
  gStyle->SetLabelFont(132,"z");
	
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gROOT->ForceStyle();
  

  std::string sipFile, sHistoName, sLegend, sXaxis, sYaxis, sSaveAs;
  int setLogZ;
  int rebinX, rebinY;

  
  sipFile = "/home/siddhesh/Work/CMS/htoaa/analysis/tmp9/analyze_htoaa_SUSY_GluGluH_01J_HToAATo4B_Pt150_mH-70_mA-12_wH-70_wA-70_TuneCP5_13TeV_madgraph_pythia8_1_0.root"; 
  /*
  sHistoName = "evt/SUSY_GluGluH_01J_HToAATo4B/hMass_GenA1_vs_GenA2_all_central";
  sLegend = "mH-70_mA-12_wH-70_wA-70";
  sXaxis = "m(Gen A1) [GeV]";
  sYaxis = "m(Gen A2) [GeV]";
  sSaveAs = "/home/siddhesh/Work/CMS/htoaa/analysis/tmp9/m_GenA1_vs_GenA2_mH-70_mA-12_wH-70_wA-70.png";
  */
  /*
  sHistoName = "evt/SUSY_GluGluH_01J_HToAATo4B/hMass_GenA1ToBBbar_vs_GenA2ToBBbar_all_central";
  sLegend = "mH-70_mA-12_wH-70_wA-70";
  sXaxis = "m(Gen A1ToBBbar) [GeV]";
  sYaxis = "m(Gen A2ToBBbar) [GeV]";
  sSaveAs = "/home/siddhesh/Work/CMS/htoaa/analysis/tmp9/m_GenA1_vs_GenA2_mH-70_mA-12_wH-70_wA-70.png"; 
  */
  /*
  sHistoName = "evt/SUSY_GluGluH_01J_HToAATo4B/hMass_GenAHeavy_vs_GenALight_all_central";
  sLegend = "mH-70_mA-12_wH-70_wA-70";
  sSaveAs = "/home/siddhesh/Work/CMS/htoaa/analysis/tmp9/m_GenAHeavy_vs_GenALight_mH-70_mA-12_wH-70_wA-70.png";
  */
  /*
  sHistoName = "evt/SUSY_GluGluH_01J_HToAATo4B/hMass_GenH_vs_GenAHeavy_all_central";
  sLegend = "mH-70_mA-12_wH-70_wA-70";
  sSaveAs = "/home/siddhesh/Work/CMS/htoaa/analysis/tmp9/m_GenH_vs_GenAHeavy_mH-70_mA-12_wH-70_wA-70.png";
  */
  /*
  sHistoName = "evt/SUSY_GluGluH_01J_HToAATo4B/hMass_GenH_vs_GenALight_all_central";
  sLegend = "mH-70_mA-12_wH-70_wA-70";
  sSaveAs = "/home/siddhesh/Work/CMS/htoaa/analysis/tmp9/m_GenH_vs_GenALight_mH-70_mA-12_wH-70_wA-70.png";
  */
  /*
  sHistoName = "evt/SUSY_GluGluH_01J_HToAATo4B/hMassGenH_vs_maxDRGenHGenB_all_central";
  sLegend = "mH-70_mA-12_wH-70_wA-70";
  sSaveAs = "/home/siddhesh/Work/CMS/htoaa/analysis/tmp9/m_GenH_vs_maxDRGenHGenB_mH-70_mA-12_wH-70_wA-70.png";
  */
  /*
  sHistoName = "evt/SUSY_GluGluH_01J_HToAATo4B/hMassGenAHeavy_vs_maxDRGenHGenB_all_central";
  sLegend = "mH-70_mA-12_wH-70_wA-70";
  sSaveAs = "/home/siddhesh/Work/CMS/htoaa/analysis/tmp9/m_GenAHeavy_vs_maxDRGenHGenB_mH-70_mA-12_wH-70_wA-70.png";
  */
  sHistoName = "evt/SUSY_GluGluH_01J_HToAATo4B/hMassGenALight_vs_maxDRGenHGenB_all_central";
  sLegend = "mH-70_mA-12_wH-70_wA-70";
  sSaveAs = "/home/siddhesh/Work/CMS/htoaa/analysis/tmp9/m_GenALight_vs_maxDRGenHGenB_mH-70_mA-12_wH-70_wA-70.png";




  rebinX = rebinY = 6;
  

  setLogZ = 0;

  
  TFile *tFIn = new TFile(sipFile.data());
  if ( ! tFIn->IsOpen()) {
    printf("File %s couldn't open \t\t\t *** ERROR **** \n",sipFile.data());
    return;
  }

  TH2D *h = nullptr;
  h = (TH2D*)tFIn->Get(sHistoName.data());
  if ( ! h) {
    printf("Couldn't fetch histogram %s from file %s  \t\t\t *** ERROR **** \n",sHistoName.data(),sipFile.data());
    return;
  }

  h->RebinX(rebinX);
  h->RebinY(rebinY);


  TAxis *axis;

  axis = h->GetXaxis();
  printf("\nX-axis binning:: %i [",axis->GetNbins());
  for (int i=1; i<=axis->GetNbins(); i++)
  {
    printf(" %g,",axis->GetBinLowEdge(i));
    if (i == axis->GetNbins()) printf(" %g] \n",axis->GetBinUpEdge(i));
  }
  
  axis = h->GetYaxis();
  /*
  printf("\nY-axis binning:: %i [",axis->GetNbins());
  for (int i=1; i<=axis->GetNbins(); i++)
  {
    printf(" %g,",axis->GetBinLowEdge(i));
    if (i == axis->GetNbins()) printf(" %g] \n",axis->GetBinUpEdge(i));
  }
  */

  
  TCanvas *c1 = new TCanvas("c1","c1", 600, 450);
  c1->cd();


  gPad->SetLogz(setLogZ);

  //h->GetXaxis()->SetTitle(sXaxis.data());
  //h->GetYaxis()->SetTitle(sYaxis.data());

  h->Draw("colz");
  
  
  TLegend *leg = new TLegend(0.2,0.89,0.8,0.99);
  leg->SetHeader(sLegend.data());

  leg->Draw();

  c1->SaveAs(Form("%s.png",sSaveAs.data()));
}
