// macro to compare f_Bg and f_NPBg for 2 cases
// fit w/o Gaussian excluding first 3 bins
// counting method excluding first 3 bins
void compfSB()
{ 
  // get f_Bg for fit
  TFile *fin_fitPR = new TFile("../bkgFits/files/bkgFrac.root");
  TH1D *fsb_fitPR = (TH1D*)fin_fitPR->Get("fbkg_unc");
  fsb_fitPR->SetDirectory(0);
  fin_fitPR->Close();
  // same for counting w/0 first bins
  TFile *fin_cPR = new TFile("/home/mariana/Documents/2023_PhD_work/CERN/07/0731_fBgComp/Psi2/bkgFits/files/bkgFrac.root");
  TH1D *fsb_cPR = (TH1D*)fin_cPR->Get("fbkg_unc");
  fsb_cPR->SetDirectory(0);
  fin_cPR->Close();

  // get f_NPBg for fit
  TFile *fin_fitNP = new TFile("../bkgFits/files/bkgFrac_NP.root");
  TH1D *fsb_fitNP = (TH1D*)fin_fitNP->Get("fbkg_unc");
  fsb_fitNP->SetDirectory(0);
  fin_fitNP->Close();
  // same for counting w/0 first bins
  TFile *fin_cNP = new TFile("/home/mariana/Documents/2023_PhD_work/CERN/07/0731_fBgComp/Psi2/bkgFits/files/bkgFrac_NP.root");
  TH1D *fsb_cNP = (TH1D*)fin_cNP->Get("fbkg_unc");
  fsb_cNP->SetDirectory(0);
  fin_cNP->Close();

  // now the plotting
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);
  c->SetRightMargin(0.015);
  c->SetTopMargin(0.015);
  c->SetLogy(0);
  
  // scaled and comparing to f_bkg from fit
  TH1F *fc = c->DrawFrame(15, 0., 105, 80);
  fc->SetXTitle("p_{T} (GeV)");
  fc->SetYTitle("f_{Bg} (%)");
  fc->GetYaxis()->SetTitleOffset(1.45);
  fc->GetYaxis()->SetLabelOffset(0.01);

  fsb_cPR->SetStats(0);
  fsb_cPR->Scale(100.);
  fsb_cPR->SetLineColor(kBlue);
  fsb_cPR->SetMarkerColor(kBlue);
  fsb_cPR->SetMarkerStyle(20);
  fsb_cPR->SetMarkerSize(.5);
  fsb_cPR->Draw("error same");

  fsb_fitPR->Scale(100.);
  fsb_fitPR->SetLineColor(kBlack);
  fsb_fitPR->SetMarkerColor(kBlack);
  fsb_fitPR->Draw("error same");

  TLegend *leg = new TLegend(0.7, 0.85, 1., 0.95);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0);
  leg->AddEntry(fsb_cPR, "yield ratio", "pl");
  leg->AddEntry(fsb_fitPR, "mass fit", "pl");
  leg->Draw();

  leg->Draw();

  c->SaveAs("plots/fitcomp_fBg.pdf");
  c->Clear();  

  // scaled and comparing to f_NPbkg from fit
  TH1F *fc2 = c->DrawFrame(15, 0., 105, 40);
  fc2->SetXTitle("p_{T} (GeV)");
  fc2->SetYTitle("f_{NPBg} (%)");
  fc2->GetYaxis()->SetTitleOffset(1.45);
  fc2->GetYaxis()->SetLabelOffset(0.01);

  fsb_cNP->SetStats(0);
  fsb_cNP->Scale(100.);
  fsb_cNP->SetLineColor(kBlue);
  fsb_cNP->SetMarkerColor(kBlue);
  fsb_cNP->SetMarkerStyle(20);
  fsb_cNP->SetMarkerSize(.5);
  fsb_cNP->Draw("error same");

  fsb_fitNP->Scale(100.);
  fsb_fitNP->SetLineColor(kBlack);
  fsb_fitNP->SetMarkerColor(kBlack);
  fsb_fitNP->Draw("error same");

  leg->Draw();

  c->SaveAs("plots/fitcomp_fNPBg.pdf");
  c->Clear();  
  c->Destructor();

}
