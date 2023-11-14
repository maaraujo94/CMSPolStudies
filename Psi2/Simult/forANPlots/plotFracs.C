// macro to plot the bkg fracs and prompt psi(2S) frac vs pT
void plotFracs()
{
  // get fit fBG and fNP
  TFile *fin1 = new TFile("../bkgFits/files/bkgFrac.root");
  TH1D *h_fSB = (TH1D*)fin1->Get("fbkg_unc"); 
  h_fSB->SetDirectory(0);
  fin1->Close();
  TFile *fin2 = new TFile("../PR_fit/files/NPFrac.root");
  TH2D *fNP2d = (TH2D*)fin2->Get("h_fnp");
  fNP2d->SetDirectory(0);
  fin2->Close();

  // get the binning
  int nBinsX = fNP2d->GetNbinsX(), nBinsY = fNP2d->GetNbinsY();
  const double *yBins = fNP2d->GetYaxis()->GetXbins()->GetArray();
  double minX = fNP2d->GetXaxis()->GetBinLowEdge(1);
  double maxX = fNP2d->GetXaxis()->GetBinUpEdge(nBinsX);
  
  TH1D *h_fNP = fNP2d->ProjectionY("fnp_1d", 1, 1);
  
  // define prompt psi(2S)
  TH1D *h_fJ = new TH1D("h_fJ", "h_fJ", nBinsY, yBins);
  for(int i = 0; i < nBinsY; i++) {
    h_fJ->SetBinContent(i+1, 1.-h_fSB->GetBinContent(i+1)-h_fNP->GetBinContent(i+1));
    h_fJ->SetBinError(i+1, sqrt( pow(h_fSB->GetBinError(i+1), 2) + pow(h_fNP->GetBinError(i+1), 2) ) );
  }

  // scale all histos before plotting
  h_fSB->Scale(100.);
  h_fNP->Scale(100.);
  h_fJ->Scale(100.);
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);

  TH1F *fr1 = c->DrawFrame(15, 0.0, 105, 100);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("fraction (%)");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  
  h_fNP->SetLineColor(kRed);
  h_fNP->SetMarkerColor(kRed);
  h_fNP->SetMarkerStyle(20);
  h_fNP->SetMarkerSize(.5);
  h_fNP->Draw("error same");

  h_fSB->SetLineColor(kGreen+1);
  h_fSB->SetMarkerColor(kGreen+1);
  h_fSB->SetMarkerStyle(20);
  h_fSB->SetMarkerSize(.5);
  h_fSB->Draw("error same");
  
  h_fJ->SetLineColor(kBlue);
  h_fJ->SetMarkerColor(kBlue);
  h_fJ->SetMarkerStyle(20);
  h_fJ->SetMarkerSize(.5);
  h_fJ->Draw("error same");

  TLatex lc;  
  lc.SetTextSize(0.03);
  lc.SetTextColor(kBlue);
  lc.DrawLatex(28, 62, "prompt #psi(2S)");
  lc.SetTextColor(kRed);
  lc.DrawLatex(30, 15, "non-prompt #psi(2S)");
  lc.SetTextColor(kGreen+1);
  lc.DrawLatex(70, 60, "continuum muon pairs");

  c->SaveAs("plots/f_comp.pdf");
  c->Clear();
  c->Destructor();
}
