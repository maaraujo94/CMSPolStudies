// macro to plot the bkg fracs and prompt psi(2S) frac over pT
void plotFracs_NP()
{
  // get fit fBG
  TFile *fin1 = new TFile("../bkgFits/files/bkgFrac_NP.root");
  TH1D *h_fSB = (TH1D*)fin1->Get("fbkg_unc"); 
  h_fSB->SetDirectory(0);
  fin1->Close();

  // this one just for the binning
  TFile *fin2 = new TFile("../PR_fit/files/NPFrac.root");
  TH2D *fNP2d = (TH2D*)fin2->Get("h_fnp");
  fNP2d->SetDirectory(0);
  fin2->Close();

  // get the binning
  int nBinsX = fNP2d->GetNbinsX(), nBinsY = fNP2d->GetNbinsY();
  const double *yBins = fNP2d->GetYaxis()->GetXbins()->GetArray();
  double minX = fNP2d->GetXaxis()->GetBinLowEdge(1);
  double maxX = fNP2d->GetXaxis()->GetBinUpEdge(nBinsX);
  
  // define non-prompt psi(2S)
  TH1D *h_fJ = new TH1D("h_fJ", "h_fJ", nBinsY, yBins);
  for(int i = 0; i < nBinsY; i++) {
    h_fJ->SetBinContent(i+1, 1.-h_fSB->GetBinContent(i+1));
    h_fJ->SetBinError(i+1, h_fSB->GetBinError(i+1));
  }

  // scale all histos before plotting
  h_fSB->Scale(100.);
  h_fJ->Scale(100.);
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.015);

  TH1F *fr1 = c->DrawFrame(15, 0.0, 105, 100);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("fraction (%)");
  fr1->GetYaxis()->SetTitleOffset(1.4);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  
  h_fJ->SetLineColor(kRed);
  h_fJ->SetMarkerColor(kRed);
  h_fJ->SetMarkerStyle(24);
  h_fJ->SetMarkerSize(.75);
  h_fJ->Draw("error same");

  h_fSB->SetLineColor(kGreen+1);
  h_fSB->SetMarkerColor(kGreen+1);
  h_fSB->SetMarkerStyle(20);
  h_fSB->SetMarkerSize(.75);
  h_fSB->Draw("error same");

  TLatex lc;
  lc.SetTextSize(0.03);
  lc.SetTextColor(kRed);
  lc.DrawLatex(28, 80, "non-prompt #psi(2S)");
  lc.SetTextColor(kGreen+1);
  lc.DrawLatex(50, 30, "non-prompt continuum muon pairs");

  c->SaveAs("plots/fNP_comp.pdf");
  c->Clear();
  
  c->Destructor();
}
