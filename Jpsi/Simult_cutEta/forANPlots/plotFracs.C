// macro to plot the bkg fracs and prompt J/psi frac over pT
void plotFracs()
{
  // get fit fBG and fNP
  TFile *fin1 = new TFile("../PR_fit/files/bkgFrac.root");
  TH1D *fSB_b = (TH1D*)fin1->Get("fbkg_unc"); // wide-pT f_bkg
  fSB_b->SetDirectory(0);
  TH2D *fSB2d = (TH2D*)fin1->Get("h_fbkg"); // fine-pT f_bkg
  fSB2d->SetDirectory(0);
  fin1->Close();
  TFile *fin2 = new TFile("../PR_fit/files/NPFrac.root");
  TH2D *fNP2d = (TH2D*)fin2->Get("h_fNPc");
  fNP2d->SetDirectory(0);
  fin2->Close();

  // get the binning
  int nBinsX = fNP2d->GetNbinsX(), nBinsY = fNP2d->GetNbinsY();
  const double *yBins = fNP2d->GetYaxis()->GetXbins()->GetArray();
  double minX = fNP2d->GetXaxis()->GetBinLowEdge(1);
  double maxX = fNP2d->GetXaxis()->GetBinUpEdge(nBinsX);
  
  TH1D *h_fSB = fSB2d->ProjectionY("fbkg_1d", 1, 1);
  TH1D *h_fNP = fNP2d->ProjectionY("fnp_1d", 1, 1);
  
  // define prompt J/psi
  TH1D *h_fJ = new TH1D("h_fJ", "h_fJ", nBinsY, yBins);
  for(int i = 0; i < nBinsY; i++) {
    h_fJ->SetBinContent(i+1, 1.-h_fSB->GetBinContent(i+1)-h_fNP->GetBinContent(i+1));
    h_fJ->SetBinError(i+1, sqrt( pow(h_fSB->GetBinError(i+1), 2) + pow(h_fNP->GetBinError(i+1), 2) ) );
  }

  // scale all histos before plotting
  h_fSB->Scale(100.);
  h_fNP->Scale(100.);
  h_fJ->Scale(100.);
  fSB_b->Scale(100.);
  
  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1F *fr1 = c->DrawFrame(20, 0.0, 125, 100);
  fr1->SetXTitle("p_{T} (GeV)");
  fr1->SetYTitle("f (%)");
  fr1->GetYaxis()->SetTitleOffset(1.3);
  fr1->GetYaxis()->SetLabelOffset(0.01);
  fr1->SetTitle("Run 2 f comparison");

  h_fNP->SetLineColor(kRed);
  h_fNP->SetMarkerColor(kRed);
  h_fNP->SetMarkerStyle(20);
  h_fNP->SetMarkerSize(.5);
  h_fNP->Draw("error same");

  h_fSB->SetLineColor(kBlack);
  h_fSB->SetLineStyle(kDashed);
  h_fSB->SetMarkerColor(kBlack);
  h_fSB->SetMarkerStyle(20);
  h_fSB->SetMarkerSize(.5);
  h_fSB->Draw("hist same c");

  fSB_b->SetMarkerColor(kGreen);
  fSB_b->SetLineColor(kGreen);
  fSB_b->SetMarkerStyle(20);
  fSB_b->SetMarkerSize(.5);
  fSB_b->Draw("error same");
  
  h_fJ->SetLineColor(kBlue);
  h_fJ->SetMarkerColor(kBlue);
  h_fJ->SetMarkerStyle(20);
  h_fJ->SetMarkerSize(.5);
  h_fJ->Draw("error same");

  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_fJ, "prompt J/#psi", "pl");
  leg->AddEntry(h_fNP, "NP", "pl");
  leg->AddEntry(fSB_b, "bkg", "pl");
  leg->Draw();

  c->SaveAs("plots/f_comp.pdf");
  c->Clear();
  c->Destructor();
}
