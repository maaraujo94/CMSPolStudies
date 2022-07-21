// macro to plot the bkg fracs and prompt psi(2S) frac over pT
void plotFracs()
{
  // get fit fBG and fNP
  TFile *fin1 = new TFile("../PR_fit/files/bkgFrac.root");
  TH2D *fSB2d = (TH2D*)fin1->Get("h_fbkg"); // fine-pT f_bkg
  fSB2d->SetDirectory(0);
  fin1->Close();
  TFile *fin2 = new TFile("../PR_fit/files/NPFrac.root");
  TH2D *fNP2d = (TH2D*)fin2->Get("h_fnp");
  TH2D *fNPc2d = (TH2D*)fin2->Get("h_fNPc");
  fNP2d->SetDirectory(0);
  fNPc2d->SetDirectory(0);
  fin2->Close();

  // get the binning
  int nBinsX = fNP2d->GetNbinsX(), nBinsY = fNP2d->GetNbinsY();
  const double *yBins = fNP2d->GetYaxis()->GetXbins()->GetArray();
  double minX = fNP2d->GetXaxis()->GetBinLowEdge(1);
  double maxX = fNP2d->GetXaxis()->GetBinUpEdge(nBinsX);
  
  TH1D *h_fSB = fSB2d->ProjectionY("fbkg_1d", 1, 1);
  TH1D *h_fNP = fNP2d->ProjectionY("fnp_1d", 1, 1);
  TH1D *h_fNPc = fNPc2d->ProjectionY("fnpc_1d", 1, 1);

  cout << "fNP = " << h_fNP->GetBinContent(1) << "; fNP_c = " << h_fNPc->GetBinContent(1) << endl;

  
  // define prompt psi(2S)
  TH1D *h_fJ = new TH1D("h_fJ", "h_fJ", nBinsY, yBins);
  TH1D *h_fJc = new TH1D("h_fJc", "h_fJc", nBinsY, yBins);
  for(int i = 0; i < nBinsY; i++) {
    h_fJ->SetBinContent(i+1, 1.-h_fSB->GetBinContent(i+1)-h_fNP->GetBinContent(i+1));
    h_fJ->SetBinError(i+1, sqrt( pow(h_fSB->GetBinError(i+1), 2) + pow(h_fNP->GetBinError(i+1), 2) ) );

    h_fJc->SetBinContent(i+1, 1.-h_fSB->GetBinContent(i+1)-h_fNPc->GetBinContent(i+1));
    h_fJc->SetBinError(i+1, sqrt( pow(h_fSB->GetBinError(i+1), 2) + pow(h_fNPc->GetBinError(i+1), 2) ) );
}

  // scale all histos before plotting
  h_fSB->Scale(100.);
  h_fNP->Scale(100.);
  h_fNPc->Scale(100.);
  h_fJ->Scale(100.);
  h_fJc->Scale(100.);
 
  cout << "fNP = " << h_fNP->GetBinContent(1) << "; fNP_c = " << h_fNPc->GetBinContent(1) << endl;
  
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

  h_fSB->SetLineColor(kGreen);
  h_fSB->SetMarkerColor(kGreen);
  h_fSB->SetMarkerStyle(20);
  h_fSB->SetMarkerSize(.5);
  h_fSB->Draw("error same");
  
  h_fJ->SetLineColor(kBlue);
  h_fJ->SetMarkerColor(kBlue);
  h_fJ->SetMarkerStyle(20);
  h_fJ->SetMarkerSize(.5);
  h_fJ->Draw("error same");

  TLegend *leg = new TLegend(0.65, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_fJ, "prompt #psi(2S)", "pl");
  leg->AddEntry(h_fNP, "NP", "pl");
  leg->AddEntry(h_fSB, "bkg", "pl");
  leg->Draw();

  c->SaveAs("plots/f_comp.pdf");
  c->Clear();

  // now with double-subtraction corrected
  TH1F *fr2 = c->DrawFrame(20, 0.0, 125, 100);
  fr2->SetXTitle("p_{T} (GeV)");
  fr2->SetYTitle("f (%)");
  fr2->GetYaxis()->SetTitleOffset(1.3);
  fr2->GetYaxis()->SetLabelOffset(0.01);
  fr2->SetTitle("Run 2 f comparison");

  h_fNPc->SetLineColor(kRed);
  h_fNPc->SetMarkerColor(kRed);
  h_fNPc->SetMarkerStyle(20);
  h_fNPc->SetMarkerSize(.5);
  h_fNPc->Draw("error same");

  h_fSB->SetLineColor(kGreen);
  h_fSB->SetMarkerColor(kGreen);
  h_fSB->SetMarkerStyle(20);
  h_fSB->SetMarkerSize(.5);
  h_fSB->Draw("error same");
  
  h_fJc->SetLineColor(kBlue);
  h_fJc->SetMarkerColor(kBlue);
  h_fJc->SetMarkerStyle(20);
  h_fJc->SetMarkerSize(.5);
  h_fJc->Draw("error same");

  leg->Draw();

  c->SaveAs("plots/f_comp_corr.pdf");
  c->Clear();
  c->Destructor();
}
