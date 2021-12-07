void effPlot()
{
  TFile *fin = new TFile("eff_mm_trailing_v5.root");

  TH1D *h_NP1 = (TH1D*)fin->Get("eff_pt");

  // get the binning
  int nBinsX = h_NP1->GetNbinsX();
  const double *xBins = h_NP1->GetXaxis()->GetXbins()->GetArray();

  TCanvas *c = new TCanvas("", "", 900, 900);

  TH1F *fl = c->DrawFrame(0, 0, 40, 1);
  fl->SetXTitle("p_{T} (GeV)");
  fl->SetYTitle("#varepsilon(#mu)");
  fl->GetYaxis()->SetTitleOffset(1.3);
  fl->GetYaxis()->SetLabelOffset(0.01);
  fl->SetTitle("Single #mu efficiency");

  h_NP1->SetStats(0);
  h_NP1->SetLineColor(kBlack);
  h_NP1->SetMarkerColor(kBlack);
  h_NP1->SetMarkerStyle(20);
  h_NP1->SetMarkerSize(0.5);

  h_NP1->Draw("error same");

  c->SaveAs(Form("mu_eff.pdf"));
  c->Clear();

  TH1F *fl2 = c->DrawFrame(0, 0.9, 40, 1);
  fl2->SetXTitle("p_{T} (GeV)");
  fl2->SetYTitle("#varepsilon(#mu)");
  fl2->GetYaxis()->SetTitleOffset(1.3);
  fl2->GetYaxis()->SetLabelOffset(0.01);
  fl2->SetTitle("Single #mu efficiency");

  h_NP1->Draw("error same");

  c->SaveAs(Form("mu_eff_zoom.pdf"));
  c->Clear();
  
  c->Destructor();
  fin->Close();
}
