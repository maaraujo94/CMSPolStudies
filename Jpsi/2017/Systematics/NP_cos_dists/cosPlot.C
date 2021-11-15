void cosPlot()
{
  TFile *fin = new TFile("cosStore.root");

  TH2D *h_NP1 = (TH2D*)fin->Get("NP1H_ab");
  TH2D *h_NP2 = (TH2D*)fin->Get("NP2H_ab");
  TH2D *h_NPr1 = (TH2D*)fin->Get("rNP1H_ab");
  TH2D *h_NPr2 = (TH2D*)fin->Get("rNP2H_ab");

  // get the binning
  int nBinsX = h_NP1->GetNbinsX(), nBinsY = h_NP1->GetNbinsY();
  const double *yBins = h_NP1->GetYaxis()->GetXbins()->GetArray();
  double minX = h_NP1->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_NP1->GetXaxis()->GetBinUpEdge(nBinsX);

  // get the base data and MC 1d projections
  TH1D ***h_NP1d = new TH1D**[2];
  h_NP1d[0] = new TH1D*[nBinsY];
  h_NP1d[1] = new TH1D*[nBinsY];
  for(int i = 0; i< nBinsY; i++) {
    h_NP1d[0][i] = h_NP1->ProjectionX(Form("h_NP1_%d", i), i+1, i+1);
    h_NP1d[1][i] = h_NP2->ProjectionX(Form("h_NP2_%d", i), i+1, i+1);
  }
  // get the base data and MC 1d projections
  TH1D ***h_NPr1d = new TH1D**[2];
  h_NPr1d[0] = new TH1D*[nBinsY];
  h_NPr1d[1] = new TH1D*[nBinsY];
  for(int i = 0; i< nBinsY; i++) {
    h_NPr1d[0][i] = h_NPr1->ProjectionX(Form("h_NPr1_%d", i), i+1, i+1);
    h_NPr1d[1][i] = h_NPr2->ProjectionX(Form("h_NPr2_%d", i), i+1, i+1);
  }

  TCanvas *c = new TCanvas("", "", 900, 900);
  for(int i = 0; i < nBinsY; i++) {
    h_NP1d[0][i]->SetTitle(Form("2017 %.0f < p_{T} < %.0f GeV", yBins[i], yBins[i+1]));
    h_NP1d[0][i]->SetStats(0);
    h_NP1d[0][i]->SetLineColor(kBlack);
    h_NP1d[0][i]->SetMarkerColor(kBlack);
    h_NP1d[0][i]->SetMarkerStyle(20);
    h_NP1d[0][i]->SetMarkerSize(0.5);
    h_NP1d[0][i]->SetMinimum(0);
    h_NP1d[0][i]->SetMaximum(h_NP1d[0][i]->GetBinContent(1)*1.8);
    h_NP1d[0][i]->Draw("error");

    h_NP1d[1][i]->SetStats(0);
    h_NP1d[1][i]->Scale(h_NP1d[0][i]->Integral()/h_NP1d[1][i]->Integral());
    h_NP1d[1][i]->SetLineColor(kBlue);
    h_NP1d[1][i]->SetMarkerColor(kBlue);
    h_NP1d[1][i]->SetMarkerStyle(20);
    h_NP1d[1][i]->SetMarkerSize(0.5);
    h_NP1d[1][i]->Draw("error same");

    TLegend *leg = new TLegend(0.65, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(h_NP1d[0][i], "[100,200] #mum", "pl");
    leg->AddEntry(h_NP1d[1][i], "[300,500] #mum", "pl");
    leg->Draw();
    
    c->SaveAs(Form("plots/NPcos_%d.pdf", i));
    c->Clear();

    h_NPr1d[0][i]->SetTitle(Form("%.0f < p_{T} < %.0f GeV", yBins[i], yBins[i+1]));
    h_NPr1d[0][i]->SetStats(0);
    h_NPr1d[0][i]->SetLineColor(kBlack);
    h_NPr1d[0][i]->SetMarkerColor(kBlack);
    h_NPr1d[0][i]->SetMinimum(0);
    h_NPr1d[0][i]->SetMaximum(h_NPr1d[0][i]->GetBinContent(1)*1.3);
    h_NPr1d[0][i]->Draw("error");

    h_NPr1d[1][i]->SetStats(0);
    h_NPr1d[1][i]->Scale(h_NPr1d[0][i]->Integral(1,12)/h_NPr1d[1][i]->Integral(1,12));
    h_NPr1d[1][i]->SetLineColor(kBlue);
    h_NPr1d[1][i]->SetMarkerColor(kBlue);
    h_NPr1d[1][i]->Draw("error same");

    c->SaveAs(Form("plots/NPrcos_%d.pdf", i));
    c->Clear();
}
  c->Destructor();
  fin->Close();
}
