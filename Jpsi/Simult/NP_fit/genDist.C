// macro to generate the sideband costh dists in the final binning, with unc

// the SB/MC model function
double fit_model(double x, double l2, double l4) {
  return 1 + l2 * pow(x,2) + l4 * pow(x,4);
}

void genDist()
{
  // get binning from the stored data histos
  TFile *infile = new TFile("../PR_fit/files/histoStore.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("dataH_ab"), hist);
  hist->SetDirectory(0);
  infile->Close();

  int nBinsX = hist->GetNbinsX(), nBinsY = hist->GetNbinsY();
  const double *yBins = hist->GetYaxis()->GetXbins()->GetArray();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/nBinsX;

  // get fit parameters from storage
  TFile *infL = new TFile("files/store_fL.root");
  double fL = ((TGraphErrors*)infL->Get("g_fL"))->GetY()[0];
  infL->Close();
  // get LSB lambda2, lambda4
  TFile *inLSB = new TFile("files/LSB2d_fitres.root");
  int nPtBins = ((TGraphErrors*)inLSB->Get("fit_l2"))->GetN();
  double L_l2 = ((TGraphErrors*)inLSB->Get("fit_l2"))->GetY()[0];
  double L_l4 = ((TGraphErrors*)inLSB->Get("fit_l4"))->GetY()[0];
  double L_el2 = ((TGraphErrors*)inLSB->Get("fit_l2"))->GetEY()[0];
  double L_el4 = ((TGraphErrors*)inLSB->Get("fit_l4"))->GetEY()[0];
  double L_cov = ((TFitResult*)inLSB->Get("fitres"))->GetCovarianceMatrix()[nPtBins][nPtBins+1];
  inLSB->Close();
  // get RSB lambda2, lambda4
  TFile *inRSB = new TFile("files/RSB2d_fitres.root");
  double R_l2 = ((TGraphErrors*)inRSB->Get("fit_l2"))->GetY()[0];
  double R_l4 = ((TGraphErrors*)inRSB->Get("fit_l4"))->GetY()[0];
  double R_el2 = ((TGraphErrors*)inRSB->Get("fit_l2"))->GetEY()[0];
  double R_el4 = ((TGraphErrors*)inRSB->Get("fit_l4"))->GetEY()[0];
  double R_cov = ((TFitResult*)inRSB->Get("fitres"))->GetCovarianceMatrix()[nPtBins][nPtBins+1];
  inRSB->Close();

  // create the histograms - SB
  TH2D *h_LSB = new TH2D(Form("h_LSB"), "Run 2 LSB/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_RSB = new TH2D(Form("h_RSB"), "Run 2 RSB/MC", nBinsX, minX, maxX, nBinsY, yBins);

  // determine the uncertainty band at each pT and cost value
  double ln = 10000;
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    for(int i_c = 0; i_c < nBinsX; i_c++) {
      double cost = minX + (i_c+0.5) * dX;
      // LSB calculations
      h_LSB->SetBinContent(i_c+1, i_pt+1, fit_model(cost, L_l2, L_l4));
      double d2 = (fit_model(cost, L_l2 + L_el2/ln, L_l4) - fit_model(cost, L_l2, L_l4))/(L_el2/ln);
      double d4 = (fit_model(cost, L_l2, L_l4 + L_el4/ln) - fit_model(cost, L_l2, L_l4))/(L_el4/ln);
      h_LSB->SetBinError(i_c+1, i_pt+1, sqrt( pow(d2 * L_el2, 2) + pow(d4 * L_el4, 2) + 2*d2*d4*L_cov ));

      // RSB calculations
      h_RSB->SetBinContent(i_c+1, i_pt+1, fit_model(cost, R_l2, R_l4));
      d2 = (fit_model(cost, R_l2 + R_el2/ln, R_l4) - fit_model(cost, R_l2, R_l4))/(R_el2/ln);
      d4 = (fit_model(cost, R_l2, R_l4 + R_el4/ln) - fit_model(cost, R_l2, R_l4))/(R_el4/ln);
      h_RSB->SetBinError(i_c+1, i_pt+1, sqrt( pow(d2 * R_el2, 2) + pow(d4 * R_el4, 2) + 2*d2*d4*R_cov ));
    }
  }

  // define the final bkg/MC dist as the weighted sum using fL
  TH2D *h_SB = new TH2D(Form("h_SB"), "Run 2 bkg/MC", nBinsX, minX, maxX, nBinsY, yBins);
  h_SB->Add(h_LSB, h_RSB, fL, 1.-fL);

  cout << "all SB histos filled" << endl;
  
  TFile *fout = new TFile("files/bkgCosModel.root", "recreate");
  h_SB->Write();
  fout->Close();

  // plot the 1d projection of the result (there's no pT dependence)
  TH1D *h_LSB1d = h_LSB->ProjectionX("h_LSB_1d", 1, 1);
  TH1D *h_RSB1d = h_RSB->ProjectionX("h_RSB_1d", 1, 1);
  TH1D *h_SB1d = h_SB->ProjectionX("h_SB_1d", 1, 1);
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  h_SB1d->SetStats(0);
  h_SB1d->SetMinimum(0);
  h_SB1d->SetMaximum(1.7);//h_LSB1d->GetMaximum()*1.1);
  h_SB1d->SetLineColor(kGreen+1);
  h_SB1d->GetXaxis()->SetTitle("|cos#theta|");
  h_SB1d->Draw();
  
  h_LSB1d->SetLineColor(kBlack);
  h_LSB1d->Draw("same");
  h_RSB1d->SetLineColor(kBlue);
  h_RSB1d->Draw("same");

  TLegend *leg = new TLegend(0.7, 0.2, 0.9, 0.45);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_SB1d, "SR", "l");
  leg->AddEntry(h_LSB1d, "LSB", "l");
  leg->AddEntry(h_RSB1d, "RSB", "l");
  leg->Draw();

  c->SaveAs("plots/SB_base.pdf");
  c->Clear();
  c->Destructor();


}
