// macro to generate the sideband costh dists in the final binning

// macro for rounding to integers
int do_round(double val)
{
  int valR = (int)val;
  if (val-valR > 0.5) return valR+1;
  else return valR;
}

void genDistL()
{
  const double nGen = 1e7;
  
  // get binning from the stored data histos
  TFile *infile = new TFile("../../PR_fit/files/histoStore.root");
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
  TFile *infL = new TFile("../../PR_fit/files/store_fL.root");
  double fL = ((TGraphErrors*)infL->Get("g_fL"))->GetY()[0];
  infL->Close();
  // get LSB lambda2, lambda4
  TFile *inLSB = new TFile("files/LSBlin_fitres.root");
  double *L_l2 = ((TGraphErrors*)inLSB->Get("fit_l2"))->GetY();
  double *L_l4 = ((TGraphErrors*)inLSB->Get("ld4_lin"))->GetY();
  inLSB->Close();
  // get RSB lambda2, lambda4
  TFile *inRSB = new TFile("files/RSBlin_fitres.root");
  double *R_l2 = ((TGraphErrors*)inRSB->Get("fit_l2"))->GetY();
  double *R_l4 = ((TGraphErrors*)inRSB->Get("ld4_lin"))->GetY();
  inRSB->Close();

  // create the histograms - SB
  TH2D *h_LSB = new TH2D(Form("h_LSB"), "2017 LSB/MC", nBinsX, minX, maxX, nBinsY, yBins);
  TH2D *h_RSB = new TH2D(Form("h_RSB"), "2017 RSB/MC", nBinsX, minX, maxX, nBinsY, yBins);
  double intv = nGen * dX;
  
  TF1 *cthLSB = new TF1("cthLSB", "(1+[0]*x*x+[1]*pow(x,4))", minX, maxX);
  TF1 *cthRSB = new TF1("cthRSB", "(1+[0]*x*x+[1]*pow(x,4))", minX, maxX);
  TF1 *f_ld2 = new TF1("lambda_2", "[0]", yBins[0], yBins[nBinsY]);
  TF1 *f_ld4 = new TF1("lambda_4", "[0]*x + [1]", yBins[0], yBins[nBinsY]);

  // cycle over each pT bin
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    double pt = 0.5 * (yBins[i_pt+1] + yBins[i_pt]);

    // define the LSB costh dist and fill histo
    f_ld2->SetParameter(0, L_l2[0]);
    f_ld4->SetParameters(L_l4[0], L_l4[1]);
    cthLSB->SetParameters(f_ld2->Eval(pt), f_ld4->Eval(pt));
    double scF = intv/cthLSB->Integral(minX,maxX);
    for(int i_c = 0; i_c < nBinsX; i_c++) {
      double cost = minX + (i_c+0.5) * dX;
      h_LSB->SetBinContent(i_c+1, i_pt+1, scF*cthLSB->Eval(cost));
      h_LSB->SetBinError(i_c+1, i_pt+1, sqrt(h_LSB->GetBinContent(i_c+1, i_pt+1)));
    }

    // define the LSB costh dist and fill histo
    f_ld2->SetParameter(0, R_l2[0]);
    f_ld4->SetParameters(R_l4[0], R_l4[1]);
    cthRSB->SetParameters(f_ld2->Eval(pt), f_ld4->Eval(pt));
    scF = intv/cthRSB->Integral(minX,maxX);
    for(int i_c = 0; i_c < nBinsX; i_c++) {
      double cost = minX + (i_c+0.5) * dX;
      h_RSB->SetBinContent(i_c+1, i_pt+1, scF*cthRSB->Eval(cost));
      h_RSB->SetBinError(i_c+1, i_pt+1, sqrt(h_RSB->GetBinContent(i_c+1, i_pt+1)));
    }
  }

  // define the final bkg/MC dist as the weighted sum using fL
  TH2D *h_SB = new TH2D(Form("h_SB"), "2017 bkg/MC", nBinsX, minX, maxX, nBinsY, yBins);
  h_SB->Add(h_LSB, h_RSB, fL, 1.-fL);
  
  cout << "bkg/MC fully filled" << endl;
  
  TFile *fout = new TFile("files/bkgCosModelL.root", "recreate");
  h_SB->Write();
  fout->Close();

  // plot the 1d projection of the result (there's no pT dependence)
  TH1D *h_LSB1d = h_LSB->ProjectionX("h_LSB_1d", 1, 1);
  TH1D *h_RSB1d = h_RSB->ProjectionX("h_RSB_1d", 1, 1);
  TH1D *h_SB1d = h_SB->ProjectionX("h_SB_1d", 1, 1);
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  h_SB1d->SetStats(0);
  h_SB1d->SetMinimum(0);
  h_SB1d->SetMaximum(h_RSB1d->GetMaximum()*1.1);
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
