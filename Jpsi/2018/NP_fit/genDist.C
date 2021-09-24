// macro to generate the sideband costh dists in the final binning

// macro for rounding to integers
int do_round(double val)
{
  int valR = (int)val;
  if (val-valR > 0.5) return valR+1;
  else return valR;
}

void genDist()
{
  const double nGen = 1e7;
  
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

  // get fit parameters from storage
  TFile *infL = new TFile("files/store_fL.root");
  double fL = ((TGraphErrors*)infL->Get("g_fL"))->GetY()[0];
  infL->Close();
  // get LSB lambda2, lambda4
  TFile *inLSB = new TFile("files/LSB2d_fitres.root");
  double L_l2 = ((TGraphErrors*)inLSB->Get("fit_l2"))->GetY()[0];
  double L_l4 = ((TGraphErrors*)inLSB->Get("fit_l4"))->GetY()[0];
  inLSB->Close();
  // get RSB lambda2, lambda4
  TFile *inRSB = new TFile("files/RSB2d_fitres.root");
  double R_l2 = ((TGraphErrors*)inRSB->Get("fit_l2"))->GetY()[0];
  double R_l4 = ((TGraphErrors*)inRSB->Get("fit_l4"))->GetY()[0];
  inRSB->Close();

  // create the histograms - SB
  TF1 *cthLSB = new TF1("cthLSB", "(1+[0]*x*x+[1]*pow(x,4))", minX, maxX);
  TF1 *cthRSB = new TF1("cthRSB", "(1+[0]*x*x+[1]*pow(x,4))", minX, maxX);
  cthLSB->SetParameters(L_l2, L_l4);
  cthRSB->SetParameters(R_l2, R_l4);
  int genL = do_round(nGen*fL);
  int genR = nGen - genL;
    TH1D *h_SB_base = new TH1D(Form("h_SB_base"), "2018 Interpolated sideband cos#theta", nBinsX, minX, maxX);
  for(int i = 0; i < genL; i++) {
    h_SB_base->Fill(cthLSB->GetRandom());
  }
  for(int i = 0; i < genR; i++) {
    h_SB_base->Fill(cthRSB->GetRandom());
  }

  TH1D **h_SB = new TH1D*[nBinsY];
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
     h_SB[i_pt] = (TH1D*)h_SB_base->Clone(Form("h_SB_%d", i_pt));
    
    h_SB[i_pt]->SetTitle(Form("SB cos#theta (%.0f < p_{T} < %.0f GeV)", yBins[i_pt], yBins[i_pt+1]));
  }

  cout << "all SB histos filled" << endl;
  
  TFile *fout = new TFile("files/bkgCosModel.root", "recreate");
  for(int i = 0; i < nBinsY; i++) {
    h_SB[i]->Write();
  }
  fout->Close();

  double norm = ( fL*cthLSB->Integral(0,1) + (1.-fL)*cthRSB->Integral(0,1) ) / h_SB_base->Integral("width");
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  h_SB_base->SetStats(0);
  h_SB_base->Scale(norm);
  h_SB_base->SetMinimum(0);
  h_SB_base->SetMaximum(cthRSB->Eval(0.95)*1.2);
  h_SB_base->GetXaxis()->SetTitle("|cos#theta|");
  h_SB_base->Draw();
  
  cthLSB->SetLineColor(kBlack);
  cthLSB->Draw("same");
  cthRSB->SetLineColor(kBlue);
  cthRSB->Draw("same");

  TLegend *leg = new TLegend(0.7, 0.2, 0.9, 0.45);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_SB_base, "SR", "l");
  leg->AddEntry(cthLSB, "LSB", "l");
  leg->AddEntry(cthRSB, "RSB", "l");
  leg->Draw();

  c->SaveAs("plots/SB_base.pdf");
  c->Clear();
  c->Destructor();

}
