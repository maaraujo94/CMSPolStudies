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
  TFile *infile = new TFile("files/histoStore.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("dataH_ab"), hist);
  hist->SetDirectory(0);
  infile->Close();

  int nBinsX = hist->GetNbinsX(), nBinsY = hist->GetNbinsY();
  const double *yBins = hist->GetYaxis()->GetXbins()->GetArray();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);

  // get fit parameters from storage
  TFile *infL = new TFile("../bkgFits/files/store_fL.root");
  double fL = ((TGraphErrors*)infL->Get("g_fL"))->GetY()[0];
  infL->Close();
  // get LSB lambda2, lambda4
  TFile *inLSB = new TFile("files/LSBlin_fitres.root");
  double *L_l2 = ((TGraphErrors*)inLSB->Get("ld2_lin"))->GetY();
  double *L_l4 = ((TGraphErrors*)inLSB->Get("ld4_lin"))->GetY();
  inLSB->Close();
  // get RSB lambda2, lambda4
  TFile *inRSB = new TFile("files/RSBlin_fitres.root");
  double *R_l2 = ((TGraphErrors*)inRSB->Get("ld2_lin"))->GetY();
  double *R_l4 = ((TGraphErrors*)inRSB->Get("ld4_lin"))->GetY();
  inRSB->Close();

  // create the histograms - define everything
  TF1 *cthLSB = new TF1("cthLSB", "(1+[0]*x*x+[1]*pow(x,4))", minX, maxX);
  TF1 *cthRSB = new TF1("cthRSB", "(1+[0]*x*x+[1]*pow(x,4))", minX, maxX);
  TF1 *f_ld2 = new TF1("lambda_2", "[0]*x + [1]", yBins[0], yBins[nBinsY]);
  TF1 *f_ld4 = new TF1("lambda_4", "[0]*x + [1]", yBins[0], yBins[nBinsY]);
  TH1D **h_SB = new TH1D*[nBinsY];
  int genL = do_round(nGen*fL);
  int genR = nGen - genL;

  // cycle over each pT bin
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    cout << Form("SB cos#theta (%.0f < p_{T} < %.0f GeV)", yBins[i_pt], yBins[i_pt+1]) << endl;
    
    h_SB[i_pt] = new TH1D(Form("h_SB_%d", i_pt), Form("SB cos#theta (%.0f < p_{T} < %.0f GeV)", yBins[i_pt], yBins[i_pt+1]), nBinsX, minX, maxX);

    double pt = 0.5 * (yBins[i_pt+1] + yBins[i_pt]);

    // define the LSB costh dist and fill histo
    f_ld2->SetParameters(L_l2[0], L_l2[1]);
    f_ld4->SetParameters(L_l4[0], L_l4[1]);
    cthLSB->SetParameters(f_ld2->Eval(pt), f_ld4->Eval(pt));
    for(int i = 0; i < genL; i++) {
      h_SB[i_pt]->Fill(cthLSB->GetRandom());
    }

    // define the RSB costh dist and fill histo
    f_ld2->SetParameters(R_l2[0], R_l2[1]);
    f_ld4->SetParameters(R_l4[0], R_l4[1]);
    cthRSB->SetParameters(f_ld2->Eval(pt), f_ld4->Eval(pt));
    for(int i = 0; i < genR; i++) {
      h_SB[i_pt]->Fill(cthRSB->GetRandom());
    }
  }
  
  cout << "all SB histos filled" << endl;
  
  TFile *fout = new TFile("files/bkgCosModelL.root", "recreate");
  for(int i = 0; i < nBinsY; i++) {
    h_SB[i]->Write();
  }
  fout->Close();


}
