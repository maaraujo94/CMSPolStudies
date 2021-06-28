// macro to generate the bkg costh dists in the fine binning

// macro for rounding to integers
int do_round(double val)
{
  int valR = (int)val;
  if (val-valR > 0.5) return valR+1;
  else return valR;
}

void genDist()
{
  const double nGen = 1e6;
  
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
  const double fL = 0.533;
  // get LSB lambda2, lambda4
  TFile *inLSB = new TFile("files/LSB_fitres.root");
  double L_l2 = ((TGraphErrors*)inLSB->Get("graph_l2"))->GetY()[0];
  double L_l4 = ((TGraphErrors*)inLSB->Get("graph_l4"))->GetY()[0];
  inLSB->Close();
  // get RSB lambda2, lambda4
  TFile *inRSB = new TFile("files/RSB_fitres.root");
  double R_l2 = ((TGraphErrors*)inRSB->Get("graph_l2"))->GetY()[0];
  double R_l4 = ((TGraphErrors*)inRSB->Get("graph_l4"))->GetY()[0];
  inRSB->Close();
  // get NP lambda
  TFile *inNP = new TFile("files/NP_fitres.root");
  double NP_l2 = ((TGraphErrors*)inNP->Get("graph_l2"))->GetY()[0];
  inNP->Close();

  // create the histograms - NP
  TF1 *cthNP = new TF1("cthNP", "(1+[0]*x*x)", minX, maxX);
  cthNP->SetParameter(0, NP_l2);
  TH1D **h_NP = new TH1D*[nBinsY];
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    h_NP[i_pt] = new TH1D(Form("h_NP_%d", i_pt), Form("NP cos#theta (%.1f < p_{T} < %.1f GeV)", yBins[i_pt], yBins[i_pt+1]), nBinsX, minX, maxX);
    for(int i = 0; i < nGen; i++) {
      h_NP[i_pt]->Fill(cthNP->GetRandom());
    }
  }
  cout << "all NP histos filled" << endl;

  // create the histograms - SB
  TF1 *cthLSB = new TF1("cthLSB", "(1+[0]*x*x+[1]*pow(x,4))", minX, maxX);
  TF1 *cthRSB = new TF1("cthRSB", "(1+[0]*x*x+[1]*pow(x,4))", minX, maxX);
  cthLSB->SetParameters(L_l2, L_l4);
  cthRSB->SetParameters(R_l2, R_l4);
  int genL = do_round(nGen*fL);
  int genR = nGen - genL;
  TH1D **h_SB = new TH1D*[nBinsY];
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    h_SB[i_pt] = new TH1D(Form("h_SB_%d", i_pt), Form("SB cos#theta (%.1f < p_{T} < %.1f GeV)", yBins[i_pt], yBins[i_pt+1]), nBinsX, minX, maxX);
    for(int i = 0; i < genL; i++) {
      h_SB[i_pt]->Fill(cthLSB->GetRandom());
    }
    for(int i = 0; i < genR; i++) {
      h_SB[i_pt]->Fill(cthRSB->GetRandom());
    }
  }

  cout << "all SB histos filled" << endl;
  
  TFile *fout = new TFile("files/bkgCosModel.root", "recreate");
  for(int i = 0; i < nBinsY; i++) {
    h_NP[i]->Write();
    h_SB[i]->Write();
  }
  fout->Close();
}
