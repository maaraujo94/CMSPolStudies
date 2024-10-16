// macro to generate the sideband costh dists in the final binning, with unc

#import "../ptbins.C"

void genDist()
{
  // get binning from the stored data histos
  TH2D **h_LSB = new TH2D*[nPtBins];
  TH2D **h_RSB = new TH2D*[nPtBins];
  TFile *infile = new TFile("files/histoStore.root");
  for(int i = 0; i < nPtBins; i++) {
    h_LSB[i] = (TH2D*)infile->Get(Form("PRLH_%d", i));
    h_RSB[i] = (TH2D*)infile->Get(Form("PRRH_%d", i));
    h_LSB[i]->SetDirectory(0);
    h_RSB[i]->SetDirectory(0);
  }
  infile->Close();

  int nBinsX = h_LSB[0]->GetNbinsX(), nBinsY = h_LSB[0]->GetNbinsY();
  double minX = h_LSB[0]->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_LSB[0]->GetXaxis()->GetBinUpEdge(nBinsX);
  double dX = (maxX-minX)/nBinsX;
  double minY = h_LSB[0]->GetYaxis()->GetBinLowEdge(1);
  double maxY = h_LSB[0]->GetYaxis()->GetBinUpEdge(nBinsY);
  double dY = (maxY-minY)/nBinsY;

  // get fit parameters from storage
  TFile *infL = new TFile("../../Simult/PR_fit/files/store_fL.root");
  double *fL = ((TGraphErrors*)infL->Get("g_fL"))->GetY();
  infL->Close();

  // scale to nr events per pt bin
  double n_s[nPtBins];
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    n_s[i_pt] = fL[i_pt] * h_LSB[i_pt]->GetEntries() + (1.-fL[i_pt])*h_RSB[i_pt]->GetEntries();
    h_LSB[i_pt]->Scale(1./h_LSB[i_pt]->GetEntries());
    h_RSB[i_pt]->Scale(1./h_RSB[i_pt]->GetEntries());
  }

  // get the sideband histos by summing with proportion fL
  TH2D **h_SB = new TH2D*[nPtBins];
  for(int i_pt = 0; i_pt < nPtBins; i_pt++) {
    h_SB[i_pt] = new TH2D(Form("h_SB_%d", i_pt), Form("Bg (|cos#theta|,#phi) (%.1f < p_{T} < %.1f GeV)", ptBins[i_pt], ptBins[i_pt+1]), nBinsX, minX, maxX, nBinsY, minY, maxY);
    
    h_SB[i_pt]->Sumw2();
    h_SB[i_pt]->Add(h_LSB[i_pt], h_RSB[i_pt], fL[i_pt], 1.-fL[i_pt]);

    h_SB[i_pt]->Scale(n_s[i_pt]);
    h_LSB[i_pt]->Scale(n_s[i_pt]);
    h_RSB[i_pt]->Scale(n_s[i_pt]);
  }

  cout << "all SB histos filled" << endl;
  
  TFile *fout = new TFile("files/bkgCosModel.root", "recreate");
  for(int i = 0; i < nPtBins; i++) {
    h_SB[i]->Write();
  }
  fout->Close();
}
