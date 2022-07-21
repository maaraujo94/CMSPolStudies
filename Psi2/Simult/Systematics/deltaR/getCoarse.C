void getCoarse()
{
  TH2D **h_data = new TH2D*[3];
  
  // get baseline evts
  TFile *finB = new TFile("../../PR_fit/files/bkgSubRes.root");
  h_data[0] = (TH2D*)finB->Get("h_J");
  h_data[0]->SetName("h_JB");
  h_data[0]->SetDirectory(0);
  finB->Close();

  // get tight cut events
  TFile *finT = new TFile("../../../Simult_dR1/PR_fit/files/bkgSubRes.root");
  h_data[2] = (TH2D*)finT->Get("h_J");
  h_data[2]->SetName("h_JT");
  h_data[2]->SetDirectory(0);
  finT->Close();

  // get loose cut events
  TFile *finL = new TFile("../../../Simult_dR2/PR_fit/files/bkgSubRes.root");
  h_data[1] = (TH2D*)finL->Get("h_J");
  h_data[1]->SetName("h_JL");
  h_data[1]->SetDirectory(0);
  finL->Close();

  //get the binning
  int nBinsX = h_data[0]->GetNbinsX();
  double minX = h_data[0]->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_data[0]->GetXaxis()->GetBinUpEdge(nBinsX);
  
  // the new binning
  const int nBinsY = 3;
  int ptBins[nBinsY+1]   = { 1,  4, 6, 8};
  double binsY[nBinsY+1] = {25, 46, 66, 120};

  // new histos
  string lbl[3] = {"coarse bins baseline",
		   "coarse bins #DeltaR>0.15",
		   "coarse bins #DeltaR>0.17"};
  string dl[3] = {"B", "L", "T"};
  TH2D **cHist = new TH2D*[3];
  for(int i = 0; i< 3; i++) {
    cHist[i] = new TH2D(Form("cHist%s", dl[i].c_str()), lbl[i].c_str(), nBinsX, minX, maxX, nBinsY, binsY);
    
    // get first 1D projections in the right binning
    TH1D *pHist[nBinsY];
    for(int j = 0; j < nBinsY; j++) {
      pHist[j] = h_data[i]->ProjectionX(Form("cbin%s_%d", dl[i].c_str(), j), ptBins[j], ptBins[j+1]-1);
      pHist[j]->SetTitle(Form("PR/MC c bin %d: [%.0f, %.0f] GeV", j, binsY[j], binsY[j+1]));
    }
  
    // then fill histo with 1D histo values and errors
    for(int iX = 0; iX < nBinsX; iX++) {
      for(int iY = 0; iY < nBinsY; iY++) {
	cHist[i]->SetBinContent(iX+1, iY+1, pHist[iY]->GetBinContent(iX+1)/(ptBins[iY+1]-ptBins[iY]));
	cHist[i]->SetBinError(iX+1, iY+1, pHist[iY]->GetBinError(iX+1)/(ptBins[iY+1]-ptBins[iY]));
      }
    }
  }

  TFile *fout = new TFile("files/chistStore.root", "recreate");
  for(int i = 0; i < 3; i++) {
    h_data[i]->Write();
    cHist[i]->Write();
  }
  fout->Close();
}
