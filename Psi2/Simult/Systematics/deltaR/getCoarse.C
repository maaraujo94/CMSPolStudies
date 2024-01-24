// macro to get coarse binning for deltaR checks
void getCoarse()
{
  TH2D **h_data = new TH2D*[3];
  TH2D **h_dataNP = new TH2D*[3];

  // get baseline evts
  TFile *finB = new TFile("../../PR_fit/files/bkgSubRes.root");
  h_data[0] = (TH2D*)finB->Get("h_J");
  h_data[0]->SetName("h_JB");
  h_data[0]->SetDirectory(0);
  finB->Close();
  TFile *finB_NP = new TFile("../../NP_fit/files/bkgSubRes.root");
  h_dataNP[0] = (TH2D*)finB_NP->Get("h_NPc");
  h_dataNP[0]->SetName("h_NPB");
  h_dataNP[0]->SetDirectory(0);
  finB_NP->Close();

  // get tight cut events
  TFile *finT = new TFile("../../../Simult_dR1/PR_fit/files/bkgSubRes.root");
  h_data[2] = (TH2D*)finT->Get("h_J");
  h_data[2]->SetName("h_JT");
  h_data[2]->SetDirectory(0);
  finT->Close();
  TFile *finT_NP = new TFile("../../../Simult_dR1/NP_fit/files/bkgSubRes.root");
  h_dataNP[2] = (TH2D*)finT_NP->Get("h_NPc");
  h_dataNP[2]->SetName("h_NPT");
  h_dataNP[2]->SetDirectory(0);
  finT_NP->Close();

  // get loose cut events
  TFile *finL = new TFile("../../../Simult_dR2/PR_fit/files/bkgSubRes.root");
  h_data[1] = (TH2D*)finL->Get("h_J");
  h_data[1]->SetName("h_JL");
  h_data[1]->SetDirectory(0);
  finL->Close();
  TFile *finL_NP = new TFile("../../../Simult_dR2/NP_fit/files/bkgSubRes.root");
  h_dataNP[1] = (TH2D*)finL_NP->Get("h_NPc");
  h_dataNP[1]->SetName("h_NPL");
  h_dataNP[1]->SetDirectory(0);
  finL_NP->Close();

  //get the binning
  int nBinsX = h_data[0]->GetNbinsX();
  double minX = h_data[0]->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_data[0]->GetXaxis()->GetBinUpEdge(nBinsX);
  
  // the new binning
  const int nBinsY = 3;
  int binMin[nBinsY]     = { 1,  9, 15};
  int binMax[nBinsY]     = { 8, 14, 19};
  double binsY[nBinsY+1] = {25, 45, 70, 120};

  // new histos
  string lbl[3] = {"coarse bins baseline",
		   "coarse bins #DeltaR>0.15",
		   "coarse bins #DeltaR>0.17"};
  string dl[3] = {"B", "L", "T"};
  TH2D **cHist = new TH2D*[3];
  TH2D **cHistNP = new TH2D*[3];
  for(int i = 0; i< 3; i++) {
    cHist[i] = new TH2D(Form("cHist%s", dl[i].c_str()), lbl[i].c_str(), nBinsX, minX, maxX, nBinsY, binsY);
    cHistNP[i] = new TH2D(Form("cHistNP%s", dl[i].c_str()), lbl[i].c_str(), nBinsX, minX, maxX, nBinsY, binsY);
    
    // get first 1D projections in the right binning
    TH1D *pHist[nBinsY];
    TH1D *pHistNP[nBinsY];
    for(int j = 0; j < nBinsY; j++) {
      pHist[j] = h_data[i]->ProjectionX(Form("cbin%s_%d", dl[i].c_str(), j), binMin[j], binMax[j]);
      pHist[j]->SetTitle(Form("PR/MC c bin %d: [%.0f, %.0f] GeV", j, binsY[j], binsY[j+1]));

      pHistNP[j] = h_dataNP[i]->ProjectionX(Form("cbinNP%s_%d", dl[i].c_str(), j), binMin[j], binMax[j]);
      pHistNP[j]->SetTitle(Form("NP/MC c bin %d: [%.0f, %.0f] GeV", j, binsY[j], binsY[j+1]));
}
  
    // then fill histo with 1D histo values and errors
    for(int iX = 0; iX < nBinsX; iX++) {
      for(int iY = 0; iY < nBinsY; iY++) {
	cHist[i]->SetBinContent(iX+1, iY+1, pHist[iY]->GetBinContent(iX+1)/(1+binMax[iY]-binMin[iY]));
	cHist[i]->SetBinError(iX+1, iY+1, pHist[iY]->GetBinError(iX+1)/(1+binMax[iY]-binMin[iY]));

	cHistNP[i]->SetBinContent(iX+1, iY+1, pHistNP[iY]->GetBinContent(iX+1)/(1+binMax[iY]-binMin[iY]));
	cHistNP[i]->SetBinError(iX+1, iY+1, pHistNP[iY]->GetBinError(iX+1)/(1+binMax[iY]-binMin[iY]));
      }
    }
  }

  TFile *fout = new TFile("files/chistStore.root", "recreate");
  for(int i = 0; i < 3; i++) {
    cHist[i]->Write();
    cHistNP[i]->Write();
  }
  fout->Close();
}
