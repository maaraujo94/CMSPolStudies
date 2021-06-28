// starting from fine binned histogram get 2D coarse binned histo

void saveCoarseBinned()
{
  // get the ratio histograms
  TFile *infile = new TFile("files/bkgSubRes.root");
  TH2D **h_fine = new TH2D*[5];
  string lbl[] = {"Data", "NP", "SB", "PR", "J"};
  for(int i = 0; i < 5; i++) {
    infile->GetObject(Form("h_%s", lbl[i].c_str()), h_fine[i]);
    h_fine[i]->SetDirectory(0);
  }
  infile->Close();
  
  // get the binning
  int nBinsX = h_fine[0]->GetNbinsX(), nBinsY = h_fine[0]->GetNbinsY();
  const double *yBins = h_fine[0]->GetYaxis()->GetXbins()->GetArray();
  double minX = h_fine[0]->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_fine[0]->GetXaxis()->GetBinUpEdge(nBinsX);

  // the new binning
  const int nPtBins = 17;
  double ptBins[nPtBins+1];
  int yBins_c[nPtBins+1];
  for(int i = 0; i < 7; i++) ptBins[i] = 25 + 3.*i;
  for(int i = 0; i < 6; i++) ptBins[i+7] = 46 + 5.*i;
  for(int i = 0; i < 3; i++) ptBins[i+13] = 76 + 8.*i;
  for(int i = 0; i < 2; i++) ptBins[i+16] = 100 + 20.*i;

  // get matching between both binnings
  int ct = 0;
  for(int i = 0; i < nBinsY+1; i++) 
    if(yBins[i] == ptBins[ct]) {
      yBins_c[ct] = i+1;
      ct++;
    }

  // cycle over 5 histos to get coarse-binning versions
  TH2D **h_coarse = new TH2D*[5];
  for(int i_t = 0; i_t < 5; i_t++) {
    h_coarse[i_t] = new TH2D(Form("%s_c", h_fine[i_t]->GetName()), h_fine[i_t]->GetTitle(), nBinsX, minX, maxX, nPtBins, ptBins); 
    
    // to go around limitations, get first 1D projections in the right binning
    TH1D *pHist[nPtBins];
    for(int i = 0; i < nPtBins; i++) 
      pHist[i] = h_fine[i_t]->ProjectionX(Form("bin%d_%d", i, i_t), yBins_c[i], yBins_c[i+1]-1);
    
    // then fill histo with 1D histo values and errors
    for(int iX = 0; iX < nBinsX; iX++) {
      for(int iY = 0; iY < nPtBins; iY++) {
	h_coarse[i_t]->SetBinContent(iX+1, iY+1, pHist[iY]->GetBinContent(iX+1)/(yBins_c[iY+1]-yBins_c[iY]));
	h_coarse[i_t]->SetBinError(iX+1, iY+1, pHist[iY]->GetBinError(iX+1)/(yBins_c[iY+1]-yBins_c[iY]));
	cout << iX << " " << iY << " " << pHist[iY]->GetBinContent(iX+1) << " " << yBins_c[iY+1]-yBins_c[iY] << endl;
      }
    }
  }
  
  // store in outfile the coarse-binned
  TFile* outfile = new TFile("files/histoStore.root", "update");
  for(int i_t = 0; i_t < 5; i_t++) {
    h_coarse[i_t]->Write(0, TObject::kOverwrite);
  }
  outfile->Close();
  
  
}
