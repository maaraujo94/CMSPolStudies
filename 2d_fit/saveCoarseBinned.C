// starting from fine binned histogram
// - get 2D coarse binned histo
// - save fine + coarse binned 1d + 2d histos

void saveCoarseBinned()
{
  TFile *infile = new TFile("files/ratioHist.root");
  TH2D *hist = new TH2D();
  infile->GetObject("ratioHist", hist);
  hist->SetDirectory(0);
  infile->Close();

  // Get the 2D coarse binned histo
  // new bins
  int nBinsX = hist->GetNbinsX();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);

  int nBinsY = 8;
  int ptBins[9]={0, 4, 8, 10, 13, 18, 23, 27, 29};
  double binsY[9] = {12, 16, 20, 24, 30, 40, 50, 60, 70};

  // new histo
  TH2D *cHist = new TH2D("cHist", "Data/MC coarse bins", nBinsX, minX, maxX, nBinsY, binsY);
  
  // to go around limitations, get first 1D projections in the right binning
  TH1D *pHist[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    pHist[i] = hist->ProjectionX(Form("coarse_bin%d_1d",i+1), ptBins[i]+1, ptBins[i+1]);
    pHist[i]->SetTitle(Form("Data/MC c bin %d: [%.0f, %.0f] GeV", i+1, binsY[i], binsY[i+1]));
  }

  // then fill histo with 1D histo values and errors
  for(int iX = 0; iX < nBinsX; iX++) {
    for(int iY = 0; iY < nBinsY; iY++) {
      cHist->SetBinContent(iX+1, iY+1, pHist[iY]->GetBinContent(iX+1)/(ptBins[iY+1]-ptBins[iY]));
      cHist->SetBinError(iX+1, iY+1, pHist[iY]->GetBinError(iX+1)/(ptBins[iY+1]-ptBins[iY]));
    }
  }

  // get fine binned 1D projections while we're at it
  int nPtBins = hist->GetNbinsY();
  TH1D *fHist[nPtBins];
  for(int i = 1; i <= nPtBins; i++) {
    fHist[i-1] = hist->ProjectionX(Form("fine_bin%d_1d",i), i, i+1);
    fHist[i-1]->SetTitle(Form("Data/MC f bin %d: [%.1f, %.1f] GeV", i, hist->GetYaxis()->GetBinLowEdge(i), hist->GetYaxis()->GetBinUpEdge(i)));
  }

  // store in outfile the fine-binned histo, the coarse-binned
  // and each of the coarse-binned 1D projections
  TFile* outfile = new TFile("files/store_hist.root", "recreate");
  hist->Write();
  cHist->Write();
  for(int i = 0; i < nBinsY; i++) pHist[i]->Write();
  for(int i = 0; i < nPtBins; i++) fHist[i]->Write();
  outfile->Close();
} 
