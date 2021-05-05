// starting from fine binned histogram (normal and |costh|)
// - get 2D coarse binned histo
// - save fine + coarse binned 2d histos

void saveCoarseBinned()
{
  double M_q = 3.686; // using J/psi mass
  
  // run over the ratio histogram (|costh|)
  TFile *infile = new TFile("files/ratioHist.root");
  TH2D *hist = new TH2D();
  infile->GetObject(Form("ratioH_ab_PS"), hist);
  hist->SetDirectory(0);
  infile->Close();
  
  // Get the 2D coarse binned histo
  // new bins
  int nBinsX = hist->GetNbinsX();
  double minX = hist->GetXaxis()->GetBinLowEdge(1);
  double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);
  
  // the new binning
  const int nBinsY = 11;
  int ptBins[nBinsY+1]   = { 0,  2,  4,  6,  8, 10, 15, 20, 25, 30, 32, 36};
  double binsY[nBinsY+1] = {25, 27, 29, 31, 33, 35, 40, 50, 60, 70, 80, 100};
  // Carlos suggested binning
  /*  const int nBinsY = 9;
  int ptBins[nBinsY+1]   = { 0,  3,  6,  9, 12, 15, 20, 25, 32, 36};
  double binsY[nBinsY+1] = {25, 28, 31, 34, 37, 40, 50, 60, 80, 100};*/
  
  // new histo
  TH2D *cHist = new TH2D(Form("cHist"), "PR/MC coarse bins", nBinsX, minX, maxX, nBinsY, binsY);
  
  // to go around limitations, get first 1D projections in the right binning
  TH1D *pHist[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    pHist[i] = hist->ProjectionX(Form("coarse_bin%d_1d", i+1), ptBins[i]+1, ptBins[i+1]);
    pHist[i]->SetTitle(Form("PR/MC c bin %d: [%.0f, %.0f] GeV", i+1, binsY[i], binsY[i+1]));
  }
  
  // then fill histo with 1D histo values and errors
  for(int iX = 0; iX < nBinsX; iX++) {
    for(int iY = 0; iY < nBinsY; iY++) {
      cHist->SetBinContent(iX+1, iY+1, pHist[iY]->GetBinContent(iX+1)/(ptBins[iY+1]-ptBins[iY]));
      cHist->SetBinError(iX+1, iY+1, pHist[iY]->GetBinError(iX+1)/(ptBins[iY+1]-ptBins[iY]));
    }
  }
  
  // store in outfile the fine-binned histo, the coarse-binned
  TFile* outfile_2 = new TFile("files/ratioHist.root", "update");
  hist->Write(0, TObject::kOverwrite);
  cHist->Write(0, TObject::kOverwrite);
  outfile_2->Close();
  
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.11);
  c->SetLogz(0);
  cHist->SetStats(0);
  cHist->SetMinimum(0);
  cHist->Draw("COLZ");
  c->SaveAs(Form("plots/ratio_2_PS_coarse.pdf"));
  c->Clear();
  c->Destructor(); 
}
