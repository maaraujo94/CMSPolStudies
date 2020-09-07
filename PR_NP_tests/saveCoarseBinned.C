// starting from fine binned histogram (normal and |costh|)
// - get 2D coarse binned histo
// - save fine + coarse binned 1d + 2d histos

void saveCoarseBinned()
{
  double M_q = 3.097; // using J/psi mass
  
  string files[2] = {"", "_ab"};

  // run over the two ratio histograms (normal and |costh|)
  for(int i_file = 0; i_file < 2; i_file++) {
    TFile *infile = new TFile("files/ratioHist.root");
    TH2D *hist = new TH2D();
    infile->GetObject(Form("ratioHist%s", files[i_file].c_str()), hist);
    hist->SetDirectory(0);
    infile->Close();

    // Get the 2D coarse binned histo
    // new bins
    int nBinsX = hist->GetNbinsX();
    double minX = hist->GetXaxis()->GetBinLowEdge(1);
    double maxX = hist->GetXaxis()->GetBinUpEdge(nBinsX);

    int nBinsY = 8;
    int ptBins[9]={0, 4, 8, 10, 13, 18, 23, 27, 29};
    double binsY[9] = {12./M_q, 16./M_q, 20./M_q, 24./M_q, 30./M_q, 40./M_q, 50./M_q, 60./M_q, 70./M_q};

    // new histo
    TH2D *cHist = new TH2D(Form("cHist%s", files[i_file].c_str()), "PR/NP coarse bins", nBinsX, minX, maxX, nBinsY, binsY);
  
    // to go around limitations, get first 1D projections in the right binning
    TH1D *pHist[nBinsY];
    for(int i = 0; i < nBinsY; i++) {
      pHist[i] = hist->ProjectionX(Form("coarse_bin%d_1d%s", i+1, files[i_file].c_str()), ptBins[i]+1, ptBins[i+1], "e");
      pHist[i]->SetTitle(Form("PR/NP c bin %d: [%.0f, %.0f] GeV", i+1, binsY[i]*M_q, binsY[i+1]*M_q));
    }

    // then fill histo with 1D histo values and errors
    for(int iX = 0; iX < nBinsX; iX++) {
      for(int iY = 0; iY < nBinsY; iY++) {
	cHist->SetBinContent(iX+1, iY+1, pHist[iY]->GetBinContent(iX+1)/(ptBins[iY+1]-ptBins[iY]));
	cHist->SetBinError(iX+1, iY+1, pHist[iY]->GetBinError(iX+1)/(ptBins[iY+1]-ptBins[iY]));
      }
    }
    
    // store in outfile the fine-binned histo, the coarse-binned
    // and each of the coarse-binned 1D projections
    TFile* outfile_2 = new TFile("files/ratioHist.root", "update");
    hist->Write();
    cHist->Write();
    outfile_2->Close();
  } 
}
