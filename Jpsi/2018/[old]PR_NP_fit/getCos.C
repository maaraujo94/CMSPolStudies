void getCos()
{
  TFile *infile = new TFile("files/ratioHist.root");

  TH2D *hist = new TH2D();
  infile->GetObject("ratioHist_ab", hist);
  int nBinsY = hist->GetNbinsY();
  const double* yBins = hist->GetYaxis()->GetXbins()->GetArray();
  hist->SetDirectory(0);
  
  infile->Close();
  TH1D *pHist[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    pHist[i] = hist->ProjectionX(Form("fine_bin%d_1d_ab", i+1), i+1, i+1);
  }
  cout << yBins[0] << " " << yBins[1] << endl;

  double cosMax[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    int nBinsX = pHist[i]->GetNbinsX();
    for(int j = 1; j < nBinsX; j++) {
      if(pHist[i]->GetBinError(j+2) > 2.*pHist[i]->GetBinError(j+1)) {
	  cosMax[i] = pHist[i]->GetBinLowEdge(j+2);
	  cout << "[" << yBins[i] << "," << yBins[i+1] << "]: cosMax = " << cosMax[i] << endl;
	  break;
      }
    }
  }

  ofstream fout;
  fout.open("text_output/cos_max.txt");
  for(int i = 0; i < nBinsY; i++) {
    if(i == 0)
      fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    else if(cosMax[i] >= cosMax[i-1])
      fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    else {
      cosMax[i] = cosMax[i-1];
      fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    }
  }
  fout.close();
}
