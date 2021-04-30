void getCos()
{
  double aux = 0;

  ofstream fout;
  fout.open("text_output/cos_max.txt");
  fout.close();  
  
  // read the coarse histos in |costh|
  TFile *infile = new TFile("files/ratioHist.root");
  TH2D *hist = new TH2D();
  //  infile->GetObject(Form("cHist"), hist);
  infile->GetObject(Form("ratioHist_ab_S"), hist);
  hist->SetDirectory(0);
  
  int nBinsY = hist->GetNbinsY();
  const double* yBins = hist->GetYaxis()->GetXbins()->GetArray();
  hist->SetDirectory(0);
  
  infile->Close();

  TH1D *pHist[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    pHist[i] = hist->ProjectionX(Form("fine_bin%d_1d_ab", i+1), i+1, i+1);
  }
    
  double cosMax[nBinsY];
  for(int i = 0; i < nBinsY; i++) {
    int nBinsX = pHist[i]->GetNbinsX();
    for(int j = 1; j < nBinsX; j++) {
      if(pHist[i]->GetBinError(j+2) > 2.*pHist[i]->GetBinError(j+1)) {
	cosMax[i] = pHist[i]->GetBinLowEdge(j+2);
	cout << i << " " << cosMax[i] << endl;
	break;
      }
    }
  }
  
  fout.open("text_output/cos_max.txt", std::ofstream::app);
  for(int i = 0; i < nBinsY; i++) {
    if(i == 0 && cosMax[i] >= aux)
      fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    else if(i == 0 && cosMax[i] < aux) {
      cosMax[i] = aux;
      fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    }
    else if(cosMax[i] >= cosMax[i-1])
      fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    else {
      cosMax[i] = cosMax[i-1];
      fout << yBins[i] << "\t" << yBins[i+1] << "\t" << cosMax[i] << endl;
    }
    aux = cosMax[i];
  }
  fout.close();
}
