// macro to get f_bkg in phi binning

void fbkgProp()
{
  // get the original f_bkg in costheta binning
  TH2D *h_fbkg = new TH2D();
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/NP_fit/files/bkgFrac.root");
  h_fbkg = (TH2D*)fin->Get("h_fbkg");
  h_fbkg->SetDirectory(0);
  fin->Close();
    
  // get the binning
  int nBinsX = h_fbkg->GetNbinsX(), nBinsY = h_fbkg->GetNbinsY();
  const double *yBins = h_fbkg->GetYaxis()->GetXbins()->GetArray();
  double minX = h_fbkg->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_fbkg->GetXaxis()->GetBinUpEdge(nBinsX);

  //define new h_fbkg with phi bins
  double minP = -180, maxP = 180;
  TH2D *h_fbkg2d = new TH2D("h_fbkg", "Run 2 f_{bkg}", nBinsX, minP, maxP, nBinsY, yBins);

  // copy the values over
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    // same result for all phi bins
    for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
      h_fbkg2d->SetBinContent(i_cos+1, i_pt+1, h_fbkg->GetBinContent(i_cos+1, i_pt+1));
      h_fbkg2d->SetBinError(i_cos+1, i_pt+1, h_fbkg->GetBinError(i_cos+1, i_pt+1));
    }
  }

  // save histo
  TFile *fout = new TFile("files/bkgFrac.root", "recreate");
  h_fbkg2d->SetName("h_fbkg");
  h_fbkg2d->Write();
  fout->Close();

}
