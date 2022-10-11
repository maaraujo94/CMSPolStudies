// macro to get f_np in phi binning

void fnpProp()
{
  // get the original f_np in costheta binning
  TH2D *h_fnp = new TH2D();
  TFile *fin = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Simult/PR_fit/files/NPFrac.root");
  h_fnp = (TH2D*)fin->Get("h_fNPc");
  h_fnp->SetDirectory(0);
  fin->Close();
    
  // get the binning
  int nBinsX = h_fnp->GetNbinsX(), nBinsY = h_fnp->GetNbinsY();
  const double *yBins = h_fnp->GetYaxis()->GetXbins()->GetArray();
  double minX = h_fnp->GetXaxis()->GetBinLowEdge(1);
  double maxX = h_fnp->GetXaxis()->GetBinUpEdge(nBinsX);

  //define new h_fnp with phi bins
  double minP = -180, maxP = 180;
  TH2D *h_fnp2d = new TH2D("h_fNPc", "Run 2 f_{np}", nBinsX, minP, maxP, nBinsY, yBins);

  // copy the values over
  for(int i_pt = 0; i_pt < nBinsY; i_pt++) {
    // same result for all phi bins
    for(int i_cos = 0; i_cos < nBinsX; i_cos++) {
      h_fnp2d->SetBinContent(i_cos+1, i_pt+1, h_fnp->GetBinContent(i_cos+1, i_pt+1));
      h_fnp2d->SetBinError(i_cos+1, i_pt+1, h_fnp->GetBinError(i_cos+1, i_pt+1));
    }
  }

  // save histo  
  TFile *fout = new TFile("files/NPFrac.root", "recreate");
  h_fnp2d->SetName("h_fNPc");
  h_fnp2d->Write();
  fout->Close();
}
