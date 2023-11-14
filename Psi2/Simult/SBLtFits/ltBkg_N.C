#import "ltPerPt_N.C"

// macro to run all the 1d lifetime fits
void ltBkg_N()
{
  // fits - all params free
  TFile *fout = new TFile("files/ltfitres_N.root", "recreate");
  fout->Close();

  for(int i = 0; i < nPtBins; i++) {
    ltPerPt_N(ptBins[i], ptBins[i+1], i);
  }  
}
