#import "ltNPPerPt.C"

// macro to run all the 1d lifetime fits
void ltBkg()
{
  // fits - all params free
  TFile *fout = new TFile("files/ltfitres.root", "recreate");
  fout->Close();

  string lbl[] = {"LSB", "RSB"};
  for(int i = 0; i < 2; i++) {
    ltPerPt(lbl[i], i);
    cout << endl;
  }
}
