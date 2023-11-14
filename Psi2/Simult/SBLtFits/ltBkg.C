#import "ltNPPerPt.C"

// macro to run all the 1d lifetime fits
void ltBkg()
{
  // fits - all params free
  TFile *fout = new TFile("files/ltfitres.root", "recreate");
  fout->Close();

  string lbl[] = {"LSB1", "LSB2", "LSB3", "LSB4", "RSB1", "RSB2", "RSB3", "RSB4"};
  for(int i = 0; i < 8; i++) {
    ltPerPt(lbl[i], i);
    cout << endl;
  }
}
