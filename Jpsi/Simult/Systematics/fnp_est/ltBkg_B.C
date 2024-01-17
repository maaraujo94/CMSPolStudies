#import "ltBPerPt.C"
#import "/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Simult/ptbins.C"

// macro to run all the 1d lifetime fits
void ltBkg_B()
{
  // fits - all params free
  TFile *fout = new TFile("files/ltfitresB.root", "recreate");
  fout->Close();

 double fracNP[nPtBins], pt[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    pt[i] = 0.5*(ptBins[i+1]+ptBins[i]);
    fracNP[i] = ltPerPt(ptBins[i], ptBins[i+1], i);
    cout << endl;
  }
  TGraph *g_fNP = new TGraph(nPtBins, pt, fracNP);
  
  TFile *fout2 = new TFile("files/ltfitresB.root", "update");
  g_fNP->SetName("g_fNP");
  g_fNP->Write(0, TObject::kOverwrite);
  fout2->Close();
 
}
