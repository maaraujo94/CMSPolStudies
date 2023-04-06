// macro to scale f_BG by 0.8

void fbkgProp()
{
  TFile *fin = new TFile("../../bkgFits/files/bkgFrac.root");
  TH1D *fbg_1d = (TH1D*)fin->Get("fbkg_unc");
  fbg_1d->SetDirectory(0);
  TH2D *fbg_2d = (TH2D*)fin->Get("h_fbkg");
  fbg_2d->SetDirectory(0);
  fin->Close();

  fbg_1d->Scale(1.2);
  fbg_2d->Scale(1.2);

  TFile *fout = new TFile("files/bkgFrac.root", "recreate");
  fbg_1d->Write();
  fbg_2d->Write();
  fout->Close();  

}
