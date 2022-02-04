void plotRellpT()
{
  TFile *fin = new TFile("files/RellpTstore_14.root");
  TH2D *r_RpT = (TH2D*)fin->Get("r_RpT");
  TH2D *t_RpT = (TH2D*)fin->Get("t_RpT");
  r_RpT->SetDirectory(0);
  t_RpT->SetDirectory(0);
  fin->Close();

  TH2D *h_RpT = (TH2D*)t_RpT->Clone("h_RpT");
  h_RpT->Sumw2();
  h_RpT->Divide(r_RpT);

  TCanvas *c = new TCanvas("", "", 600, 600);

  gStyle->SetPalette(kRainBow);
  h_RpT->SetStats(0);
  h_RpT->SetTitle("trig/reco, p_ {T} > 50");
  h_RpT->Draw("colz");

  TLine *rlin = new TLine(0.12, 0, 0.08, 120);
  rlin->SetLineColor(kBlue);
  //rlin->Draw("same");

  /*TFile *fout = new TFile("RpTstore.root", "update");
  h_RpT->Write(0, TObject::kOverwrite);
  fout->Close();*/

  c->SaveAs("plots/RellpT_14.pdf");
  c->Destructor();
}
