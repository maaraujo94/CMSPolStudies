void plotRpT()
{
  TFile *fin = new TFile("files/RpTstore.root");
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
  h_RpT->SetXTitle("#DeltaR");
  h_RpT->SetYTitle("|#Deltap_{T}| (GeV)");
  h_RpT->GetYaxis()->SetTitleOffset(1.3);
  h_RpT->SetTitle("trig/reco, p_{T} > 50");
  h_RpT->Draw("colz");

  TLine *rlin = new TLine(0.128, 0, 0.128, 120);
  rlin->SetLineColor(kBlue);
  //rlin->Draw("same");
  TLine *rlin2 = new TLine(0.096, 0, 0.096, 120);
  rlin2->SetLineColor(kBlue);
  //rlin2->Draw("same");

  TF1 *fc = new TF1("fc", "[0]*exp(-x*[1])", 0, 0.3);
  fc->SetParameters(5600, 45);
  fc->SetLineColor(kBlack);
  fc->SetLineStyle(kDashed);
  fc->Draw("lsame");

  /*TFile *fout = new TFile("RpTstore.root", "update");
  h_RpT->Write(0, TObject::kOverwrite);
  fout->Close();*/

  c->SaveAs("plots/RpT.pdf");
  c->Destructor();
}
