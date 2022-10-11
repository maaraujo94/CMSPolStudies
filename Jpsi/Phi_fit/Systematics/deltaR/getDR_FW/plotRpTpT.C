void plotRpTpT()
{
  TFile *fin = new TFile("files/RpTpTstore.root");
  TH2D *r_RpTpT = (TH2D*)fin->Get("r_RpTpT");
  TH2D *t_RpTpT = (TH2D*)fin->Get("t_RpTpT");
  r_RpTpT->SetDirectory(0);
  t_RpTpT->SetDirectory(0);
  fin->Close();

  TH2D *h_RpTpT = (TH2D*)t_RpTpT->Clone("h_RpTpT");
  h_RpTpT->Sumw2();
  h_RpTpT->Divide(r_RpTpT);

  TCanvas *c = new TCanvas("", "", 600, 600);

  gStyle->SetPalette(kRainBow);
  h_RpTpT->SetStats(0);
  h_RpTpT->SetMaximum(1);
  h_RpTpT->SetXTitle("#DeltaR_{p_{T}}");
  h_RpTpT->SetYTitle("p_{T} (GeV)");
  h_RpTpT->GetYaxis()->SetTitleOffset(1.3);
  h_RpTpT->SetTitle("trig/reco");
  h_RpTpT->Draw("colz");

  TLine *rlin = new TLine(0.17, 25, 0.17, 120);
  rlin->SetLineColor(kBlack);
  rlin->SetLineStyle(kDashed);
  rlin->Draw("same");
  TLine *rlin2 = new TLine(0.15, 25, 0.15, 120);
  rlin2->SetLineColor(kBlack);
  rlin2->SetLineStyle(kDashed);
  rlin2->Draw("same");
  TLine *rlin3 = new TLine(0.12, 25, 0.12, 120);
  rlin3->SetLineColor(kBlue);
  rlin3->SetLineStyle(kDashed);
  //rlin3->Draw("same");

  int vc = kGreen+4;
  // vertical lines - MC sample separation
  TLine *vlin = new TLine(0, 46, 0.4, 46);
  vlin->SetLineColor(vc);
  vlin->SetLineStyle(kDashed);
  vlin->SetLineWidth(2);
  vlin->Draw("same");
  TLine *vlin2 = new TLine(0, 66, 0.4, 66);
  vlin2->SetLineColor(vc);
  vlin2->SetLineWidth(2);
  vlin2->SetLineStyle(kDashed);
  vlin2->Draw("same");

  c->SaveAs("plots/RpTpT.pdf");
  c->Clear();
  c->SetLogz();
  
  r_RpTpT->SetStats(0);
  r_RpTpT->SetTitle("reco");
  r_RpTpT->SetXTitle("#DeltaR_{p_{T}}");
  r_RpTpT->SetYTitle("p_{T}| (GeV)");
  r_RpTpT->GetYaxis()->SetTitleOffset(1.3);
  r_RpTpT->Draw("colz");

  rlin->Draw("same");
  rlin2->Draw("same");
  //rlin3->Draw("same");
  
  c->SaveAs("plots/RpTpT_r.pdf");
  c->Clear();
  c->Destructor();
}
