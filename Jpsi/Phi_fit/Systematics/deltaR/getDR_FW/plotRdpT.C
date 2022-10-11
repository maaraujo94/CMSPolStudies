void plotRdpT()
{
  TFile *fin = new TFile("files/RdpTstore.root");
  TH2D *r_RdpT = (TH2D*)fin->Get("r_RdpT");
  TH2D *t_RdpT = (TH2D*)fin->Get("t_RdpT");
  r_RdpT->SetDirectory(0);
  t_RdpT->SetDirectory(0);
  fin->Close();

  TH2D *h_RdpT = (TH2D*)t_RdpT->Clone("h_RdpT");
  h_RdpT->Sumw2();
  h_RdpT->Divide(r_RdpT);

  TCanvas *c = new TCanvas("", "", 600, 600);

  gStyle->SetPalette(kRainBow);
  h_RdpT->SetStats(0);
  h_RdpT->SetXTitle("#DeltaR_{p_{T}}");
  h_RdpT->SetYTitle("|#Deltap_{T}| (GeV)");
  h_RdpT->GetYaxis()->SetTitleOffset(1.3);
  h_RdpT->SetTitle("trig/reco, p_{T} > 50");
  h_RdpT->Draw("colz");

  TLine *rlin = new TLine(0.192, 0, 0.192, 120);
  rlin->SetLineColor(kBlack);
  rlin->SetLineStyle(kDashed);
  rlin->Draw("same");
  TLine *rlin2 = new TLine(0.096, 0, 0.096, 120);
  rlin2->SetLineColor(kBlue);
  //rlin2->Draw("same");

  c->SaveAs("plots/RdpT.pdf");
  c->Destructor();
}
