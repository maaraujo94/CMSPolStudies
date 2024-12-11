// code to plot the Delta R_pT vs pT maps, using MC reco and trig samples
void plotRpTpT()
{
  // read the stored dists
  TFile *fin = new TFile("files/RpTpTstore.root");
  TH2D *r_RpTpT = (TH2D*)fin->Get("r_RpTpT");
  TH2D *t_RpTpT = (TH2D*)fin->Get("t_RpTpT");
  r_RpTpT->SetDirectory(0);
  t_RpTpT->SetDirectory(0);
  fin->Close();

  TH2D *h_RpTpT = (TH2D*)t_RpTpT->Clone("h_RpTpT");
  h_RpTpT->Sumw2();
  h_RpTpT->Divide(r_RpTpT);

  // start plotting
  TCanvas *c = new TCanvas("", "", 600, 600);
  c->SetRightMargin(0.11);
  c->SetLeftMargin(0.11);
  c->SetTopMargin(0.03);

  // first the trig/reco plot
  gStyle->SetPalette(kRainBow);
  h_RpTpT->SetStats(0);
  h_RpTpT->SetMaximum(1);
  h_RpTpT->SetXTitle("#DeltaR_{#it{p}_{T}}");
  h_RpTpT->SetYTitle("#it{p}_{T} (GeV)");
  h_RpTpT->GetYaxis()->SetTitleOffset(1.6);
  h_RpTpT->GetYaxis()->SetLabelOffset(0.01);
  h_RpTpT->GetXaxis()->SetTitleOffset(1.1);
  h_RpTpT->GetYaxis()->SetLabelOffset(0.01);
  h_RpTpT->GetXaxis()->CenterTitle(true);
  h_RpTpT->SetTitle("");
  h_RpTpT->Draw("colz");

  gPad->Update();
  auto palette = (TPaletteAxis*)h_RpTpT->GetListOfFunctions()->FindObject("palette");
  palette->SetY2(120);

  // auxiliary lines
  TLine *rlin = new TLine(0.17, 25, 0.17, 120);
  rlin->SetLineColor(kBlack);
  rlin->SetLineStyle(kDashed);
  rlin->Draw("same");
  TLine *rlin2 = new TLine(0.15, 25, 0.15, 120);
  rlin2->SetLineColor(kBlack);
  rlin2->SetLineStyle(kDashed);
  rlin2->Draw("same");

  int vc = kGreen+4;
  // vertical lines - MC sample separation
  // not currently drawing
  TLine *vlin = new TLine(0, 46, 0.4, 46);
  vlin->SetLineColor(vc);
  vlin->SetLineStyle(kDashed);
  vlin->SetLineWidth(2);
  //  vlin->Draw("same");
  TLine *vlin2 = new TLine(0, 66, 0.4, 66);
  vlin2->SetLineColor(vc);
  vlin2->SetLineWidth(2);
  vlin2->SetLineStyle(kDashed);
  //vlin2->Draw("same");

  // draw text on plot
  TLatex lc;
  lc.SetTextSize(0.03);
  double xp = 0.3, yp = 105;
  lc.DrawLatex(xp, yp, "#bf{trig / reco}"); 

  lc.SetTextSize(0.04);
  yp = 112.5;
  lc.DrawLatex(xp,yp,"#bf{J/#psi}");


  c->SaveAs("plots/RpTpT.pdf");
  c->Clear();
  c->SetLogz();

  // now plot just reco
  r_RpTpT->SetStats(0);
  r_RpTpT->SetTitle("");
  r_RpTpT->SetXTitle("#DeltaR_{#it{p}_{T}}");
  r_RpTpT->SetYTitle("#it{p}_{T} (GeV)");
  r_RpTpT->GetYaxis()->SetTitleOffset(1.6);
  r_RpTpT->GetYaxis()->SetLabelOffset(0.01);
  r_RpTpT->GetXaxis()->SetTitleOffset(1.1);
  r_RpTpT->GetYaxis()->SetLabelOffset(0.01);
  r_RpTpT->GetXaxis()->CenterTitle(true);
  r_RpTpT->Draw("colz");

  gPad->Update();
  auto palette2 = (TPaletteAxis*)r_RpTpT->GetListOfFunctions()->FindObject("palette");
  palette2->SetY2(120);


  rlin->Draw("same");
  rlin2->Draw("same");
  //rlin3->Draw("same");

  lc.SetTextSize(0.03);
  yp = 105;
  lc.DrawLatex(xp, yp, "#bf{reco}"); 

  lc.SetTextSize(0.04);
  yp = 112.5;
  lc.DrawLatex(xp,yp,"#bf{J/#psi}");

  
  c->SaveAs("plots/RpTpT_r.pdf");
  c->Clear();
  c->Destructor();
}
