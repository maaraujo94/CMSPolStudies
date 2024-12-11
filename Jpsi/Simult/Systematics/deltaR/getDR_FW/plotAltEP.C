// code to plot Delta Eta vs Delta Phi maps, using MC reco and trig

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

void plotAltEP()
{
  // read plots (from store code)
  TFile *fin = new TFile("files/AltEPstore.root");
  TH2D *ra_EtaPhi = (TH2D*)fin->Get("r_EtaPhi_0");
  TH2D *ta_EtaPhi = (TH2D*)fin->Get("t_EtaPhi_0");
  TH2D *r_EtaPhi = (TH2D*)fin->Get("rC_EtaPhi");
  TH2D *t_EtaPhi = (TH2D*)fin->Get("tC_EtaPhi");
  r_EtaPhi->SetDirectory(0);
  t_EtaPhi->SetDirectory(0);
  ra_EtaPhi->SetDirectory(0);
  ta_EtaPhi->SetDirectory(0);
  fin->Close();

  // make canvas
  TCanvas *c = new TCanvas("", "", 600, 600);
  gStyle->SetPalette(kRainBow);

  // plot reco map
  r_EtaPhi->SetStats(0);
  r_EtaPhi->Draw("colz");

  c->SaveAs("plots/AEtaPhi_r.pdf");
  c->Clear();

  // plot trig map
  t_EtaPhi->SetStats(0);
  t_EtaPhi->Draw("colz");

  c->SaveAs("plots/AEtaPhi_t.pdf");
  c->Clear();

  // plot ratio
  c->SetTopMargin(0.03);
  TH2D *h_EtaPhi = new TH2D();
  h_EtaPhi = (TH2D*)t_EtaPhi->Clone("h_EtaPhi");
  h_EtaPhi->Sumw2();
  h_EtaPhi->Divide(r_EtaPhi);

  h_EtaPhi->SetStats(0);
  h_EtaPhi->SetXTitle("#Delta#eta");
  h_EtaPhi->SetYTitle("#Delta#phi");
  h_EtaPhi->GetXaxis()->SetRangeUser(-0.29,0.29);
  h_EtaPhi->GetYaxis()->SetRangeUser(-0.29,0.29);
  h_EtaPhi->GetYaxis()->SetTitleOffset(1.4);
  h_EtaPhi->GetYaxis()->SetLabelOffset(0.01);
  h_EtaPhi->GetXaxis()->SetTitleOffset(1.1);
  h_EtaPhi->GetYaxis()->SetLabelOffset(0.01);
  h_EtaPhi->GetXaxis()->CenterTitle(true);
  h_EtaPhi->SetTitle("");
  h_EtaPhi->Draw("colz");

  gPad->Update();
  auto palette = (TPaletteAxis*)h_EtaPhi->GetListOfFunctions()->FindObject("palette");
  palette->SetY2(0.29);

  TEllipse *el2 = new TEllipse(0, 0, 0.12, 0.12, 0, 360);
  el2->SetLineColor(kBlack);
  el2->SetLineStyle(kDashed);
  el2->SetFillColorAlpha(kBlack, 0);
  el2->Draw("same");

  // draw text on plot
  TLatex lc;
  lc.SetTextSize(0.03);
  double xp = 0.15, yp = 0.25;
  lc.DrawLatex(xp, yp, "#bf{trig / reco}"); 
  yp = 0.225;
  lc.DrawLatex(xp, yp, "#bf{#it{p}_{T} > 50 GeV}");

  lc.SetTextSize(0.04);
  xp = -0.25;
  yp = 0.25;
  lc.DrawLatex(xp,yp,"#bf{J/#psi}");
  
  c->SaveAs("plots/AEtaPhi_full.pdf");
  c->Clear();

  // |Delta| map
  TH2D *ha_EtaPhi = (TH2D*)ta_EtaPhi->Clone("ha_EtaPhi");
  ha_EtaPhi->Sumw2();
  ha_EtaPhi->Divide(ra_EtaPhi);

  ha_EtaPhi->SetStats(0);
  ha_EtaPhi->SetTitle("trig/reco, p_{T} > 50");
  ha_EtaPhi->Draw("colz");

  /*  TArc *arc1 = new TArc(0, 0, 0.196, 0, 90);
  arc1->SetLineColor(kBlue);
  arc1->SetFillColorAlpha(kBlack, 0);
  arc1->Draw("same");*/

  TArc *arc2 = new TArc(0, 0, 0.122, 0, 90);
  arc2->SetLineColor(kBlue);
  arc2->SetFillColorAlpha(kBlack, 0);
  //arc2->Draw("same");

  //cout << 0.122/0.116 << endl;
  
  /* TEllipse *el2 = new TEllipse(0, 0, 0.122, 0.116, 0, 90);
  el2->SetLineColor(kViolet);
  el2->SetFillColorAlpha(kBlack, 0);
  el2->Draw("same");*/

  TLine *lin1 = new TLine(0, 0.12, 0.06, 0.12);
  //lin1->Draw("same");
  TLine *lin2 = new TLine(0.06, 0.12, 0.12, 0.03);
  //lin2->Draw("same");
  TLine *lin3 = new TLine(0.12, 0.03, 0.12, 0.0);
  //lin3->Draw("same");

  c->SaveAs("plots/AEtaPhi_abs.pdf");
  c->Clear();

  c->Destructor();
}
