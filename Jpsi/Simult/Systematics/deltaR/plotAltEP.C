void plotAltEP()
{
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

  TCanvas *c = new TCanvas("", "", 600, 600);
  gStyle->SetPalette(kRainBow);

  r_EtaPhi->SetStats(0);
  r_EtaPhi->Draw("colz");

  c->SaveAs("plots/AEtaPhi_r.pdf");
  c->Clear();

  t_EtaPhi->SetStats(0);
  t_EtaPhi->Draw("colz");

  c->SaveAs("plots/AEtaPhi_t.pdf");
  c->Clear();

  TH2D *h_EtaPhi = (TH2D*)t_EtaPhi->Clone("h_EtaPhi");
  h_EtaPhi->Sumw2();
  h_EtaPhi->Divide(r_EtaPhi);

  h_EtaPhi->SetStats(0);
  h_EtaPhi->SetTitle("trig/reco, p_ {T} > 50");
  h_EtaPhi->Draw("colz");

  TArc *arc1 = new TArc(0, 0, 0.122, 0, 360);
  arc1->SetLineColor(kBlue);
  arc1->SetFillColorAlpha(kBlack, 0);
  //arc1->Draw("same");

  /* TLine *flin1 = new TLine(-0.06, 0.12, 0.06, 0.12);
  flin1->Draw("same");
  TLine *flin4 = new TLine(-0.06, -0.12, 0.06, -0.12);
  flin4->Draw("same");

  TLine *flin3 = new TLine(0.12, 0.03, 0.12, -0.03);
  flin3->Draw("same");
  TLine *flin6 = new TLine(-0.12, 0.03, -0.12, -0.03);
  flin6->Draw("same");
  
  TLine *flin2 = new TLine(0.06, 0.12, 0.12, 0.03);
  flin2->Draw("same");
  TLine *flin5 = new TLine(-0.06, 0.12, -0.12, 0.03);
  flin5->Draw("same");
  TLine *flin7 = new TLine(-0.06, -0.12, -0.12, -0.03);
  flin7->Draw("same");
  TLine *flin8 = new TLine(0.06, -0.12, 0.12, -0.03);
  flin8->Draw("same");*/

  TEllipse *el2 = new TEllipse(0, 0, 0.12, 0.124, 0, 360);
  el2->SetLineColor(kViolet);
  el2->SetFillColorAlpha(kBlack, 0);
  el2->Draw("same");
  
  c->SaveAs("plots/AEtaPhi_full.pdf");
  c->Clear();

  TH2D *ha_EtaPhi = (TH2D*)ta_EtaPhi->Clone("ha_EtaPhi");
  ha_EtaPhi->Sumw2();
  ha_EtaPhi->Divide(ra_EtaPhi);

  ha_EtaPhi->SetStats(0);
  ha_EtaPhi->SetTitle("trig/reco, p_ {T} > 50");
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
