void plotEP()
{
  TFile *fin = new TFile("files/EPstore.root");
  TH2D *r_EtaPhi = (TH2D*)fin->Get("r_EtaPhi_0");
  TH2D *t_EtaPhi = (TH2D*)fin->Get("t_EtaPhi_0");
  //TH2D *r_EtaPhi = (TH2D*)fin->Get("rC_EtaPhi");
  //TH2D *t_EtaPhi = (TH2D*)fin->Get("tC_EtaPhi");
  r_EtaPhi->SetDirectory(0);
  t_EtaPhi->SetDirectory(0);
  fin->Close();

  TCanvas *c = new TCanvas("", "", 600, 600);
  gStyle->SetPalette(kRainBow);

  r_EtaPhi->SetStats(0);
  r_EtaPhi->Draw("colz");

  c->SaveAs("plots/EtaPhi_r.pdf");
  c->Clear();

  t_EtaPhi->SetStats(0);
  t_EtaPhi->Draw("colz");

  c->SaveAs("plots/EtaPhi_t.pdf");
  c->Clear();

  TH2D *h_EtaPhi = (TH2D*)t_EtaPhi->Clone("h_EtaPhi");
  h_EtaPhi->Sumw2();
  h_EtaPhi->Divide(r_EtaPhi);

  h_EtaPhi->SetStats(0);
  h_EtaPhi->SetTitle("trig/reco, p_ {T} > 50");
  h_EtaPhi->Draw("colz");

  TArc *arc1 = new TArc(0, 0, 0.08, 0, 90);
  arc1->SetLineColor(kBlue);
  arc1->SetFillColorAlpha(kBlack, 0);
  arc1->Draw("same");

  TArc *arc2 = new TArc(0, 0, 0.12, 0, 90);
  arc2->SetLineColor(kBlue);
  arc2->SetFillColorAlpha(kBlack, 0);
  arc2->Draw("same"); 

  //cout << 0.122/0.116 << endl;
  
  TEllipse *el2 = new TEllipse(0, 0, 0.12, 0.1, 0, 90);
  el2->SetLineColor(kCyan);
  el2->SetFillColorAlpha(kBlack, 0);
  el2->Draw("same");

  c->SaveAs("plots/EtaPhi_neg.pdf");
  c->Destructor();
}
