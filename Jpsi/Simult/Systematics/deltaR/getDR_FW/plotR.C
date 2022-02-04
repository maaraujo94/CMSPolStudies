void plotR()
{
  TFile *fin = new TFile("files/Rstore.root");
  TH2D *r_R = (TH2D*)fin->Get("r_R");
  TH2D *t_R = (TH2D*)fin->Get("t_R");
  r_R->SetDirectory(0);
  t_R->SetDirectory(0);
  fin->Close();

  TH2D *h_R = (TH2D*)t_R->Clone("h_R");
  h_R->Sumw2();
  h_R->Divide(r_R);

  TCanvas *c = new TCanvas("", "", 600, 600);

  h_R->SetStats(0);
  h_R->SetMaximum(1);
  h_R->SetTitle("trig/reco, p_ {T} > 50");
  h_R->Draw();

  TLine *rlin = new TLine(0.12, 0, 0.12, h_R->GetMaximum());
  rlin->SetLineColor(kBlue);
  //rlin->Draw("same");

  c->SaveAs("plots/R.pdf");
  c->Destructor();
}
