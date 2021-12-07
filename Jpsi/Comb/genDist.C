// macro to comp the sideband costh dists between years

void genDist()
{
  TFile *fin7 = new TFile("../2017/PR_fit/files/bkgCosModel.root");
  TH2D *h_SB7 = (TH2D*)fin7->Get("h_SB");
  h_SB7->SetDirectory(0);
  fin7->Close();

  TFile *fin8 = new TFile("../2018/PR_fit/files/bkgCosModel.root");
  TH2D *h_SB8 = (TH2D*)fin8->Get("h_SB");
  h_SB8->SetDirectory(0);
  fin8->Close();

  // plot the 1d projection of the result (there's no pT dependence)
  TH1D *h_SB1d7 = h_SB7->ProjectionX("h_SB_1d7", 1, 1);
  TH1D *h_SB1d8 = h_SB8->ProjectionX("h_SB_1d8", 1, 1);
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  h_SB1d7->SetStats(0);
  h_SB1d7->SetMinimum(0);
  h_SB1d7->SetMaximum(2.5);
  h_SB1d7->SetLineColor(kBlack);
  h_SB1d7->SetTitle("bkg/MC comp");
  h_SB1d7->GetXaxis()->SetTitle("|cos#theta|");
  h_SB1d7->Draw();

  h_SB1d8->SetStats(0);
  h_SB1d8->SetLineColor(kBlue);
  h_SB1d8->Draw("same");

  TLegend *leg = new TLegend(0.7, 0.2, 0.9, 0.4);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_SB1d7, "2017", "l");
  leg->AddEntry(h_SB1d8, "2018", "l");
  leg->Draw();

  c->SaveAs("SB_base.pdf");
  c->Clear();
  c->Destructor();

}
