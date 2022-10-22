// macro to comp the sideband costh dists between years

void compSB()
{
  const int n_pt = 19;
  
  TFile *fin7 = new TFile("../2017/PR_fit/files/bkgCosModel.root");
  TH1D **h_SB7 = new TH1D*[n_pt];
  for(int i = 0; i < n_pt; i++){
    h_SB7[i] = (TH1D*)fin7->Get(Form("h_SB_%d",i));
    h_SB7[i]->SetDirectory(0);
  }
  fin7->Close();

  TFile *fin8 = new TFile("../2018/PR_fit/files/bkgCosModel.root");
  TH1D **h_SB8 = new TH1D*[n_pt];
  for(int i = 0; i < n_pt; i++){
    h_SB8[i] = (TH1D*)fin8->Get(Form("h_SB_%d",i));
    h_SB8[i]->SetDirectory(0);
  }
  fin8->Close();

  // plot the 1d projection of the result (there's no pT dependence)
  TCanvas *c = new TCanvas("", "", 900, 900);
  h_SB7[0]->SetStats(0);
  h_SB7[0]->SetMinimum(0);
  h_SB7[0]->SetLineColor(kBlack);
  h_SB7[0]->SetTitle("bkg/MC comp");
  h_SB7[0]->GetXaxis()->SetTitle("|cos#theta|");
  h_SB7[0]->Draw();

  h_SB8[0]->SetStats(0);
  h_SB8[0]->SetLineColor(kBlue);
  h_SB8[0]->Draw("same");

  TLegend *leg = new TLegend(0.7, 0.2, 0.9, 0.4);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_SB7[0], "2017", "l");
  leg->AddEntry(h_SB8[0], "2018", "l");
  leg->Draw();

  c->SaveAs("SB_base_0.pdf");
  c->Clear();

  h_SB7[10]->SetStats(0);
  h_SB7[10]->SetMinimum(0);
  h_SB7[10]->SetLineColor(kBlack);
  h_SB7[10]->SetTitle("bkg/MC comp");
  h_SB7[10]->GetXaxis()->SetTitle("|cos#theta|");
  h_SB7[10]->Draw();

  h_SB8[10]->SetStats(0);
  h_SB8[10]->SetLineColor(kBlue);
  h_SB8[10]->Draw("same");

  leg->Draw();

  c->SaveAs("SB_base_10.pdf");
  c->Clear();
  c->Destructor();

}
