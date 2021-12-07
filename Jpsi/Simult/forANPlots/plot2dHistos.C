// macro to plot the 2d |costh|:pT maps
// peak data, mc, peak data/MC, NP data/MC, SB data/MC
void plot2dHistos()
{
  // first get the ratio plots
  TFile *fin = new TFile("../PR_fit/files/bkgSubRes.root");
  TH2D *h_Data = (TH2D*)fin->Get("h_Data");
  TH2D *h_NP = (TH2D*)fin->Get("h_NP");
  TH2D *h_SB = (TH2D*)fin->Get("h_SB");

  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);
  c->SetRightMargin(0.13);
  
  h_Data->SetStats(0);
  h_Data->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_Data->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_Data->SetTitle("Full Peak Data/MC");
  h_Data->Draw("COLZ");
  c->SaveAs("plots/2dMaps/ratio_Peak.pdf");
  c->Clear();

  h_NP->SetStats(0);
  h_NP->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_NP->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_NP->SetTitle("Full NP Data/MC");
  h_NP->Draw("COLZ");
  c->SaveAs("plots/2dMaps/ratio_NP.pdf");
  c->Clear();

  h_SB->SetStats(0);
  h_SB->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_SB->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_SB->SetTitle("Full SB Data/MC");
  h_SB->Draw("COLZ");
  c->SaveAs("plots/2dMaps/ratio_SB.pdf");
  c->Clear();

  fin->Close();

  // then get the data and mc plots
  TFile *fin2 = new TFile("../PR_fit/files/histoStore.root");
  TH2D *h_Data2 = (TH2D*)fin2->Get("dataH_ab");
  TH2D *h_NP2 = (TH2D*)fin2->Get("NPH_ab");
  TH2D *h_MC2 = (TH2D*)fin2->Get("mcH_ab");

  c->SetLogz();
 
  h_Data2->SetStats(0);
  h_Data2->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_Data2->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_Data2->SetTitle("Full Peak Data");
  h_Data2->Draw("COLZ");
  c->SaveAs("plots/2dMaps/data_2d_plot.pdf");
  c->Clear();

  h_NP2->SetStats(0);
  h_NP2->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_NP2->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_NP2->SetTitle("Full NP Data");
  h_NP2->Draw("COLZ");
  c->SaveAs("plots/2dMaps/np_2d_plot.pdf");
  c->Clear();

  h_MC2->SetStats(0);
  h_MC2->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_MC2->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_MC2->SetTitle("Full Peak MC");
  h_MC2->Draw("COLZ");
  c->SaveAs("plots/2dMaps/mc_2d_plot.pdf");
  c->Clear();

  c->Destructor();  
  fin2->Close();
}
