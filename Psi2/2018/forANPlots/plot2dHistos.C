// macro to plot the 2d |costh|:pT maps
// peak data, mc, peak data/MC, NP data/MC, SB data/MC
void plot2dHistos()
{
  // get everything from histoStore
  TFile *fin2 = new TFile("../PR_fit/files/histoStore.root");
  // first the ratios
  TH2D *h_Data = (TH2D*)fin2->Get("rPRH");
  TH2D *h_NP = (TH2D*)fin2->Get("rNPH"); // ORIGINAL RATIO: not including bkg sub
  // then the base dists
  TH2D *h_Data2 = (TH2D*)fin2->Get("PRH");
  TH2D *h_NP2 = (TH2D*)fin2->Get("NPH");
  TH2D *h_MC2 = (TH2D*)fin2->Get("MCH");

  // first draw the ratios
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);
  c->SetRightMargin(0.12);
  c->SetTopMargin(0.02);

  h_Data->SetStats(0);
  h_Data->GetXaxis()->SetTitle("|cos #theta_{HX}|");
  h_Data->GetYaxis()->SetTitle("#it{p}_{T} (GeV)");
  h_Data->SetTitle("");
  h_Data->GetYaxis()->SetTitleOffset(1.6);
  h_Data->GetYaxis()->SetLabelOffset(0.01);
  h_Data->GetXaxis()->SetLabelOffset(0.01);
  h_Data->GetXaxis()->SetTitleOffset(1.25);
  h_Data->GetXaxis()->CenterTitle(true);
  h_Data->Draw("COLZ");
  c->SaveAs("plots/2dMaps/ratio_PRS.pdf");
  c->Clear();

  h_NP->SetStats(0);
  h_NP->GetXaxis()->SetTitle("|cos #theta_{HX}|");
  h_NP->GetYaxis()->SetTitle("#it{p}_{T} (GeV)");
  h_NP->SetTitle("");
  h_NP->Draw("COLZ");
  c->SaveAs("plots/2dMaps/ratio_NPS.pdf");
  c->Clear();

  // then draw the base dists
  c->SetLogz();
 
  h_Data2->SetStats(0);
  h_Data2->GetXaxis()->SetTitle("|cos #theta_{HX}|");
  h_Data2->GetYaxis()->SetTitle("#it{p}_{T} (GeV)");
  h_Data2->SetTitle("");
  h_Data2->GetYaxis()->SetTitleOffset(1.6);
  h_Data2->GetYaxis()->SetLabelOffset(0.01);
  h_Data2->GetXaxis()->SetLabelOffset(0.01);
  h_Data2->GetXaxis()->SetTitleOffset(1.25);
  h_Data2->GetXaxis()->CenterTitle(true);
  h_Data2->Draw("COLZ");
  c->SaveAs("plots/2dMaps/data_2d_plot.pdf");
  c->Clear();

  h_NP2->SetStats(0);
  h_NP2->GetXaxis()->SetTitle("|cos #theta_{HX}|");
  h_NP2->GetYaxis()->SetTitle("#it{p}_{T} (GeV)");
  h_NP2->SetTitle("");
  h_NP2->Draw("COLZ");
  c->SaveAs("plots/2dMaps/np_2d_plot.pdf");
  c->Clear();

  h_MC2->SetStats(0);
  h_MC2->GetXaxis()->SetTitle("|cos #theta_{HX}|");
  h_MC2->GetYaxis()->SetTitle("#it{p}_{T} (GeV)");
  h_MC2->SetTitle("");
  h_MC2->GetYaxis()->SetTitleOffset(1.6);
  h_MC2->GetYaxis()->SetLabelOffset(0.01);
  h_MC2->GetXaxis()->SetLabelOffset(0.01);
  h_MC2->GetXaxis()->SetTitleOffset(1.25);
  h_MC2->GetXaxis()->CenterTitle(true);
  h_MC2->Draw("COLZ");
  c->SaveAs("plots/2dMaps/mc_2d_plot.pdf");
  c->Clear();

  c->Destructor();  
  fin2->Close();
}
