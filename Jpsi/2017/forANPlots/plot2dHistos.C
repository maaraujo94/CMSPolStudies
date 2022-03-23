// macro to plot the 2d |costh|:pT maps
// peak data, np, mc, peak data/MC, NP data/MC
void plot2dHistos()
{
  // get everything from histoStore
  TFile *fin2 = new TFile("../PR_fit/files/histoStore.root");
  // first the ratios
  TH2D *h_Data = (TH2D*)fin2->Get("ratioH_ab");
  TH2D *h_NP = (TH2D*)fin2->Get("ratNPH_ab"); // ORIGINAL RATIO: not including bkg sub
  // then the base dists
  TH2D *h_Data2 = (TH2D*)fin2->Get("dataH_ab");
  TH2D *h_NP2 = (TH2D*)fin2->Get("NPH_ab");
  TH2D *h_MC2 = (TH2D*)fin2->Get("mcH_ab");

  // first draw the ratios
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);
  c->SetRightMargin(0.13);

  // scaling Peak ratio
  double maxR[3];
  int bin_min[] = {1,8,12};
  int bin_max[] = {7,11,17};
  for(int i = 0; i < 3; i++) {
    h_Data->GetYaxis()->SetRange(bin_min[i], bin_max[i]);
    maxR[i] = h_Data->GetMaximum();
  }
  h_Data->GetYaxis()->SetRange(bin_min[0], bin_max[2]);

  double scalF;
  for(int i_x = 0; i_x < h_Data->GetNbinsX(); i_x++) {
    for(int i_y = 0; i_y < h_Data->GetNbinsY(); i_y++){
      // get scaling factor
      for(int i = 0; i < 3; i++) 
	if(i_y+1 >= bin_min[i] && i_y+1 <= bin_max[i]) 
	  scalF = maxR[0]/maxR[i];
      h_Data->SetBinContent(i_x+1, i_y+1, h_Data->GetBinContent(i_x+1, i_y+1)*scalF);
    }
  }
  
  h_Data->SetStats(0);
  h_Data->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_Data->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_Data->SetTitle("2017 Peak Data/MC");
  h_Data->Draw("COLZ");
  c->SaveAs("plots/2dMaps/ratio_Peak.pdf");
  c->Clear();

  // scaling NP ratio
  for(int i = 0; i < 3; i++) {
    h_NP->GetYaxis()->SetRange(bin_min[i], bin_max[i]);
    maxR[i] = h_NP->GetMaximum();
  }
  h_NP->GetYaxis()->SetRange(bin_min[0], bin_max[2]);

  for(int i_x = 0; i_x < h_NP->GetNbinsX(); i_x++) {
    for(int i_y = 0; i_y < h_NP->GetNbinsY(); i_y++){
      // get scaling factor
      for(int i = 0; i < 3; i++) 
	if(i_y+1 >= bin_min[i] && i_y+1 <= bin_max[i]) 
	  scalF = maxR[0]/maxR[i];
      h_NP->SetBinContent(i_x+1, i_y+1, h_NP->GetBinContent(i_x+1, i_y+1)*scalF);
    }
  }

  h_NP->SetStats(0);
  h_NP->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_NP->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_NP->SetTitle("2017 NP Data/MC");
  h_NP->Draw("COLZ");
  c->SaveAs("plots/2dMaps/ratio_NP.pdf");
  c->Clear();

  // then draw the base dists
  c->SetLogz();
 
  h_Data2->SetStats(0);
  h_Data2->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_Data2->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_Data2->SetTitle("2017 Peak Data");
  h_Data2->Draw("COLZ");
  c->SaveAs("plots/2dMaps/data_2d_plot.pdf");
  c->Clear();

  h_NP2->SetStats(0);
  h_NP2->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_NP2->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_NP2->SetTitle("2017 NP Data");
  h_NP2->Draw("COLZ");
  c->SaveAs("plots/2dMaps/np_2d_plot.pdf");
  c->Clear();

  h_MC2->SetStats(0);
  h_MC2->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_MC2->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_MC2->SetTitle("2017 Peak MC");
  h_MC2->Draw("COLZ");
  c->SaveAs("plots/2dMaps/mc_2d_plot.pdf");
  c->Clear();

  c->Destructor();  
  fin2->Close();
}
