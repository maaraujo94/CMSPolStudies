// macro to plot the 2d |costh|:pT maps
// peak data, mc, peak data/MC, NP data/MC, SB data/MC
void plot2dHistos()
{
  // first get the data for 2018
  TFile *fin18 = new TFile("../2018/PR_fit/files/histoStore.root");
  TH2D *h_Data18 = (TH2D*)fin18->Get("dataH_ab");
  h_Data18->SetDirectory(0);
  fin18->Close();

  // then get the data for 2017
  TFile *fin17 = new TFile("../2017/PR_fit/files/histoStore.root");
  TH2D *h_Data17 = (TH2D*)fin17->Get("dataH_ab");
  h_Data17->SetDirectory(0);
  fin17->Close();

  // then add both together
  h_Data17->Add( h_Data18);
  
  TCanvas *c = new TCanvas("", "", 900, 900);
  c->SetLeftMargin(0.11);
  c->SetRightMargin(0.13);
  
  c->SetLogz();
 
  h_Data17->SetStats(0);
  h_Data17->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  h_Data17->GetYaxis()->SetTitle("p_{T} (GeV)");
  h_Data17->SetTitle("Peak Data");
  h_Data17->Draw("COLZ");
  c->SaveAs("data_2d_plot.pdf");
  c->Clear();

  c->Destructor();  
}
