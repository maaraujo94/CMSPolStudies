void plotCos()
{
  // open files and read TTrees
  TFile *finD = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/data18_cos.root");
  TTree *treeD = (TTree*)finD->Get("data_cos");
  TFile *finM = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/MC18_cos.root");
  TTree *treeM = (TTree*)finM->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int mEvt = treeM->GetEntries();

  TH1D **dataHist = new TH1D*[2];
  dataHist[0] = new TH1D("dataH_PR", "2018 Data (PR)", 20, 0, 1.);
  dataHist[1] = new TH1D("dataH_NP", "2018 Data (NP)", 20, 0, 1.);
  TH1D *mcHist = new TH1D("mcH", "2018 MC", 20, 0, 1.);

  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m, data_y;
  Double_t mc_th, mc_pt, mc_lt, mc_m, mc_y;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  treeM->SetBranchAddress("theta", &mc_th);
  treeM->SetBranchAddress("dimPt", &mc_pt);
  treeM->SetBranchAddress("Rap", &mc_y);
  treeM->SetBranchAddress("Mass", &mc_m);
  treeM->SetBranchAddress("lt", &mc_lt);
  
  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(data_pt > 39 && data_pt < 46) {
	// NP SR
	if(data_lt > 0.014 && data_lt < 0.05 && data_m < 3.2 && data_m > 3.0) {
	  dataHist[1]->Fill(abs(cos(data_th)));
	}
	// PR SR
	if(abs(data_lt) < 0.01 && data_m < 3.2 && data_m > 3.0) {
	  dataHist[0]->Fill(abs(cos(data_th)));
	}
      }
    }
  
  for(int i = 0; i < mEvt; i++)
    {
      treeM->GetEntry(i);
      if(mc_pt > 39 && mc_pt < 46 && abs(mc_lt) < 0.01 && mc_m < 3.2 && mc_m > 3.0) {
	mcHist->Fill(abs(cos(mc_th)), mc_pt);
      }
    }

  // plot costh histograms
  TCanvas *c = new TCanvas("", "", 700, 700);

  for(int i = 0; i < 2; i++) {
    dataHist[i]->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  }
  mcHist->GetXaxis()->SetTitle("|cos#theta_{HX}|");
    
  // get the ratio histogram for each bin
  TH1D **ratioHist = new TH1D*[2];

  ratioHist[0] = new TH1D("ratioH_PR", "2018 Data/MC (PR)", 20, 0, 1.);
  ratioHist[0] = (TH1D*)dataHist[0]->Clone("ratioH_PR");
  ratioHist[0]->Sumw2();
  ratioHist[0]->Divide(mcHist);
  ratioHist[0]->SetTitle("2018 Data/MC (PR)");
  
  ratioHist[1] = new TH1D("ratioH_NP", "2018 Data/MC (NP)", 20, 0, 1.);
  ratioHist[1] = (TH1D*)dataHist[1]->Clone("ratioH_NP");
  ratioHist[1]->Sumw2();
  ratioHist[1]->Divide(mcHist);
  ratioHist[1]->SetTitle("2018 Data/MC (NP)");

  dataHist[0]->SetTitle("Data |cos#theta| (39 < p_{T} < 46 GeV)");
  dataHist[0]->SetStats(0);
  dataHist[0]->SetLineColor(kBlack);
  dataHist[0]->SetMarkerColor(kBlack);
  dataHist[0]->SetMinimum(0);
  dataHist[0]->SetMaximum(dataHist[0]->GetBinContent(1)*1.5);
  dataHist[0]->Draw("error");

  dataHist[1]->SetLineColor(kBlue);
  dataHist[1]->SetMarkerColor(kBlue);
  dataHist[1]->Draw("error same");

  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(dataHist[0], "PR SR", "pl");
  leg->AddEntry(dataHist[1], "NP", "pl");
  leg->Draw();

  c->SaveAs("plots/costh_data.pdf");
  c->Clear();

  ratioHist[0]->SetTitle("Data/MC |cos#theta| (39 < p_{T} < 46 GeV)");
  ratioHist[0]->SetStats(0);
  ratioHist[0]->SetLineColor(kBlack);
  ratioHist[0]->SetMarkerColor(kBlack);
  ratioHist[0]->SetMinimum(0);
  ratioHist[0]->SetMaximum(ratioHist[0]->GetBinContent(1)*1.5);
  ratioHist[0]->Draw("error");

  ratioHist[1]->SetLineColor(kBlue);
  ratioHist[1]->SetMarkerColor(kBlue);
  ratioHist[1]->Draw("same");

  leg->Draw();

  c->SaveAs("plots/costh_ratio.pdf");
  c->Clear();

}
