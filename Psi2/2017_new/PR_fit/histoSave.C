// code to get the 2d data/mc ratio hist (pr, np)
// saves data, mc, ratio (normal and |costh|)

void histoSave()
{
  // fine pT binning - to be rebinned after background subtraction
  const int nPtBins = 7;
  double ptBins[nPtBins+1];
  for(int i=0; i<3; i++) ptBins[i] = 7.*i+25.;
  for(int i=0; i<4; i++) ptBins[i+3] = 46.+10.*i;
  ptBins[7] = 120;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  double m_min = 3.57, m_max = 3.81;
 
  // histograms for data (peak and NP) and MC
  TH2D *dataHist = new TH2D("dataH", "2017 Data (PR)", 40, -1., 1., nPtBins, ptBins);
  TH2D *NPHist = new TH2D("NPH", "2017 Data (NP)", 40, -1., 1., nPtBins, ptBins);
  TH2D *mcHist = new TH2D("mcH", "2017 MC", 40, -1., 1., nPtBins, ptBins);
  
  TH2D *dataHist_ab = new TH2D("dataH_ab", "2017 Data (PR)", 20, 0, 1., nPtBins, ptBins);
  TH2D *NPHist_ab = new TH2D("NPH_ab", "2017 Data (NP)", 20, 0, 1., nPtBins, ptBins);
  TH2D *mcHist_ab = new TH2D("mcH_ab", "2017 MC", 20, 0, 1., nPtBins, ptBins);

  // open files and read TTrees
  TFile *fin = new TFile("../../Store_data_codes/data17_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/MC17_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin4 = new TFile("../../Store_data_codes/MCvh17_cos.root");
  TTree *treeM3 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m3Evt = treeM3->GetEntries();
  
  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m, data_y;
  Double_t mc_th, mc_pt, mc_lt, mc_m, mc_y;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  treeM1->SetBranchAddress("theta", &mc_th);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("lt", &mc_lt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && data_m > m_min && data_m < m_max && abs(data_y) < 1.2) {
	if(abs(data_lt) < 0.005 ) {
	  dataHist->Fill(cos(data_th), data_pt);
	  dataHist_ab->Fill(abs(cos(data_th)), data_pt);
	}
	else if(data_lt > 0.01 && data_lt < 0.05 ) {
	  NPHist->Fill(cos(data_th), data_pt);
	  NPHist_ab->Fill(abs(cos(data_th)), data_pt);
	}
      }
    }
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if(mc_pt > ptBins[0] && mc_pt < 46 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > m_min && mc_m < m_max) {

	mcHist->Fill(cos(mc_th), mc_pt);
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
      }
    }

  treeM3->SetBranchAddress("theta", &mc_th);
  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Rap", &mc_y);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("lt", &mc_lt);
  
  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);
      if(mc_pt > 46 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > m_min && mc_m < m_max) {

	mcHist->Fill(cos(mc_th), mc_pt);
	mcHist_ab->Fill(abs(cos(mc_th)), mc_pt);
	
      }
    }
  
  dataHist_ab->SetStats(0);
  dataHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  dataHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  NPHist_ab->SetStats(0);
  NPHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  NPHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  mcHist_ab->SetStats(0);
  mcHist_ab->GetXaxis()->SetTitle("|cos#theta_{HX}|");
  mcHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  TH2D *ratioHist_ab = new TH2D("ratioH_ab", "2017 Data/MC", 20, 0, 1., nPtBins, ptBins);
  ratioHist_ab = (TH2D*)dataHist_ab->Clone("ratioH_ab");
  ratioHist_ab->Sumw2();
  ratioHist_ab->Divide(mcHist_ab);
  ratioHist_ab->SetTitle("2017 Data/MC");

  TH2D *ratNPHist_ab = new TH2D("ratNPH_ab", "2017 NP/MC", 20, 0, 1., nPtBins, ptBins);
  ratNPHist_ab = (TH2D*)NPHist_ab->Clone("ratNPH_ab");
  ratNPHist_ab->Sumw2();
  ratNPHist_ab->Divide(mcHist_ab);
  ratNPHist_ab->SetTitle("2017 NP/MC");

  TFile *outfile = new TFile("files/histoStore.root", "recreate");
  dataHist->Write();
  dataHist_ab->Write();
  NPHist->Write();
  NPHist_ab->Write();
  mcHist->Write();
  mcHist_ab->Write();
  ratioHist_ab->Write();
  ratNPHist_ab->Write();
  outfile->Close();
  
  cout << Form("%.0f data (PR) events, %.0f data (NP) events and %.0f MC events", dataHist->GetEntries(), NPHist->GetEntries(), mcHist->GetEntries()) << endl;

  // split between pT ranges - must set by hand
  
  double pt_min[] = {dataHist->GetYaxis()->GetBinLowEdge(1), dataHist->GetYaxis()->GetBinLowEdge(4)};
  double pt_max[] = {dataHist->GetYaxis()->GetBinUpEdge(3), dataHist->GetYaxis()->GetBinUpEdge(7)};

  cout << pt_min[0] << " - " << pt_max[0] << endl;
  cout << pt_min[1] << " - " << pt_max[1] << endl;
  
  
  cout << Form("%.0f data (PR) events, %.0f data (NP) events and %.0f MC events in pT range 1", dataHist->Integral(1, 40, 1, 3), NPHist->Integral(1, 40, 1, 3), mcHist->Integral(1, 40, 1, 3)) << endl;
  cout << Form("%.0f data (PR) events, %.0f data (NP) events and %.0f MC events in pT range 2", dataHist->Integral(1, 40, 4, 7), NPHist->Integral(1, 40, 4, 7), mcHist->Integral(1, 40, 4, 7)) << endl;
}
