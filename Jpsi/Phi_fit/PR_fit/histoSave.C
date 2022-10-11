// code to get the 2d data/mc ratio hist (pr, np)
// saves data, mc, ratio (normal and |costh|)

void histoSave()
{
  // fine pT binning - to be rebinned after background subtraction
  const int nPtBins = 17;
  double ptBins[nPtBins+1];
  for(int i = 0; i < 7; i++) ptBins[i] = 25 + 3.*i;
  for(int i = 0; i < 6; i++) ptBins[i+7] = 46 + 5.*i;
  for(int i = 0; i < 3; i++) ptBins[i+13] = 76 + 8.*i;
  for(int i = 0; i < 2; i++) ptBins[i+16] = 100 + 20.*i;
  for(int i=0; i<nPtBins+1; i++) cout << ptBins[i] << ",";
  cout << endl;

  // histograms for data (peak and NP) and MC
  TH2D *dataHist_ab = new TH2D("dataH_ab", "Full Data (PR)", 20, -180, 180, nPtBins, ptBins);
  TH2D *NPHist_ab = new TH2D("NPH_ab", "Full Data (NP)", 20, -180, 180, nPtBins, ptBins);
  TH2D *mcHist_ab = new TH2D("mcH_ab", "Full MC", 20, -180, 180, nPtBins, ptBins);

  // open files and read TTrees
  TFile *fin = new TFile("../../Store_data_codes/dataS_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/MCS_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../../Store_data_codes/MChS_cos.root");
  TTree *treeM2 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("../../Store_data_codes/MCvhS_cos.root");
  TTree *treeM3 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries();
  
  // definitions to store data and MC events
  Double_t data_phi, data_pt, data_lt, data_m, data_y;
  Double_t mc_phi, mc_pt, mc_lt, mc_m, mc_y;
  
  treeD->SetBranchAddress("phi", &data_phi);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  treeM1->SetBranchAddress("phi", &mc_phi);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("lt", &mc_lt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && data_m > 3.0 && data_m < 3.2 && abs(data_y) < 1.2) {
	if(abs(data_lt) < 0.005 ) {
	  dataHist_ab->Fill(data_phi, data_pt);
	}
	else if(data_lt > 0.01 && data_lt < 0.05 ) {
	  NPHist_ab->Fill(data_phi, data_pt);
	}
      }
    }
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if(mc_pt > ptBins[0] && mc_pt < 46 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {

	mcHist_ab->Fill(mc_phi, mc_pt);
      }
    }

  treeM2->SetBranchAddress("phi", &mc_phi);
  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Rap", &mc_y);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("lt", &mc_lt);

  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
      if(mc_pt > 46 && mc_pt < 66 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {

	mcHist_ab->Fill(mc_phi, mc_pt);
	
      }
    }

  treeM3->SetBranchAddress("phi", &mc_phi);
  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Rap", &mc_y);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("lt", &mc_lt);
  
  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);
      if(mc_pt > 66 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {

	mcHist_ab->Fill(mc_phi, mc_pt);
	
      }
    }
  
  dataHist_ab->SetStats(0);
  dataHist_ab->GetXaxis()->SetTitle("#phi_{HX}");
  dataHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  NPHist_ab->SetStats(0);
  NPHist_ab->GetXaxis()->SetTitle("#phi_{HX}");
  NPHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  mcHist_ab->SetStats(0);
  mcHist_ab->GetXaxis()->SetTitle("#phi_{HX}");
  mcHist_ab->GetYaxis()->SetTitle("p_{T} (GeV)");

  TH2D *ratioHist_ab = new TH2D("ratioH_ab", "Full Data/MC", 20, -180, 180, nPtBins, ptBins);
  ratioHist_ab = (TH2D*)dataHist_ab->Clone("ratioH_ab");
  ratioHist_ab->Sumw2();
  ratioHist_ab->Divide(mcHist_ab);
  ratioHist_ab->SetTitle("Full Data/MC");

  TH2D *ratNPHist_ab = new TH2D("ratNPH_ab", "Full NP/MC", 20, -180, 180, nPtBins, ptBins);
  ratNPHist_ab = (TH2D*)NPHist_ab->Clone("ratNPH_ab");
  ratNPHist_ab->Sumw2();
  ratNPHist_ab->Divide(mcHist_ab);
  ratNPHist_ab->SetTitle("Full NP/MC");

  TFile *outfile = new TFile("files/histoStore.root", "recreate");
  dataHist_ab->Write();
  NPHist_ab->Write();
  mcHist_ab->Write();
  ratioHist_ab->Write();
  ratNPHist_ab->Write();
  outfile->Close();
}
