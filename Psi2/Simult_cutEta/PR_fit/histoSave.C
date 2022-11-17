// code to get the 2d data/mc ratio hist (pr, np)
// saves data, mc, ratio (normal and |costh|)

#import "../ptbins.C"
#import "../etacut.C"

void histoSave()
{
  // histograms for data (PR and NP, peak and sideband) and MC
  TH2D *PRHist  = new TH2D("PRH",  "Data (PR Peak)", 20, 0, 1., nPtBins, ptBins);
  TH2D *PRLHist = new TH2D("PRLH", "Data (PR LSB)",  20, 0, 1., nPtBins, ptBins);
  TH2D *PRRHist = new TH2D("PRRH", "Data (PR RSB)",  20, 0, 1., nPtBins, ptBins);
  TH2D *NPHist  = new TH2D("NPH",  "Data (NP Peak)", 20, 0, 1., nPtBins, ptBins);
  TH2D *NPLHist = new TH2D("NPLH", "Data (NP LSB)",  20, 0, 1., nPtBins, ptBins);
  TH2D *NPRHist = new TH2D("NPRH", "Data (NP RSB)",  20, 0, 1., nPtBins, ptBins);
  TH2D *MCHist  = new TH2D("MCH",  "MC",             20, 0, 1., nPtBins, ptBins);

  // open files and read TTrees
  TFile *fin = new TFile("../../Store_data_codes/dataS_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/MCmS_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();

  double m_min[] = {3.4, 3.57, 3.82};
  double m_max[] = {3.52, 3.81, 4.0};

  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m, data_y;
  Double_t mc_th, mc_pt, mc_lt, mc_m, mc_y;
  double mPEta, mMEta;
    
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  treeD->SetBranchAddress("muonPEta", &mPEta);
  treeD->SetBranchAddress("muonMEta", &mMEta);

  treeM1->SetBranchAddress("theta", &mc_th);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("lt", &mc_lt);
  treeM1->SetBranchAddress("muonPEta", &mPEta);
  treeM1->SetBranchAddress("muonMEta", &mMEta);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      // pt and y conditions
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_y) < 1.2) {
	// eta cut
	if((abs(mPEta) < eta_lo || abs(mPEta) > eta_hi) && (abs(mMEta) < eta_lo || abs(mMEta) > eta_hi)) {

	  // PR peak and sidebands
	  if(abs(data_lt) < 0.005) {
	    if(data_m > m_min[1] && data_m < m_max[1])
	      PRHist->Fill(abs(cos(data_th)), data_pt);
	    else if(data_m < m_max[0] && data_m > m_min[0])
	      PRLHist->Fill(abs(cos(data_th)), data_pt);
	    else if(data_m < m_max[2] && data_m > m_min[2])
	      PRRHist->Fill(abs(cos(data_th)), data_pt);
	  }
	  // NP peak and sidebands
	  else if(data_lt > 0.01 && data_lt < 0.05 ) {
	    if(data_m > m_min[1] && data_m < m_max[1])
	      NPHist->Fill(abs(cos(data_th)), data_pt);
	    else if(data_m < m_max[0] && data_m > m_min[0])
	      NPLHist->Fill(abs(cos(data_th)), data_pt);
	    else if(data_m < m_max[2] && data_m > m_min[2])
	      NPRHist->Fill(abs(cos(data_th)), data_pt);
	  }
	}
      }
    }
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if((abs(mPEta) < eta_lo || abs(mPEta) > eta_hi) && (abs(mMEta) < eta_lo || abs(mMEta) > eta_hi))
	if(mc_pt > ptBins[0] && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > m_min[1] && mc_m < m_max[1]) {
	  MCHist->Fill(abs(cos(mc_th)), mc_pt);
	}
    }
  
  // get ratios
  TH2D *rPRHist  = new TH2D();
  TH2D *rPRLHist = new TH2D();
  TH2D *rPRRHist = new TH2D();
  TH2D *rNPHist  = new TH2D();
  TH2D *rNPLHist = new TH2D();
  TH2D *rNPRHist = new TH2D();

  rPRHist = (TH2D*)PRHist->Clone("rPRH");
  rPRHist->Sumw2();
  rPRHist->Divide(MCHist);
  rPRHist->SetTitle("Data/MC (PR Peak)");

  rPRLHist = (TH2D*)PRLHist->Clone("rPRLH");
  rPRLHist->Sumw2();
  rPRLHist->Divide(MCHist);
  rPRLHist->SetTitle("Data/MC (PR LSB)");

  rPRRHist = (TH2D*)PRRHist->Clone("rPRRH");
  rPRRHist->Sumw2();
  rPRRHist->Divide(MCHist);
  rPRRHist->SetTitle("Data/MC (PR RSB)");
  
  rNPHist = (TH2D*)NPHist->Clone("rNPH");
  rNPHist->Sumw2();
  rNPHist->Divide(MCHist);
  rNPHist->SetTitle("Data/MC (NP Peak)");

  rNPLHist = (TH2D*)NPLHist->Clone("rNPLH");
  rNPLHist->Sumw2();
  rNPLHist->Divide(MCHist);
  rNPLHist->SetTitle("Data/MC (NP LSB)");

  rNPRHist = (TH2D*)NPRHist->Clone("rNPRH");
  rNPRHist->Sumw2();
  rNPRHist->Divide(MCHist);
  rNPRHist->SetTitle("Data/MC (NP RSB)");

  // store data, MC and data/MC ratios
  TFile *outfile = new TFile("files/histoStore.root", "recreate");
  PRHist->Write();
  PRLHist->Write();
  PRRHist->Write();
  NPHist->Write();
  NPLHist->Write();
  NPRHist->Write();
  MCHist->Write();
  rPRHist->Write();
  rPRLHist->Write();
  rPRRHist->Write();
  rNPHist->Write();
  rNPLHist->Write();
  rNPRHist->Write();
  outfile->Close();
  
  cout << Form("%.0f data (PR) events, %.0f data (NP) events and %.0f MC events", PRHist->GetEntries(), NPHist->GetEntries(), MCHist->GetEntries()) << endl;
}
