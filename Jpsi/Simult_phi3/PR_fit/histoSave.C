#import "../phiCode.C"

// code to get the 2d data/mc ratio hist (pr, np, peak and sidebands)
// saves data, mc, ratio

#import "../ptbins.C"

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
  TFile *fin1 = new TFile("../../Store_data_codes/MCS_cos.root");
  TTree *treeM1 = (TTree*)fin1->Get("MC_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/MCmS_cos.root");
  TTree *treeM2 = (TTree*)fin2->Get("MC_cos");
  TFile *fin3 = new TFile("../../Store_data_codes/MChS_cos.root");
  TTree *treeM3 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("../../Store_data_codes/MCvhS_cos.root");
  TTree *treeM4 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries();
  int m4Evt = treeM4->GetEntries(); 
 
  // definitions to store data and MC events
  Double_t data_th, data_pt, data_lt, data_m, data_y;
double phi, w_phi;
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
  treeM1->SetBranchAddress("phi", &phi);

  treeM2->SetBranchAddress("theta", &mc_th);
  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Rap", &mc_y);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("lt", &mc_lt);
  treeM2->SetBranchAddress("phi", &phi);

  treeM3->SetBranchAddress("theta", &mc_th);
  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Rap", &mc_y);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("lt", &mc_lt);
  treeM3->SetBranchAddress("phi", &phi);

  treeM4->SetBranchAddress("theta", &mc_th);
  treeM4->SetBranchAddress("dimPt", &mc_pt);
  treeM4->SetBranchAddress("Rap", &mc_y);
  treeM4->SetBranchAddress("Mass", &mc_m);
  treeM4->SetBranchAddress("lt", &mc_lt);
  treeM4->SetBranchAddress("phi", &phi);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      // pt and y conditions
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_y) < 1.2) {
	// PR peak and sidebands
	if(abs(data_lt) < 0.005) {
	  if(data_m > 3.0 && data_m < 3.2)
	    PRHist->Fill(abs(cos(data_th)), data_pt);
	  else if(data_m < 2.95 && data_m > 2.92)
	    PRLHist->Fill(abs(cos(data_th)), data_pt);
	  else if(data_m < 3.28 && data_m > 3.21)
	    PRRHist->Fill(abs(cos(data_th)), data_pt);
	}
	// NP peak and sidebands
	else if(data_lt > 0.01 && data_lt < 0.08 ) {
	  if(data_m > 3.0 && data_m < 3.2)
	    NPHist->Fill(abs(cos(data_th)), data_pt);
	  else if(data_m < 2.95 && data_m > 2.92)
	    NPLHist->Fill(abs(cos(data_th)), data_pt);
	  else if(data_m < 3.28 && data_m > 3.21)
	    NPRHist->Fill(abs(cos(data_th)), data_pt);
	}
      }
    }

  // MC sample 1
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
	  w_phi = f_phi(phi);
      if(mc_pt > ptBins[0] && mc_pt < 45 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	MCHist->Fill(abs(cos(mc_th)), mc_pt, w_phi);
      }
    }

  // MC sample 2
  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
	  w_phi = f_phi(phi);
      if(mc_pt > 45 && mc_pt < 50 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	MCHist->Fill(abs(cos(mc_th)), mc_pt, w_phi);
      }
    }

  // MC sample 3
  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);
	  w_phi = f_phi(phi);
      if(mc_pt > 50 && mc_pt < 70 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	MCHist->Fill(abs(cos(mc_th)), mc_pt, w_phi);
      }
    }

  // MC sample 4
  for(int i = 0; i < m4Evt; i++)
    {
      treeM4->GetEntry(i);
	  w_phi = f_phi(phi);
      if(mc_pt > 70 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	MCHist->Fill(abs(cos(mc_th)), mc_pt, w_phi);
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
