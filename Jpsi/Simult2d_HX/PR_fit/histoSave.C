// code to get the 2d data/mc ratio hist (pr, np, peak and sidebands)
// saves data, mc, ratio

#import "../ptbins.C"

double *ang_fold(double costh, double phi)
{
  static double ang[2];
  if(phi < -90 && phi > -180) {
    ang[0] = -costh;
    ang[1] = phi+180;
  }
  if(phi < 0 && phi > -90) {
    ang[0] = costh;
    ang[1] = -phi;
  }
  if(phi < 90 && phi > 0) {
    ang[0] = costh;
    ang[1] = phi;
  }
  if(phi < 180 && phi > 90) {
    ang[0] = -costh;
    ang[1] = 180-phi;
  }
  return ang;
}

void histoSave()
{
  // histograms for data (PR and NP, peak and sideband) and MC
  // one per pT so they can be 2d angular
  
  TH2D **PRHist  = new TH2D*[nPtBins];
  TH2D **PRLHist  = new TH2D*[nPtBins];
  TH2D **PRRHist  = new TH2D*[nPtBins];
  TH2D **NPHist  = new TH2D*[nPtBins];
  TH2D **NPLHist  = new TH2D*[nPtBins];
  TH2D **NPRHist  = new TH2D*[nPtBins];
  TH2D **MCHist  = new TH2D*[nPtBins];
  for(int i = 0; i < nPtBins; i++) {
    PRHist[i] = new TH2D(Form("PRH_%d", i),  "Data (PR Peak)", 20, 0, 1, 20, 0, 90);
    PRLHist[i] = new TH2D(Form("PRLH_%d", i), "Data (PR LSB)",  20, 0, 1, 20, 0, 90);
    PRRHist[i] = new TH2D(Form("PRRH_%d", i), "Data (PR RSB)",  20, 0, 1, 20, 0, 90);
    NPHist[i] = new TH2D(Form("NPH_%d", i),  "Data (NP Peak)", 20, 0, 1, 20, 0, 90);
    NPLHist[i] = new TH2D(Form("NPLH_%d", i), "Data (NP LSB)",  20, 0, 1, 20, 0, 90);
    NPRHist[i] = new TH2D(Form("NPRH_%d", i), "Data (NP RSB)",  20, 0, 1, 20, 0, 90);
    MCHist[i] = new TH2D(Form("MCH_%d", i), "MC",  20, 0, 1, 20, 0, 90);
  }

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
  Double_t data_th, data_ph, data_pt, data_lt, data_m, data_y;
  Double_t mc_th, mc_ph, mc_pt, mc_lt, mc_m, mc_y;
  double *ang;
  
  treeD->SetBranchAddress("theta", &data_th);
  treeD->SetBranchAddress("phi", &data_ph);
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  
  treeM1->SetBranchAddress("theta", &mc_th);
  treeM1->SetBranchAddress("phi", &mc_ph);
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("lt", &mc_lt);

  treeM2->SetBranchAddress("theta", &mc_th);
  treeM2->SetBranchAddress("phi", &mc_ph);
  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Rap", &mc_y);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("lt", &mc_lt);

  treeM3->SetBranchAddress("theta", &mc_th);
  treeM3->SetBranchAddress("phi", &mc_ph);
  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Rap", &mc_y);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("lt", &mc_lt);

  treeM4->SetBranchAddress("theta", &mc_th);
  treeM4->SetBranchAddress("phi", &mc_ph);
  treeM4->SetBranchAddress("dimPt", &mc_pt);
  treeM4->SetBranchAddress("Rap", &mc_y);
  treeM4->SetBranchAddress("Mass", &mc_m);
  treeM4->SetBranchAddress("lt", &mc_lt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);
      // pt and y conditions
      if(data_pt > ptBins[0] && data_pt < ptBins[nPtBins] && abs(data_y) < 1.2) {
	for(int j = 0; j < nPtBins; j++) {
	  if(data_pt > ptBins[j] && data_pt < ptBins[j+1]) {
	    // PR peak and sidebands
	    if(abs(data_lt) < 0.005) {
	      ang = ang_fold(cos(data_th), data_ph);
	      if(data_m > 3.0 && data_m < 3.2)
		PRHist[j]->Fill(ang[0], ang[1]);
	      else if(data_m < 2.95 && data_m > 2.92)
		PRLHist[j]->Fill(ang[0], ang[1]);
	      else if(data_m < 3.28 && data_m > 3.21)
		PRRHist[j]->Fill(ang[0], ang[1]);
	    }
	    // NP peak and sidebands
	    else if(data_lt > 0.01 && data_lt < 0.08 ) {
	      ang = ang_fold(cos(data_th), data_ph);
	      if(data_m > 3.0 && data_m < 3.2)
		NPHist[j]->Fill(ang[0], ang[1]);
	      else if(data_m < 2.95 && data_m > 2.92)
		NPLHist[j]->Fill(ang[0], ang[1]);
	      else if(data_m < 3.28 && data_m > 3.21)
		NPRHist[j]->Fill(ang[0], ang[1]);
	    }
	  }
	}
      }
    }

  // MC sample 1
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      if(mc_pt > ptBins[0] && mc_pt < 45 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	for(int j = 0; j < nPtBins; j++) {
	  if(mc_pt > ptBins[j] && mc_pt < ptBins[j+1]) {
	    ang = ang_fold(cos(mc_th), mc_ph);
	    MCHist[j]->Fill(ang[0], ang[1]);
	  }
	}
      }
    }

  // MC sample 2
  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);
      if(mc_pt > 45 && mc_pt < 50 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	for(int j = 0; j < nPtBins; j++) {
	  if(mc_pt > ptBins[j] && mc_pt < ptBins[j+1]) {
	    ang = ang_fold(cos(mc_th), mc_ph);
	    MCHist[j]->Fill(ang[0], ang[1]);
	  }
	}
      }
    }

  // MC sample 3
  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);
      if(mc_pt > 50 && mc_pt < 70 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	for(int j = 0; j < nPtBins; j++) {
	  if(mc_pt > ptBins[j] && mc_pt < ptBins[j+1]) {
	    ang = ang_fold(cos(mc_th), mc_ph);
	    MCHist[j]->Fill(ang[0], ang[1]);
	  }
	}
      }
    }

  // MC sample 4
  for(int i = 0; i < m4Evt; i++)
    {
      treeM4->GetEntry(i);
      if(mc_pt > 70 && mc_pt < ptBins[nPtBins] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_m > 3.0 && mc_m < 3.2) {
	for(int j = 0; j < nPtBins; j++) {
	  if(mc_pt > ptBins[j] && mc_pt < ptBins[j+1]) {
	    ang = ang_fold(cos(mc_th), mc_ph);
	    MCHist[j]->Fill(ang[0], ang[1]);
	  }
	}
      }
    }

  // get ratios
  TH2D **rPRHist  = new TH2D*[nPtBins];
  TH2D **rPRLHist = new TH2D*[nPtBins];
  TH2D **rPRRHist = new TH2D*[nPtBins];
  TH2D **rNPHist  = new TH2D*[nPtBins];
  TH2D **rNPLHist = new TH2D*[nPtBins];
  TH2D **rNPRHist = new TH2D*[nPtBins];

  for(int i = 0; i < nPtBins; i++) {
    rPRHist[i] = (TH2D*)PRHist[i]->Clone(Form("rPRH_%d", i));
    rPRHist[i]->Sumw2();
    rPRHist[i]->Divide(MCHist[i]);
    rPRHist[i]->SetTitle("Data/MC (PR Peak)");

    rPRLHist[i] = (TH2D*)PRLHist[i]->Clone(Form("rPRLH_%d", i));
    rPRLHist[i]->Sumw2();
    rPRLHist[i]->Divide(MCHist[i]);
    rPRLHist[i]->SetTitle("Data/MC (PR LSB)");

    rPRRHist[i] = (TH2D*)PRRHist[i]->Clone(Form("rPRRH_%d", i));
    rPRRHist[i]->Sumw2();
    rPRRHist[i]->Divide(MCHist[i]);
    rPRRHist[i]->SetTitle("Data/MC (PR RSB)");
  
    rNPHist[i] = (TH2D*)NPHist[i]->Clone(Form("rNPH_%d", i));
    rNPHist[i]->Sumw2();
    rNPHist[i]->Divide(MCHist[i]);
    rNPHist[i]->SetTitle("Data/MC (NP Peak)");

    rNPLHist[i] = (TH2D*)NPLHist[i]->Clone(Form("rNPLH_%d", i));
    rNPLHist[i]->Sumw2();
    rNPLHist[i]->Divide(MCHist[i]);
    rNPLHist[i]->SetTitle("Data/MC (NP LSB)");

    rNPRHist[i] = (TH2D*)NPRHist[i]->Clone(Form("rNPRH_%d", i));
    rNPRHist[i]->Sumw2();
    rNPRHist[i]->Divide(MCHist[i]);
    rNPRHist[i]->SetTitle("Data/MC (NP RSB)");
  }
  
  // store data, MC and data/MC ratios
  TFile *outfile = new TFile("files/histoStore.root", "recreate");
  for(int i = 0; i < nPtBins; i++) {
    PRHist[i]->Write();
    PRLHist[i]->Write();
    PRRHist[i]->Write();
    NPHist[i]->Write();
    NPLHist[i]->Write();
    NPRHist[i]->Write();
    MCHist[i]->Write();
    rPRHist[i]->Write();
    rPRLHist[i]->Write();
    rPRRHist[i]->Write();
    rNPHist[i]->Write();
    rNPLHist[i]->Write();
    rNPRHist[i]->Write();
  }
  outfile->Close();
}
