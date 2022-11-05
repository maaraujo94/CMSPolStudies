// macro to plot the mass, rapidity and pT dists of data and MC
// also plots the lifetime of data

void store_ANdists()
{
  // PART 1 : creating the histograms

  // four pT dists: PRSR data + the 4 MC
  TH1D **h_pT = new TH1D*[5]; 
  for(int i = 0; i < 5; i++)
    h_pT[i] = new TH1D(Form("h_pT%d", i), Form("p_{T} distributions"), 100, 0, 200);

  // 8 y dists: (PRSR data + MC) over 4 pT regions
  TH1D **h_y = new TH1D*[8]; 
  for(int i = 0; i < 8; i++)
    h_y[i] = new TH1D(Form("h_y%d", i), Form("y distributions"), 60, -1.5, 1.5);

  // 9 M dists: (PR data + MC) over 4 pT regions, full pT data
  TH1D **h_m = new TH1D*[9]; 
  for(int i = 0; i < 9; i++)
    h_m[i] = new TH1D(Form("h_m%d", i), Form("M distributions"), 80, 2.9, 3.3);

  // 4 lifetime dists: data in 5 pT regions
  // in micron!
  TH1D **h_lt = new TH1D*[5];
  for(int i = 0; i < 5; i++) 
    h_lt[i] = new TH1D(Form("h_lt%d", i), "c#tau distribution", 75, -150, 600);

  // 12 costh dists: PR data + NP data + MC over 5 pT regions
  TH1D **h_cos = new TH1D*[15]; 
  for(int i = 0; i < 15; i++)
    h_cos[i] = new TH1D(Form("h_cos%d", i), Form("cos#theta distributions"), 20, 0, 1.);
  
  // PART 2 : open files and read TTrees
  TFile *fin = new TFile("../../Store_data_codes/dataS_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("../../Store_data_codes/MCS_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin2a = new TFile("../../Store_data_codes/MCmS_cos.root");
  TTree *treeM1a = (TTree*)fin2a->Get("MC_cos");
  TFile *fin3 = new TFile("../../Store_data_codes/MChS_cos.root");
  TTree *treeM2 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("../../Store_data_codes/MCvhS_cos.root");
  TTree *treeM3 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m1aEvt = treeM1a->GetEntries();
  int m2Evt = treeM2->GetEntries();
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

  treeM1a->SetBranchAddress("theta", &mc_th);
  treeM1a->SetBranchAddress("dimPt", &mc_pt);
  treeM1a->SetBranchAddress("Rap", &mc_y);
  treeM1a->SetBranchAddress("Mass", &mc_m);
  treeM1a->SetBranchAddress("lt", &mc_lt);

  treeM2->SetBranchAddress("theta", &mc_th);
  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Rap", &mc_y);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("lt", &mc_lt);
  
  treeM3->SetBranchAddress("theta", &mc_th);
  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Rap", &mc_y);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("lt", &mc_lt);

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);

      // fill the pT histo
      if(data_m > 3.0 && data_m < 3.2 && abs(data_lt) < 0.005 && abs(data_y) < 1.2) {
	h_pT[0]->Fill(data_pt);
      }

      // fill the y histo
      if(data_m > 3.0 && data_m < 3.2 && abs(data_lt) < 0.005) {
	if(data_pt > 25 && data_pt < 45) 
	  h_y[0]->Fill(data_y);
	else if(data_pt > 45 && data_pt < 50) 
	  h_y[1]->Fill(data_y);
	else if(data_pt > 50 && data_pt < 70) 
	  h_y[2]->Fill(data_y);
	else if(data_pt > 70 && data_pt < 120) 
	  h_y[3]->Fill(data_y);
      }
      
      // fill the M histo
      if(abs(data_lt) < 0.005 && abs(data_y) < 1.2) {
	if(data_pt > 25 && data_pt < 120)
	  h_m[8]->Fill(data_m);
	if(data_pt > 25 && data_pt < 45) 
	  h_m[0]->Fill(data_m);
	else if(data_pt > 45 && data_pt < 50) 
	  h_m[1]->Fill(data_m);
	else if(data_pt > 50 && data_pt < 70) 
	  h_m[2]->Fill(data_m);
	else if(data_pt > 70 && data_pt < 120) 
	  h_m[3]->Fill(data_m);
      }
      
      // fill the lt histo - in MICRON
      if(data_m > 3.0 && data_m < 3.2 && abs(data_y) < 1.2) {
	if(data_pt > 25 && data_pt < 45) 
	  h_lt[0]->Fill(data_lt*1e4);
	else if(data_pt > 45 && data_pt < 50) 
	  h_lt[1]->Fill(data_lt*1e4);
	else if(data_pt > 50 && data_pt < 70) 
	  h_lt[2]->Fill(data_lt*1e4);
	else if(data_pt > 70 && data_pt < 120) 
	  h_lt[3]->Fill(data_lt*1e4);
	if(data_pt > 25 && data_pt < 120) 
	  h_lt[4]->Fill(data_lt*1e4);
      }

      // fill the costh histo - PRSR
      if(data_m > 3.0 && data_m < 3.2 && abs(data_y) < 1.2 && abs(data_lt) < 0.005) {
	if(data_pt > 25 && data_pt < 45) 
	  h_cos[0]->Fill(abs(cos(data_th)));
	else if(data_pt > 45 && data_pt < 50) 
	  h_cos[1]->Fill(abs(cos(data_th)));
	else if(data_pt > 50 && data_pt < 70) 
	  h_cos[2]->Fill(abs(cos(data_th)));
	else if(data_pt > 70 && data_pt < 120) 
	  h_cos[3]->Fill(abs(cos(data_th)));
	if(data_pt > 25 && data_pt < 120) 
	  h_cos[4]->Fill(abs(cos(data_th)));
      }
      // fill the costh histo - NPSR
      if(data_m > 3.0 && data_m < 3.2 && abs(data_y) < 1.2 && data_lt > 0.01 && data_lt < 0.05) {
	if(data_pt > 25 && data_pt < 45) 
	  h_cos[5]->Fill(abs(cos(data_th)));
	else if(data_pt > 45 && data_pt < 50) 
	  h_cos[6]->Fill(abs(cos(data_th)));
	else if(data_pt > 50 && data_pt < 70) 
	  h_cos[7]->Fill(abs(cos(data_th)));
	else if(data_pt > 70 && data_pt < 120) 
	  h_cos[8]->Fill(abs(cos(data_th)));
	if(data_pt > 25 && data_pt < 120) 
	  h_cos[9]->Fill(abs(cos(data_th)));
      }

    }

  cout << "data filled" << endl;
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);

      // fill the pT histo
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	h_pT[1]->Fill(mc_pt);
      }

      // fill the y histo
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005) {
	if(mc_pt > 25 && mc_pt < 45) 
	  h_y[4]->Fill(mc_y);
      }

      // fill the m histo
      if(abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	if(mc_pt > 25 && mc_pt < 45) 
	  h_m[4]->Fill(mc_m);
      }

      // fill the costh histo
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_y) < 1.2 && abs(mc_lt) < 0.005) {
	if(mc_pt > 25 && mc_pt < 45) {
	  h_cos[10]->Fill(abs(cos(mc_th)));
	  h_cos[14]->Fill(abs(cos(mc_th)));
	}
      }

    }

  cout << "low-pT MC filled" << endl;

  for(int i = 0; i < m1aEvt; i++)
    {
      treeM1a->GetEntry(i);

      // fill the pT histo
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	h_pT[2]->Fill(mc_pt);
      }
      
      // fill the y histo
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005) {
	if(mc_pt > 45 && mc_pt < 50) 
	  h_y[5]->Fill(mc_y);
      }
      
      // fill the m histo
      if(abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	if(mc_pt > 45 && mc_pt < 50) 
	  h_m[5]->Fill(mc_m);
      }

      // fill the costh histo
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_y) < 1.2 && abs(mc_lt) < 0.005) {
	if(mc_pt > 45 && mc_pt < 50) {
	  h_cos[11]->Fill(abs(cos(mc_th)));
	  h_cos[14]->Fill(abs(cos(mc_th)));
	}
      }

    }

  cout << "mid-pT MC filled" << endl;

  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);

      // fill the pT histo
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	h_pT[3]->Fill(mc_pt);
      }

      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005) {
	if(mc_pt > 50 && mc_pt < 70) 
	  h_y[6]->Fill(mc_y);
      }
      
      // fill the m histo
      if(abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	if(mc_pt > 50 && mc_pt < 70) 
	  h_m[6]->Fill(mc_m);
      }
      
      // fill the costh histo
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_y) < 1.2 && abs(mc_lt) < 0.005) {
	if(mc_pt > 50 && mc_pt < 70) {
	  h_cos[12]->Fill(abs(cos(mc_th)));
	  h_cos[14]->Fill(abs(cos(mc_th)));
	}
      }

    }

  cout << "high-pT MC filled" << endl;

  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);

      // fill the pT histo
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	h_pT[4]->Fill(mc_pt);
      }

      // fill the y histo
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005) {
	if(mc_pt > 70 && mc_pt < 120) 
	  h_y[7]->Fill(mc_y);
      }

      // fill the m histo
      if(abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	if(mc_pt > 70 && mc_pt < 120) 
	  h_m[7]->Fill(mc_m);
      }

      // fill the costh histo
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_y) < 1.2 && abs(mc_lt) < 0.005) {
	if(mc_pt > 70 && mc_pt < 120) {
	  h_cos[13]->Fill(abs(cos(mc_th)));
	  h_cos[14]->Fill(abs(cos(mc_th)));
	}
      }

    }

  cout << "highest-pT MC filled" << endl;
  
  fin->Close();
  fin2->Close();
  fin2a->Close();
  fin3->Close();
  fin4->Close();

  TFile *fout = new TFile("files/store_ANdists.root", "recreate");

  // store the pT dists
  string lbl_pt[] = {"Data", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC"};
  for(int i = 0; i < 5; i++) {
    h_pT[i]->Write(Form("h_pT_%s", lbl_pt[i].c_str()));
  }

  // store the y dists
  string lbl_y[] = {"lowPtData", "midPtData", "highPtData", "highestPtData", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC"};
  for(int i = 0; i < 8; i++) {
    h_y[i]->Write(Form("h_y_%s", lbl_y[i].c_str()));
  }

  // store the M dists
  string lbl_m[] = {"lowPtData", "midPtData", "highPtData", "highestPtData", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC", "fullData"};
  for(int i = 0; i < 9; i++) {
    h_m[i]->Write(Form("h_m_%s", lbl_m[i].c_str()));
  }

  // store the lifetime dists
  string lbl_lt[] = {"lowPt", "midPt", "highPt", "highestPt", "full"};
  for(int i = 0; i < 5; i++) {
    h_lt[i]->Write(Form("h_lt_%s", lbl_lt[i].c_str()));
  }

  // store the cos dists
  string lbl_cos[] = {"lowPtPR", "midPtPR", "highPtPR", "highestPtPR", "fullPR", "lowPtNP", "midPtNP", "highPtNP", "highestPtNP", "fullNP", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC", "fullMC"};
  for(int i = 0; i < 15; i++) {
    h_cos[i]->Write(Form("h_cos_%s", lbl_cos[i].c_str()));
  }

  cout << "all histos stored" << endl;
  
  fout->Close();
}
