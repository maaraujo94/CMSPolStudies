// macro to plot the mass, rapidity and pT dists of data and MC
// also plots the lifetime of data

void MC_comp_store()
{
  // the pT intervals
  double pT_i[] = {41, 47, 67};
  double pT_f[] = {45, 51, 120};
  
  // PART 1 : creating the histograms
  // for all variables, there are 6 histos covering 3 intervals: low-pT, mid-pT, high-pT
  
  // dimuon pT
  TH1D **h_pT = new TH1D*[6]; 
  for(int i = 0; i < 6; i++)
    h_pT[i] = new TH1D(Form("h_pT%d", i), Form("p_{T} distributions"), 130, 0, 130);

  //dimuon y
  TH1D **h_y = new TH1D*[6]; 
  for(int i = 0; i < 6; i++)
    h_y[i] = new TH1D(Form("h_y%d", i), Form("y distributions"), 60, -1.5, 1.5);

  // dimuon mass
  TH1D **h_m = new TH1D*[6]; 
  for(int i = 0; i < 6; i++)
    h_m[i] = new TH1D(Form("h_m%d", i), Form("M distributions"), 80, 2.9, 3.3);

  // single muon pT
  TH1D **hP_pT = new TH1D*[6]; 
  TH1D **hM_pT = new TH1D*[6]; 
  for(int i = 0; i < 6; i++) {
    hP_pT[i] = new TH1D(Form("hP_pT%d", i), Form("p_{T}(#mu^{+}) distributions"), 130, 0, 130);
    hM_pT[i] = new TH1D(Form("hM_pT%d", i), Form("p_{T}(#mu^{-}) distributions"), 130, 0, 130);
  }
  
  //single muon eta
  TH1D **hP_eta = new TH1D*[6]; 
  TH1D **hM_eta = new TH1D*[6]; 
  for(int i = 0; i < 6; i++) {
    hP_eta[i] = new TH1D(Form("hP_eta%d", i), Form("#eta(#mu^{+}) distributions"), 60, -1.5, 1.5);
    hM_eta[i] = new TH1D(Form("hM_eta%d", i), Form("#eta(#mu^{-}) distributions"), 60, -1.5, 1.5);
  }
  
  // PART 2 : open files and read TTrees
  TFile *fin2 = new TFile("../../Store_data_codes/MC18_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin2a = new TFile("../../Store_data_codes/MCm18_cos.root");
  TTree *treeM1a = (TTree*)fin2a->Get("MC_cos");
  TFile *fin3 = new TFile("../../Store_data_codes/MCh18_cos.root");
  TTree *treeM2 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("../../Store_data_codes/MCvh18_cos.root");
  TTree *treeM3 = (TTree*)fin4->Get("MC_cos");
  
  int m1Evt = treeM1->GetEntries();
  int m1aEvt = treeM1a->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries();
  
  // definitions to store data and MC events
  Double_t mc_pt, mc_lt, mc_m, mc_y;
  double mMPt, mPPt, mMEta, mPEta;
  
  treeM1->SetBranchAddress("dimPt", &mc_pt);
  treeM1->SetBranchAddress("Rap", &mc_y);
  treeM1->SetBranchAddress("Mass", &mc_m);
  treeM1->SetBranchAddress("lt", &mc_lt);
  treeM1->SetBranchAddress("muonPPt", &mPPt);
  treeM1->SetBranchAddress("muonMPt", &mMPt);
  treeM1->SetBranchAddress("muonPEta", &mPEta);
  treeM1->SetBranchAddress("muonMEta", &mMEta);

  treeM1a->SetBranchAddress("dimPt", &mc_pt);
  treeM1a->SetBranchAddress("Rap", &mc_y);
  treeM1a->SetBranchAddress("Mass", &mc_m);
  treeM1a->SetBranchAddress("lt", &mc_lt);
  treeM1a->SetBranchAddress("muonPPt", &mPPt);
  treeM1a->SetBranchAddress("muonMPt", &mMPt);
  treeM1a->SetBranchAddress("muonPEta", &mPEta);
  treeM1a->SetBranchAddress("muonMEta", &mMEta);

  treeM2->SetBranchAddress("dimPt", &mc_pt);
  treeM2->SetBranchAddress("Rap", &mc_y);
  treeM2->SetBranchAddress("Mass", &mc_m);
  treeM2->SetBranchAddress("lt", &mc_lt);
  treeM2->SetBranchAddress("muonPPt", &mPPt);
  treeM2->SetBranchAddress("muonMPt", &mMPt);
  treeM2->SetBranchAddress("muonPEta", &mPEta);
  treeM2->SetBranchAddress("muonMEta", &mMEta);

  treeM3->SetBranchAddress("dimPt", &mc_pt);
  treeM3->SetBranchAddress("Rap", &mc_y);
  treeM3->SetBranchAddress("Mass", &mc_m);
  treeM3->SetBranchAddress("lt", &mc_lt);
  treeM3->SetBranchAddress("muonPPt", &mPPt);
  treeM3->SetBranchAddress("muonMPt", &mMPt);
  treeM3->SetBranchAddress("muonPEta", &mPEta);
  treeM3->SetBranchAddress("muonMEta", &mMEta);

  // cycle over MC, fill the var histograms
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);

      // fill the histos
      if(abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_pt > pT_i[0] && mc_pt < pT_f[0])  {
	h_m[0]->Fill(mc_m);
	if(mc_m > 3.0 && mc_m < 3.2) {
	  h_pT[0]->Fill(mc_pt);
	  h_y[0]->Fill(mc_y);
	  hP_pT[0]->Fill(mPPt);
	  hM_pT[0]->Fill(mMPt);
	  hP_eta[0]->Fill(mPEta);
	  hM_eta[0]->Fill(mMEta);
	}
      }

    }
  
  cout << "low-pT MC filled" << endl;

  for(int i = 0; i < m1aEvt; i++)
    {
      treeM1a->GetEntry(i);

      // fill the histos
      if(abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	if(mc_pt > pT_i[0] && mc_pt < pT_f[0])  {
	  h_m[1]->Fill(mc_m);
	  if(mc_m > 3.0 && mc_m < 3.2) {
	    h_pT[1]->Fill(mc_pt);
	    h_y[1]->Fill(mc_y);
	    hP_pT[1]->Fill(mPPt);
	    hM_pT[1]->Fill(mMPt);
	    hP_eta[1]->Fill(mPEta);
	    hM_eta[1]->Fill(mMEta);
	  }
	}
	
	if(mc_pt > pT_i[1] && mc_pt < pT_f[1])  {
	  h_m[2]->Fill(mc_m);
	  if(mc_m > 3.0 && mc_m < 3.2) {
	    
	    h_pT[2]->Fill(mc_pt);
	    h_y[2]->Fill(mc_y);
	    hP_pT[2]->Fill(mPPt);
	    hM_pT[2]->Fill(mMPt);
	    hP_eta[2]->Fill(mPEta);
	    hM_eta[2]->Fill(mMEta);
	  }
	}
      }
    }
	
  cout << "mid-pT MC filled" << endl;

  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);

      // fill the histos
      if(abs(mc_lt) < 0.005 && abs(mc_y) < 1.2) {
	if(mc_pt > pT_i[1] && mc_pt < pT_f[1])  {
	  h_m[3]->Fill(mc_m);
	  if(mc_m > 3.0 && mc_m < 3.2) {
	    h_pT[3]->Fill(mc_pt);
	    h_y[3]->Fill(mc_y);
	    hP_pT[3]->Fill(mPPt);
	    hM_pT[3]->Fill(mMPt);
	    hP_eta[3]->Fill(mPEta);
	    hM_eta[3]->Fill(mMEta);
	  }
	}
	
	if(mc_pt > pT_i[2] && mc_pt < pT_f[2])  {
	  h_m[4]->Fill(mc_m);
	  if(mc_m > 3.0 && mc_m < 3.2) {
	    
	    h_pT[4]->Fill(mc_pt);
	    h_y[4]->Fill(mc_y);
	    hP_pT[4]->Fill(mPPt);
	    hM_pT[4]->Fill(mMPt);
	    hP_eta[4]->Fill(mPEta);
	    hM_eta[4]->Fill(mMEta);
	  }
	}
      }
    }

  cout << "high-pT MC filled" << endl;

  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);

      // fill the histos
      if(abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_pt > pT_i[2] && mc_pt < pT_f[2])  {
	h_m[5]->Fill(mc_m);
	if(mc_m > 3.0 && mc_m < 3.2) {
	  h_pT[5]->Fill(mc_pt);
	  h_y[5]->Fill(mc_y);
	  hP_pT[5]->Fill(mPPt);
	  hM_pT[5]->Fill(mMPt);
	  hP_eta[5]->Fill(mPEta);
	  hM_eta[5]->Fill(mMEta);
	}
      }
    }

  cout << "highest-pT MC filled" << endl;
  
  fin2->Close();
  fin2a->Close();
  fin3->Close();
  fin4->Close();

  TFile *fout = new TFile("files/store_MCdists.root", "recreate");

  // store the dists
  string lbl_s[] = {"lowPtMC0", "midPtMC0",
		    "midPtMC1", "highpTMC1",
		    "highpTMC2", "vhighpTMC2"};
  for(int i = 0; i < 6; i++) {
    h_pT[i]->Write(Form("h_pT_%s", lbl_s[i].c_str()));
    h_y[i]->Write(Form("h_y_%s", lbl_s[i].c_str()));
    h_m[i]->Write(Form("h_m_%s", lbl_s[i].c_str()));
    hP_pT[i]->Write(Form("hP_pT_%s", lbl_s[i].c_str()));
    hM_pT[i]->Write(Form("hM_pT_%s", lbl_s[i].c_str()));
    hP_eta[i]->Write(Form("hP_eta_%s", lbl_s[i].c_str()));
    hM_eta[i]->Write(Form("hM_eta_%s", lbl_s[i].c_str()));
  }

  cout << "all histos stored" << endl;
  
  fout->Close();
}
