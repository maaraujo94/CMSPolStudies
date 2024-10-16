// macro to plot the single muon pT, eta dists of data and MC

void store_muDists()
{
  // PART 1 : creating the histograms

  // mu(+)
  // 8 pT dists: PRSR data + MC over 4 pT regions
  TH1D **hP_pT = new TH1D*[8]; 
  for(int i = 0; i < 8; i++)
    hP_pT[i] = new TH1D(Form("hP_pT%d", i), Form("p_{T}(#mu^{+}) distributions"), 65, 0, 130);

  // 8 y dists: (PRSR data + MC) over 4 pT regions
  TH1D **hP_eta = new TH1D*[8]; 
  for(int i = 0; i < 8; i++)
    hP_eta[i] = new TH1D(Form("hP_eta%d", i), Form("#eta(mu^{+}) distributions"), 60, -1.5, 1.5);

  // mu(-)
  // 8 pT dists: PRSR data + MC over 4 pT regions
  TH1D **hM_pT = new TH1D*[8]; 
  for(int i = 0; i < 8; i++)
    hM_pT[i] = new TH1D(Form("hM_pT%d", i), Form("p_{T}(#mu^{-}) distributions"), 65, 0, 130);

  // 8 y dists: (PRSR data + MC) over 4 pT regions
  TH1D **hM_eta = new TH1D*[8]; 
  for(int i = 0; i < 8; i++)
    hM_eta[i] = new TH1D(Form("hM_eta%d", i), Form("#eta(mu^{-}) distributions"), 60, -1.5, 1.5);

  // PART 2 : open files and read TTrees
  TFile *fin = new TFile("files/muData_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  TFile *fin2 = new TFile("files/muMC_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");
  TFile *fin2a = new TFile("files/muMCm_cos.root");
  TTree *treeM1a = (TTree*)fin2a->Get("MC_cos");
  TFile *fin3 = new TFile("files/muMCh_cos.root");
  TTree *treeM2 = (TTree*)fin3->Get("MC_cos");
  TFile *fin4 = new TFile("files/muMCvh_cos.root");
  TTree *treeM3 = (TTree*)fin4->Get("MC_cos");
  
  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  int m1aEvt = treeM1a->GetEntries();
  int m2Evt = treeM2->GetEntries();
  int m3Evt = treeM3->GetEntries();
  
  // definitions to store data and MC events
  Double_t data_pt, data_lt, data_m, data_y;
  Double_t mc_pt, mc_lt, mc_m, mc_y;
  double mMPt, mPPt, mMEta, mPEta;
  
  treeD->SetBranchAddress("dimPt", &data_pt);
  treeD->SetBranchAddress("Rap", &data_y);
  treeD->SetBranchAddress("Mass", &data_m);
  treeD->SetBranchAddress("lt", &data_lt);
  treeD->SetBranchAddress("muonPPt", &mPPt);
  treeD->SetBranchAddress("muonMPt", &mMPt);
  treeD->SetBranchAddress("muonPEta", &mPEta);
  treeD->SetBranchAddress("muonMEta", &mMEta);
  
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

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);

      // fill the histos - split by pT interval
      if(data_m > 3.0 && data_m < 3.2 && abs(data_lt) < 0.005 && abs(data_y) < 1.2 && data_pt > 25 && data_pt < 120) {
	if(data_pt < 45) {
	  hP_pT[0]->Fill(mPPt);
	  hM_pT[0]->Fill(mMPt);
	  hP_eta[0]->Fill(mPEta);
	  hM_eta[0]->Fill(mMEta);
	}
	if(data_pt > 45 && data_pt < 50) {
	  hP_pT[1]->Fill(mPPt);
	  hM_pT[1]->Fill(mMPt);
	  hP_eta[1]->Fill(mPEta);
	  hM_eta[1]->Fill(mMEta);
	}
	if(data_pt > 50 && data_pt < 70) {
	  hP_pT[2]->Fill(mPPt);
	  hM_pT[2]->Fill(mMPt);
	  hP_eta[2]->Fill(mPEta);
	  hM_eta[2]->Fill(mMEta);
	}
	if(data_pt > 70) {
	  hP_pT[3]->Fill(mPPt);
	  hM_pT[3]->Fill(mMPt);
	  hP_eta[3]->Fill(mPEta);
	  hM_eta[3]->Fill(mMEta);
	}
      }
    }

  cout << "data filled" << endl;
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);

      
      // fill the histos 
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_pt > 25 && mc_pt < 45) {
	
	hP_pT[4]->Fill(mPPt);
	hM_pT[4]->Fill(mMPt);
	hP_eta[4]->Fill(mPEta);
	hM_eta[4]->Fill(mMEta);
      }
    }
  
  cout << "low-pT MC filled" << endl;

  for(int i = 0; i < m1aEvt; i++)
    {
      treeM1a->GetEntry(i);

      // fill the histos 
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_pt > 45 && mc_pt < 50) {
	hP_pT[5]->Fill(mPPt);
	hM_pT[5]->Fill(mMPt);
	hP_eta[5]->Fill(mPEta);
	hM_eta[5]->Fill(mMEta);
      }
    }

  cout << "mid-pT MC filled" << endl;

  for(int i = 0; i < m2Evt; i++)
    {
      treeM2->GetEntry(i);

      // fill the histos 
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_pt > 50 && mc_pt < 70) {
	hP_pT[6]->Fill(mPPt);
	hM_pT[6]->Fill(mMPt);
	hP_eta[6]->Fill(mPEta);
	hM_eta[6]->Fill(mMEta);
      }
    }

  cout << "high-pT MC filled" << endl;

  for(int i = 0; i < m3Evt; i++)
    {
      treeM3->GetEntry(i);

      // fill the histos 
      if(mc_m > 3.0 && mc_m < 3.2 && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_pt > 70 && mc_pt < 120) {
	hP_pT[7]->Fill(mPPt);
	hM_pT[7]->Fill(mMPt);
	hP_eta[7]->Fill(mPEta);
	hM_eta[7]->Fill(mMEta);
      }
    }

  cout << "highest-pT MC filled" << endl;
  
  fin->Close();
  fin2->Close();
  fin2a->Close();
  fin3->Close();
  fin4->Close();

  TFile *fout = new TFile("files/store_muDists.root", "recreate");

  // store the dists
  string lbl_pT[] = {"lowPtData", "midPtData", "highPtData", "highestPtData", "lowPtMC", "midPtMC", "highPtMC", "highestPtMC"};
  for(int i = 0; i < 8; i++) {
    hP_pT[i]->Write(Form("hP_pT_%s", lbl_pT[i].c_str()));
    hM_pT[i]->Write(Form("hM_pT_%s", lbl_pT[i].c_str()));
    hP_eta[i]->Write(Form("hP_eta_%s", lbl_pT[i].c_str()));
    hM_eta[i]->Write(Form("hM_eta_%s", lbl_pT[i].c_str()));
  }

  cout << "all histos stored" << endl;
  
  fout->Close();
}
