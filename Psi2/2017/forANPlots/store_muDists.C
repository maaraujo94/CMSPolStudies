// macro to plot the single muon pT, eta dists of data and MC

void store_muDists()
{
  // PART 1 : creating the histograms

  // mu(+)
  // 2 pT dists: PRSR data + MC 
  TH1D **hP_pT = new TH1D*[2]; 
  for(int i = 0; i < 2; i++)
    hP_pT[i] = new TH1D(Form("hP_pT%d", i), Form("p_{T}(#mu^{+}) distributions"), 55, 0, 110);

  // 2 y dists: (PRSR data + MC) 
  TH1D **hP_eta = new TH1D*[2]; 
  for(int i = 0; i < 2; i++)
    hP_eta[i] = new TH1D(Form("hP_eta%d", i), Form("#eta(mu^{+}) distributions"), 60, -1.5, 1.5);

  // mu(-)
  // 2 pT dists: PRSR data + MC 
  TH1D **hM_pT = new TH1D*[2]; 
  for(int i = 0; i < 2; i++)
    hM_pT[i] = new TH1D(Form("hM_pT%d", i), Form("p_{T}(#mu^{-}) distributions"), 55, 0, 110);

  // 2 y dists: (PRSR data + MC) 
  TH1D **hM_eta = new TH1D*[2]; 
  for(int i = 0; i < 2; i++)
    hM_eta[i] = new TH1D(Form("hM_eta%d", i), Form("#eta(mu^{-}) distributions"), 60, -1.5, 1.5);

  double m_min[] = {3.4, 3.57, 3.82};
  double m_max[] = {3.52, 3.81, 4.0};

  // PART 2 : open files and read TTrees
  TFile *fin = new TFile("files/muData_cos.root");
  TTree *treeD = (TTree*)fin->Get("data_cos");
  //TFile *fin2 = new TFile("files/muMC_cos.root");
  TFile *fin2 = new TFile("../../Store_data_codes/MCm17_cos.root");
  TTree *treeM1 = (TTree*)fin2->Get("MC_cos");

  int dEvt = treeD->GetEntries();
  int m1Evt = treeM1->GetEntries();
  
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

  // cycle over data and MC, fill the costh histogram acc to binning
  for(int i = 0; i < dEvt; i++)
    {
      treeD->GetEntry(i);

      // fill the histos - split by pT interval
      if(data_m > m_min[1] && data_m < m_max[1] && abs(data_lt) < 0.005 && abs(data_y) < 1.2 && data_pt > 20 && data_pt < 100) {
	hP_pT[0]->Fill(mPPt);
	hM_pT[0]->Fill(mMPt);
	hP_eta[0]->Fill(mPEta);
	hM_eta[0]->Fill(mMEta);
      }
    }
  
  cout << "data filled" << endl;
  
  for(int i = 0; i < m1Evt; i++)
    {
      treeM1->GetEntry(i);
      
      // fill the histos 
      if(mc_m > m_min[1] && mc_m < m_max[1] && abs(mc_lt) < 0.005 && abs(mc_y) < 1.2 && mc_pt > 20 && mc_pt < 100) {
	
	hP_pT[1]->Fill(mPPt);
	hM_pT[1]->Fill(mMPt);
	hP_eta[1]->Fill(mPEta);
	hM_eta[1]->Fill(mMEta);
      }
    }
  
  cout << "MC filled" << endl;

  
  fin->Close();
  fin2->Close();

  TFile *fout = new TFile("files/store_muDists.root", "recreate");

  // store the dists
  string lbl_pT[] = {"Data", "MC"};
  for(int i = 0; i < 2; i++) {
    hP_pT[i]->Write(Form("hP_pT_%s", lbl_pT[i].c_str()));
    hM_pT[i]->Write(Form("hM_pT_%s", lbl_pT[i].c_str()));
    hP_eta[i]->Write(Form("hP_eta_%s", lbl_pT[i].c_str()));
    hM_eta[i]->Write(Form("hM_eta_%s", lbl_pT[i].c_str()));
  }

  cout << "all histos stored" << endl;
  
  fout->Close();
}
