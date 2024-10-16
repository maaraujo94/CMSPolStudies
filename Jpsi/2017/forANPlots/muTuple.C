// code that runs over 2017 data and MC to store cut variables
// current cuts: based on trigger only

void muTuple()
{
  Double_t mmPt, ct, ctErr, vProb,  lt, lterr, mass, rap, mPPt, mMPt, mPEta, mMEta;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
  UInt_t trigger;

  // data tree
  /*  TFile *finD = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/orig_files/Data/filtered-all-psi-UL17.root");
  TTree *treeD = (TTree*)finD->Get("jpsitree");
  
  // setting branch address  
  treeD->SetBranchAddress("muP_p4", &muP_p4);
  treeD->SetBranchAddress("muM_p4", &muM_p4);
  treeD->SetBranchAddress("mumu_p4", &mumu_p4);
  treeD->SetBranchAddress("vProb", &vProb);
  treeD->SetBranchAddress("trigger", &trigger);
  treeD->SetBranchAddress("ctpv", &ct);
  treeD->SetBranchAddress("ctpv_error", &ctErr);
  
  int dEvt = treeD->GetEntries();
  int perc = dEvt / 100;

  TFile *foutD = new TFile("files/muData_cos.root", "recreate");
  TTree *newtreeD = new TTree("data_cos", "");
  
  newtreeD->Branch("dimPt", &mmPt);
  newtreeD->Branch("lt", &lt);
  newtreeD->Branch("lterr", &lterr);
  newtreeD->Branch("Mass", &mass);
  newtreeD->Branch("Rap", &rap);
  newtreeD->Branch("muonPPt", &mPPt);
  newtreeD->Branch("muonPEta", &mPEta);
  newtreeD->Branch("muonMPt", &mMPt);
  newtreeD->Branch("muonMEta", &mMEta);
  
  for(int i = 0; i < dEvt; i++) {
    treeD->GetEntry(i);
    if( vProb > 0.01 &&	(trigger&16) == 16)
      {
	mPPt = muP_p4->Pt();
	mMPt = muM_p4->Pt();
	mPEta = muP_p4->Eta();
	mMEta = muM_p4->Eta();

	mmPt = mumu_p4->Pt();
	rap = mumu_p4->Rapidity();
	mass = mumu_p4->M();
	lt = ct;
	lterr = ctErr;
	
	newtreeD->Fill();
	
      }
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2017 data" << endl; 
  }

  foutD->Write();
  foutD->Close();
  finD->Close();*/
  
  // MC (low pT) tree
  Float_t ct_MC, ct_MCErr;
  TFile *finM1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/orig_files/MC/filtered-all-psi-mc-UL17-25_46v1.root");
  TTree *treeM1 = (TTree*)finM1->Get("jpsitree");
  
  treeM1->SetBranchAddress("muP_p4", &muP_p4);
  treeM1->SetBranchAddress("muM_p4", &muM_p4);
  treeM1->SetBranchAddress("mumu_p4", &mumu_p4);
  treeM1->SetBranchAddress("ctpv", &ct_MC);
  treeM1->SetBranchAddress("ctpv_error", &ct_MCErr);

  int mEvt = treeM1->GetEntries();
  // perc = mEvt / 100;
  int perc = mEvt / 100;

  TFile *foutM1 = new TFile("files/muMC_cos.root", "recreate");
  TTree *newtreeM1 = new TTree("MC_cos", "");

  newtreeM1->Branch("dimPt", &mmPt);
  newtreeM1->Branch("lt", &lt);
  newtreeM1->Branch("lterr", &lterr);
  newtreeM1->Branch("Mass", &mass);
  newtreeM1->Branch("Rap", &rap);
  newtreeM1->Branch("muonPPt", &mPPt);
  newtreeM1->Branch("muonPEta", &mPEta);
  newtreeM1->Branch("muonMPt", &mMPt);
  newtreeM1->Branch("muonMEta", &mMEta);

  for(int i = 0; i < mEvt; i++) {
    treeM1->GetEntry(i);

    mPPt = muP_p4->Pt();
    mMPt = muM_p4->Pt();
    mPEta = muP_p4->Eta();
    mMEta = muM_p4->Eta();

    mmPt = mumu_p4->Pt();
    rap = mumu_p4->Rapidity();
    mass = mumu_p4->M();
    lt = ct_MC;
    lterr = ct_MCErr;

    newtreeM1->Fill();
    
    if((i+1)%perc == 0) {
      cout << (i+1)/perc << "% done with low-pT MC" << endl;
    }
  }
  
  foutM1->Write();
  foutM1->Close();
  finM1->Close();

  // MC (mid pT) tree
  TFile *finM2 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/orig_files/MC/filtered-all-psi-mc-LOCAL17-40_52v1.root");
  TTree *treeM2 = (TTree*)finM2->Get("jpsitree");
  
  treeM2->SetBranchAddress("muP_p4", &muP_p4);
  treeM2->SetBranchAddress("muM_p4", &muM_p4);
  treeM2->SetBranchAddress("mumu_p4", &mumu_p4);
  treeM2->SetBranchAddress("ctpv", &ct_MC);
  treeM2->SetBranchAddress("ctpv_error", &ct_MCErr);

  mEvt = treeM2->GetEntries();
  perc = mEvt / 100;

  TFile *foutM2 = new TFile("files/muMCm_cos.root", "recreate");
  TTree *newtreeM2 = new TTree("MC_cos", "");

  newtreeM2->Branch("dimPt", &mmPt);
  newtreeM2->Branch("lt", &lt);
  newtreeM2->Branch("lterr", &lterr);
  newtreeM2->Branch("Mass", &mass);
  newtreeM2->Branch("Rap", &rap);
  newtreeM2->Branch("muonPPt", &mPPt);
  newtreeM2->Branch("muonPEta", &mPEta);
  newtreeM2->Branch("muonMPt", &mMPt);
  newtreeM2->Branch("muonMEta", &mMEta);

  for(int i = 0; i < mEvt; i++) {
    treeM2->GetEntry(i);

    mPPt = muP_p4->Pt();
    mMPt = muM_p4->Pt();
    mPEta = muP_p4->Eta();
    mMEta = muM_p4->Eta();

    mmPt = mumu_p4->Pt();
    rap = mumu_p4->Rapidity();
    mass = mumu_p4->M();
    lt = ct_MC;
    lterr = ct_MCErr;

    newtreeM2->Fill();
    
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with mid-pT MC" << endl; 
  }
  
  foutM2->Write();
  foutM2->Close();
  finM2->Close();
  
  // MC (high pT) tree
  TFile *finM3 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/orig_files/MC/filtered-all-psi-mc-LOCAL17-highpt.root");
  TTree *treeM3 = (TTree*)finM3->Get("jpsitree");
  
  treeM3->SetBranchAddress("muP_p4", &muP_p4);
  treeM3->SetBranchAddress("muM_p4", &muM_p4);
  treeM3->SetBranchAddress("mumu_p4", &mumu_p4);
  treeM3->SetBranchAddress("ctpv", &ct_MC);
  treeM3->SetBranchAddress("ctpv_error", &ct_MCErr);

  mEvt = treeM3->GetEntries();
  perc = mEvt / 100;

  TFile *foutM3 = new TFile("files/muMCh_cos.root", "recreate");
  TTree *newtreeM3 = new TTree("MC_cos", "");

  newtreeM3->Branch("dimPt", &mmPt);
  newtreeM3->Branch("lt", &lt);
  newtreeM3->Branch("lterr", &lterr);
  newtreeM3->Branch("Mass", &mass);
  newtreeM3->Branch("Rap", &rap);
  newtreeM3->Branch("muonPPt", &mPPt);
  newtreeM3->Branch("muonPEta", &mPEta);
  newtreeM3->Branch("muonMPt", &mMPt);
  newtreeM3->Branch("muonMEta", &mMEta);

  for(int i = 0; i < mEvt; i++) {
    treeM3->GetEntry(i);

    mPPt = muP_p4->Pt();
    mMPt = muM_p4->Pt();
    mPEta = muP_p4->Eta();
    mMEta = muM_p4->Eta();

    mmPt = mumu_p4->Pt();
    rap = mumu_p4->Rapidity();
    mass = mumu_p4->M();
    lt = ct_MC;
    lterr = ct_MCErr;

    newtreeM3->Fill();
    
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with high-pT MC" << endl; 
  }
  
  foutM3->Write();
  foutM3->Close();
  finM3->Close();

  // MC (very high pT) tree
  TFile *finM4 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Jpsi/Store_data_codes/orig_files/MC/filtered-all-psi-mc-LOCAL17-veryhighpt.root");
  TTree *treeM4 = (TTree*)finM4->Get("jpsitree");
  
  treeM4->SetBranchAddress("muP_p4", &muP_p4);
  treeM4->SetBranchAddress("muM_p4", &muM_p4);
  treeM4->SetBranchAddress("mumu_p4", &mumu_p4);
  treeM4->SetBranchAddress("ctpv", &ct_MC);
  treeM4->SetBranchAddress("ctpv_error", &ct_MCErr);

  mEvt = treeM4->GetEntries();
  perc = mEvt / 100;

  TFile *foutM4 = new TFile("files/muMCvh_cos.root", "recreate");
  TTree *newtreeM4 = new TTree("MC_cos", "");

  newtreeM4->Branch("dimPt", &mmPt);
  newtreeM4->Branch("lt", &lt);
  newtreeM4->Branch("lterr", &lterr);
  newtreeM4->Branch("Mass", &mass);
  newtreeM4->Branch("Rap", &rap);
  newtreeM4->Branch("muonPPt", &mPPt);
  newtreeM4->Branch("muonPEta", &mPEta);
  newtreeM4->Branch("muonMPt", &mMPt);
  newtreeM4->Branch("muonMEta", &mMEta);

  for(int i = 0; i < mEvt; i++) {
    treeM4->GetEntry(i);

    mPPt = muP_p4->Pt();
    mMPt = muM_p4->Pt();
    mPEta = muP_p4->Eta();
    mMEta = muM_p4->Eta();

    mmPt = mumu_p4->Pt();
    rap = mumu_p4->Rapidity();
    mass = mumu_p4->M();
    lt = ct_MC;
    lterr = ct_MCErr;

    newtreeM4->Fill();
    
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with very-high-pT MC" << endl; 
  }
  
  foutM4->Write();
  foutM4->Close();
  finM4->Close();

}
