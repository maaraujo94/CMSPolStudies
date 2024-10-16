// code that runs over 2018 data and MC to store cut variables
// current cuts: based on trigger only

void muTuple()
{
  Double_t mmPt, ct, ctErr, vProb,  lt, lterr, mass, rap, mPPt, mMPt, mPEta, mMEta;
  TLorentzVector *mumu_p4 = 0, *muM_p4 = 0, *muP_p4 = 0;
  UInt_t trigger;

  // data tree
  TFile *finD = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/orig_files/filtered-all-psi2s-UL18.root");
  TTree *treeD = (TTree*)finD->Get("psi2stree");
  
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
    if( vProb > 0.01 &&	(trigger&8) == 8)
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
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2018 data" << endl; 
  }

  foutD->Write();
  foutD->Close();
  finD->Close();
  
  // MC  tree
  Float_t ct_MC, ct_MCErr;
  TFile *finM1 = new TFile("/home/mariana/Documents/2020_PhD_work/CERN/CMSPolStudies/Psi2/Store_data_codes/orig_files/filtered-all-psi2s-mc-LOCAL18-midpt-v1.root");
  TTree *treeM1 = (TTree*)finM1->Get("psi2stree");
  
  treeM1->SetBranchAddress("muP_p4", &muP_p4);
  treeM1->SetBranchAddress("muM_p4", &muM_p4);
  treeM1->SetBranchAddress("mumu_p4", &mumu_p4);
  treeM1->SetBranchAddress("ctpv", &ct_MC);
  treeM1->SetBranchAddress("ctpv_error", &ct_MCErr);

  int mEvt = treeM1->GetEntries();
  perc = mEvt / 100;
  
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
      cout << (i+1)/perc << "% done with MC" << endl;
    }
  }
  
  foutM1->Write();
  foutM1->Close();
  finM1->Close();

}
