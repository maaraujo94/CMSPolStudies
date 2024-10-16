// macro to run through both year datafiles and combine

void comb_mid()
{
  Double_t mmPt, th, phi, lt, lterr, mass, rap, dR;
  Double_t mPPt, mMPt, mPEta, mMEta;
  int tEvt, perc;
  
  TFile *fMC = new TFile("../MCmS_cosCS.root", "recreate");
  TTree *tMC = new TTree("MC_cos", "");
  
  tMC->Branch("theta", &th);
  tMC->Branch("phi", &phi);
  tMC->Branch("dimPt", &mmPt);
  tMC->Branch("lt", &lt);
  tMC->Branch("lterr", &lterr);
  tMC->Branch("Mass", &mass);
  tMC->Branch("Rap", &rap);
  tMC->Branch("muonPPt", &mPPt);
  tMC->Branch("muonPEta", &mPEta);
  tMC->Branch("muonMPt", &mMPt);
  tMC->Branch("muonMEta", &mMEta);
  tMC->Branch("DeltaR", &dR);
  
  // open files and read TTrees - 2017
  TFile *finM7 = new TFile("../MCm17_cosCS.root");
  TTree *treeM7 = (TTree*)finM7->Get("MC_cos");

  tEvt = treeM7->GetEntries();
  perc = tEvt / 100;
  treeM7->SetBranchAddress("dimPt", &mmPt);
  treeM7->SetBranchAddress("theta", &th);
  treeM7->SetBranchAddress("phi", &phi);
  treeM7->SetBranchAddress("lt", &lt);
  treeM7->SetBranchAddress("lterr", &lterr);
  treeM7->SetBranchAddress("Mass", &mass);
  treeM7->SetBranchAddress("Rap", &rap);
  treeM7->SetBranchAddress("muonPPt", &mPPt);
  treeM7->SetBranchAddress("muonPEta", &mPEta);
  treeM7->SetBranchAddress("muonMPt", &mMPt);
  treeM7->SetBranchAddress("muonMEta", &mMEta);
  treeM7->SetBranchAddress("DeltaR", &dR);

  for(int i = 0; i < tEvt; i++) {
    treeM7->GetEntry(i);
    tMC->Fill();
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2017 MC (mid-pT)" << endl; 
  }
  finM7->Close();

  // open files and read TTrees - 2018
  TFile *finM8 = new TFile("../MCm18_cosCS.root");
  TTree *treeM8 = (TTree*)finM8->Get("MC_cos");

  tEvt = treeM8->GetEntries();
  perc = tEvt / 100;
  treeM8->SetBranchAddress("dimPt", &mmPt);
  treeM8->SetBranchAddress("theta", &th);
  treeM8->SetBranchAddress("phi", &phi);
  treeM8->SetBranchAddress("lt", &lt);
  treeM8->SetBranchAddress("lterr", &lterr);
  treeM8->SetBranchAddress("Mass", &mass);
  treeM8->SetBranchAddress("Rap", &rap);
  treeM8->SetBranchAddress("muonPPt", &mPPt);
  treeM8->SetBranchAddress("muonPEta", &mPEta);
  treeM8->SetBranchAddress("muonMPt", &mMPt);
  treeM8->SetBranchAddress("muonMEta", &mMEta);
  treeM8->SetBranchAddress("DeltaR", &dR);

  for(int i = 0; i < tEvt; i++) {
    treeM8->GetEntry(i);
    tMC->Fill();
    if((i+1)%perc == 0) cout << (i+1)/perc << "% done with 2018 MC (mid-pT)" << endl; 
  }
  finM8->Close();

  fMC->Write();
  fMC->Close();
}
